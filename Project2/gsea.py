import sys
import pandas as pd
import math
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

class GSEA:
    """
    Class to perform GSEA analysis.
    Attributes:
        expfile: file path to the expression data
        sampfile: file path to the sample data
        keggfile: file path to the gene set data
        exp_df: dataframe containing the expression data
        samp_df: dataframe containing the sample data
        kegg_sets: dictionary containing the gene sets
        enrichment_scores: dictionary containing the enrichment scores

    Methods:
        load_data: load the data from the files
        get_gene_rank_order: get the gene rank order
        get_enrichment_score: get the enrichment score
        background_distribution: create the background distribution
        get_sig_sets: get the significant gene sets
        run_gsea: run the GSEA analysis
    """
    def __init__(self):
        self.expfile = None
        self.sampfile = None
        self.keggfile = None
        self.exp_df = None
        self.samp_df = None
        self.kegg_sets = None 
        self.enrichment_scores = {}

    def load_data(self, expfile, sampfile, genesets):
        """
        Load the data from the files.
        Should take the file paths to the expression, sample and gene set data, read them in and store within the GSEA instance.
        The expression file is a tab-separated file with the first column being the gene names and the rest of the columns being the expression values
        """
        self.exp_df = pd.read_csv(expfile, sep='\t')
        # normalize gene symbols to uppercase/stripped for consistent overlap
        if 'SYMBOL' in self.exp_df.columns:
            self.exp_df['SYMBOL'] = (
                self.exp_df['SYMBOL'].astype(str).str.strip().str.upper()
            )

        # if the sample file has no header; enforce names and drop blank lines
        self.samp_df = pd.read_csv(sampfile, sep='\t', header=None, names=['sample', 'label']).dropna()

        # parse GMT manually, dont know how many columns per row
        kegg_sets = {}
        with open(genesets, 'r') as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue  # skip malformed lines
                set_name = parts[0]
                genes = [g.strip().upper() for g in parts[2:] if g]
                kegg_sets[set_name] = set(genes)
        self.kegg_sets = kegg_sets

    def get_fc_for_gene(self, gene):
        """
        Return the log fold-change (patients - controls) for a single gene.
        If the gene is not present in the expression table, return None.
        """
        # 1) Identify patient/control columns 
        samp = self.samp_df.copy()
        samp['label'] = samp['label'].astype(int)
        patient_cols = [s for s in samp.loc[samp['label'] == 1, 'sample'] if s in self.exp_df.columns]
        control_cols = [s for s in samp.loc[samp['label'] == 0, 'sample'] if s in self.exp_df.columns]

        # 2) Index expression by gene symbol
        exp = self.exp_df.copy()
        gene_col = 'SYMBOL' if 'SYMBOL' in exp.columns else exp.columns[0]
        exp = exp.set_index(gene_col)

        # 3) Vectorized logFC for all genes, then select the one requested
        patient_mean = exp[patient_cols].mean(axis=1)
        control_mean = exp[control_cols].mean(axis=1)
        logfc = (patient_mean - control_mean)

        # 4) Normalize the lookup key to match symbols
        key = str(gene).strip().upper()
        return float(logfc.loc[key]) if key in logfc.index else None


    def get_gene_rank_order(self):
        """
        should return a list of all genes (as strings) ranked by their logFC between patient and control, 
        with the gene with the highest logFC ordered first.
        """
        # get patient and control samples
        samp = self.samp_df.copy()
        samp['label'] = samp['label'].astype(int)
        patient_cols = [s for s in samp.loc[samp['label'] == 1, 'sample'] if s in self.exp_df.columns]
        control_cols = [s for s in samp.loc[samp['label'] == 0, 'sample'] if s in self.exp_df.columns]
        # Choose gene identifier column and compute logFC
        exp = self.exp_df.copy()
        
        gene_col = 'SYMBOL' if 'SYMBOL' in exp.columns else exp.columns[0]
        exp = exp.set_index(gene_col)

        # compute the mean expression for the patient and control samples
        # EFFICIENCY FIX:
        patient_mean = exp[patient_cols].mean(axis=1)
        control_mean = exp[control_cols].mean(axis=1)
        # compute FC leverage pandas dataframe operations as opposed to looping through each gene
        logfc = patient_mean - control_mean

        # rank the genes by their logFC
        ranked_genes = logfc.sort_values(ascending=False).index.astype(str).tolist()
        return ranked_genes

    def get_enrichment_score(self, geneset):
        """
        should return the enrichment score, a float correct to two decimal places for a given gene set, 
        such as ‘KEGG_CITRATE_CYCLE_TCA_CYCLE’ (which is the string for the gene set name corresponding to the gene set).
        This method should run get_gene_rank_order at some point to initiate enrichment calculations.
        """

        # S is the set of genes (normalized to uppercase during load)
        S = self.kegg_sets[geneset]
        # L is ranked list of all genes by degree of correlation with differential expression 
        L = self.get_gene_rank_order()
        # N is the number of genes in L 
        N = len(L)

        # S' used in the ES walk, only the overlapping genes are considered
        S_prime = S & set(L)
        G = len(S_prime)
        if G == 0 or G == N:
            return 0.0

        # 1) walk down from L computing sum at each step 
        SUM = 0
        scores = []
        for gene in L: 
            if gene in S_prime:
                # sqrt weighting uses the size of geneset that actually overlaps with the ranked list
                score = math.sqrt((N - G) / G)
            else:
                # negative sqrt of (G / (N - G))
                score = -math.sqrt(G / (N - G))
            SUM += score
            scores.append(SUM)
        
        # find the maximum enrichment score
        supremum = max(scores)
        return round(supremum, 2)

    def background_distribution(self, geneset):
        """
        Create a background distribution of enrichment scores by permuting the sample labels
        and re-calculating the ranked gene list 100 times 
        and calculate a p-value for each gene set by counting the number of times in the 
        permuted iterations a gene set had an equal or higher score than its actual score 
        over the total number of iterations.
        """
        # make copy of original sample_df to reset later 
        sample_df_orig= self.samp_df.copy()

        # randomize the sample_df that is accessed by the get_gene_rank_order method
        scores = []
        for i in range(100):
            # randomize the labels
            self.samp_df['label'] = np.random.randint(0, 2, len(self.samp_df))
            # get the ranked gene list
            enrichment_score = self.get_enrichment_score(geneset)
            scores.append(enrichment_score)
    
        self.samp_df = sample_df_orig
        return scores

    def get_sig_sets(self, p):
        """
        Should return the list of significant gene sets (as strings), at a corrected threshold of p, by name. 
        If no gene sets are significant, return an empty list. This method should run get_gene_rank_order and/or get_enrichment_score 
        at some point to initiate enrichment calculations and then identify significant gene sets. 
        Make sure you write an efficient implementation of this method, 
        as inefficient solutions will time out on the autograder.
        """
        significant_sets = []
        num_tests = len(self.enrichment_scores) if self.enrichment_scores is not None else 0
        if num_tests == 0:
            return significant_sets

        # calculate the p-value for each gene set using the distribution of enrichment scores and bonferroni correction
        for gene_set, data in self.enrichment_scores.items():
            es_obs = data['enrichment_score']
            perm_scores = data.get('scores', [])
            n = len(perm_scores)
            if n == 0:
                continue

            # smoothed permutation p-value
            k = sum(1 for s in perm_scores if s >= es_obs)
            p_raw = (k + 1) / (n + 1)

            # Bonferroni correction
            p_adj = min(1.0, p_raw * num_tests)
            if p_adj <= p:
                significant_sets.append(gene_set)

        return significant_sets

    def count_unique_kegg_genes(self) -> int:
        """
        Return the number of unique genes across all KEGG sets.
        Requires load_data(...) to have been called.
        """
        if self.kegg_sets is None:
            raise ValueError("Data not loaded. Call load_data(expfile, sampfile, keggfile) first.")
        unique = set()
        for genes in self.kegg_sets.values():
            unique.update(genes)
        return len(unique)

    def plot_kegg_gene_occurrence_histogram(self, out_file: str = "kegg_gene_occurrences_hist.png", show: bool = False) -> None:
        """
        Plot a histogram where the x-axis is the number of KEGG sets a gene appears in,
        and the y-axis is the number of genes with that occurrence count.
        Saves the figure to 'out_file' by default; set show=True to display.
        """
        if self.kegg_sets is None:
            raise ValueError("Data not loaded. Call load_data(expfile, sampfile, keggfile) first.")

        # Count occurrences per gene across all sets
        from collections import Counter
        gene_counts = Counter()
        for genes in self.kegg_sets.values():
            gene_counts.update(genes)

        # Frequency of occurrence counts (how many genes have count c)
        count_freq = Counter(gene_counts.values())
        xs = sorted(count_freq.keys())
        ys = [count_freq[x] for x in xs]

        plt.figure(figsize=(8, 5))
        plt.bar(xs, ys, color="steelblue")
        plt.xlabel("number of occurrences")
        plt.ylabel("number of genes")
        plt.title("Gene occurrences across KEGG sets")
        # Limit tick clutter for large ranges
        if len(xs) <= 30:
            plt.xticks(xs)
        plt.tight_layout()
        if out_file:
            plt.savefig(out_file, dpi=150)
        if show:
            plt.show()
        plt.close()

    def most_frequent_kegg_genes(self):
        """
        Return (genes, count) where genes is a sorted list of gene symbols that
        appear in the largest number of KEGG sets, and count is that maximum
        occurrence count. Requires load_data to have been called.
        """
        if self.kegg_sets is None:
            raise ValueError("Data not loaded. Call load_data(expfile, sampfile, keggfile) first.")
        from collections import Counter
        gene_counts = Counter()
        for genes in self.kegg_sets.values():
            gene_counts.update(genes)
        if not gene_counts:
            return [], 0
        max_count = max(gene_counts.values())
        top_genes = sorted([g for g, c in gene_counts.items() if c == max_count])
        return top_genes, max_count
        
    def run_gsea(self, expfile, sampfile, keggfile, threshold):
        """
        should run get_gene_rank_order, get_enrichment_score and get_sig_sets, and then return the list of significant gene sets.
        """
        self.load_data(expfile, sampfile, keggfile)
        all_gene_sets = list(self.kegg_sets.keys())

        # plot the histogram of gene occurrences across KEGG sets
        self.plot_kegg_gene_occurrence_histogram(out_file="kegg_gene_occurrences_hist.png", show=False)

        # print the gene(s) that appear most often across KEGG sets
        top_genes, max_count = self.most_frequent_kegg_genes()
        print(f"Most frequent KEGG gene occurrences: {max_count}")
        print(f"Gene(s): {', '.join(top_genes)}")

        # get BMP4 logFC for the quiz 
        bmp4_logfc = self.get_fc_for_gene('BMP4')
        print(f"BMP4 logFC: {bmp4_logfc}")

        # get the number of unique genes across all KEGG sets
        num_unique_genes = self.count_unique_kegg_genes()
        print(f"Number of unique genes across all KEGG sets: {num_unique_genes}")

        max_enrichment_score = 0
        max_enrichment_score_gene_set = None
        # loop through all gene sets and calculate the enrichment score and background distribution
        for gene_set in tqdm(all_gene_sets, desc="Running GSEA for all gene sets", total=len(all_gene_sets)):
            self.enrichment_scores[gene_set] = {}
            enrichment_score = self.get_enrichment_score(gene_set)
            print(f"Enrichment score for {gene_set}: {enrichment_score}")
            if enrichment_score > max_enrichment_score:
                max_enrichment_score = enrichment_score
                max_enrichment_score_gene_set = gene_set
            self.enrichment_scores[gene_set]['enrichment_score'] = enrichment_score
            scores = self.background_distribution(gene_set)
            self.enrichment_scores[gene_set]['scores'] = scores
        
        # get the significant gene sets at a corrected threshold of threshodl
        significant_sets = self.get_sig_sets(threshold)
        print(f"Max enrichment score: {max_enrichment_score} for gene set: {max_enrichment_score_gene_set}")
        return significant_sets


def main():


    """
    Main method to run the GSEA analysis.
    Loads in the expression, sample and gene set data from the command line arguments.
    Runs the GSEA analysis with a given threshold
    Outputs the significant gene sets to a file.
    """

    if len(sys.argv) != 4:
        print("Usage: python gsea.py <expfile> <sampfile> <keggfile>")
        sys.exit(1) 
    
    # read in command line arguments 
    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    keggfile = sys.argv[3]

    # set inclusion threshold for significant gene sets
    threshold = 0.2

    

    # run GSEA methd
    gsea = GSEA()
    significant_sets = gsea.run_gsea(expfile, sampfile, keggfile, threshold)
    print(significant_sets, sep='\t')
    print(f"Number of significant gene sets: {len(significant_sets)}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python gsea.py <expfile> <sampfile> <keggfile>")
        sys.exit(1)
    main()