import sys
import pandas as pd
import math
import numpy as np
from tqdm import tqdm

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
        patient_mean = exp[patient_cols].mean(axis=1)
        control_mean = exp[control_cols].mean(axis=1)
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
        
    def run_gsea(self, expfile, sampfile, keggfile, threshold):
        """
        should run get_gene_rank_order, get_enrichment_score and get_sig_sets, and then return the list of significant gene sets.
        """
        self.load_data(expfile, sampfile, keggfile)
        all_gene_sets = list(self.kegg_sets.keys())

        # loop through all gene sets and calculate the enrichment score and background distribution
        for gene_set in tqdm(all_gene_sets, desc="Running GSEA for all gene sets", total=len(all_gene_sets)):
            self.enrichment_scores[gene_set] = {}
            enrichment_score = self.get_enrichment_score(gene_set)
            self.enrichment_scores[gene_set]['enrichment_score'] = enrichment_score
            scores = self.background_distribution(gene_set)
            self.enrichment_scores[gene_set]['scores'] = scores
        
        # get the significant gene sets at a corrected threshold of threshodl
        return self.get_sig_sets(threshold)

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
    threshold = 0.5

    # run GSEA methd
    gsea = GSEA()
    significant_sets = gsea.run_gsea(expfile, sampfile, keggfile, threshold)
    print(significant_sets, sep='\t')


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python gsea.py <expfile> <sampfile> <keggfile>")
        sys.exit(1)
    main()