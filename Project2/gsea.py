import sys
import pandas as pd

class GSEA:
    def __init__(self, expfile, sampfile, keggfile):
        self.expfile = expfile
        self.sampfile = sampfile
        self.keggfile = keggfile
        self.exp_df = None
        self.samp_df = None
        self.kegg_df = None
        self.kegg_sets = None  # dict: set_name -> set(list_of_genes)

    def load_data(self, expfile, sampfile, genesets):
        # should take the file paths to the expression, sample and gene set data, read them in and store within the GSEA instance.
        # the expression file is a tab-separated file with the first column being the gene names and the rest of the columns being the expression values
        self.exp_df = pd.read_csv(expfile, sep='\t')
        # sample file has no header; enforce names and drop blank lines
        self.samp_df = pd.read_csv(sampfile, sep='\t', header=None, names=['sample', 'label']).dropna()
        # parse GMT manually (variable columns per row)
        kegg_sets = {}
        with open(genesets, 'r') as fh:
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue  # skip malformed lines
                set_name = parts[0]
                genes = [g for g in parts[2:] if g]
                kegg_sets[set_name] = set(genes)
        self.kegg_sets = kegg_sets

    def get_gene_rank_order(self):
        """
        should return a list of all genes (as strings) ranked by their logFC between patient and control, 
        with the gene with the highest logFC ordered first.
        """
        # Harmonize sample columns
        samp = self.samp_df.copy()
        samp['label'] = samp['label'].astype(int)
        patient_cols = [s for s in samp.loc[samp['label'] == 1, 'sample'] if s in self.exp_df.columns]
        control_cols = [s for s in samp.loc[samp['label'] == 0, 'sample'] if s in self.exp_df.columns]
        # Choose gene identifier column and compute logFC
        exp = self.exp_df.copy()

        # print(exp.columns)
        
        gene_col = 'SYMBOL' if 'SYMBOL' in exp.columns else exp.columns[0]
        exp = exp.set_index(gene_col)

        print(exp.head())

        patient_mean = exp[patient_cols].mean(axis=1)
        control_mean = exp[control_cols].mean(axis=1)
        logfc = patient_mean - control_mean
        print(logfc.head())
        ranked_genes = logfc.sort_values(ascending=False).index.astype(str).tolist()
        print(ranked_genes[:10])
        return ranked_genes


    def get_enrichment_score(geneset):
        """
        should return the enrichment score, a float correct to two decimal places for a given gene set, 
        such as ‘KEGG_CITRATE_CYCLE_TCA_CYCLE’ (which is the string for the gene set name corresponding to the gene set).
        This method should run get_gene_rank_order at some point to initiate enrichment calculations.
        """
    
    def get_sig_sets(p):
        """
        Should return the list of significant gene sets (as strings), at a corrected threshold of p, by name. 
        If no gene sets are significant, return an empty list. This method should run get_gene_rank_order and/or get_enrichment_score 
        at some point to initiate enrichment calculations and then identify significant gene sets. 
        Make sure you write an efficient implementation of this method, 
        as inefficient solutions will time out on the autograder.
        """
        
    def run_gsea(self, expfile, sampfile, keggfile):
        """
        should run get_gene_rank_order, get_enrichment_score and get_sig_sets, and then return the list of significant gene sets.
        """
        self.load_data(expfile, sampfile, keggfile)
        ranked_genes =self.get_gene_rank_order()

def main():

    if len(sys.argv) != 4:
        print("Usage: python gsea.py <expfile> <sampfile> <keggfile>")
        sys.exit(1)
    
    # read in command line arguments 
    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    keggfile = sys.argv[3]

    # run GSEA methd
    gsea = GSEA(expfile, sampfile, keggfile)
    gsea.run_gsea(expfile, sampfile, keggfile)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python gsea.py <expfile> <sampfile> <keggfile>")
        sys.exit(1)
    main()