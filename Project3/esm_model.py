import pandas as pd
import sys

class ESM_model:
    def __init__(self, fasta_file, output_dir, pro_fam_datafile):
        self.fasta_file = fasta_file
        self.output_dir = output_dir
        self.prot_fam_datafile = prot_fam_datafile
        self.prot_fam_data = pd.read_csv(prot_fam_datafile, sep="\t")

    def run(self):
        pass 

    def predict(self):
        pass 

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_dir = sys.argv[2]
    prot_fam_datafile = sys.argv[3]



