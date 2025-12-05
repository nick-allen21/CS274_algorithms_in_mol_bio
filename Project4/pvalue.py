import sys
import argparse
import random
import pandas as pd
from chemoUtils import load_data

class PValue:
    def __init__(self, drugs_csv, targets_csv, protein_a, protein_b):
        self.drugs_csv = drugs_csv
        self.targets_csv = targets_csv
        self.protein_a = protein_a
        self.protein_b = protein_b
        self.drugs_df, self.targets_df, self.drug_mapping_dict, _ = load_data(drugs_csv, targets_csv)

    def get_ligand_sets(self, protein_a, protein_b):
        self.ligand_set_a = {key for key, value in self.drug_mapping_dict.items() if protein_a in value}
        self.ligand_set_b = {key for key, value in self.drug_mapping_dict.items() if protein_b in value}
        self.size_a = len(self.ligand_set_a)
        self.size_b = len(self.ligand_set_b)

    def calculate_pvalue(self):
        for key, value in self.drug_mapping_dict.items():
            print(key, value)

    # add tanimoto score calcualtion to the chemutils functiosn
    def calculate_tsummary(self):

    def run_pvalue(self):
        self.calculate_pvalue()
        self.get_ligand_sets(self.protein_a, self.protein_b)

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, default=500)
    parser.add_argument('-r', type=int, default=214)
    parser.add_argument('drugs_csv', type=str)
    parser.add_argument('targets_csv', type=str)
    parser.add_argument('protein_a', type=str)
    parser.add_argument('protein_b', type=str)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    random.seed(args.r)

    pvalue = PValue(args.drugs_csv, args.targets_csv, args.protein_a, args.protein_b)
    pvalue.run_pvalue()

