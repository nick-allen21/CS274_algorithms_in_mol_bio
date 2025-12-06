"""
Bootstrap P-value Calculator for Protein Ligand Set Similarity

This program calculates a bootstrap p-value for comparing the ligand sets
of two proteins using Tanimoto coefficient similarity.

Author: Nicholas Allen
Course: CS274 - Algorithms in Molecular Biology
"""

import sys
import argparse
import random
import pandas as pd
from chemoUtils import load_data, tanimoto_coeff


class PValue:
    """
    Class to calculate bootstrap p-values for protein ligand set comparisons.
    """
    
    def __init__(self, drugs_csv, targets_csv):
        """
        Initialize PValue calculator with data files and protein IDs.
        
        Inputs:
            drugs_csv: string, path to CSV file containing drug fingerprints
            targets_csv: string, path to CSV file containing drug-protein targets
            protein_a: string, uniprot_accession ID for first protein
            protein_b: string, uniprot_accession ID for second protein
        
        Returns:
            None
        """
        self.drugs_csv = drugs_csv
        self.targets_csv = targets_csv
        self.protein_a = None
        self.protein_b = None
        # Load data and create lookup dictionaries
        self.drugs_df, self.targets_df, self.drug_mapping_dict, self.fingerprint_dict = load_data(drugs_csv, targets_csv)
        
    def get_ligand_sets(self):
        # Get ligand sets for each protein (drugs that bind to each protein)
        self.ligand_set_a = {key for key, value in self.drug_mapping_dict.items() if self.protein_a in value}
        self.ligand_set_b = {key for key, value in self.drug_mapping_dict.items() if self.protein_b in value}
        
        # Store sizes for bootstrap sampling
        self.size_a = len(self.ligand_set_a)
        self.size_b = len(self.ligand_set_b)

    def calculate_tsummary(self, drug_set_a, drug_set_b, fingerprint_dict):
        """
        Calculate Tsummary: sum of pairwise Tanimoto coefficients above cutoff.
        
        Only includes Tanimoto coefficients > 0.5 in the sum, as per the
        Similarity Ensemble Approach (SEA) method.
        
        Inputs:
            drug_set_a: list or set of drug IDs for first set
            drug_set_b: list or set of drug IDs for second set
            fingerprint_dict: dictionary mapping drug IDs to fingerprint sets
        
        Returns:
            float, sum of all pairwise Tanimoto coefficients > 0.5
        """
        tsummary = 0
        # Calculate all pairwise Tanimoto coefficients
        for drug_a in drug_set_a:
            for drug_b in drug_set_b:
                tc = tanimoto_coeff(drug_a, drug_b, fingerprint_dict)
                # Only sum values above cutoff threshold
                if tc > 0.5:
                    tsummary += tc
        return tsummary

    def run_pvalue(self, n_iterations, random_seed, protein_a, protein_b):
        """
        Calculate bootstrap p-value for protein ligand set similarity.
        
        Uses bootstrap resampling to determine if the Tsummary for the two
        protein ligand sets is significantly higher than random expectation.
        
        Algorithm:
        1. Calculate Tsummary for actual protein pair
        2. For each bootstrap iteration:
           - Sample random ligand sets of same sizes
           - Calculate Ti_summary for random sets
           - Count if Ti_summary >= Tsummary
        3. p-value = count / n_iterations
        
        Inputs:
            n_iterations: int, number of bootstrap sampling iterations
            random_seed: int, random seed for reproducibility
        Returns:
            None (stores p-value in self.p_value)
        """
        # Set random seed for reproducibility
        random.seed(random_seed)

        # get ligand sets for the proteins
        self.protein_a = protein_a
        self.protein_b = protein_b
        # get ligand sets for the proteins
        self.get_ligand_sets()

        # Calculate Tsummary for actual protein ligand sets
        tsummary = self.calculate_tsummary(self.ligand_set_a, self.ligand_set_b, self.fingerprint_dict)
        
        # Get all ligands for bootstrap sampling
        all_ligands = list(self.fingerprint_dict.keys())
        
        # Bootstrap sampling iterations
        count_greater_equal = 0
        for i in range(n_iterations):
            # Sample with replacement: na random ligands for set A
            random_set_a = random.choices(all_ligands, k=self.size_a)
            # Sample with replacement: nb random ligands for set B
            random_set_b = random.choices(all_ligands, k=self.size_b)
            
            # Calculate Ti_summary for random sets
            ti_summary = self.calculate_tsummary(random_set_a, random_set_b, self.fingerprint_dict)
            
            # Count if Ti_summary >= Tsummary (indicates random is as good or better)
            if ti_summary >= tsummary:
                count_greater_equal += 1
        
        # Calculate bootstrap p-value
        p_bootstrap = count_greater_equal / n_iterations
        self.p_value = p_bootstrap


def parse_arguments():
    """
    Parse command-line arguments for pvalue calculation.
    
    Inputs:
        None 
    
    Returns:
        argparse.Namespace containing parsed arguments:
            - n: int, number of bootstrap iterations (default 500)
            - r: int, random seed (default 214)
            - drugs_csv: string, path to drugs CSV file
            - targets_csv: string, path to targets CSV file
            - protein_a: string, first protein uniprot_accession ID
            - protein_b: string, second protein uniprot_accession ID
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, default=500)
    parser.add_argument('-r', type=int, default=214)
    parser.add_argument('drugs_csv', type=str)
    parser.add_argument('targets_csv', type=str)
    parser.add_argument('protein_a', type=str)
    parser.add_argument('protein_b', type=str)
    return parser.parse_args()

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Create PValue calculator and run bootstrap analysis
    pvalue_calc = PValue(args.drugs_csv, args.targets_csv)
    pvalue_calc.run_pvalue(args.n, args.r, args.protein_a, args.protein_b)
    
    # Print the p-value
    print(pvalue_calc.p_value)

