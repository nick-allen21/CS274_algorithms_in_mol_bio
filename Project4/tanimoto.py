"""
Tanimoto Coefficient Calculator for Drug Comparison

This program calculates pairwise Tanimoto similarity coefficients for drugs
based on their molecular fingerprints. It also determines if drug pairs share
any protein targets.

Author: Nicholas Allen
Course: CS274 - Algorithms in Molecular Biology
"""

import sys
import pandas as pd
import os

class Tanimoto:
    """
    Class to calculate Tanimoto coefficients between drug pairs and identify
    shared protein targets.
    """
    
    def __init__(self, drugs_data_path, targets_data_path, output_path):
        """
        Initialize Tanimoto calculator with input and output file paths.
        
        Inputs:
            drugs_data_path: string, path to CSV file containing drug fingerprints
            targets_data_path: string, path to CSV file containing drug-protein targets
            output_path: string, path where output CSV will be written
        
        Returns:
            None
        """
        self.drugs_data_path = drugs_data_path
        self.targets_data_path = targets_data_path
        self.output_path = output_path

    def load_data(self):
        """
        Load drug and target data from CSV files and create lookup dictionaries.
        
        Creates two dictionaries:
        1. drug_mapping_dict: maps each drug to its set of protein targets
        2. fingerprint_dict: maps each drug to its molecular fingerprint set
        
        Inputs:
            None (uses paths stored in instance variables)
        
        Returns:
            None (stores data in instance variables)
        """
        # Load CSV files into pandas DataFrames
        self.drugs_df = pd.read_csv(self.drugs_data_path)
        self.targets_df = pd.read_csv(self.targets_data_path)
        
        # Create dictionary mapping each drug to set of protein targets it binds
        # Key: db_id (drug), Value: set of uniprot_accession (proteins)
        self.drug_mapping_dict = {}
        for _, row in self.targets_df.iterrows():
            drug_id = row['db_id']
            protein = row['uniprot_accession']
            
            if drug_id not in self.drug_mapping_dict:
                self.drug_mapping_dict[drug_id] = set()
            self.drug_mapping_dict[drug_id].add(protein)
        
        # Create dictionary mapping each drug to its fingerprint set
        # Key: db_id (drug), Value: set of fingerprint features (integers)
        self.fingerprint_dict = {}
        for _, row in self.drugs_df.iterrows():
            drug_id = row['db_id']
            # Parse space-separated fingerprint string into set of integers
            fingerprint = set(map(int, row['maccs'].split()))
            self.fingerprint_dict[drug_id] = fingerprint

    def shared_targets(self, drug_a, drug_b):
        """
        Check if two drugs share any protein targets.
        
        Inputs:
            drug_a: string, db_id of first drug
            drug_b: string, db_id of second drug
        
        Returns:
            int, 1 if drugs share at least one protein target, 0 otherwise
        """
        # Use set intersection to find common targets
        return 1 if self.drug_mapping_dict[drug_a] & self.drug_mapping_dict[drug_b] else 0

    def tanimoto_coeff(self, drug_a, drug_b):
        """
        Calculate Tanimoto coefficient between two drugs based on fingerprints.
        
        The Tanimoto coefficient (Jaccard Index) measures chemical similarity:
        Tc = |fpt(A) ∩ fpt(B)| / |fpt(A) ∪ fpt(B)|
        
        Where fpt(A) and fpt(B) are the fingerprint feature sets.
        
        Inputs:
            drug_a: string, db_id of first drug
            drug_b: string, db_id of second drug
        
        Returns:
            float, Tanimoto coefficient between 0.0 (no similarity) and 1.0 (identical)
        """
        # Get fingerprint sets for both drugs
        fpt_a = self.fingerprint_dict[drug_a]
        fpt_b = self.fingerprint_dict[drug_b]
        
        # Calculate size of intersection (common features)
        intersection = len(fpt_a & fpt_b)
        # Calculate size of union (all unique features)
        union = len(fpt_a | fpt_b)
        
        # Return ratio, handling edge case of empty union
        return intersection / union if union > 0 else 0.0

    def calculate_tanimoto(self):
        """
        Calculate pairwise Tanimoto coefficients for all unique drug pairs.
        
        For each pair of drugs, calculates:
        1. Tanimoto coefficient (chemical similarity)
        2. Whether they share any protein targets
        
        Inputs:
            None (uses data loaded in instance variables)
        
        Returns:
            None (stores results in self.tanimoto_output_df)
        
        Output format: drug_a, drug_b, tanimoto_score (6 decimals), shared_target (0 or 1)
        Note: Only includes each pair once (A,B but not B,A)
        """
        results = []
        drug_ids = list(self.fingerprint_dict.keys())
        
        # Generate all unique pairs using nested loops
        # i starts from 0, j starts from i+1 to avoid duplicates and self-comparisons
        for i in range(len(drug_ids)):
            for j in range(i + 1, len(drug_ids)):
                drug_a = drug_ids[i]
                drug_b = drug_ids[j]
                
                # Calculate Tanimoto coefficient for chemical similarity
                tc_score = self.tanimoto_coeff(drug_a, drug_b)
                
                # Check if they share protein targets
                # Handle case where drug may not have any known targets
                if drug_a in self.drug_mapping_dict and drug_b in self.drug_mapping_dict:
                    shared = self.shared_targets(drug_a, drug_b)
                else:
                    shared = 0
                
                # Store results with Tanimoto rounded to 6 decimal places
                results.append({
                    'drug_a': drug_a,
                    'drug_b': drug_b,
                    'tanimoto': f'{tc_score:.6f}',
                    'shared_target': shared
                })
        
        # Convert results to DataFrame for easy CSV writing
        self.tanimoto_output_df = pd.DataFrame(results)
    
    def write_output(self):
        """
        Write calculated Tanimoto results to CSV file without header row.
        
        Creates output directory if it doesn't exist.
        Output format: drug_a,drug_b,tanimoto_score,shared_target (no header)
        
        Inputs:
            None (uses data in self.tanimoto_output_df and self.output_path)
        
        Returns:
            None (writes file to disk)
        """
        # Check if the output directory exists, if not create it
        output_dir = os.path.dirname(self.output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Write to CSV without header row or index column
        self.tanimoto_output_df.to_csv(self.output_path, index=False, header=False)

    def run_tanimoto(self):
        """
        Execute the complete Tanimoto analysis pipeline.
        
        Performs three steps:
        1. Load drug and target data
        2. Calculate all pairwise Tanimoto coefficients
        3. Write results to output file
        
        Inputs:
            None (uses data stored in instance variables)
        
        Returns:
            None (writes output file and prints progress messages)
        """
        print("Running Tanimoto...")
        print("Loading data...")
        self.load_data()
        print("Calculating Tanimoto coefficients...")
        self.calculate_tanimoto()
        print("Writing output...")
        self.write_output()

if __name__ == "__main__":
    """
    Main execution block for command-line usage.
    
    Usage: python tanimoto.py <drugs.csv> <targets.csv> <output.csv>
    """
    print("starting tanimoto.py")
    
    # Check for correct number of command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python tanimoto.py <drugs_data_path> <targets_data_path> <output_path>")
        sys.exit(1)
    
    # Parse command-line arguments
    drugs_data_path = sys.argv[1]
    targets_data_path = sys.argv[2]
    output_path = sys.argv[3]

    # Create Tanimoto object and run analysis
    tanimoto = Tanimoto(drugs_data_path, targets_data_path, output_path)
    tanimoto.run_tanimoto()
