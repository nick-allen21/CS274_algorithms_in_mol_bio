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
import matplotlib.pyplot as plt
from chemoUtils import load_data, tanimoto_coeff

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
        self.drugs_df, self.targets_df, self.drug_mapping_dict, self.fingerprint_dict = load_data(drugs_data_path, targets_data_path)


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


    def calculate_tanimoto(self):
        """
        Calculate pairwise Tanimoto coefficients for all unique drug pairs.
        
        For each pair of drugs, calculates:
        1. Tanimoto coefficient
        2. Whether they share any protein targets
        
        Inputs:
            None (uses self stored vars)
        
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
                tc_score = tanimoto_coeff(drug_a, drug_b, self.fingerprint_dict)
                
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
        Output format: drug_a,drug_b,tanimoto_score,shared_target
        
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

    def create_histogram(self, tanimoto_values, title, filename, output_dir='results'):
        """
        Create and save a histogram of Tanimoto coefficient values.
        
        Inputs:
            tanimoto_values: list or array of float, Tanimoto coefficient values to plot
            title: string, title for the histogram figure
            filename: string, filename to save the histogram
            output_dir: string, directory where histogram will be saved
        
        Returns:
            None (saves histogram to disk)
        """
        # Create new figure
        plt.figure(figsize=(10, 6))
        
        # Create histogram with appropriate bins
        # Using 50 bins for good resolution 
        plt.hist(tanimoto_values, bins=50, edgecolor='black', alpha=0.7)
        
        # Set labels and title
        plt.xlabel('Tanimoto Coefficient', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.title(title, fontsize=14)
        
        # Add grid for better readability
        plt.grid(axis='y', alpha=0.3)
        
        # Ensure output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Save figure
        output_path = os.path.join(output_dir, filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Saved histogram: {output_path}")

    def generate_all_histograms(self, sunet_id):
        """
        Generate all three required histograms for Tanimoto analysis.
        
        Creates:
        1. Histogram of all Tanimoto values
        2. Histogram of Tanimotos for drug pairs that share targets
        3. Histogram of Tanimotos for drug pairs that don't share targets
        
        Inputs:
            sunet_id: string, SUNet ID for labeling histograms
        
        Returns:
            None (saves three histogram files to disk)
        """
        # Convert Tanimoto scores from string to float for plotting
        self.tanimoto_output_df['tanimoto_float'] = self.tanimoto_output_df['tanimoto'].astype(float)
        
        # Get output directory 
        output_dir = 'results'
        
        # 1. All Tanimoto values
        all_values = self.tanimoto_output_df['tanimoto_float']
        self.create_histogram(
            all_values,
            f"{sunet_id} All",
            'all_tanimoto.png',
            output_dir
        )
        
        # 2. Tanimotos for drug pairs that share a target (shared_target == 1)
        shared_values = self.tanimoto_output_df[
            self.tanimoto_output_df['shared_target'] == 1
        ]['tanimoto_float']
        self.create_histogram(
            shared_values,
            f"{sunet_id} Shared",
            'shared_tanimoto.png',
            output_dir
        )
        
        # 3. Tanimotos for drug pairs that don't share a target (shared_target == 0)
        not_shared_values = self.tanimoto_output_df[
            self.tanimoto_output_df['shared_target'] == 0
        ]['tanimoto_float']
        self.create_histogram(
            not_shared_values,
            f"{sunet_id} Not Shared",
            'notshared_tanimoto.png',
            output_dir
        )

    def run_tanimoto(self, sunet_id='nallen21'):
        """
        Execute the complete Tanimoto analysis pipeline.
        
        Performs five steps:
        1. Load drug and target data
        2. Calculate all pairwise Tanimoto coefficients
        3. Write results to output file
        4. Generate histograms of Tanimoto distributions
        
        Inputs:
            sunet_id: string, SUNet ID for labeling histograms
        
        Returns:
            None (writes output file, generates histograms, and prints progress messages)
        """
        print("Running Tanimoto...")
        print("Calculating Tanimoto coefficients...")
        self.calculate_tanimoto()
        print("Writing output...")
        self.write_output()
        print("Generating histograms...")
        self.generate_all_histograms(sunet_id)

if __name__ == "__main__":
    """
    Main execution block for command-line usage.
    
    Usage: python tanimoto.py <drugs.csv> <targets.csv> <output.csv> 
    """
    print("starting tanimoto.py")
    
    # Check for correct number of command-line arguments (3 required, 1 optional)
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print("Usage: python tanimoto.py <drugs_data_path> <targets_data_path> <output_path> ")
        sys.exit(1)
    
    # Parse command-line arguments
    drugs_data_path = sys.argv[1]
    targets_data_path = sys.argv[2]
    output_path = sys.argv[3]
    
    if len(sys.argv) != 4:
        print("Usage: python tanimoto.py <drugs_data_path> <targets_data_path> <output_path> ")
        sys.exit(1)
        
    # Create Tanimoto object and run analysis
    tanimoto = Tanimoto(drugs_data_path, targets_data_path, output_path)
    tanimoto.run_tanimoto()
