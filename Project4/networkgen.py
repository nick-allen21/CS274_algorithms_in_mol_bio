"""
Network Generation for Protein Ligand Set Similarity

Generates a network of proteins connected by significant ligand set similarity
based on bootstrap p-values.

Author: Nicholas Allen
Course: CS274 - Algorithms in Molecular Biology
"""

from pvalue import PValue
from tanimoto import Tanimoto
import sys
import pandas as pd
import networkx as nx


class NetworkGen:
    """
    Class to generate protein similarity networks based on ligand set comparisons.
    """
    
    def __init__(self, drugs_csv, targets_csv, protein_nodes_csv):
        """
        Initialize NetworkGen with data file paths.
        
        Inputs:
            drugs_csv: string, path to drugs CSV file
            targets_csv: string, path to targets CSV file
            protein_nodes_csv: string, path to protein nodes CSV file
        
        Returns:
            None
        """
        self.drugs_csv_path = drugs_csv
        self.targets_csv_path = targets_csv
        self.protein_nodes_csv_path = protein_nodes_csv
        self.protein_nodes_df = pd.read_csv(self.protein_nodes_csv_path)
        self.pvalue_calc = PValue(self.drugs_csv_path, self.targets_csv_path)

    def generate_network(self, n_iterations=500, random_seed=214, p_threshold=0.05, output_file='network_edgelist.txt'):
        """
        Generate network edgelist for protein pairs with significant ligand similarity.
        
        For each unique pair of proteins in protein_nodes.csv:
        1. Calculate bootstrap p-value for ligand set similarity
        2. If p-value <= threshold, add edge to network
        3. Write edgelist to output file
        
        Inputs:
            n_iterations: int, number of bootstrap iterations (default 500)
            random_seed: int, random seed for reproducibility (default 214)
            p_threshold: float, p-value threshold for significance (default 0.05)
            output_file: string, output filename for edgelist 
        
        Returns:
            None (writes edgelist file to disk)
        """
        # Get list of protein uniprot_accession IDs
        protein_ids = self.protein_nodes_df['uniprot_accession'].tolist()
        n_proteins = len(protein_ids)
        total_pairs = n_proteins * (n_proteins - 1) // 2
        
        print(f"Starting network generation...")
        print(f"Proteins: {n_proteins}")
        print(f"Total pairs to evaluate: {total_pairs}")
        print(f"Bootstrap iterations per pair: {n_iterations}")
        print(f"Significance threshold: p <= {p_threshold}\n")
        
        # Store edges that meet significance threshold
        edges = []
        current_pair = 0
        
        # Generate all unique pairs (loop structure ensures no duplicates or self-loops)
        # i < j guarantees each pair processed exactly once
        for i in range(len(protein_ids)):
            for j in range(i + 1, len(protein_ids)):
                protein_a = protein_ids[i]
                protein_b = protein_ids[j]
                
                # Calculate bootstrap p-value for this protein pair
                self.pvalue_calc.run_pvalue(n_iterations, random_seed, protein_a, protein_b)
                p_value = self.pvalue_calc.p_value
                
                # If significant, add edge to network
                if p_value <= p_threshold:
                    edges.append((protein_a, protein_b))
                
                # Progress update every 200 pairs
                current_pair += 1
                if (current_pair % 200) == 0 or (current_pair == total_pairs):
                    print(f"Progress: {current_pair}/{total_pairs} pairs evaluated ({100*current_pair/total_pairs:.1f}%)")
        
        # Write edgelist to file
        print(f"\nWriting edges to file...")
        with open(output_file, 'w') as f:
            for edge in edges:
                f.write(f"{edge[0]} {edge[1]}\n")
        
        print(f"\nNetwork generation complete!")
        print(f"Evaluated {total_pairs} protein pairs")
        print(f"Found {len(edges)} significant edges (p <= {p_threshold})")
        print(f"Output: {output_file}")

if __name__ == "__main__":
    
    # Parse command-line arguments
    drugs_csv = sys.argv[1]
    targets_csv = sys.argv[2]
    protein_nodes_csv = sys.argv[3]
    
    # Create NetworkGen instance and generate network
    networkgen = NetworkGen(drugs_csv, targets_csv, protein_nodes_csv)
    networkgen.generate_network()