"""
Network Visualization for Protein Ligand Set Similarity

Creates a network graph visualization with nodes colored by disease indication.

Author: Nicholas Allen
Course: CS274 - Algorithms in Molecular Biology
"""

import sys
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def load_protein_info(protein_nodes_csv):
    """
    Load protein information from CSV file.
    
    Inputs:
        protein_nodes_csv: string, path to protein nodes CSV file
    
    Returns:
        tuple of two dictionaries:
            - label_mapping: dict mapping uniprot_accession to uniprot_id
            - color_mapping: dict mapping uniprot_accession to color based on indication
    """
    # Read protein nodes data
    protein_df = pd.read_csv(protein_nodes_csv)
    
    # Create mapping from uniprot_accession to uniprot_id for node labels
    label_mapping = dict(zip(protein_df['uniprot_accession'], protein_df['uniprot_id']))
    
    # Color mapping based on indications
    indication_colors = {
        "bp": "red",
        "bp;cholesterol": "green",
        "bp;cholesterol;diabetes": "blue",
        "bp;diabetes": "purple"
    }
    
    # Create color mapping for each protein
    color_mapping = {}
    for _, row in protein_df.iterrows():
        protein_id = row['uniprot_accession']
        indication = row['indications']
        color_mapping[protein_id] = indication_colors.get(indication, "gray")
    
    return label_mapping, color_mapping


def create_network_graph(edgelist_file, label_mapping, color_mapping, output_file):
    """
    Create and save network graph visualization.
    
    Inputs:
        edgelist_file: string, path to network edgelist file
        label_mapping: dict, mapping uniprot_accession to uniprot_id
        color_mapping: dict, mapping uniprot_accession to color
        output_file: string, path to save output PNG file
    
    Returns:
        None (saves figure to disk)
    """
    # Read edgelist and create graph
    G = nx.read_edgelist(edgelist_file)
    
    # Create node colors list in same order as nodes
    node_colors = [color_mapping.get(node, "gray") for node in G.nodes()]
    
    # Create node labels using uniprot_id instead of uniprot_accession
    node_labels = {node: label_mapping.get(node, node) for node in G.nodes()}
    
    # Create figure with specified size (8x8 inches)
    plt.figure(figsize=(8, 8))
    
    # Use spring layout for node positioning
    pos = nx.spring_layout(G, seed=214, k=0.5, iterations=50)
    
    # Draw network
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=300, alpha=0.9)
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=1)
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=6)
    
    # Add legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
                   markersize=10, label='bp'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', 
                   markersize=10, label='bp;cholesterol'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', 
                   markersize=10, label='bp;cholesterol;diabetes'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='purple', 
                   markersize=10, label='bp;diabetes')
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    plt.axis('off')
    
    # Save with exact specifications: 8x8 inches at 150 DPI
    plt.savefig(output_file, dpi=150)
    plt.close()
    
    print(f"Network visualization saved to: {output_file}")
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges: {G.number_of_edges()}")


if __name__ == "__main__":
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python3 plot_graph.py <network_edgelist.txt> <protein_nodes.csv> <output_figure.png>")
        sys.exit(1)
    
    edgelist_file = sys.argv[1]
    protein_nodes_csv = sys.argv[2]
    output_file = sys.argv[3]
    
    # Load protein information
    print("Loading protein information...")
    label_mapping, color_mapping = load_protein_info(protein_nodes_csv)
    
    # Create and save network visualization
    print("Creating network visualization...")
    create_network_graph(edgelist_file, label_mapping, color_mapping, output_file)
