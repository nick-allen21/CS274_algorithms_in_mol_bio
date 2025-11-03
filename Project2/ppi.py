"""
This script is used to calculate PPI for given sequences. 
It is part of an assingment for CS274
Written by: Nick Allen 
"""
import sys
import os
import re
import logging
import pandas as pd
import networkx as nx
from gensim.models import Word2Vec
from node2vec import Node2Vec
import numpy as np
from sklearn.metrics.pairwise import cosine_distances


class PPI:
    """
    Class to calculate PPI for given sequences.
    Attributes:
        disease_gene_file: file path to the disease gene list
        interaction_file: file path to the interaction network
        disease_gene_df: dataframe containing the disease genes
        interaction_df: dataframe containing the interaction network
        interaction_graph: graph containing the interaction network
        gene_nodes_dict: dictionary containing the gene nodes and their embeddings

    Methods:
        run_ppi: run the PPI pipeline
        load_data: load the data from the files
        calculate_embedding: calculate the embedding for the genes
        get_close_genes: get the close genes for the given set of disease genes
        output_similar_genes: output the similar genes to a file
    """
    def __init__(self):
        self.disease_gene_file = None
        self.interaction_file = None
        self.disease_gene_df = None
        self.interaction_df = None
        self.interaction_graph = None
        self.gene_nodes_dict = None

    def run_ppi(self, disease_gene_file, interaction_file, threshold):
        """
        Run the PPI pipeline.
        """

        self.load_data(disease_gene_file, interaction_file)
        gene_nodes, gene_embeddings = self.calculate_embedding()
        similar_genes = self.get_close_genes(gene_nodes, gene_embeddings, threshold)
        self.output_similar_genes(similar_genes)

    def load_data(self, disease_gene_file, interaction_file):

        # disease genes: one gene per line
        disease_gene_df = pd.read_csv(disease_gene_file, header=None, names=["gene"]) 

        # interaction network: whitespace-separated: geneA geneB weight
        interaction_df = pd.read_csv(
            interaction_file, sep=r"\s+", header=None, names=["gene1", "gene2", "weight"], engine="python"
        )
        self.disease_gene_df = disease_gene_df
        self.interaction_df = interaction_df

        # store interaction network as an UNWEIGHTED networkx graph (ignore weights completely)
        self.interaction_graph = nx.from_pandas_edgelist(
            interaction_df, source="gene1", target="gene2", create_using=nx.Graph()
        )
        
    def calculate_embedding(self):
        """
        Should return a list of the nodes in the graph and a list of their vector embeddings. 
        The order of the nodes should correspond to the order of the embeddings.
        Train the embedding model using the following parameters: window=3, min_count=1, batch_words=4
        """
        if os.path.exists("node2vec_pretrained"):
            model = Word2Vec.load("node2vec_pretrained")

        else:
            # create and fit the model here if it doesn't already exist
            # When calling Node2Vec use the following parameters: 
            # dimensions=64, walk_length=30, num_walks=100, workers=1, seed=42 
            node2vec = Node2Vec(
                self.interaction_graph,
                dimensions=64,
                walk_length=30,
                num_walks=100,
                workers=1,
                seed=42,
            )

            # train the embedding model using the following parameters: 
            # window=3, min_count=1, batch_words=4
            model = node2vec.fit(window=3, min_count=1, batch_words=4)
            model.save("node2vec_pretrained")
            
        # use the model to create the embeddings for the genes
        gene_nodes = list(self.interaction_graph.nodes())
        gene_nodes_str = [str(n) for n in gene_nodes]
        gene_embeddings = [model.wv[n] for n in gene_nodes_str]

        # save the gene nodes and embeddings to a dictionary
        gene_nodes_dict = dict(zip(gene_nodes_str, gene_embeddings))
        self.gene_nodes_dict = gene_nodes_dict

        return gene_nodes, gene_embeddings

    def get_close_genes(self, gene_nodes, gene_embeddings, threshold):
        """
        where gene_nodes is the list of nodes produced and gene_embeddings is the list of node vectors produced by the embedding 
        (they should be in corresponding order). The threshold is the maximum distance two genes can be to be considered similar. 
        This function should return a set of genes that are predicted similar to the given set of disease genes, 
        including the set of known disease genes. When multiple disease genes are provided, you will take a union of their similar genes.
        If no genes are predicted, return simply the known disease genes.
        """
        # create a index for each node in order to building an embedding matrix 
        # each column of the matrix is a node embedding
        # EFFICIENCY FIX:
        node_to_index = {node: i for i, node in enumerate(gene_nodes)}

        # take embedding list and build a matrix where each column is a node embedding as opposed to looping through each gene
        embeddings_matrix = np.vstack(gene_embeddings)

        # disease genes from input file clean and restrict to nodes present in graph
        disease_genes = [g.strip() for g in self.disease_gene_df["gene"].dropna().astype(str) if g and str(g).strip()]
        disease_genes_in_graph = [g for g in disease_genes if g in node_to_index]

        # if none of the disease genes are in the embedding, return known disease genes
        if not disease_genes_in_graph:
            return set(disease_genes)

        similar_genes = set()
        for dg in disease_genes_in_graph:
            idx = node_to_index[dg]

            # get the embedding for the disease gene matrix 
            dg_vec = embeddings_matrix[idx:idx+1, :]

            # flatten to 1D array to get distances between disease gene and all other known disease genes
            # far more efficient than looping through each gene and computing the distance for each gene
            dists = cosine_distances(dg_vec, embeddings_matrix).ravel()
            close_indices = np.where(dists <= float(threshold))[0]
            for ci in close_indices:
                similar_genes.add(gene_nodes[ci])

        # always include original disease genes that are present in the graph
        similar_genes.update(disease_genes_in_graph)
        return similar_genes

    def output_similar_genes(self, similar_genes):
        """
        Output the similar genes to a .txt file
        format: <pathway_name> <\t> <description> <\t> <gene1> <\t> <gene2> <\t>…
        """
        pathway_name = "PredictedSet"
        description = "PPI Node2Vec similar genes"
        genes_sorted = sorted(map(str, similar_genes))
        line = "\t".join([pathway_name, description] + genes_sorted)
        with open("similar_genes.txt", "w") as f:
            f.write(line + "\n")
        print(f"Number of similar genes: {len(similar_genes)}")

    def count_interaction_pairs_outside_threshold(self, gene_nodes, gene_embeddings, threshold: float) -> None:
        """
        Count interacting gene pairs (edges in the interaction file) whose cosine distance
        in the embedding space is greater than the given threshold. Prints the count and a
        brief rationale for why this can happen.
        """
        # Build quick lookup and normalized embedding matrix for fast cosine distance via dot product
        node_to_index = {g: i for i, g in enumerate(gene_nodes)}
        emb = np.vstack(gene_embeddings)
        norms = np.linalg.norm(emb, axis=1, keepdims=True)
        emb_norm = emb / (norms + 1e-12)

        total = 0
        outside = 0
        for _, row in self.interaction_df[["gene1", "gene2"]].iterrows():
            g1, g2 = row["gene1"], row["gene2"]
            if g1 in node_to_index and g2 in node_to_index:
                total += 1
                i, j = node_to_index[g1], node_to_index[g2]
                # Cosine distance = 1 - cosine similarity
                cos_sim = float(np.dot(emb_norm[i], emb_norm[j]))
                cos_dist = 1.0 - cos_sim
                if cos_dist > threshold:
                    outside += 1

        print(f"At distance threshold {threshold}, {outside} of {total} documented interacting pairs exceed the threshold.")
        print("These pairs can appear distant because Node2Vec captures broader random-walk co-occurrence, "
              "not just single edges. Using an unweighted graph, noise/indirect associations, and specific "
              "hyperparameters can place some directly interacting genes farther apart in the embedding space.")


def main():
    """
    Main method to run the PPI pipeline.
    Loads in the disease gene file and interaction network file from the command line arguments.
    Runs the PPI pipeline with a given threshold
    Outputs the similar genes to a file.
    """

    if len(sys.argv) != 3:
        sys.exit(1)

    disease_gene_file = sys.argv[1]
    interaction_file = sys.argv[2]
    ppi = PPI()
    threshold = .2
    ppi.run_ppi(disease_gene_file, interaction_file, threshold)

    # # Separate analysis: count interacting pairs beyond a small threshold using current embeddings
    # gene_nodes, gene_embeddings = ppi.calculate_embedding()
    # ppi.count_interaction_pairs_outside_threshold(gene_nodes, gene_embeddings, 0.1)


if __name__ == "__main__":
    main()