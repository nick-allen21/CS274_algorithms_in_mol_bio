"""
This script is used to calculate PPI for given sequences. 
It is part of an assingment for CS274
Written by: Nick Allen 
"""
# 0) Import necessary libraries
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
# 1) Read in the interaction file, interaction_network.txt, and disease gene file, disease_gene_list.txt. The file names for this input should be command line arguments.

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s - %(message)s"
)
logger = logging.getLogger(__name__)

# helper to convert camelCase/kebab-case to snake_case
def camel_to_snake(name):
    name = name.replace('-', '_')
    name = re.sub(r'(.)([A-Z][a-z]+)', r'\1_\2', name)
    name = re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', name)
    return name.lower()

# 2) Define and implement PPI class 
class PPI:
    def __init__(self):
        self.disease_gene_file = None
        self.interaction_file = None
        self.disease_gene_df = None
        self.interaction_df = None
        self.interaction_graph = None
        self.gene_nodes_dict = None

    def run_ppi(self, disease_gene_file, interaction_file, threshold):

        logger.info("Starting PPI pipeline with threshold=%s", threshold)
        self.load_data(disease_gene_file, interaction_file)
        gene_nodes, gene_embeddings = self.calculate_embedding()
        similar_genes = self.get_close_genes(gene_nodes, gene_embeddings, threshold)
        self.output_similar_genes(similar_genes)
        logger.info("Completed PPI pipeline. Similar genes written to similar_genes.txt")

    def load_data(self, disease_gene_file, interaction_file):

        logger.info("Loading disease genes from %s", disease_gene_file)
        # disease genes: one gene per line
        disease_gene_df = pd.read_csv(disease_gene_file, header=None, names=["gene"]) 

        logger.info("Loading interaction network from %s", interaction_file)
        # interaction network: whitespace-separated: geneA geneB weight
        interaction_df = pd.read_csv(
            interaction_file, sep=r"\s+", header=None, names=["gene1", "gene2", "weight"], engine="python"
        )
        self.disease_gene_df = disease_gene_df
        self.interaction_df = interaction_df

        # store interaction network as a networkx graph
        self.interaction_graph = nx.from_pandas_edgelist(
            interaction_df, source="gene1", target="gene2", edge_attr="weight", create_using=nx.Graph()
        )
        logger.info(
            "Loaded %d disease genes, %d interactions -> graph with %d nodes, %d edges",
            len(self.disease_gene_df), len(self.interaction_df),
            self.interaction_graph.number_of_nodes(), self.interaction_graph.number_of_edges()
        )

    def calculate_embedding(self):
        """
        Should return a list of the nodes in the graph and a list of their vector embeddings. 
        The order of the nodes should correspond to the order of the embeddings.
        Train the embedding model using the following parameters: window=3, min_count=1, batch_words=4
        """
        if os.path.exists("node2vec_pretrained"):

            logger.info("Loading pretrained Node2Vec model from node2vec_pretrained")
            model = Word2Vec.load("node2vec_pretrained")
            logger.info("Loaded pretrained model successfully")

        else:
            # create and fit the model here if it doesn't already exist
            # When calling Node2Vec use the following parameters: 
            # dimensions=64, walk_length=30, num_walks=100, workers=1, seed=42 
            # create and fit the model here if it doesn't already exist
            logger.info(
                "Training Node2Vec: dimensions=64, walk_length=30, num_walks=100, workers=1, seed=42"
            )
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
            logger.info("Fitting Word2Vec on generated walks: window=3, min_count=1, batch_words=4")
            model = node2vec.fit(window=3, min_count=1, batch_words=4)
            model.save("node2vec_pretrained")
            logger.info("Saved pretrained model to node2vec_pretrained")
            
        # use the model to create the embeddings for the genes
        gene_nodes = list(self.interaction_graph.nodes())
        gene_nodes_str = [str(n) for n in gene_nodes]
        gene_embeddings = [model.wv[n] for n in gene_nodes_str]
        logger.info("Created embeddings for %d genes", len(gene_nodes))

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
        node_to_index = {node: i for i, node in enumerate(gene_nodes)}

        # take embedding list and build a matrix where each column is a node embedding
        embeddings_matrix = np.vstack(gene_embeddings)

        # disease genes from input file (clean: drop NaNs/empties) and restrict to nodes present in graph
        disease_genes = [g.strip() for g in self.disease_gene_df["gene"].dropna().astype(str) if g and str(g).strip()]
        disease_genes_in_graph = [g for g in disease_genes if g in node_to_index]
        logger.info(
            "Computing close genes: %d disease genes provided, %d present in graph",
            len(disease_genes), len(disease_genes_in_graph)
        )

        # if none of the disease genes are in the embedding, return known disease genes
        if not disease_genes_in_graph:
            logger.warning("No disease genes found in graph; returning input disease genes only")
            return set(disease_genes)

        similar_genes = set()
        for dg in disease_genes_in_graph:
            idx = node_to_index[dg]

            # get the embedding for the disease gene matrix 
            dg_vec = embeddings_matrix[idx:idx+1, :]  # shape (1, d)

            # flatten to 1D array to get distances between disease gene and all other known disease genes
            dists = cosine_distances(dg_vec, embeddings_matrix).ravel()
            close_indices = np.where(dists <= float(threshold))[0]
            for ci in close_indices:
                similar_genes.add(gene_nodes[ci])

        # always include original disease genes that are present in the graph
        similar_genes.update(disease_genes_in_graph)
        logger.info("Selected %d similar genes (including disease genes)", len(similar_genes))
        return similar_genes

    def output_similar_genes(self, similar_genes):
        """
        Output the similar genes to a .txt file
        format:
        <pathway_name> <\t> <description> <\t> <gene1> <\t> <gene2> <\t>…
        """
        pathway_name = "PredictedSet"
        description = "PPI Node2Vec similar genes"
        genes_sorted = sorted(map(str, similar_genes))
        line = "\t".join([pathway_name, description] + genes_sorted)
        with open("similar_genes.txt", "w") as f:
            f.write(line + "\n")
        logger.info("Wrote %d genes to similar_genes.txt in pathway-format line", len(genes_sorted))


# Main method
def main():

    if len(sys.argv) != 3:
        logger.error("Usage: python3 ppi.py diseaseGeneFile interactionNetworkFile")
        sys.exit(1)

    disease_gene_file = sys.argv[1]
    interaction_file = sys.argv[2]
    logger.info("Invoked with disease_gene_file=%s, interaction_file=%s", disease_gene_file, interaction_file)
    ppi = PPI()
    threshold = .5
    logger.info("Using default threshold=%s", threshold)
    ppi.run_ppi(disease_gene_file, interaction_file, threshold)


if __name__ == "__main__":
    main()