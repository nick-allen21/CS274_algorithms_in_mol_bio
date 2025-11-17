import pandas as pd
import sys
from Bio import SeqIO
import numpy as np
import torch
import esm
import os
import pickle
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples
from tqdm import tqdm

class ESM_model:
    def __init__(self, fasta_file, output_dir, pro_fam_datafile):
        self.fasta_file = fasta_file
        self.output_dir = output_dir
        self.prot_fam_datafile = pro_fam_datafile
        self.prot_fam_data = pd.read_csv(pro_fam_datafile, sep="\t")
        self.mean_pooled = None
        self.max_pooled = None

    def load_fasta(self, fasta_file):
        """
        Takes in a fasta file and extracts the UniProt protein ID and sequence from each entry. 
        This function should return a list of tuples where the first element of each tuple is a UniProt protein id and the second element is the corresponding sequence.
        Use the BioPython package’s SeqIO functions to do this (https://biopython.org/wiki/SeqIOLinks to an external site.)
        The fasta file we are working with is downloaded from UniProt. The header takes the format: >sp|<<UniProt Id>>|<<Protein Name>>_<<organism>> 
        """
        print(f"[ESM_model] Loading FASTA: {fasta_file}")
        id_seq_list = []
        for record in tqdm(SeqIO.parse(fasta_file, "fasta"), desc="Reading FASTA records"):
            header = record.description if record.description else record.id
            if "|" in header:
                parts = header.split("|")
                uniprot_id = parts[1] if len(parts) > 1 else parts[0]
            else:
                uniprot_id = header.split()[0]
            sequence = str(record.seq)
            id_seq_list.append((uniprot_id, sequence))
        self.id_seq_list = id_seq_list
        print(f"[ESM_model] Parsed {len(id_seq_list)} sequences.")
        return id_seq_list

    def get_vectors(self, list_tuples_protId_seq):
        """
        Takes in the list of tuples from load_fasta and returns a dictionary with the UniProt protein ids as the keys and 
        returns their associated representation vectors from the ESM model as the values.
        We will use a small version of the ESM2 model for this project due to limited compute availability: esm2_t33_650M_UR50D
        Pass only one sequence to the model at a time. When passing more than one sequence together to the model
        ESM2 instructions:https://github.com/facebookresearch/esm/blob/main/README.md
        """
        # Load model once
        print("[ESM_model] Loading ESM2 model (esm2_t33_650M_UR50D)...")
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        model.eval()
        batch_converter = alphabet.get_batch_converter()
        embeddings_dict = {}

        print(f"[ESM_model] Generating per-residue embeddings for {len(list_tuples_protId_seq)} proteins...")
        for prot_id, seq in tqdm(list_tuples_protId_seq, desc="Embedding proteins"):
            data = [(prot_id, seq)]
            _, _, batch_tokens = batch_converter(data)
            with torch.no_grad():
                outputs = model(batch_tokens, repr_layers=[33], return_contacts=False)
            token_reps = outputs["representations"][33]  # shape: (1, L+2, 1280)
            token_reps = token_reps[0, 1:-1, :].cpu().numpy()  # remove BOS/EOS -> (L, 1280)
            embeddings_dict[prot_id] = token_reps

        self.embeddings_dict = embeddings_dict
        print(f"[ESM_model] Finished embeddings. Stored embeddings for {len(embeddings_dict)} proteins.")
        return embeddings_dict

    def pool_representations(self, dict_protId_embedding, mean_max_param):
        """
        Takes in your dictionary of per token vectors and a string parameter denoting the pooling regime ("mean" or "max") 
        and outputs a dictionary with the UniProt protein ids as keys and the pooled embedding vectors (as a numpy object) as the values (dimension = [1,1280]).
        The output embedding vectors from the ESM model are per-residue (each vector corresponds to one amino acid in a protein), 
        but some downstream tasks require a single vector representing an entire protein. 
        Per-residue vectors can be pooled together to create one representation vector for a protein.
        Two common pooling strategies are mean pooling and max pooling. 
        Mean pooling averages the embeddings across all residues for each sequence, while max pooling takes the maximum value for each embedding feature across all residues for each sequence.
        """
        if mean_max_param is None:
            raise ValueError("mean_max_param must be 'mean' or 'max'.")
        mode = str(mean_max_param).strip().lower()
        if mode not in ("mean", "max"):
            raise ValueError("mean_max_param must be 'mean' or 'max'.")

        pooled = {}
        print(f"[ESM_model] Pooling embeddings with mode='{mode}' for {len(dict_protId_embedding)} proteins...")
        for prot_id, per_residue in tqdm(dict_protId_embedding.items(), desc=f"Pooling ({mode})"):
            if per_residue is None:
                continue
            arr = np.asarray(per_residue)
            if arr.ndim != 2:
                raise ValueError(f"Embedding for {prot_id} must have shape (L, 1280); got {arr.shape}")
            if mode == "mean":
                pooled_vec = np.mean(arr, axis=0, keepdims=True)
            else:
                pooled_vec = np.max(arr, axis=0, keepdims=True)
            pooled[prot_id] = pooled_vec
        if mean_max_param == "mean":
            self.mean_pooled = pooled
        else:
            self.max_pooled = pooled
        print(f"[ESM_model] Completed pooling: {len(pooled)} vectors.")
        return pooled

    def create_t_SNE(self, mean_pooled, max_pooled):
        """
        Use these two dictionaries to create a t-SNE plot for each pooling regime (the t-SNE function is from the scikit-learn package). 
        t-Distributed Stochastic Neighbor Embedding (t-SNE) is a powerful dimensionality reduction technique primarily used for visualizing high-dimensional data in a lower-dimensional space,
        typically two or three dimensions. It works by modeling the similarities between data points as probabilities and aims to preserve these similarities when mapping the data to a lower-dimensional representation. 
        The algorithm emphasizes retaining local structures, making it particularly effective for revealing clusters in data. 

        Use n_components = 2, random_state=42 when creating the t-SNE
        Save the two figures in the outputs folder as maxPooled_viz.png and meanPooled_viz.png.
        """
        def _stack_vectors(emb_dict, desc):
            ids = []
            vectors = []
            for prot_id, vec in tqdm(emb_dict.items(), desc=desc):
                arr = np.asarray(vec)
                if arr.ndim == 2 and arr.shape[0] == 1:
                    arr = arr[0]
                ids.append(prot_id)
                vectors.append(arr)
            X = np.vstack(vectors) if vectors else np.empty((0, 1280))
            return ids, X

        # Mean pooled t-SNE
        mean_ids, mean_X = _stack_vectors(mean_pooled, desc="Stacking mean-pooled vectors")
        if mean_X.size > 0:
            print(f"[ESM_model] Running t-SNE (mean pooled) on {mean_X.shape[0]} proteins...")
            tsne = TSNE(n_components=2, random_state=42, init="pca", learning_rate="auto")
            mean_emb = tsne.fit_transform(mean_X)
            plt.figure(figsize=(6, 5))
            plt.scatter(mean_emb[:, 0], mean_emb[:, 1], s=12, alpha=0.8)
            plt.xlabel("t-SNE 1")
            plt.ylabel("t-SNE 2")
            plt.title("Mean Pooled t-SNE")
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "meanPooled_viz.png"), dpi=200)
            plt.close()
            print(f"[ESM_model] Saved meanPooled_viz.png")

        # Max pooled t-SNE
        max_ids, max_X = _stack_vectors(max_pooled, desc="Stacking max-pooled vectors")
        if max_X.size > 0:
            print(f"[ESM_model] Running t-SNE (max pooled) on {max_X.shape[0]} proteins...")
            tsne = TSNE(n_components=2, random_state=42, init="pca", learning_rate="auto")
            max_emb = tsne.fit_transform(max_X)
            plt.figure(figsize=(6, 5))
            plt.scatter(max_emb[:, 0], max_emb[:, 1], s=12, alpha=0.8)
            plt.xlabel("t-SNE 1")
            plt.ylabel("t-SNE 2")
            plt.title("Max Pooled t-SNE")
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, "maxPooled_viz.png"), dpi=200)
            plt.close()
            print(f"[ESM_model] Saved maxPooled_viz.png")

    def calculate_silhouette_score (self, family_data_file, pooled_vectors_dict):
        """
        Takes in the protein family metadata file and a pooled vector dictionary, calculates the silhouette score for each family, and returns a dataframe with the family name and corresponding score per row. 
        You may want to use this function: https://scikit-learn.org/dev/modules/generated/sklearn.metrics.silhouette_score.htmlLinks to an external site.
        We are considering each protein family to be a cluster for the purposes of this question
        """
        print(f"[ESM_model] Loading family metadata from: {family_data_file}")
        fam_df = pd.read_csv(family_data_file)
        
        # Handle duplicates by keeping the first occurrence
        before = len(fam_df)
        fam_df = fam_df.drop_duplicates(subset=["protid"], keep="first").reset_index(drop=True)
        dropped = before - len(fam_df)
        if dropped > 0:
            print(f"[ESM_model] Warning: {dropped} duplicate protid rows dropped from family metadata.")

        fam_map = dict(zip(fam_df["protid"], fam_df["family"]))

        # Intersect proteins present in pooled vectors and family file
        common_ids = [pid for pid in pooled_vectors_dict.keys() if pid in fam_map]
        print(f"[ESM_model] Computing silhouettes for {len(common_ids)} proteins (intersection).")
        if len(common_ids) < 2:
            print("[ESM_model] Not enough proteins to compute silhouette. Returning empty dataframe.")
            return pd.DataFrame(columns=["family", "silhouette_score"])

        # Stack vectors and build labels
        X_list = []
        y_list = []
        for pid in tqdm(common_ids, desc="Stacking vectors for silhouette"):
            vec = pooled_vectors_dict[pid]
            arr = np.asarray(vec)
            if arr.ndim == 2 and arr.shape[0] == 1:
                arr = arr[0]
            elif arr.ndim != 1:
                raise ValueError(f"Vector for {pid} must be shape (1,1280) or (1280,), got {arr.shape}")
            X_list.append(arr)
            y_list.append(fam_map[pid])

        # Need at least 2 unique labels
        if len(set(y_list)) < 2:
            print("[ESM_model] Only one family present after intersection; silhouette undefined. Returning empty dataframe.")
            return pd.DataFrame(columns=["family", "silhouette_score"])

        X = np.vstack(X_list)
        s_samples = silhouette_samples(X, y_list)
        per_sample_df = pd.DataFrame({"protid": common_ids, "family": y_list, "silhouette": s_samples})
        by_family = per_sample_df.groupby("family", as_index=False)["silhouette"].mean()
        by_family = by_family.rename(columns={"silhouette": "silhouette_score"})
        print(f"[ESM_model] Computed silhouette scores for {len(by_family)} families.")
        return by_family

    def create_family_t_SNE(self, family_data_file, mean_pooled, max_pooled):
        """
        Create two t-SNE plots (one for each pooling regime) such that the points are colored by their family membership. 
        Save these figures in the outputs folder as meanPooled_viz_family.png and maxPooled_viz_family.png
        """
        print(f"[ESM_model] Loading family metadata from: {family_data_file}")
        fam_df = pd.read_csv(family_data_file)
        fam_df = fam_df.drop_duplicates(subset=["protid"], keep="first").reset_index(drop=True)
        fam_map = dict(zip(fam_df["protid"], fam_df["family"]))

        def _prep_XY(pooled_dict, desc):
            ids = []
            X_list = []
            y_list = []
            for pid, vec in tqdm(pooled_dict.items(), desc=desc):
                if pid not in fam_map:
                    continue
                arr = np.asarray(vec)
                if arr.ndim == 2 and arr.shape[0] == 1:
                    arr = arr[0]
                ids.append(pid)
                X_list.append(arr)
                y_list.append(fam_map[pid])
            if not X_list:
                return [], np.empty((0, 1280)), []
            X = np.vstack(X_list)
            return ids, X, y_list

        def _plot_tsne(X, labels, title, out_name):
            if X.size == 0 or len(set(labels)) == 0:
                print(f"[ESM_model] Skipping {title}: no data.")
                return
            print(f"[ESM_model] Running t-SNE for {title} on {X.shape[0]} proteins...")
            tsne = TSNE(n_components=2, random_state=42, init="pca", learning_rate="auto")
            emb = tsne.fit_transform(X)
            # Map families to integers/colors
            families = sorted(set(labels))
            fam_to_idx = {fam: i for i, fam in enumerate(families)}
            colors = plt.cm.tab20(np.linspace(0, 1, max(20, len(families))))
            plt.figure(figsize=(7, 6))
            for fam in families:
                idxs = [i for i, l in enumerate(labels) if l == fam]
                plt.scatter(emb[idxs, 0], emb[idxs, 1], s=14, alpha=0.8,
                            color=colors[fam_to_idx[fam] % len(colors)], label=fam)
            plt.xlabel("t-SNE 1")
            plt.ylabel("t-SNE 2")
            plt.title(title)
            plt.legend(loc="best", fontsize=8, markerscale=1.0, frameon=False)
            plt.tight_layout()
            out_path = os.path.join(self.output_dir, out_name)
            plt.savefig(out_path, dpi=200)
            plt.close()
            print(f"[ESM_model] Saved {out_name}")

        # Mean pooled family-colored t-SNE
        mean_ids, mean_X, mean_labels = _prep_XY(mean_pooled, desc="Preparing mean-pooled (family)")
        _plot_tsne(mean_X, mean_labels, "Mean Pooled t-SNE (colored by family)", "meanPooled_viz_family.png")

        # Max pooled family-colored t-SNE
        max_ids, max_X, max_labels = _prep_XY(max_pooled, desc="Preparing max-pooled (family)")
        _plot_tsne(max_X, max_labels, "Max Pooled t-SNE (colored by family)", "maxPooled_viz_family.png")

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    output_dir = sys.argv[2]
    prot_fam_datafile = sys.argv[3]
    ESMmodel = ESM_model(fasta_file, output_dir, prot_fam_datafile)
    ESMmodel.load_fasta(fasta_file)
    embeddings_dict = ESMmodel.get_vectors(ESMmodel.id_seq_list)

    # Use pickle to write the dictionary from get_vectors to a file in the output folder called: embeddings.pkl. 
    with open(os.path.join(output_dir, "embeddings.pkl"), "wb") as f:
        pickle.dump(embeddings_dict, f)

    # use the pool_representations function to perform both mean pooling and max pooling. 
    mean_pooled = ESMmodel.pool_representations(embeddings_dict, "mean")
    max_pooled = ESMmodel.pool_representations(embeddings_dict, "max")
    # Save both the dictionaries in the output folder directory as meanPooled_embeddings.pkl and maxPooled_embeddings.pkl. 
    with open(os.path.join(output_dir, "meanPooled_embeddings.pkl"), "wb") as f:
        pickle.dump(mean_pooled, f)
    with open(os.path.join(output_dir, "maxPooled_embeddings.pkl"), "wb") as f:
        pickle.dump(max_pooled, f)

    # create t SNE plots for both mean and max pooling
    ESMmodel.create_t_SNE(mean_pooled, max_pooled)

    # Run the calculate_silhouette_score function for both the mean pooled and max pooled dictionaries 
    mean_silhouette = ESMmodel.calculate_silhouette_score(prot_fam_datafile, mean_pooled)
    max_silhouette = ESMmodel.calculate_silhouette_score(prot_fam_datafile, max_pooled)
    # Save each dataframe as a .csv in the outputs folder called meanPooled_silhouette.csv or maxPooled_silhouette.csv
    mean_silhouette.to_csv(os.path.join(output_dir, "meanPooled_silhouette.csv"), index=False)
    max_silhouette.to_csv(os.path.join(output_dir, "maxPooled_silhouette.csv"), index=False)

    # Print the mean and max silhouette scores for each family
    print(f"Mean Pooled Silhouette Scores: {mean_silhouette}")
    print(f"Max Pooled Silhouette Scores: {max_silhouette}")

    # create family t SNE plots for both mean and max pooling   
    ESMmodel.create_family_t_SNE(prot_fam_datafile, mean_pooled, max_pooled)
    print(f"[ESM_model] Completed family t-SNE plots.")



