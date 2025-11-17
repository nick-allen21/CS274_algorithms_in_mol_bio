import os
import pandas as pd
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt

class AnalyzePredictions:
    def __init__(self):
        self.predictions_file = None
        self.family_data_file = None

    def load_data(self, predicted_results_dir, cutoff):
        """
        Write a function load_data(predicted_results_dir, cutoff) that takes in the directory with the predicted results and a cutoff value.
        This function should iterate through the files in the predict_results folder (use the function os.listdir(<folder path>)) and read in each tsv file.
        Save any rows where the predicted probability is greater than or equal to the cutoff in a pandas dictionary that has the columns protid, Index, AA, Prob and return this dataframe.
        """
        print(f"Loading data from {predicted_results_dir} with cutoff {cutoff}")
        records = []
        files = [f for f in os.listdir(predicted_results_dir) if f.lower().endswith(".tsv")]
        for fname in tqdm(files, desc="Reading prediction TSVs"):
            protid = os.path.splitext(fname)[0]
            fpath = os.path.join(predicted_results_dir, fname)
           
            df = pd.read_csv(fpath, sep="\t")
            df = df[df["Prob"] >= float(cutoff)].copy()
            df.insert(0, "protid", protid)
            records.append(df[["protid", "Index", "AA", "Prob"]])
    
        return pd.concat(records, ignore_index=True)

    def get_bindingSite_labels(self, protein_bindingSite_file):
        """
        Reads the file protein_bindingSite.csv into a dataframe and returns this dataframe, such that each row has a protein id and one binding site. 
        A protein will have multiple binding sites, so a protein id will appear in multiple rows of the dataframe.
        """
  
        print(f"Loading binding site labels from {protein_bindingSite_file}")
        df = pd.read_csv(protein_bindingSite_file)
        # Normalize expected columns
        expected_cols = {"protid", "binding_site"}
        # Keep required columns and coerce types
        if not expected_cols.issubset(df.columns):
            raise ValueError("Binding site file must have columns 'protid' and 'binding_site'.")
        df = df[["protid", "binding_site"]].copy()
        df["binding_site"] = pd.to_numeric(df["binding_site"], errors="coerce").astype("Int64")
        df = df.dropna(subset=["protid", "binding_site"]).reset_index(drop=True)
        return df

    def calculate_metrics(self, predicted_results_dir, cutoff, bindingSite_dataframe):
        """
        Takes in a cutoff value, the predicted results directory, and the known binding sites dataframe. 
        This function should return the accuracy, precision, true positive rate, and false positive rate for your predictions (in that order). 
        - It may be helpful to call load_data within this function.
        - There are some predictions (proteins in the predicted_results directory) that are not in the known binding sites dataframe. 
        Be sure to include only proteins that are in the known binding sites dataframe in these metric calculations.
        """
        # Load predicted positives at cutoff
        preds_df = self.load_data(predicted_results_dir, cutoff)
        if preds_df.empty:
            return 0.0, 0.0, 0.0, 0.0

        # Prepare known binding sites
        known_df = bindingSite_dataframe.copy()
        if "protid" not in known_df.columns or "binding_site" not in known_df.columns:
            raise ValueError("bindingSite_dataframe must have columns 'protid' and 'binding_site'.")

        # Restrict to proteins present in both predictions and labels
        pred_proteins = set(preds_df["protid"].unique().tolist())
        known_proteins = set(known_df["protid"].unique().tolist())
        use_proteins = pred_proteins.intersection(known_proteins)
        if not use_proteins:
            return 0.0, 0.0, 0.0, 0.0

        preds_df = preds_df[preds_df["protid"].isin(use_proteins)].copy()
        known_df = known_df[known_df["protid"].isin(use_proteins)].copy()

        # Build maps for quick lookup
        known_map = {}
        for pid, g in known_df.groupby("protid"):
            known_map[pid] = set(g["binding_site"].astype(int).tolist())

        # Sequence length map if available; else will infer later
        seq_len_map = {}
        if "seq_length" in known_df.columns:
            seq_len_map = {pid: int(g["seq_length"].iloc[0]) for pid, g in known_df.groupby("protid")}

        # Accumulate confusions globally across proteins
        TP = 0
        FP = 0
        FN = 0
        TN = 0

        for pid, g in preds_df.groupby("protid"):
            pred_pos = set(pd.to_numeric(g["Index"], errors="coerce").dropna().astype(int).tolist())
            true_pos = known_map.get(pid, set())
            # Infer sequence length if not provided
            if pid in seq_len_map:
                L = seq_len_map[pid]
            else:
                max_idx = -1
                if pred_pos:
                    max_idx = max(max_idx, max(pred_pos))
                if true_pos:
                    max_idx = max(max_idx, max(true_pos))
                L = max_idx + 1 if max_idx >= 0 else 0
            # Confusion counts
            tp = len(pred_pos & true_pos)
            fp = len(pred_pos - true_pos)
            fn = len(true_pos - pred_pos)
            tn = max(L - tp - fp - fn, 0)
            TP += tp
            FP += fp
            FN += fn
            TN += tn

        denom_acc = TP + TN + FP + FN
        accuracy = (TP + TN) / denom_acc if denom_acc > 0 else 0.0
        precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
        tpr = TP / (TP + FN) if (TP + FN) > 0 else 0.0
        fpr = FP / (FP + TN) if (FP + TN) > 0 else 0.0
        return accuracy, precision, tpr, fpr

    def plot_histogram(self, predicted_results_dir, cutoff):
        """
        Plot a histogram of the frequency that each amino acid is predicted to be an ATP binding site across all proteins. Save this plot in the output directory.
        """
        df = self.load_data(predicted_results_dir, cutoff)
        if df.empty:
            print("No predictions above cutoff; skipping histogram.")
            return
        counts = df["AA"].value_counts().sort_index()
        
        outputs_dir = "outputs"
        os.makedirs(outputs_dir, exist_ok=True)
        plt.figure(figsize=(8, 4))
        counts.plot(kind="bar", color="#4C72B0")
        plt.xlabel("Amino Acid")
        plt.ylabel("Count predicted as binding (>= cutoff)")
        plt.title(f"Predicted ATP-binding residues by amino acid (cutoff={cutoff})")
        plt.tight_layout()
        out_path = os.path.join(outputs_dir, "AA_binding_histogram.png")
        plt.savefig(out_path, dpi=200)
        plt.close()
        print(f"Saved histogram to {out_path}")

    def plot_roc_curve(self, predicted_results_dir, cutoff, bindingSite_dataframe):
        # Generate at least 10 cutoff points in [0,1]
        thresholds = np.linspace(0.0, 1.0, 11)
        fprs = []
        tprs = []
        for thr in thresholds:
            acc, prec, tpr, fpr = self.calculate_metrics(predicted_results_dir, thr, bindingSite_dataframe)
            fprs.append(fpr)
            tprs.append(tpr)
        # Determine outputs directory
        parent = os.path.dirname(predicted_results_dir.rstrip("/"))
        root = os.path.dirname(parent)
        outputs_dir = os.path.join(root, "outputs")
        os.makedirs(outputs_dir, exist_ok=True)
        plt.figure(figsize=(5, 5))
        plt.plot(fprs, tprs, marker="o", label="ROC")
        plt.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Random")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC Curve across probability cutoffs")
        plt.legend()
        plt.tight_layout()
        out_path = os.path.join(outputs_dir, "ROC.png")
        plt.savefig(out_path, dpi=200)
        plt.close()
        print(f"Saved ROC curve to {out_path}")


    def run_calculate_metrics(self, predicted_results_dir, cutoff, bindingSite_dataframe):
        # Compute metrics over multiple thresholds and return as DataFrame
        thresholds = np.linspace(0.0, 1.0, 21)
        rows = []
        for thr in thresholds:
            acc, prec, tpr, fpr = self.calculate_metrics(predicted_results_dir, thr, bindingSite_dataframe)
            rows.append({
                "cutoff": thr,
                "accuracy": acc,
                "precision": prec,
                "tpr": tpr,
                "fpr": fpr
            })
        return pd.DataFrame(rows)


if __name__ == "__main__":

    # In main: using cutoff = 0.5, plot a histogram of the frequency that each amino acid is predicted to be an ATP binding site across all proteins. Save this plot in the output directory. 
    # Run calculate metrics for at least 10 cutoff points within the range 0 to 1 and plot an ROC curve. Save this curve as ROC.png in the outputs folder.
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    predicted_results_dir = os.path.join(base_dir, "predict_results")
    outputs_dir = os.path.join(base_dir, "outputs")
    binding_sites_file = os.path.join(base_dir, "binding_sites.csv")

    analyzer = AnalyzePredictions()
    cutoff = 0.5
    print(f"[AnalyzePredictions] Using predicted results from: {predicted_results_dir}")
    print(f"[AnalyzePredictions] Outputs directory: {outputs_dir}")

    # Histogram of amino acids predicted as binding at cutoff
    analyzer.plot_histogram(predicted_results_dir, cutoff)

    # Load binding site labels
    binding_df = analyzer.get_bindingSite_labels(binding_sites_file)

    # Plot ROC curve (sweeps thresholds internally)
    analyzer.plot_roc_curve(predicted_results_dir, cutoff, binding_df)

    # Also dump metrics across thresholds to CSV for reference
    metrics_df = analyzer.run_calculate_metrics(predicted_results_dir, cutoff, binding_df)
    metrics_csv = os.path.join(outputs_dir, "metrics_by_cutoff.csv")
    metrics_df.to_csv(metrics_csv, index=False)
    print(f"Saved metrics across cutoffs to {metrics_csv}")