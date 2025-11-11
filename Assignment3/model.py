# 0) import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
import sys
import warnings
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import BernoulliNB
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
from sklearn.feature_selection import mutual_info_classif
warnings.filterwarnings("ignore")

# 1) Read in the data and convert the features and outcome variables to binary (1 or 0) instead of y/n and site/non-site.

def _get_features_csv_path() -> str:
    """
    Resolve the path to features.csv.
    Prefer a CLI-provided path; otherwise use the absolute workspace path.
    """
    if len(sys.argv) > 1 and os.path.isfile(sys.argv[1]):
        return os.path.abspath(sys.argv[1])
    # Absolute path per workspace layout
    return "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Assignment3/features.csv"


# Read CSV
csv_path = _get_features_csv_path()
df = pd.read_csv(csv_path)

# Split features and target
feature_columns = [c for c in df.columns if c != "SITE"]

# Convert y/n -> 1/0 for features
X = df[feature_columns].replace({"y": 1, "n": 0}).astype(int)

# Convert site/nonsite -> 1/0 for target
y = df["SITE"].replace({"site": 1, "nonsite": 0}).astype(int)

# Optional combined binary DataFrame for convenience
df_binary = X.copy()
df_binary["SITE"] = y

print(f"Loaded {len(df_binary)} rows from {os.path.basename(csv_path)}. "
      f"X shape: {X.shape}, positives in y: {int(y.sum())}")

#2) Divide the data into 80% training and 20% test.
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=1115)


#3) Build randoms bayes model using BernoulliNB.
nb = BernoulliNB()
nb.fit(X_train, y_train)

# 4) Predict on test set and produce confusion matrix
y_pred = nb.predict(X_test)
cm = confusion_matrix(y_test, y_pred, labels=[0, 1])
cm_df = pd.DataFrame(cm, index=["actual_nonsite", "actual_site"], columns=["pred_nonsite", "pred_site"])
acc = accuracy_score(y_test, y_pred)

print("Confusion matrix (rows: actual, cols: predicted):")
print(cm_df)
print(f"Accuracy: {acc:.4f}")
print("Classification report:")
print(classification_report(y_test, y_pred, target_names=["nonsite", "site"]))
tn, fp, fn, tp = cm[0, 0], cm[0, 1], cm[1, 0], cm[1, 1]
tpr = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
fpr = fp / (fp + tn) if (fp + tn) > 0 else float("nan")
print(f"True Positive Rate (Recall for 'site'): {tpr:.4f}")
print(f"False Positive Rate: {fpr:.4f}")

# 5) Feature evaluation via Mutual Information (features are discrete/binary)
mi_scores = mutual_info_classif(X_train, y_train, discrete_features=True, random_state=1115)
mi_series = pd.Series(mi_scores, index=feature_columns).sort_values(ascending=False)

print("Top 15 features by mutual information:")
print(mi_series.head(15))

# Optionally save full ranking
mi_out_path = os.path.join(os.path.dirname(csv_path), "feature_mi_ranking.csv")
mi_series.to_csv(mi_out_path, header=["mutual_information"])
print(f"Saved full MI ranking to {mi_out_path}")

# Additional answers required by the assignment:
# Top 5 MI features
top5 = mi_series.head(5)
print("Top 5 MI features:")
print(top5)

# Compare shell 1 vs shell 4 (aggregate MI across all amino acids in each shell)
def shell_mi(shell_idx: int) -> float:
    shell_cols = [c for c in feature_columns if c.endswith(str(shell_idx))]
    return float(mi_series.loc[shell_cols].mean()) if len(shell_cols) > 0 else float("nan")

shell1_mi = shell_mi(1)
shell4_mi = shell_mi(4)
more_predictive_shell = "shell 1" if shell1_mi >= shell4_mi else "shell 4"
print(f"Mean MI for shell 1: {shell1_mi:.6f}")
print(f"Mean MI for shell 4: {shell4_mi:.6f}")
print(f"More predictive: {more_predictive_shell}")

# Compare amino acids ASP vs CYS across all shells
def aa_total_mi(aa: str) -> float:
    cols = [c for c in feature_columns if c.startswith(aa)]
    return float(mi_series.loc[cols].sum()) if len(cols) > 0 else float("nan")

asp_mi = aa_total_mi("ASP")
cys_mi = aa_total_mi("CYS")
more_predictive_aa = "ASP" if asp_mi >= cys_mi else "CYS"
print(f"Total MI for ASP: {asp_mi:.6f}")
print(f"Total MI for CYS: {cys_mi:.6f}")
print(f"More predictive amino acid: {more_predictive_aa}")