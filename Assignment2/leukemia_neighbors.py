"""
Now we will perform classification using the K-Nearest Neighbors
(KNN) algorithm.

Train a KNN classifier with KNeighborsClassifier function from scikit-learn with default settings and K=5.  Perform leave-one-out cross-validation to get an estimate of test error and generate a confusion matrix for the classifier you built.
"""

import pandas as pd
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.naive_bayes import GaussianNB

# Load the data
leukemia_csv_path = "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Assignment2/leukemia.csv"
df = pd.read_csv(leukemia_csv_path)
# number of neighbors
k = 5

# Split features/labels
y = df["leukemia_type"].astype(str)
x = df.drop(columns=["leukemia_type"])

# Build pipeline: scale -> KNN
pipeline = Pipeline([
    ("knn", KNeighborsClassifier(n_neighbors=k))
])

# Leave-One-Out CV predictions with progress bar
loo = LeaveOneOut()

# use cross_val_predict
y_pred = cross_val_predict(pipeline, x, y, cv=loo)

# Q1) What is the percentage of correctly classified instances? In other words, what is the accuracy of this classifier?# Metrics
acc = accuracy_score(y, y_pred)
labels = sorted(np.unique(y))
cm = confusion_matrix(y, y_pred, labels=labels)

# Report
print("\n===============================[K-Nearest Neighbors - K={k}]===============================")
print(f"LOOCV Accuracy (K={k}): {acc:.2%}")
print(f"Labels order: {labels}")
print("Confusion matrix (rows=true, cols=pred):")
print(cm)
print(f"Percentage of correctly classified instances: {acc:.2%}")


# Q2) compute specificity, sensitivity, and accuracy for each class
def compute_binary_rates(y_true: pd.Series, y_pred: np.ndarray, positive_label: str):
    """Compute TP, FP, TN, FN and derived rates for a given positive label."""
    y_true_pos = (y_true == positive_label).to_numpy()
    y_pred_pos = (y_pred == positive_label)

    tp = int(np.sum(y_true_pos & y_pred_pos))
    fp = int(np.sum(~y_true_pos & y_pred_pos))
    tn = int(np.sum(~y_true_pos & ~y_pred_pos))
    fn = int(np.sum(y_true_pos & ~y_pred_pos))

    tp_denom = tp + fn
    fp_denom = fp + tn

    tpr = (tp / tp_denom) if tp_denom > 0 else float("nan")  # sensitivity/recall
    fpr = (fp / fp_denom) if fp_denom > 0 else float("nan")
    specificity = 1.0 - fpr if fp_denom > 0 else float("nan")

    return {
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "tpr": tpr,
        "fpr": fpr,
        "specificity": specificity,
    }


rates_all = compute_binary_rates(y, y_pred, positive_label="ALL")
rates_aml = compute_binary_rates(y, y_pred, positive_label="AML")

print("\nPer-class rates (treating the class as positive):")
print(f"ALL TP rate (sensitivity): {rates_all['tpr']:.2f}")
print(f"ALL FP rate: {rates_all['fpr']:.2f}")
print(f"ALL Specificity: {rates_all['specificity']:.2f}")
print(f"AML TP rate (sensitivity): {rates_aml['tpr']:.2f}")
print(f"AML FP rate: {rates_aml['fpr']:.2f}")
print(f"AML Specificity: {rates_aml['specificity']:.2f}")

# Naive Bayes Classifier
# Leave-One-Out CV with GaussianNB (default settings)
nb = GaussianNB()
y_pred_nb = cross_val_predict(nb, x, y, cv=loo)

acc_nb = accuracy_score(y, y_pred_nb)
labels_nb = sorted(np.unique(y))
cm_nb = confusion_matrix(y, y_pred_nb, labels=labels_nb)

print("\n===============================[Naive Bayes - GaussianNB]===============================")
print(f"LOOCV Accuracy (GaussianNB): {acc_nb:.2%}")
print(f"Labels order: {labels_nb}")
print("Confusion matrix (rows=true, cols=pred):")
print(cm_nb)
print(f"Percentage of correctly classified instances: {acc_nb:.2%}")

rates_all_nb = compute_binary_rates(y, y_pred_nb, positive_label="ALL")
rates_aml_nb = compute_binary_rates(y, y_pred_nb, positive_label="AML")

print("\nPer-class rates (treating the class as positive):")
print(f"ALL TP rate (sensitivity): {rates_all_nb['tpr']:.2f}")
print(f"ALL FP rate: {rates_all_nb['fpr']:.2f}")
print(f"ALL Specificity: {rates_all_nb['specificity']:.2f}")
print(f"AML TP rate (sensitivity): {rates_aml_nb['tpr']:.2f}")
print(f"AML FP rate: {rates_aml_nb['fpr']:.2f}")
print(f"AML Specificity: {rates_aml_nb['specificity']:.2f}")
