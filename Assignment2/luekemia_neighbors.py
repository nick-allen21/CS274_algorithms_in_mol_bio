"""
Now we will perform classification using the K-Nearest Neighbors
(KNN) algorithm.

Train a KNN classifier with KNeighborsClassifier function from scikit-learn with default settings and K=5.  Perform leave-one-out cross-validation to get an estimate of test error and generate a confusion matrix for the classifier you built.
"""

import pandas as pd
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import confusion_matrix, accuracy_score, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

# Load the data
CSV_PATH = "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Assignment2/leukemia.csv"
df = pd.read_csv(CSV_PATH)
# number of neighbors
k = 5

# Split features/labels
y = df["leukemia_type"].astype(str)
x = df.drop(columns=["leukemia_type"])

# Build pipeline: scale -> KNN
pipeline = Pipeline([
    ("scaler", StandardScaler()),
    ("knn", KNeighborsClassifier(n_neighbors=k))
])

# Leave-One-Out CV predictions with progress bar
loo = LeaveOneOut()
y_pred = []
for train_idx, test_idx in tqdm(loo.split(x), total=len(x), desc="LOOCV"):
    x_train, x_test = x.iloc[train_idx], x.iloc[test_idx]
    y_train = y.iloc[train_idx]
    pipeline.fit(x_train, y_train)
    pred = pipeline.predict(x_test)[0]
    y_pred.append(pred)

y_pred = np.array(y_pred)

# alternatively, we can use cross_val_predict
# y_pred = cross_val_predict(pipeline, x, y, cv=loo)
# acc = accuracy_score(y, y_pred)
# labels = sorted(np.unique(y))
# cm = confusion_matrix(y, y_pred, labels=labels)

# Q1) What is the percentage of correctly classified instances? In other words, what is the accuracy of this classifier?# Metrics
acc = accuracy_score(y, y_pred)
labels = sorted(np.unique(y))
cm = confusion_matrix(y, y_pred, labels=labels)

# Report
print(f"LOOCV Accuracy (K={k}): {acc:.2%}")
print(f"Labels order: {labels}")
print("Confusion matrix (rows=true, cols=pred):")
print(cm)
print(f"Percentage of correctly classified instances: {acc:.2%}")