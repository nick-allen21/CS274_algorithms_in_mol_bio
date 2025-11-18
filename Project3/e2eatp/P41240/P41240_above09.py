# /Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project3/e2eatp/P41240/P41240_above09.py
import os
import pandas as pd

path = "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project3/predict_results/P41240.tsv"

df = pd.read_csv(path, sep="\t")
df["Prob"] = pd.to_numeric(df["Prob"], errors="coerce")

hi = df[df["Prob"] >= 0.9].copy()
idxs = hi["Index"].astype(int).tolist()
print("Indexes with Prob >= 0.9:", idxs)

out_dir = os.path.dirname(path)
out_path = os.path.join(out_dir, "P41240_above09.tsv")
hi.to_csv(out_path, sep="\t", index=False)
print(f"Saved {len(hi)} rows to {out_path}")