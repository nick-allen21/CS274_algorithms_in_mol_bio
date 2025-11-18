import pandas as pd
from Bio import SeqIO

fasta = "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project3/kinase.fasta"
labels = "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project3/binding_sites.csv"

# All proteins from FASTA (UniProt IDs)
all_ids = set()
for rec in SeqIO.parse(fasta, "fasta"):
    hdr = rec.description or rec.id
    parts = hdr.split("|")
    uid = parts[1] if len(parts) > 1 else hdr.split()[0]
    all_ids.add(uid)

# Proteins with known ATP sites
known_ids = set(pd.read_csv(labels)["protid"].astype(str).unique())

known = all_ids & known_ids
unknown = all_ids - known_ids

print(f"Total proteins: {len(all_ids)}")
print(f"Known ATP-binding annotated: {len(known)}")
print(f"Unknown: {len(unknown)}")
print(f"Ratio known:unknown = {len(known)}:{len(unknown)} ({len(known)/max(1,len(unknown)):.3f})")
print(f"Fraction known = {len(known)/len(all_ids):.3%}")