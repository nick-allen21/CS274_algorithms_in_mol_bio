# Which amino acids are most common in ATP binding sites (Project3)
import pandas as pd
from collections import Counter
from Bio import SeqIO

FASTA = "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project3/kinase.fasta"
BINDING = "/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project3/binding_sites.csv"

# Load sequences keyed by UniProt ID from FASTA headers: >sp|<ID>|...
seqs = {}
for rec in SeqIO.parse(FASTA, "fasta"):
    hdr = rec.description or rec.id
    uid = hdr.split("|")[1] if "|" in hdr else hdr.split()[0]
    seqs[uid] = str(rec.seq)

# Load binding sites
df = pd.read_csv(BINDING)  # columns: protid, binding_site, seq_length
cnt = Counter()
n_used = 0
n_skipped = 0

for _, row in df.iterrows():
    pid = str(row["protid"])
    idx = int(row["binding_site"])
    s = seqs.get(pid)
    if not s:
        n_skipped += 1
        continue
    # Assume 0-based, but fall back to 1-based if out of range
    if 0 <= idx < len(s):
        aa = s[idx]
    elif 1 <= idx <= len(s):
        aa = s[idx - 1]
    else:
        n_skipped += 1
        continue
    cnt[aa] += 1
    n_used += 1

print(f"Counted {n_used} binding residues (skipped {n_skipped}).")
print("AA frequencies (most common first):")
total = sum(cnt.values()) or 1
for aa, c in cnt.most_common():
    print(f"{aa}\t{c}\t{c/total:.2%}")