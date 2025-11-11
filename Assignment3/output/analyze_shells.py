#!/usr/bin/env python3
import os

# Map 3-letter residue codes to full names
RESIDUE_FULL_NAMES = {
    "ALA": "Alanine",
    "ARG": "Arginine",
    "ASN": "Asparagine",
    "ASP": "Aspartic acid",
    "CYS": "Cysteine",
    "GLN": "Glutamine",
    "GLU": "Glutamic acid",
    "GLY": "Glycine",
    "HIS": "Histidine",
    "ILE": "Isoleucine",
    "LEU": "Leucine",
    "LYS": "Lysine",
    "MET": "Methionine",
    "PHE": "Phenylalanine",
    "PRO": "Proline",
    "SER": "Serine",
    "THR": "Threonine",
    "TRP": "Tryptophan",
    "TYR": "Tyrosine",
    "VAL": "Valine",
}


def top_n_nonsite_residues_for_shell(shell_index: int = 4, n: int = 4):
    """
    Read shell_{shell_index}.txt and return the top-n residues by NonSites frequency.
    Returns a list of tuples: (aa_code, full_name, nonsites_frequency).
    """
    here = os.path.dirname(os.path.abspath(__file__))
    shell_file = os.path.join(here, f"shell_{shell_index}.txt")

    residues = []
    with open(shell_file, "r") as f:
        for i, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            if i == 0:
                # header: "AA\tSites\tNonSites"
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            aa, sites_str, nonsites_str = parts[0], parts[1], parts[2]
            try:
                nonsites = float(nonsites_str)
            except ValueError:
                continue
            residues.append((aa, RESIDUE_FULL_NAMES.get(aa, aa), nonsites))

    # Sort by NonSites frequency descending and take top 4
    residues.sort(key=lambda x: x[2], reverse=True)
    return residues[:n]


def _read_scores(filepath: str):
    """
    Read a *scores.txt file where each line is:
    PDB\tX\tY\tZ\tScore
    Returns a list of tuples: (pdb, x, y, z, score_float).
    """
    records = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            pdb = parts[0]
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                score = float(parts[4])
            except ValueError:
                continue
            records.append((pdb, x, y, z, score))
    return records


def summarize_score_extrema():
    """
    Answer two similar questions by reading the score files:
    - Worst calcium binding site (lowest score) among CAsites_scores.txt
    - Best (highest) score among CAnonsites_scores.txt
    Prints the PDB filename and associated score (and coordinates).
    """
    here = os.path.dirname(os.path.abspath(__file__))
    ca_sites_path = os.path.join(here, "CAsites_scores.txt")
    ca_nonsites_path = os.path.join(here, "CAnonsites_scores.txt")

    site_records = _read_scores(ca_sites_path)
    nonsite_records = _read_scores(ca_nonsites_path)

    if site_records:
        worst_site = min(site_records, key=lambda r: r[4])
        print("Worst (lowest) calcium-site score in CAsites_scores.txt:")
        print(f"- {worst_site[0]} at ({worst_site[1]:.3f}, {worst_site[2]:.3f}, {worst_site[3]:.3f}) -> {worst_site[4]:.3f}")
    else:
        print("No site records found.")

    if nonsite_records:
        best_nonsite = max(nonsite_records, key=lambda r: r[4])
        print("Best (highest) nonsite score in CAnonsites_scores.txt:")
        print(f"- {best_nonsite[0]} at ({best_nonsite[1]:.3f}, {best_nonsite[2]:.3f}, {best_nonsite[3]:.3f}) -> {best_nonsite[4]:.3f}")
    else:
        print("No nonsite records found.")


def main():
    # Q1: Top-4 residues in shell 4 by NonSites frequency
    top4 = top_n_nonsite_residues_for_shell(shell_index=4, n=4)
    print("Top 4 residues in shell 4 by NonSites frequency:")
    for rank, (aa, full, freq) in enumerate(top4, start=1):
        print(f"{rank}. {aa} ({full}): {freq:.3f}")

    print()

    # Q2 & Q3: Extrema from score files
    summarize_score_extrema()


if __name__ == "__main__":
    main()

