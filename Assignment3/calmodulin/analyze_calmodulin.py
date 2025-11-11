#!/usr/bin/env python3
import os
import math

HERE = os.path.dirname(os.path.abspath(__file__))

def _resolve_pdb_path():
    cand1 = os.path.join(HERE, "Calmodulin.pdb")
    cand2 = os.path.join(os.path.dirname(HERE), "pdb", "Calmodulin.pdb")
    return cand1 if os.path.isfile(cand1) else cand2

def parse_ca_sites_from_pdb(pdb_path):
    sites = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not line.startswith("HETATM"):
                continue
            # Robust calcium detection:
            # - standard PDB fields:
            #   atom name  (cols 13-16)   -> line[12:16]
            #   residue    (cols 18-20)   -> line[17:20]
            #   element    (cols 77-78)   -> line[76:78]
            resname = line[17:20]           # residue name (often "CA ")
            atomname = line[12:16]          # atom name (often "CA")
            element = line[76:78] if len(line) >= 78 else ""  # element symbol
            is_ca = (
                resname.strip() == "CA" or
                atomname.strip() == "CA" or
                element.strip() == "CA"
            )
            if is_ca:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                sites.append((x, y, z))
    return sites

def parse_top_points(txt_path, k=10):
    pts = []
    with open(txt_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                score = float(parts[4])
            except ValueError:
                continue
            pts.append((x, y, z, score))
    pts.sort(key=lambda t: t[3], reverse=True)
    return pts[:k]

def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def count_represented_sites(sites, points, radius=2.0):
    represented = 0
    details = []
    for i, s in enumerate(sites):
        dmin = min(dist(s, (p[0], p[1], p[2])) for p in points) if points else float("inf")
        hit = dmin <= radius
        represented += 1 if hit else 0
        details.append((i, s, dmin, hit))
    return represented, details

def main():
    pdb_path = _resolve_pdb_path()
    top_file = os.path.join(HERE, "Calmodulin_100.txt")

    sites = parse_ca_sites_from_pdb(pdb_path)
    points = parse_top_points(top_file, k=10)
    n_repr, details = count_represented_sites(sites, points, radius=2.0)

    print(f"Found {len(sites)} annotated Ca2+ sites in PDB: {pdb_path}")
    print(f"Evaluating top {len(points)} predicted points from: {top_file}")
    print(f"Radius threshold: 2.0 Å")
    print()
    for i, (idx, s, dmin, hit) in enumerate(details):
        flag = "HIT" if hit else "MISS"
        print(f"Site {idx+1}: ({s[0]:.3f}, {s[1]:.3f}, {s[2]:.3f})  min_dist_to_top10 = {dmin:.3f} Å  -> {flag}")
    print()
    print(f"Number of real Ca2+ sites represented by the top 10 predictions: {n_repr}")

if __name__ == "__main__":
    main()