#!/usr/bin/env python3
import os


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    top100_path = os.path.join(here, "mystery_protein_100.txt")

    max_score = None
    with open(top100_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            try:
                score = float(parts[4])
            except ValueError:
                continue
            if max_score is None or score > max_score:
                max_score = score

    if max_score is None:
        print("No scores found.")
    else:
        print(f"{max_score:.3f}")


if __name__ == "__main__":
    main()

