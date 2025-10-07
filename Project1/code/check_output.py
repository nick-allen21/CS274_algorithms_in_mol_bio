import os
from collections import Counter


def _fuzzy_equal(a: float, b: float, eps: float = 1e-6) -> bool:
    return abs(a - b) < eps


def parse_alignment_output(path: str):
    """
    Parse an alignment output file.

    Returns:
      (score: float, pairs: Counter[(str, str)])
    """
    with open(path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f.readlines()]

    # first non-empty line is the score
    i = 0
    while i < len(lines) and lines[i].strip() == "":
        i += 1
    if i == len(lines):
        raise ValueError(f"No score line found in {path}")
    try:
        score = float(lines[i].strip())
    except ValueError:
        raise ValueError(f"Score line is not a float in {path}: '{lines[i]}'")

    # advance past the blank after score if present
    i += 1
    if i < len(lines) and lines[i].strip() == "":
        i += 1

    pairs = []
    while i < len(lines):
        # skip any extra blank separators
        while i < len(lines) and lines[i].strip() == "":
            i += 1
        if i >= len(lines):
            break
        a = lines[i].strip()
        i += 1
        if i >= len(lines):
            raise ValueError(f"Unpaired alignment line at end of file {path}")
        b = lines[i].strip()
        i += 1
        pairs.append((a, b))
        # optional blank line after each pair
        if i < len(lines) and lines[i].strip() == "":
            i += 1

    return score, Counter(pairs)


def compare_outputs(answer_path: str, example_path: str) -> bool:
    ans_score, ans_pairs = parse_alignment_output(answer_path)
    ex_score, ex_pairs = parse_alignment_output(example_path)

    ok = True
    if not _fuzzy_equal(ans_score, ex_score):
        print(f"Score mismatch: answer={ans_score} example={ex_score}")
        ok = False

    if ans_pairs != ex_pairs:
        missing = ex_pairs - ans_pairs  # in example but not in answer
        extra = ans_pairs - ex_pairs    # in answer but not in example
        if missing:
            print("Missing pairs (expected but not found):")
            for (a, b), cnt in missing.items():
                print(f"  x{cnt}: {a}\n      {b}")
        if extra:
            print("Extra pairs (found but not expected):")
            for (a, b), cnt in extra.items():
                print(f"  x{cnt}: {a}\n      {b}")
        ok = False

    return ok


def main():
    here = os.path.dirname(__file__)
    examples_dir = os.path.normpath(os.path.join(here, "..", "examples"))

    # case suffixes: "" (0) and 1..8
    suffixes = [""] + [str(i) for i in range(1, 8 + 1)]

    all_ok = True
    for sfx in suffixes:
        ans = os.path.join(examples_dir, f"alignment_answer{sfx}.output" if sfx else "alignment_answer.output")
        ex = os.path.join(examples_dir, f"alignment_example{sfx}.output" if sfx else "alignment_example.output")
        if not (os.path.exists(ans) and os.path.exists(ex)):
            print(f"Skipping case '{sfx or '0'}' (files not found)")
            all_ok = False
            continue
        print(f"Checking case '{sfx or '0'}':")
        ok = compare_outputs(ans, ex)
        print("Match" if ok else "Mismatch", "\n")
        all_ok = all_ok and ok

    if not all_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()


