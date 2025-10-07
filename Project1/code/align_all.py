import os
import sys
import subprocess

def run_all():
    base_dir = os.path.dirname(__file__)
    examples_dir = os.path.normpath(os.path.join(base_dir, "..", "examples"))
    align_py = os.path.join(base_dir, "align.py")

    suffixes = [""] + [str(i) for i in range(1, 9)]
    results = []

    for sfx in suffixes:
        in_name = f"alignment_example{sfx}.input" if sfx else "alignment_example.input"
        out_name = f"alignment_answer{sfx}.output" if sfx else "alignment_answer.output"
        in_path = os.path.join(examples_dir, in_name)
        out_path = os.path.join(examples_dir, out_name)

        if not os.path.exists(in_path):
            results.append((sfx or "0", False, f"missing input {in_path}"))
            continue

        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        subprocess.run([sys.executable, align_py, in_path, out_path], check=True)
        results.append((sfx or "0", True, out_path))

    return results

def main():
    runs = run_all()
    for tag, ok, msg in runs:
        print(f"[{tag}] {'ok' if ok else 'skip'}: {msg}")

    # Run the comparer
    from check_output import main as check_main
    check_main()

if __name__ == "__main__":
    main()