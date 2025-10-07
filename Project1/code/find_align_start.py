# /Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project1/code/find_align_start.py
import sys
from align_quiz_functions import find_align_location

def main():
    if len(sys.argv) != 3:
        print("Usage: python find_align_start.py <input_file> <output_file>")
        return
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    loc = find_align_location(input_file, output_file)
    print("Final alignment start indices (seq_a, seq_b):", loc)

if __name__ == "__main__":
    main()