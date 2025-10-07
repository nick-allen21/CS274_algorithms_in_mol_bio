import sys
from align_quiz_functions import count_mismatches, count_gaps


def main():
    if len(sys.argv) != 2:
        print("Usage: python count_gaps_n_mismatches.py <alignment_output_file>")
        return
    output_file = sys.argv[1]

    count_mismatches(output_file)
    count_gaps(output_file)
    
if __name__ == "__main__":
    main()


