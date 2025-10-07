import sys
from align_quiz_functions import compare_alignments


def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_alignments.py <alignment_output_file1> <alignment_output_file2>")
        return
    output_file1 = sys.argv[1]
    output_file2 = sys.argv[2]

    compare_alignments(output_file1, output_file2)
    
if __name__ == "__main__":
    main()


