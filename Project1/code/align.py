
"""

This file provides skeleton code for align.py. 

Locations with "FILL IN" in comments are where you need to add code.

Note - you MUST follow this structure or else the autograder will not run properly

Usage: python align.py input_file output_file

"""


import glob
import sys

class pointer(object):
    """
    Object to store a pointer in the score matrix.
    """
    def __init__(self, row, col, name):
        self.row = row
        self.col = col
        self.name = name

    def __repr__(self):
        return f"{self.name}({self.row},{self.col})"

    # define equality for pointer objects
    def __eq__(self, other):
        return isinstance(other, pointer) and \
               self.row == other.row and self.col == other.col and self.name == other.name


class ScoreEntry(object):
    """
    Object to store a score entry in the score matrix.
    """
    def __init__(self, row, col, score, matrix_name):
        self.row = row
        self.col = col
        self.score = score
        self.matrix_name = matrix_name
        # (row, col, name). Important that we set what matrix we came from
        self.pointers = []


#### ------ USEFUL FUNCTIONS ------- ####

def fuzzy_equals(a, b):
    """
    Checks if two floating point numbers are equivalent.
    """
    epsilon = 10**(-6) 
    return (abs(a - b) < epsilon)
    

def get_maxes(score_entry_list : list[ScoreEntry]):
    """
    This function takes in a list of score entries and returns the max scores, adjusting for equal values.
    The scores in the list will already be adjusted for the addition of the score and gap penalties. The 
    row and column in the list will be the row and column that we just came from 
    
    This function will return a score and the array of pointers that give us the max score
    """

    # if the list is empty, return 0 and an empty list
    if not score_entry_list:
        return 0, []

    # Initialize the max score and pointer array
    max_score = max(score_entry.score for score_entry in score_entry_list)
    max_pointers = []

    # Iterate through the list of score entries using special floating point logic 
    for score_entry in score_entry_list:
        if fuzzy_equals(score_entry.score, max_score):
            max_pointers.append(pointer(score_entry.row, score_entry.col, score_entry.matrix_name))
    return max_score, max_pointers


#### ------- CLASSES ------- ####

class MatchMatrix(object):
    """
    Match matrix class stores the scores of matches in a data structure
    """
    def __init__(self):
        # Init as dict mapping tuples of characters to scores
        self.match_matrix_scores = {}

    def set_score(self, a, b, score):
        """
        Updates or adds a score for a specified match

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
           score = the score to set it for
        """
        # key will be a tuple of the two characters
        self.match_matrix_scores[(a, b)] = score

    def get_score(self, a, b):
        """
        Returns the score for a particular match, where a is the
        character from sequence a and b is from sequence b.

        Inputs:
           a = the character from sequence A
           b = the character from sequence B
        Returns:
           the score of that match
        """
        # return value form the hash table
        return self.match_matrix_scores[(a, b)]


class ScoreMatrix(object):
    """
    Object to store a score matrix, which generated during the alignment process. The score matrix consists of a 2-D array of
    ScoreEntries that are updated during alignment and used to output the maximum alignment.
    """

    def __init__(self, name, nrow, ncol):
        self.name = name # identifier for the score matrix - Ix, Iy, or M
        self.nrow = nrow
        self.ncol = ncol

        # initialize the score matrix at zeroes for all entries, we don't penalize the end/start gaps
        # down the line, we won't recompute these boundaries, first entry we compute is (1,1)
        self.score_matrix = [[ScoreEntry(row, col, float(0), self.name) for col in range(ncol)] for row in range(nrow)]
    
    def get_score(self, row, col):
        # return the score entry obj for the given row and column
        return self.score_matrix[row][col].score

    def get_score_obj(self, row, col):
        # return the score value for the given row and column
        return self.score_matrix[row][col]

    def set_score(self, row, col, score):
        # set the score for the given row and column
        self.score_matrix[row][col] = ScoreEntry(row, col, score, self.name)

    def get_pointers(self, row, col):
        """
        Returns the indices of the entries that are pointed to
        This should be formatted as a list of tuples:
         ex. [(1,1, "M"), (1,0, "Ix")]
        """
        return self.score_matrix[row][col].pointers

    def set_pointers(self, row, col, pointers: list[pointer]):
        # set the pointers for the given row and column for score entry object
        self.score_matrix[row][col].pointers = pointers


    def print_scores(self):
        """
        Returns a nicely formatted string containing the scores in the score matrix. Use this for debugging!

        Example:
        M=
            0.0, 0.0, 0.0, 0.0, 0.0
            0.0, 1.0, 0.0, 0.0, 0.0
            0.0, 1.0, 1.0, 1.0, 1.0
            0.0, 0.0, 1.0, 1.0, 1.0
            0.0, 0.0, 2.0, 2.0, 1.0
            0.0, 0.0, 1.0, 2.0, 3.0

        """
        # format all cells to 2 decimals
        cells = [[f"{self.score_matrix[r][c].score:.2f}" for c in range(self.ncol)] for r in range(self.nrow)]
        # column widths
        col_w = [max(len(cells[r][c]) for r in range(self.nrow)) for c in range(self.ncol)]
        # header
        print(f"{self.name}:")
        # rows
        for r in range(self.nrow):
            line = "  ".join(cells[r][c].rjust(col_w[c]) for c in range(self.ncol))
            print("  " + line)


    def print_pointers(self):
        # Build display strings per cell
        cells = []
        for r in range(self.nrow):
            row_cells = []
            for c in range(self.ncol):
                pts = getattr(self.score_matrix[r][c], "pointers", [])
                if not pts:
                    s = "∅"
                else:
                    s = " ".join(f"{p.name}({p.row},{p.col})" for p in pts)
                row_cells.append(s)
            cells.append(row_cells)

        # Column widths
        col_w = [max(len(cells[r][c]) for r in range(self.nrow)) for c in range(self.ncol)]

        print(f"{self.name} pointers:")
        for r in range(self.nrow):
            line = "  ".join(cells[r][c].ljust(col_w[c]) for c in range(self.ncol))
            print("  " + line)


class AlignmentParameters(object):
    """
    Object to hold a set of alignment parameters from an input file.
    """

    def __init__(self):
        # default values for variables that are filled in by reading
        # the input alignment file
        self.seq_a = ""
        self.seq_b = ""
        self.global_alignment = False 
        self.local_alignment = False
        self.dx = 0 # open gap penalty in A (I_x, skip row)
        self.ex = 0 # extend gap penalty in A (I_x, skip row)
        self.dy = 0 # open gap penalty in B (I_y, skip col)
        self.ey = 0 # extend gap penalty in B (I_y, skip col)
        self.alphabet_a = "" 
        self.alphabet_b = ""
        self.len_alphabet_a = 0
        self.len_alphabet_b = 0
        self.match_matrix = MatchMatrix()

    def load_params_from_file(self, input_file): 
        """
        Reads the parameters from an input file and stores in the object

        Input:
           input_file = specially formatted alignment input file
        """
        # load the alignment parameters into the align_params object
        with open(input_file, 'r') as f:
            # keep blank-line tolerant; strip trailing newlines
            lines = [ln.strip() for ln in f if ln.strip() != '']
        it = iter(lines)

        # sequences
        self.seq_a = next(it)
        self.seq_b = next(it)

        # global vs local
        self.global_alignment = (next(it) == '0')
        self.local_alignment = not self.global_alignment

        # gap penalties
        self.dx, self.ex, self.dy, self.ey = map(float, next(it).split())

        # alphabets
        self.len_alphabet_a = int(next(it)); self.alphabet_a = next(it)
        self.len_alphabet_b = int(next(it)); self.alphabet_b = next(it)

        # create match matrix
        for row in range(1, self.len_alphabet_a + 1):
            for col in range(1, self.len_alphabet_b + 1):

                # getting indices and scores from next line in input file
                i, j, a, b, s = next(it).split()

                # checking to make sure the indices are correct
                assert int(i)== row and int(j) == col
                self.match_matrix.set_score(a, b, float(s))


class Align(object):
    """
    Object to hold and run an alignment; running is accomplished by using "align()"
    """

    def __init__(self, input_file, output_file):
        """
        Input:
            input_file = file with the input for running an alignment
            output_file = file to write the output alignments to
        """
        self.input_file = input_file
        self.output_file = output_file
        self.align_params = AlignmentParameters() 

        # Note the below three lines is ensure the autograder runs properly.
        # You should leave the below three lines as is but then
        # Initialize m_matrix, ix_matrix, and iy_matrix in populate_score_matrices
        self.m_matrix = None
        self.ix_matrix = None
        self.iy_matrix = None

        self.max_score = 0
        self.max_loc = set()
        self.paths = []

    def align(self):
        """
        Main method for running alignment.
        """

        # load the alignment parameters into the align_params object
        self.align_params.load_params_from_file(self.input_file)

        # populate the score matrices based on the input parameters
        self.populate_score_matrices()
        # final score matrixes 
        print("\n2) Final Score Matrices:")
        self.m_matrix.print_scores()
        self.ix_matrix.print_scores()
        self.iy_matrix.print_scores()

        self.m_matrix.print_pointers()
        self.ix_matrix.print_pointers()
        self.iy_matrix.print_pointers()

        # perform a traceback and write the output to an output file

        ### TO-DO! FILL IN ###
        self.max_score, self.max_loc = self.find_traceback_start()
        print("Max score: ", self.max_score, "\nMax loc: ", self.max_loc)
        self.paths = self.traceback()
        # self.print_paths()
        self.write_output()

    def populate_score_matrices(self):
        """
        Method to populate the score matrices based on the data in align_params.
        Should call update(i,j) for each entry in the score matrices
        Note: You MUST initialize M, Ix, Iy in this function rather than elsewhere
        """

        # initialize the score matrices, will all start at 0, we won't recompute the boundaries
        # col / row 0 represents starting w a gap
        nrow = len(self.align_params.seq_a) + 1
        ncol = len(self.align_params.seq_b) + 1
        self.m_matrix = ScoreMatrix("M", nrow, ncol)
        self.ix_matrix = ScoreMatrix("Ix", nrow, ncol)
        self.iy_matrix = ScoreMatrix("Iy", nrow, ncol)

        # forbid starting in gap states
        neg_inf = float('-inf')

        # forbid starting in gap states
        for j in range(self.ix_matrix.ncol):
            self.ix_matrix.set_score(0, j, neg_inf)
        for i in range(self.ix_matrix.nrow):
            self.ix_matrix.set_score(i, 0, neg_inf)

        for j in range(self.iy_matrix.ncol):
            self.iy_matrix.set_score(0, j, neg_inf)
        for i in range(self.iy_matrix.nrow):
            self.iy_matrix.set_score(i, 0, neg_inf)

        # score matrixes 
        print("\n1) Intitial Score Matrices:")
        print("\nM: \n")
        self.m_matrix.print_scores()
        print("\nIx: \n")
        self.ix_matrix.print_scores()
        print("\nIy: \n")
        self.iy_matrix.print_scores()

        # update the score matrices 
        for row in range(1, nrow):
            for col in range(1, ncol):
                self.update(row, col)

    def update(self, row, col):
        """
        Method to update the matrices at a given row and column index.

        Input:
           row = the row index to update
           col = the column index to update
        """
        # Update all three score matrices concurrently
        # With boundaries set diagonally back, vertically back, and horizontally back should always have been populated
        # across all three matrices
        self.update_m(row, col)
        self.update_ix(row, col)
        self.update_iy(row, col)

    def update_m(self, row, col):
        
        # create list of score entries
        candidate_score_entries = []
        prev = (row - 1, col - 1)

        # get score for match from hash table
        # correct score lives at index - 1 because we are 1-indexed on the
        S_ij = self.align_params.match_matrix.get_score(self.align_params.seq_a[row-1], self.align_params.seq_b[col-1])
            
        # loop through the three score matrices, create candidate score to represent score if we came from this potential space, and append to the list of score 
        for matrix in [self.m_matrix, self.ix_matrix, self.iy_matrix]:
            prev_score_obj = matrix.get_score_obj(prev[0], prev[1])
            update_score = prev_score_obj.score + S_ij

            # zero out negative scores for local alignment
            if self.align_params.local_alignment:
                update_score = max(float(0), update_score)

            cand = ScoreEntry(prev[0], prev[1], update_score, prev_score_obj.matrix_name)
            candidate_score_entries.append(cand)

        # get the max score and pointers from list of score entries``
        max_score, max_pointers = get_maxes(candidate_score_entries)
        self.m_matrix.set_score(row, col, max_score)
        self.m_matrix.set_pointers(row, col, max_pointers)

    def update_ix(self, row, col):

        # initialize list of score entry objects
        candidate_score_entries = []
        
        # declare prev and make sure it is in bound for matrices
        prev = (row - 1, col)
        # loop through the two score matrices, adjust score to represent score if we came from this potential space, and append to the list of scores
        for matrix, penalty in zip([self.m_matrix, self.ix_matrix], [self.align_params.dy, self.align_params.ey]):
            prev_score_obj = matrix.get_score_obj(prev[0], prev[1])
            update_score = prev_score_obj.score - penalty

            # zero out negative scores for local alignment
            if self.align_params.local_alignment:
                update_score = max(float(0), update_score)

            cand = ScoreEntry(prev[0], prev[1], update_score, prev_score_obj.matrix_name)
            candidate_score_entries.append(cand)

        # get the max score and pointers
        max_score, max_pointers = get_maxes(candidate_score_entries)
        self.ix_matrix.set_score(row, col, max_score)
        self.ix_matrix.set_pointers(row, col, max_pointers)
         

    def update_iy(self, row, col):
        
        # initialize list of score entry objects
        candidate_score_entries = []

        # declare prev
        prev = (row, col - 1)

        # loop through the two score matrices, adjust score to represent score if we came from this potential space, and append to the list of scores
        for matrix, penalty in zip([self.m_matrix, self.iy_matrix], [self.align_params.dx, self.align_params.ex]):
            prev_score_obj = matrix.get_score_obj(prev[0], prev[1])
            update_score = prev_score_obj.score - penalty

            # zero out negative scores for local alignment
            if self.align_params.local_alignment:
                update_score = max(float(0), update_score)

            cand = ScoreEntry(prev[0], prev[1], update_score, prev_score_obj.matrix_name)
            candidate_score_entries.append(cand)
        
        # get the max score and pointers
        max_score, max_pointers = get_maxes(candidate_score_entries)
        self.iy_matrix.set_score(row, col, max_score)
        self.iy_matrix.set_pointers(row, col, max_pointers)

    def find_traceback_start(self):
        """
        Finds the location to start the traceback..
        Think carefully about how to set this up for local 

        Returns:
            (max_val, max_loc) where max_val is the best score
            max_loc is a set() containing tuples with the (i,j) location(s) to start the traceback
             (ex. [(1,2), (3,4)])
        """
        scores = set()
        for row in range(1, self.m_matrix.nrow):
            for col in range(1, self.m_matrix.ncol):
                score = self.m_matrix.get_score(row, col)
                scores.add(score)
        
        max_score = max(scores)
        max_loc = set()
        for row in range(1, self.m_matrix.nrow):
            for col in range(1, self.m_matrix.ncol):
                if fuzzy_equals(self.m_matrix.get_score(row, col), max_score):
                    max_loc.add((row, col))
        return max_score, max_loc
                
    def traceback(self): ### TO-DO! FILL IN additional arguments ###
        """
        Performs a traceback.
        Hint: include a way to printing the traceback path. This will be helpful for debugging!
           ex. M(5,4)->Iy(4,3)->M(4,2)->Ix(3,1)->Ix(2,1)->M(1,1)->M(0,0)
        """
        name_to_matrix = {"M": self.m_matrix, "Ix": self.ix_matrix, "Iy": self.iy_matrix}
        paths = []

        def traceback_depth_first_search(ptr, path):

            name, row, col = ptr.name, ptr.row, ptr.col
            matrix = name_to_matrix[name]
            score_obj = matrix.get_score_obj(row, col)
            # base case that we hit a boundary, add path to paths and return
            if row == 0 or col == 0:
                paths.append(path)
                return

            pointers = score_obj.pointers
            # recuse down path with each pointer
            for ptr in pointers:
                traceback_depth_first_search(ptr, path + [ptr])
        
        # for each max location, perform a depth first search
        max_loc = self.max_loc
        for loc in max_loc:
            score_entry_loc = self.m_matrix.get_score_obj(loc[0], loc[1])
            pointers = score_entry_loc.pointers
            # init with a pointer to max loc 
            path = [pointer(loc[0], loc[1], "M")]
            for ptr in pointers:
                traceback_depth_first_search(ptr, path + [ptr])

        return paths

    def print_paths(self):
        for i, path in enumerate(self.paths):
            print(f"Path {i+1}: ", end="")
            for pointer in path:
                print(pointer, end=" -> ")
            print("\n")


    def write_output(self):
        alignments = []
        # complete the emission process for each path 
        # each path represents a potential alignment
        for path in self.paths:
            alignment = {"a": "", "b": ""}
            
            # we want to emit from the start to the end because that is how sequences are represented
            # while in path, we have pointers that go from the end to the start
            for ptr in reversed(path): 
                row, col, name = ptr.row, ptr.col, ptr.name
                # skip boundary pointers (no emission)
                if row == 0 or col == 0:
                    continue

                # if diagonal, emit one letter from a, on letter from b
                if name == "M":
                    alignment["a"] += self.align_params.seq_a[row - 1]
                    alignment["b"] += self.align_params.seq_b[col - 1]

                # if horizontal, emit B[i] aligned with a gap 
                elif name == "Ix":
                    alignment["a"] += self.align_params.seq_a[row - 1]
                    alignment["b"] += "_"

                # if vertical, emit A[j] aligned with a gap 
                elif name == "Iy":
                    alignment["a"] += "_"
                    alignment["b"] += self.align_params.seq_b[col - 1]
            alignments.append(alignment)
        
        # write the alignments to the output file
        with open(self.output_file, "w") as f:
            f.write(str(self.max_score) + "\n\n")
            for alignment in alignments:
                f.write(alignment["a"] + "\n")
                f.write(alignment["b"] + "\n")
                f.write("\n")


def main():

    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # create an align object and run
    align = Align(input_file, output_file)
    align.align()


if __name__=="__main__":
    main()
