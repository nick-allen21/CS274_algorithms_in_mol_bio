"""
Unit tests provide a way of testing individual components of your program.

This will help you with debugging and making sure you don't break your code!


Here, we provide some unit tests and a skeleton for a few others.
Note that you do not have to follow these exactly, they are designed to help you.


What other unit tests might you want to write?
 - Think about the traceback and writing the output file. 
 - Try to write at least one or two additional tests as you think of them.


To run:
  python align_test.py

Make sure align.py is located in the same directory, and the test_example.input file is present!
"""

import unittest

from align_skeleton import *
TEST_INPUT_FILE="/Users/nickallen/Documents/GitHub/-CS274-Algorithms-in-Molecular-Biology/Project1/code/test_example.input"

class TestAlignmentClasses(unittest.TestCase):

    def test_match_matrix(self):
        """
        Tests match matrix object
        """
        match_matrix = MatchMatrix()
        match_matrix.set_score("A", "C", 5)
        self.assertEqual(match_matrix.get_score("A", "C"), 5)

    def test_score_matrix_score(self):
        """
        Tests score matrix object score set + get methods
        """
        score_matrix = ScoreMatrix("S", 5, 4)
        score_matrix.set_score(2, 2, 5)
        # Act
        got_entry = score_matrix.get_score(2, 2)
        print("Got score:", got_entry.score, "Expected:", 5)
        # Assert
        self.assertAlmostEqual(got_entry.score, 5)

    def test_score_matrix_pointers(self):
        """
        Tests score matrix object pointer set + get methods
        """
        score_matrix = ScoreMatrix("S", 5, 4)
        score_matrix.set_score(2, 2, 5)
        expected = [pointer(2, 2, "S")]
        score_matrix.set_pointers(2, 2, expected)

        got = score_matrix.get_pointers(2, 2)
        print("Got pointers:", got, "Expected:", expected)

        # With __eq__ on pointer, this should pass
        self.assertEqual(got, expected)

    def test_param_loading(self):
        """
        Tests AlignmentParameters "load_params_from_file()" function
        """
        align_params = AlignmentParameters()
        align_params.load_params_from_file(TEST_INPUT_FILE)
        self.assertEqual(align_params.seq_a, "AATGC")
        self.assertEqual(align_params.seq_b, "AGGC")
        self.assertTrue(align_params.global_alignment)
        self.assertEqual(align_params.dx, 0.1)
        self.assertEqual(align_params.ex, 0.5)
        self.assertEqual(align_params.dy, 0.6)
        self.assertEqual(align_params.ey, 0.3)
        self.assertEqual(align_params.alphabet_a, "ATGC")
        self.assertEqual(align_params.alphabet_b, "ATGCX")
        self.assertEqual(align_params.len_alphabet_a, 4)
        self.assertEqual(align_params.len_alphabet_b, 5)

        # test that the match match is set up correctly
        #  if this fails, make sure you are loading the asymmetric matrix properly!
        match_mat = align_params.match_matrix
        self.assertEqual(match_mat.get_score("A", "X"), 0.3)
        self.assertEqual(match_mat.get_score("C", "G"), -0.3)
        self.assertEqual(match_mat.get_score("G", "C"), 0)


    def test_update_ix(self):
        """
        Test AlignmentAlgorithm's update Ix
        """

        # configure alignment params
        align_params = AlignmentParameters()
        align_params.dy = 1
        align_params.ey = 0.5
        # added these to make sure the update_ix method works
        align_params.dx = 1
        align_params.ex = 0.5

        # ive added this to make sure the update_ix method works
        align_params.alphabet_a = "ATGC"
        align_params.alphabet_b = "ATGCX"
        align_params.len_alphabet_a = 4
        align_params.len_alphabet_b = 5

        # create an alignment object
        align = Align("", "")
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.m_matrix.set_score(2,2, 3)
        align.ix_matrix.set_score(2,2, 2.5)

        # run the method!
        align.update_ix(3, 2)

        got_entry = align.ix_matrix.get_score(3, 2)
        got_pointers = align.ix_matrix.get_pointers(3, 2)

        # From M: 3 - dx = 2.0, from Ix: 2.5 - ex = 2.0, tie -> we want to include both pointers
        self.assertTrue(fuzzy_equals(got_entry.score, 2.0))
        self.assertEqual(set((p.name, p.row, p.col) for p in got_pointers),
                        {("M", 2, 2), ("Ix", 2, 2)})

    def test_update_m(self):
        align_params = AlignmentParameters()
        align_params.load_params_from_file(TEST_INPUT_FILE)

        align = Align("", "")
        align.align_params = align_params
        
        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.ix_matrix = ScoreMatrix("Ix", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2, 2, 2.5)
        align.ix_matrix.set_score(2, 2, 3)
        align.iy_matrix.set_score(2, 2, 3)

        align_params.match_matrix.set_score(align_params.seq_a[3], align_params.seq_b[3], 5)

        align.update_m(3,3)

        got_entry = align.m_matrix.get_score(3, 3)
        got_pointers = align.m_matrix.get_pointers(3, 3)

        print("Got score:", got_entry.score, "Expected:", 8)
        print("Got pointers:", got_pointers, "Expected:", [pointer(2, 2, "Ix"), pointer(2, 2, "Iy")])

        # From M: 2.5 + 5 = 7.5, from Ix: 3 + 5 = 8, from Iy: 3 + 3 = 8, will get points to Ix and Iy
        self.assertAlmostEqual(got_entry.score, 8)
        self.assertEqual(set((p.name, p.row, p.col) for p in got_pointers),
                        {("Ix", 2, 2), ("Iy", 2, 2)})
    
    def test_update_iy(self):
        """
        Test AlignmentAlgorithm's update Iy
        """
        align_params = AlignmentParameters()
        align_params.load_params_from_file(TEST_INPUT_FILE)

        align = Align("", "")
        # dy = .6, # ey = .3
        align.align_params = align_params

        align.m_matrix = ScoreMatrix("M", 5, 4)
        align.iy_matrix = ScoreMatrix("Iy", 5, 4)
        align.m_matrix.set_score(2,2, 3.6)
        align.iy_matrix.set_score(2,2, 3.2)

        print("Starting update_iy")
        align.update_iy(2,3)

        got_entry = align.iy_matrix.get_score(2, 3)
        got_pointers = align.iy_matrix.get_pointers(2, 3)
    
        
        print("Got score:", got_entry.score, "Expected:", 3.0)
        print("Got pointers:", got_pointers, "Expected:", [pointer(2, 2, "M")])

        # From M: 3.6 - .6 = 3.0, from Iy: 3.2 - .3 = 2.0, should get max score 3 and points to M
        self.assertAlmostEqual(got_entry.score, 3.0)
        self.assertEqual(set((p.name, p.row, p.col) for p in got_pointers),
                        {("M", 2, 2)})

    def test_traceback_start(self):
        """
        Tests that the traceback finds the correct start
        Should test local and global alignment!
        """
        ### FILL IN ###
        return


if __name__=='__main__':
    unittest.main(verbosity=3)