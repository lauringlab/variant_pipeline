import unittest

from scripts.trim_to_coding import get_regions, trim

class test_loci(unittest.TestCase):
    def setUp(self):
        pass

    def test_trimming_ends(self):
        solution = ["AAAAA",[[2,7]]]
        seqs = ["--AAAAA--","TTAAAAATT"]#TTAATTAAA-AGC"
        output = trim(seqs)
        self.assertEqual(output,solution)
    
    def test_trimming_middles(self):
        solution = ["AAAAAAA",[[2,5],[6,8],[9,11]]]
        seqs = ["--AAA-AA-AA-","TTAAATAATAAT"]#TTAATTAAA-AGC"
        output = trim(seqs)
        self.assertEqual(output,solution)

    def test_trimming_middles_noend(self):
        solution = ["AAAAAAA",[[2,5],[6,8],[9,11]]]
        seqs = ["--AAA-AA-AA","TTAAATAATAA"]#TTAATTAAA-AGC"
        output = trim(seqs)
        self.assertEqual(output,solution)

    def test_no_gap(self):
        solution = ["AAAA",[[0,4]]]
        seqs = ["AAAA","AAAA"]#TTAATTAAA-AGC"
        output = trim(seqs)
        self.assertEqual(output,solution) 

    def test_no_gap_at_end(self):
        solution = ["AAAA",[[0,3],[4,5]]]
        seqs = ["AAA-A","AAATA"]#TTAATTAAA-AGC"
        output = trim(seqs)
        self.assertEqual(output,solution)


if __name__ == '__main__':
    unittest.main()