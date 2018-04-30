import unittest

from scripts.seq_classes import locus

class test_loci(unittest.TestCase):
    def setUp(self):
        pass
    def test_adding_bases(self):
        solution = {"A":11,"T":4,"G":1,"C":1,"-":1}
        pileup = "AAAAATTAATTAAA-AGC"
        l = locus(chr="HA",pos =1)
        for b in pileup:
            l.update(b)
        self.asserEqual(l.counts,solution)

if __name__ == '__main__':
    unittest.main()