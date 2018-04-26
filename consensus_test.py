import unittest
import pysam
import copy

from scripts.seq_classes import locus, segment

class test_loci(unittest.TestCase):
    def setUp(self):
        
        self.solution = {"A":11,"T":4,"G":1,"C":1,"-":1}
        self.pileup = "AAAAATTAATTAAA-AGC"
        self.l = locus(chr="HA",pos =1)
        for b in self.pileup:
            self.l.update(b)
    
    def tearDown(self):
        """
        This method is called after each test
        """
        pass
    def test_adding_bases(self):
        """
        Test that we can add bases and count correctly
        """
        self.assertEqual(self.l.counts,self.solution)
        
    def test_consensus(self):
        """
        Test that consensus works with precent cutoff
        """
        self.assertEqual(self.l.consensus(0.5),"A")
    
    def test_consensus_N(self):
        """
        Test that consensus works when precent cutoff is not met
        """
        self.assertEqual(self.l.consensus(1.0),"N")
    def test_consensus_common(self):
        """
        Test that consensus works when no cutoff is given
        """
        self.assertEqual(self.l.consensus(),"A")

class test_segment(unittest.TestCase):
    
    def setUp(self):
        
        self.PB1 = segment("PB1")
        self.pileup = "AAAAATTAATTAAA-AGC"
        self.l = locus(chr="PB1",pos =1)
        for b in self.pileup:
            self.l.update(b)

    def tearDown(self):
        """
        This method is called after each test
        """
        pass
        
    def test_append_not_loci(self):
        self.assertRaises(ValueError,lambda : self.PB1.append_loci("A"))
    
    def test_append_wrong_seg(self):
        self.assertRaises(ValueError,lambda : self.PB1.append_loci(locus(chr="HA",pos=1)))
    
    def test_append_wrong_pos(self):
        self.assertRaises(ValueError,lambda : self.PB1.append_loci(locus(chr="HA",pos=2)))
   
        
    def test_append(self):
        self.PB1.append_loci(self.l)
        self.assertEqual(self.PB1.seq[0],self.l)
    
    def test_loci_return(self):
        self.PB1.append_loci(self.l)
        self.assertEqual(self.PB1.locus(1),self.l)
        
    def test_consensus(self):
        self.PB1.append_loci(self.l)
        self.l.pos=2
        self.PB1.append_loci(self.l)
        self.assertEqual(self.PB1.consensus(0.5),"AA")
    def test_coverage(self):
        self.PB1.append_loci(self.l)
        self.l.pos=2
        self.PB1.append_loci(self.l)
        self.assertEqual(self.PB1.calc_coverage(),[len(self.pileup),len(self.pileup)])
    

if __name__ == '__main__':
    unittest.main()