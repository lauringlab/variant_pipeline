import unittest
import pysam
import copy

from scripts.seq_classes import locus, segment, allele

class test_loci(unittest.TestCase):
    def setUp(self):
        
        self.solution ={'A':allele("A"),
                        'T':allele("T"),
                        'C':allele("C"),
                        'G':allele("G"),
                        '-':allele("-")}
        self.solution["A"].count = 11
        self.solution["A"].freq = 11/18.0
        self.solution["T"].count = 4
        self.solution["T"].freq = 4/18.0
        self.solution["C"].count = 1
        self.solution["C"].freq = 1/18.0
        self.solution["G"].count = 1
        self.solution["G"].freq = 1/18.0
        self.solution["-"].count = 1
        self.solution["-"].freq = 1/18.0
        self.pileup = "AAAAATTAATTAAA-AGC"
        self.l = locus(pos =0)
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
        bases = ['A','T','C','G','-']
        p = True
        for b in bases:
            if self.l.alleles[b].freq==self.solution[b].freq:
                p==False
            if self.l.alleles[b].count==self.solution[b].count:
                p==False  
        self.assertTrue(p)
        
    def test_consensus(self):
        """
        Test that consensus works with precent cutoff
        """
        self.assertEqual(self.l.calc_consensus(0.5),"A")
    
    def test_consensus_N(self):
        """
        Test that consensus works when precent cutoff is not met
        """
        self.assertEqual(self.l.calc_consensus(1.0),"N")
    def test_consensus_common(self):
        """
        Test that consensus works when no cutoff is given
        """
        self.assertEqual(self.l.consensus,"A")

class test_segment(unittest.TestCase):
    
    def setUp(self):
        
        self.PB1 = segment("PB1")
        self.pileup = "AAAAATTAATTAAA-AGC"
        self.l = locus(pos =0)
        for b in self.pileup:
            self.l.update(b)

    def tearDown(self):
        """
        This method is called after each test
        """
        pass
        
    def test_append_not_loci(self):
        self.assertRaises(ValueError,lambda : self.PB1.append_loci("A"))
    
    # def test_append_wrong_seg(self):
    #     self.assertRaises(ValueError,lambda : self.PB1.append_loci(locus(chr="HA",pos=0)))
    
    def test_append_wrong_pos(self):
        self.PB1.append_loci(self.l)
        self.assertRaises(ValueError,lambda : self.PB1.append_loci(locus(pos=4)))
   
    def test_append(self):
        self.PB1.append_loci(self.l)
        self.assertEqual(self.PB1.seq[0],self.l)
    
    def test_loci_return(self):
        self.PB1.append_loci(self.l)
        self.assertEqual(self.PB1.locus(0),self.l)
        
    def test_consensus(self):
        self.PB1.append_loci(self.l)
        self.l2 = copy.deepcopy(self.l)
        self.l2.pos = 1
        self.PB1.append_loci(self.l2)
        self.assertEqual(self.PB1.consensus(0.5),"AA")
        
    def test_coverage(self):
        self.PB1.append_loci(self.l)
        self.l2 = copy.deepcopy(self.l)
        self.l2.pos = 1
        self.PB1.append_loci(self.l2)
        self.assertEqual(self.PB1.calc_coverage(),[len(self.pileup),len(self.pileup)])
    

if __name__ == '__main__':
    unittest.main()