import unittest
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from scripts.seq_classes import allele, checkORF

class test_checkORF(unittest.TestCase):
        
    def tearDown(self):
        """
        This method is called after each test
        """
        pass

    def test_good(self):
        sequence = Seq("ATGATGTAA")
        self.assertTrue(checkORF(sequence))

    def test_noStart(self):
        sequence = Seq("AGGTAA")
        with self.assertRaises(ValueError):
            checkORF(sequence)

    def test_missplacedStart(self):
        sequence = Seq("AGGATGTAA")
        with self.assertRaises(ValueError):
            checkORF(sequence)

    def test_notfullcodon(self):
        sequence = Seq("ATGAATAA")
        with self.assertRaises(ValueError):
            checkORF(sequence)

    def test_noStop(self):
        sequence = Seq("ATGCAA")
        with self.assertRaises(ValueError):
            checkORF(sequence)

    def test_internalStop(self):
        sequence = Seq("ATGTAACAA")
        with self.assertRaises(ValueError):
            checkORF(sequence)
"""             
class test_classifier(unnitTest.TestCase):
    def setUp(self):
        self.sequence = SeqRecord("---ATGGCGT--GCGATGA-----ATAA---",id = "testGene", description = "")  
        self.primary = "MACDE*"
        self.codingRegion = {
					"name": "myFavoriteProtien",
					"regions": [
						{
							"start": 3,
							"stop": 10
						},
						{
							"start": 12,
							"stop": 19
						},{
                        "start":24,
                        "stop":28
                        }
					]
				}
        self """