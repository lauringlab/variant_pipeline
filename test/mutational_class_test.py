import unittest
import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from scripts.seq_classes import allele, checkORF,classify

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
          
class test_classifier(unittest.TestCase):
    def setUp(self):
        self.sequence = SeqRecord("---ATGGCGT--GCGATGA-----ATAA---",id = "testGene", description = "")  
        self.primary = "MACDE*"
        self.codingRegion = [{
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
				},{
                    "name": "myFavoriteProtien_alternate_splice",
					"regions": [
						{
							"start": 3,
							"stop": 7
						},
						{
                        "start":15,
                        "stop":16
                        },
                        {
                            "start":24,
                            "stop":28
                        }
					]
				}]
        # The sequences are ATG GCG TGC GAT GAA TAA = MACDE*

        # and               ATG GAA TAA  = MA*
        
    def test_outside_5prime(self):
        output = {
            "ORF": "myFavoriteProtien",
            "codonPos": None,
            "codingPos": None,
            "aminoAcidPos": None,
            "consensusAA": None,
            "varAA": None,
            "classification": "Noncoding"

        }
        a = allele("T")
        a.classifyVar(self.sequence,self.codingRegion[0],20)
        self.assertEqual(output,a.mutationalClass[0])

    def test_outside_splice(self):
        output = {
            "ORF": "myFavoriteProtien",
            "codonPos": None,
            "codingPos": None,
            "aminoAcidPos": None,
            "consensusAA": None,
            "varAA": None,
            "classification": "Noncoding"

        }
        a = allele("T")
        a.classifyVar(self.sequence,self.codingRegion[0],1)
        self.assertEqual(output,a.mutationalClass[0])

    def test_outside_3prime(self):
        output = {
            "ORF": "myFavoriteProtien",
            "codonPos": None,
            "codingPos": None,
            "aminoAcidPos": None,
            "consensusAA": None,
            "varAA": None,
            "classification": "Noncoding"

        }
        a = allele("T")
        a.classifyVar(self.sequence,self.codingRegion[0],30)
        self.assertEqual(output,a.mutationalClass[0])

    def test_indel(self):
        output = {
            "ORF": "myFavoriteProtien",
            "codonPos": 1,
            "codingPos": 1,
            "aminoAcidPos": 0,
            "consensusAA": "M",
            "varAA": None,
            "classification": "Indel"

        }
        a = allele("-")
        a.classifyVar(self.sequence,self.codingRegion[0],4)
        self.assertEqual(output,a.mutationalClass[0])

    def test_NS(self):
        output = {
            "ORF": "myFavoriteProtien",
            "codonPos": 0,
            "codingPos": 3,
            "aminoAcidPos": 1,
            "consensusAA": "A",
            "varAA": "P" ,
            "classification": "Nonsynonymous"

        }
        a = allele("C")
        a.classifyVar(self.sequence,self.codingRegion[0],6)
        self.assertEqual(output,a.mutationalClass[0])           

    def test_S(self):
        output = {
            "ORF": "myFavoriteProtien",
            "codonPos": 2,
            "codingPos": 14,
            "aminoAcidPos": 4,
            "consensusAA": "E",
            "varAA": "E" ,
            "classification": "Synonymous"

        }
        a = allele("G")
        a.classifyVar(self.sequence,self.codingRegion[0],24)
        self.assertEqual(output,a.mutationalClass[0])   

    def test_readthrough(self):
        output = {
            "ORF": "myFavoriteProtien",
            "codonPos": 0,
            "codingPos": 15,
            "aminoAcidPos": 5,
            "consensusAA": "*",
            "varAA": "E" ,
            "classification": "Readthrough"

        }
        a = allele("G")
        a.classifyVar(self.sequence,self.codingRegion[0],25)
        self.assertEqual(output,a.mutationalClass[0])   

    def test_stop(self):
        output = {
            "ORF": "myFavoriteProtien_alternate_splice",
            "codonPos": 0,
            "codingPos": 3,
            "aminoAcidPos": 1,
            "consensusAA": "E",
            "varAA": "*" ,
            "classification": "Stop"

        }
        a = allele("T")
        a.classifyVar(self.sequence,self.codingRegion[1],6)
        self.assertEqual(output,a.mutationalClass[0])   
