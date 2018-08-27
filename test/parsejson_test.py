import unittest
import copy
import json

from scripts.variantJSONtocsv import parseJson

class parsejsonfile(unittest.TestCase):
    def setUp(self):
        self.maxDiff=2000
        self.jsonString='{"Sample":"5_5","genome":[{ "chr":"PA", "seq":[ { "alleles":{ "A":{ "count":90, "freq":0.90, "mutationalClass":[ { "ORF":"PA", "aminoAcidPos":6, "classification":"Synonymous", "codingPos":20, "codonPos":2, "consensusAA":"Q", "varAA":"Q" } ], "nucleotide":"A" }, "G":{ "count":10, "freq":0.10, "mutationalClass":[ { "ORF":"PA", "aminoAcidPos":6, "classification":"Synonymous", "codingPos":20, "codonPos":2, "consensusAA":"Q", "varAA":"Q" } ], "nucleotide":"G" } }, "concat_pos":4728, "consensus":"A", "coverage":100, "pos":44 }]}]}'
        self.solution =[{
    "Sample": u"5_5",
    "chr": u"PA",
    "count": 90,
    "freq": 0.90,
    "concat_pos": 4728,
    "consensus": u"A",
    "coverage": 100,
    "pos": 44,
    "mutationalClass": [
                                {
                                    u"ORF": u"PA",
                                    u"aminoAcidPos": 6,
                                    u"classification": u"Synonymous",
                                    u"codingPos": 20,
                                    u"codonPos": 2,
                                    u"consensusAA": u"Q",
                                    u"varAA": u"Q"
                                }
                            ],
    "nucleotide": u"A"
    },{
    "Sample": u"5_5",
    "chr": u"PA",
    "count": 10,
    "freq": 0.10,
    "mutationalClass": [
                                {
                                    u"ORF": u"PA",
                                    u"aminoAcidPos": 6,
                                    u"classification": u"Synonymous",
                                    u"codingPos": 20,
                                    u"codonPos": 2,
                                    u"consensusAA": u"Q",
                                    u"varAA": u"Q"
                                }
                            ],
    "nucleotide": u"G",
    "concat_pos": 4728,
    "consensus": u"A",
    "coverage": 100,
    "pos": 44
    }]
        
    def tearDown(self):
        """
        This method is called after each test
        """
        pass
    def test_parse_correctly(self):
        """
        Test that we can add bases and count correctly
        """
        x = parseJson(json.loads(self.jsonString))
        self.assertListEqual(self.solution,parseJson(json.loads(self.jsonString)))


if __name__ == '__main__':
    unittest.main()