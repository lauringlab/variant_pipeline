import unittest
import pysam

from scripts.seq_classes import allele

class test_classifier(unittest.TestCase):
    def setUp(self):
        self.codingRegion = {
					"name": "PB2",
					"regions": [
						{
							"start": 27,
							"stop": 2307
						}
					]
				}
    
    
    def tearDown(self):
        """
        This method is called after each test
        """
        pass
    def test_classifier(self):
        """
        Test that we can add bases and count correctly
        """
        self.assertEqual()