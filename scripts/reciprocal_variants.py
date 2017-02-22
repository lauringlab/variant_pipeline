import argparse
import pandas as pd
import pysam

parser = argparse.ArgumentParser(description='This script is designed to read in a csv of putative variant calls from deepSNV and query the bam file for the number of reads matching the reference base at that position.',usage ="python reciprocal_variants.py sample.bam sample.csv")

parser.add_argument('bam', metavar='bam', nargs='+',
                    help='The bam file of the sample')
parser.add_argument('csv', metavar='csv', nargs='+',
                    help='The deepSNV csv file')


class test_args: # for debugging and testing
    csv = ["~/Desktop/scratch/1139.removed.csv"]
    bam = ["/Users/jt/Desktop/scratch/1139.removed.bam"]
   
#args=parser.parse_args()    
args=test_args()
variants=pd.read_csv(args.csv[0],index_col=0)
bam= pysam.AlignmentFile(args.bam[0], "rb") 

def get_reference = function(row,bam) :
    chr=row["chr"]
    pos=int(row["pos"])
    py_pos=pos-1 # pythonify the position
    ref=row["ref"]
    