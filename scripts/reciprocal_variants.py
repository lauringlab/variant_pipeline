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
variants=pd.read_csv(args.csv[0],index_col=False)
bam= pysam.AlignmentFile(args.bam[0], "rb") 

def get_reference(row,bam):
    #print row
    chr=row["chr"]
    pos=int(row["pos"])
    py_pos=pos-1 # pythonify the position
    ref=row["ref"]
    fw=0
    bw=0
    cov_fw=0
    cov_bw=0
    for pileupcolumn in bam.pileup(chr,py_pos,py_pos+1,truncate=True,stepper="all",max_depth=1E6): # look at the range py_pos to py_pos+1 and all the bases that align to those positions. Think of a column on each position like the samtools output of mpileup
        if pileupcolumn.pos==py_pos: #If the position is the one we want
            for pileupread in pileupcolumn.pileups: #Now for each read in that column
                if  not pileupread.is_del and not pileupread.is_refskip:
                    if pileupread.alignment.is_reverse:
                        called_base=pileupread.alignment.query_sequence[pileupread.query_position] # what is the called base at the position
                        called_phred=pileupread.alignment.query_qualities[pileupread.query_position] # what is the phred of that base
                        if called_phred>=30: # change this if you change the phred cut off in deepSNV. deepSNV only looks a phred>30. and we only want those calles that match the variant.
                            cov_bw=cov_bw+1
                            if called_base==ref:
                                bw=bw+1
                    if not pileupread.alignment.is_reverse:
                        called_base=pileupread.alignment.query_sequence[pileupread.query_position] # what is the called base at the position
                        called_phred=pileupread.alignment.query_qualities[pileupread.query_position] # what is the phred of that base
                        if called_phred>=30: # change this if you change the phred cut off in deepSNV. deepSNV only looks a phred>30. and we only want those calles that match the variant.
                            cov_fw=cov_fw+1
                            if called_base==ref:
                                fw=fw+1
    print(fw)

for index, row in variants.iterrows():
    get_reference(row,bam)                                                    