from __future__ import division
import pysam
import numpy as np
import yaml 
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from seq_classes import locus, segment, tally
import argparse


def main():
    parser = argparse.ArgumentParser(description='This scipts takes a bam file \
    and identifies the consensus file sequence of the sample.',
    usage ="python consensus.py reference.fa sample.bam sample.fa -maxDepth 1000")
    
    parser.add_argument('bed', metavar='bed', nargs='+',
                        help='a bed file with regions to compare. base 0 left inclusive right exclusive')
    parser.add_argument('bam', metavar='bam', nargs='+',
                        help='The bam file of the sample')
    parser.add_argument('sample_fa', metavar='sample_fa', nargs='+',
                        help='the consensus sequence output file')
    parser.add_argument('--maxDepth', metavar='maxDepth', type=int,
                        help='the max depth to use for pileup default is 1000')
    parser.add_argument('--cutoff', metavar='cutoff', type=float,
                        help='The frequency threshold for calling a consensus base. If it is not used then the most common base is selected')
    # for each chr in bam make a segement
        # for each position in bam make a loci and add it to the sement
        # calcuate the concensus
        # calculate the coverage
    args = parser.parse_args()
    if args.maxDepth==None:
        maxDepth = 1000
    else:
        maxDepth=args.maxDepth
    
    # get bam file
    bam = pysam.AlignmentFile(args.bam[0],"rb")
    # set up reference dictions with key for each segment and value of [0,length]
    
    ref_genome={}
    with open(args.bed[0],"r") as regions:
        for record in regions:
            record = record.split("\t")
            ref_genome.update({record[0]: [int(record[1]),int(record[2])]})
    
    # tally up base counts for each segement
    sample_genome={}
    for seg in ref_genome:
        sample_genome.update({seg: tally(bamfile=bam,chr=seg,\
        start = ref_genome[seg][0],stop = ref_genome[seg][1],maxDepth=maxDepth) })
    
    # get consensus for each segment and save as seq in list
    consensus=[]
    for seg in sample_genome:
        seq = SeqIO.SeqRecord(id =seg,
            seq = Seq(sample_genome[seg].consensus(args.cutoff),generic_dna),
            description = "Made using arguments: " + str(args))
        consensus.append(seq)
    
    # write to file
    
    SeqIO.write(consensus, args.sample_fa[0], "fasta")

if __name__ == '__main__':
   main()



