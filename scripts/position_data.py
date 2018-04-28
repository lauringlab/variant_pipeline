from __future__ import division
import pysam
import numpy as np
import yaml 
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from seq_classes import locus, segment, tally
import argparse
import os
import json

def main():
    parser = argparse.ArgumentParser(description='This scipts takes a bam file \
    and identifies sequence accuracy metric.',
    usage ="python consensus.py reference.fa sample.bam -maxDepth 1000")
    
    parser.add_argument('bed', metavar='bed', nargs='+',
                        help='a bed file with regions to compare')
    parser.add_argument('bam', metavar='bam', nargs='+',
                        help='The bam file of the sample')
    # parser.add_argument('basename', metavar='basename', nargs='+',
    #                     help='The basename of output files')
    parser.add_argument('output', metavar='output', nargs='+',
                        help='The json file to hold the segmens')
                                            
    parser.add_argument('--maxDepth', metavar='maxDepth', type=int,
                        help='the max depth to use for pileup default is 1000')
    parser.add_argument('-c','--coverage',action= 'store_true',dest = 'coverage',default = False)
    parser.add_argument('-f','--frequencies',action= 'store_true',dest = 'frequencies',default = False)
    parser.add_argument('-cdf','--cdf',action= 'store_true',dest = 'cdf',default = False)

    args = parser.parse_args()
    if args.maxDepth==None:
        maxDepth = 1000
    else:
        maxDepth=args.maxDepth

    # get bam file
    bam = pysam.AlignmentFile(args.bam[0],"rb")
    # set up reference dictions with key for each segment and value of [0,length]
    ref_genome={}
    # this is to maintain the order for concatenated pos
    chr_order  = []
    chr_length = []
    with open(args.bed[0],"r") as regions:
        for record in regions:
            record = record.split("\t")
            ref_genome.update({record[0]: [int(record[1]),int(record[2])]})
            chr_order.append(record[0])
            chr_length.append(int(record[2])-int(record[1])) # I have to subtract the start as well
    chr_cumsum = [0] + list(np.cumsum(chr_length))
    # tally up base counts for each segement
    sample_genome={}
    for seg in ref_genome:
        sample_genome.update({seg: tally(bamfile=bam,chr=seg,\
        start = ref_genome[seg][0],stop = ref_genome[seg][1],maxDepth=maxDepth) })
    #makes sure the frequencies are up to date
    for seg in sample_genome:       
        sample_genome[seg].consensus()
    
    # if args.coverage:
    #     concat_cov=[]
    #     concat_pos = []
    #     i = 1
    #     for seg in sample_genome:
    #         cov = sample_genome[seg].calc_coverage()
    #         for loci in cov:
    #             concat_cov.append(loci)
    #             concat_pos.append(i)
    #             i+=1
    #     with open(args.basename[0]+"_mqc.cov.csv","w") as covfile:
    #         i = 0
    #         while i<len(concat_cov):
    #             covfile.write("%d,%d\n" %(concat_pos[i],concat_cov[i]))
    #             i+=1
        
    # if args.frequencies:
    #     concat_pos = []
    #     max_pos = 0
    #     freqs = []
    #     for seg in sample_genome:
    #         seg_freq=[]
    #         pos = []
    #         for loci in sample_genome[seg].seq:
    #             for k,v in loci.freqs.items():
    #                 if v >0 and k!=loci.consensus():
    #                     freqs.append(v)
    #                     seg_freq.append(v)
    #                     pos.append(loci.pos)
    #         for p in pos:
    #             p = p +max_pos
    #         max_pos=max(pos)
    #         for p in pos:
    #             concat_pos.append(p)
    #         if args.cdf:
    #             seg_freq = np.sort(seg_freq)
    #             cdf = np.array(range(len(seg_freq)))/float(len(seg_freq))
    #             with open(args.basename[0]+ "-"+seg+"_mqc.cdf.csv","w") as cdfile:
    #                 i = 0
    #                 while i<len(seg_freq):
    #                     cdfile.write("%f,%f\n" %(np.log10(seg_freq[i]),cdf[i]))
    #                     i+=1
                   
    #     with open(args.basename[0]+"_mqc.frequencies.csv","w") as freqfile:
    #         i = 0
    #         while i<len(freqs):
    #             freqfile.write("%d,%f\n" %(concat_pos[i],np.log10(freqs[i])))
    #             i+=1
    i=0 # this is to make sure we handel the lase one correctly
    total = 0
    for seg in ref_genome:
        total = total + ref_genome[seg][1]-ref_genome[seg][0]
    with open(args.output[0],'w') as out:
        out.write("[")
        for seg in sample_genome:
            for pos in sample_genome[seg].seq:
                # set concatpos
                pos.concat_pos = pos.pos + chr_cumsum[chr_order.index(pos.chr)]
                if i <(total-1):
                    out.write(json.dumps(vars(pos),sort_keys=True, indent=4)+",")
                else:
                    out.write(json.dumps(vars(pos),sort_keys=True, indent=4)+"]")
                i+=1



        
              
if __name__ == '__main__':
   main()
        
        