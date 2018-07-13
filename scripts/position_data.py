from __future__ import division
import pysam
import numpy as np
import yaml 
from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from scripts.seq_classes import locus, segment, tally, allele $ %%
#from seq_classes import locus, segment, tally, allele
import argparse
import os
import json

def main():
    parser = argparse.ArgumentParser(description='This scipts takes a bam file \
    and identifies variants and according to a consensus file.',
    usage ="python position_data.py sample.bed reference.fa sample.bam  sample.json -maxDepth 1000 -mqc")
    
    parser.add_argument('bed_json', metavar='bed_json', nargs='+',
                        help='a json bed like file with regions to compare')

    parser.add_argument('reference.fa', metavar='ref',nargs='+',
                            help = 'The sample consensus file which will be used to call nonsynonymous and synonymous mutations')

    parser.add_argument('bam', metavar='bam', nargs='+',
                        help='The bam file of the sample')
                        
    parser.add_argument('output', metavar='output', nargs='+',
                        help='The json file to hold the output')
                                            
    parser.add_argument('--maxDepth', metavar='maxDepth', type=int,
                        help='the max depth to use for pileup default is 1000')
    
    parser.add_argument('-mqc','--quality_metrics',action= 'store_true',dest = 'mqc',default = False)


    args = parser.parse_args()
    if args.maxDepth==None:
        maxDepth = 1000
    else:
        maxDepth=args.maxDepth

    # get bam file
    #bam = pysam.AlignmentFile(args.bam[0],"rb")
    bam = pysam.AlignmentFile("./test/accuracy/data/5_5.removed.bam") ###%%%%
    # set up reference dictions with key for each segment and value of [0,length]
    ref_genome_main={}
    # this is to maintain the order for concatenated pos
    chr_order  = []
    chr_length = []
    chr_starts = []

    # This needs to be changed to account for the new bed file format 
    # it should be from the  min of all start codons for each ORF to the max end 
    
    #with open(args.bed_json[0],"r") as f: 
    with open("./test/accuracy/data/reference/or.bed.json") as f: # %%
        regions=json.load(f)
        for segment in regions["genome"]:
            start = []
            stop = []
            chr_order.append(segment["seg"])
            chr_length.append(segment["size"])
            for orf in segment["ORF"]:
                for reading in orf["regions"]:
                    start.append(reading["start"])
                    stop.append(reading["stop"])
            ref_genome_main.update({segment["seg"]: [min(start),max(stop)]})
            
    chr_cumsum = [0] + list(np.cumsum(chr_length))

    # tally up base counts for each segement
    sample_genome={}
    for seg in ref_genome_main:
        sample_genome.update({seg: tally(bamfile=bam,chr=seg,\
        start = ref_genome_main[seg][0],stop = ref_genome_main[seg][1],maxDepth=maxDepth)})
    #makes sure the frequencies are up to date
    # probably don't need it now
    for seg in sample_genome:       
        sample_genome[seg].consensus()
    

# Here we will classify the variants 

    for seg in sample_genome:
        for ref_seg in regions["genome"]:
            if seg == ref_seg["seg"]:
                consensus_sequence = 
                for orf in ref_seg["ORF"]:
                    print orf



    i=0 # this is to make sure we handel the lase one correctly
    total = 0

# output the results here will take some decoding
    for seg in ref_genome_main:
        total = total + ref_genome_main[seg][1]-ref_genome_main[seg][0]
    with open(args.output[0],'w') as out:
        out.write("{\n\"loci\":[")
        for seg in sample_genome:
            for pos in sample_genome[seg].seq:
                # set concatpos
                pos.concat_pos = pos.pos + chr_cumsum[chr_order.index(pos.chr)]-chr_starts[chr_order.index(pos.chr)]
                if i <(total-1):
                    out.write(json.dumps(vars(pos),sort_keys=True, indent=4)+",")
                else:
                    out.write(json.dumps(vars(pos),sort_keys=True, indent=4)+"]}")
                i+=1



    if args.mqc: 
        # check if mqc dir exists if not make it
        if not os.path.exists("./mqc_position_stats"):
            os.makedirs("./mqc_position_stats")
        # get sample name
        basename = "./mqc_position_stats/"+os.path.splitext(os.path.basename(args.bam[0]))[0]
        concat_cov=[] 
        concat_pos = [] 
        i = 1 
        for loci in sample_genome[seg].seq:
            concat_cov.append(loci.coverage) 
            concat_pos.append(loci.concat_pos) 
            i+=1 
        with open(basename+"_mqc.cov.csv","w") as covfile: 
            i = 0 
            while i<len(concat_cov): 
                covfile.write("%d,%d\n" %(concat_pos[i],concat_cov[i])) 
                i+=1 
         
    # Frequencies
        concat_pos = [] 
        max_pos = 0 
        freqs = [] 
        for seg in sample_genome: 
            seg_freq=[] 
            pos = [] 
            for loci in sample_genome[seg].seq: 
                for k,v in loci.freqs.items(): 
                    if v >0 and k!=loci.consensus: 
                        freqs.append(v) 
                        seg_freq.append(v) 
                        concat_pos.append(loci.concat_pos) 


            seg_freq = np.sort(seg_freq) 
            cdf = np.array(range(len(seg_freq)))/float(len(seg_freq)) 
            with open(basename+ "-"+seg+"_mqc.cdf.csv","w") as cdfile: 
                i = 0 
                while i<len(seg_freq): 
                    cdfile.write("%f,%f\n" %(np.log10(seg_freq[i]),cdf[i])) 
                    i+=1 
                    
        with open(basename+"_mqc.frequencies.csv","w") as freqfile: 
            i = 0 
            while i<len(freqs): 
                freqfile.write("%d,%f\n" %(concat_pos[i],np.log10(freqs[i]))) 
                i+=1 
              
if __name__ == '__main__':
   main()
        
        
