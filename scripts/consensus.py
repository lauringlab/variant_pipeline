from __future__ import division
import pysam
import numpy as np
import yaml 
from variant_pipeline.scripts.seq_classes import locus, segment, tally



#parser = argparse.ArgumentParser(description='This scipts takes a bam file \
#and identifies the consensus file sequence of the sample.',
#usage ="python consensus.py sample.bam sample.fa")

#parser.add_argument('bam', metavar='bam', nargs='+',
#                    help='The bam file of the sample')
#parser.add_argument('sample.fa', metavar='fasta', nargs='+',
#                    help='the consensus sequence output file')

# for each chr in bam make a segement
    # for each position in bam make a loci and add it to the sement
    # calcuate the concensus
    # calculate the coverage

def tally(bamfile,chr,start,stop,max_depth = 1000,phred = 30):
    """
        This function takes in a pysam alignment file and tallies base calls 
        across a defined region. Bases must pass a phred score cutoff.
    """
    seg = segment(chr)
    for pos in range(start,stop):
        pileup = bamfile.pileup(chr,pos,pos+1,
        stepper='all',truncate=True,max_depth = max_depth)
        l = locus(chr,pos=pos+1) # loci positions are base 1
        for pileupColumn in pileup:
            for pileupRead in pileupColumn.pileups:
            # what is the called base at the position
                if pileupRead.query_position==pos:
                    called_base= \
                    pileupRead.alignment.query_sequence[pileupRead.query_position] 
                    called_phred= \
                    pileupRead.alignment.query_qualities[pileupRead.query_position]
                    if called_phred>=phred: 
                        l.update(called_base)
        
        seg.append_loci(l)
    return seg


bam = pysam.AlignmentFile("./chasing_excellence/test/Plasmid_control.removed.bam","rb")







