from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import argparse
import pandas as pd
from ast import literal_eval
import copy


parser = argparse.ArgumentParser(description='This script takes in a fasta file and a csv containing the open reading frames for each sequence in the fasta file. The script trims the sequences and outputs a new fasta file containing just the open reading frames.')

parser.add_argument('file', metavar='file', nargs='+',
                    help='The input fasta file')

parser.add_argument('out_fa', metavar='out_fa', nargs='+',
                    help='The output fasta file')

parser.add_argument('csv', metavar='csv', nargs='+',
                    help='The csv with each the open reading frames as a list of regions [[start,stop][splice_start,splice_stop]]. The columns should be labled "chr" and "coding"')

args = parser.parse_args()

#class test_args:
#    csv = ["../../Desktop/flu_OR_work/wsn33.regions.csv"]
#    file = ["../../Desktop/flu_OR_work/wsn33.fa"]
#    out_fa = ["/Users/jt/Desktop/flu_OR_work/wsn33.OR.fa"]
#args=test_args()

def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    'fastafile' should specify the name of a FASTA file.

    This function reads all sequences from the FASTA file.  It returns the
        list 'headers_seqs'.  This list is composed of a seq_record objects.
    """
    seqs =[]
    header = None
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seq_record.seq.alphabet=IUPAC.unambiguous_dna
        seqs.append(seq_record)

    return seqs


fasta=ReadFASTA(args.file[0])
regions = pd.read_csv(args.csv[0],converters={"coding": literal_eval},index_col="chr") # the chromosomes are the indexes here - that's fine they should be unique

# cycle through fasta filte
OR=[]
for gene in fasta:
    seg=gene.id
    coding_regions = regions.loc[seg,"coding"]
    open_reading=copy.deepcopy(gene)
    open_reading.seq=Seq('',alphabet=IUPAC.IUPACUnambiguousDNA())
    for splice in coding_regions:
        open_reading.seq=open_reading.seq+gene.seq[splice[0]-1:splice[1]] # for python numbering and slicing
    
    OR.append(open_reading)  
    
output_handle = open(args.out_fa[0], "w")
SeqIO.write(OR, output_handle, "fasta")
output_handle.close()
