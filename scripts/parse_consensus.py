#!/Users/jt/.virtualenvs/sci-py2.7/bin/python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import argparse
import copy


parser = argparse.ArgumentParser(description='This script takes in a deepSNV consensus sequence (concatenated string of all the segments) and the reference file used for alignment. It then deconcatenates the deepSNV file and outputs a vaild fasta file. I parses the string using the lengths of the segments in the reference file.')
# parser.add_argument('aligner_path', metavar='aligner_path', nargs='+',
#                     help='The path to the muscle executable - assuming the executable is name muscle')

parser.add_argument('ref_fa', metavar='ref', nargs='+',
                    help='The reference fasta to which the sequences will be trimmed.')
                    
parser.add_argument('in_fa', metavar='in_fa', nargs='+',
                    help='The input (sample) fa')

parser.add_argument('out_fa',metavar='out_fa',nargs="+",
                    help=' The trimmed fasta file')


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


def main():
    """Main body of script."""
    # parse arguments
    args = parser.parse_args()

    with open(args.in_fa[0],'r') as deepsnv:
        sample=deepsnv.readline().strip()
    
    ref=ReadFASTA(args.ref_fa[0])

    samp_fasta=copy.deepcopy(ref) 
    samp_seqname=[]
    ref_seqname=[]
    ref_length=[]
    for seq in ref:
        ref_seqname.append(seq.id)
        samp_seqname.append(seq.id)
        ref_length.append(len(seq))
    #print ref_length	
    j=0
    start=[]
    stop=[]
    for i in ref_length:
        start.append(sum(ref_length[0:j]))
        stop.append(start[j]+ref_length[j])
        j=j+1
    print start
    print stop
    i=0
    for seq in samp_fasta:
        seq.seq=Seq(sample[start[i]:stop[i]])
        i=i+1
    #print samp_fasta

    print "writing output to %s"  % args.out_fa
    SeqIO.write(samp_fasta, args.out_fa[0], "fasta")


main()
