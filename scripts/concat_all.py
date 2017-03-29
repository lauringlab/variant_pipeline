from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import glob
import pandas as pd
import os

parser = argparse.ArgumentParser(description='This script is designed to concatenate all the segements in the input fa file')


parser.add_argument('in_fa', metavar='in)fa', nargs='+',
                    help='The input file')

parser.add_argument('out_fa', metavar='out_fa', nargs='+',
                    help='The output file')


args = parser.parse_args()


def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    'fastafile' should specify the name of a FASTA file.

    This function reads all sequences from the FASTA file.  It returns the
        list 'headers_seqs'.  This list is composed of a seq_record objects.
    """
    seqs =[]
    header = None
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seqs.append(seq_record)

    return seqs

# Read in the fasta file a list of each sequence record
fasta=ReadFASTA(args.in_fa[0])

#cycle through the list and add the sequence to a concat string
concat_seq=""
for i in fasta:
    concat_seq=concat_seq+str(i.seq)

#make the concat string a seq object and then a seqRecord
concat_seq=Seq(concat_seq)
concat_record=SeqRecord(concat_seq)

identifier = os.path.basename(args.in_fa[0])
identifier=identifier.split(".")[0]

concat_record.id=identifier
concat_record.description="All seqeunces from %s concatenated" %(args.in_fa[0])




SeqIO.write(concat_record, args.out_fa[0], "fasta")
