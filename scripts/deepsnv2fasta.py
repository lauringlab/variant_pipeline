from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import argparse


parser = argparse.ArgumentParser(description='This script takes in a file that contains just one line of a DNA sequence and converts it to a fasta file')

parser.add_argument('file', metavar='file', nargs='+',
                    help='The input file')

parser.add_argument('out_fa', metavar='out_fa', nargs='+',
                    help='The output file')

parser.add_argument('name', metavar='name', nargs='+',
                    help='The identifier to be added to the fasta file. It will follow the \'>\' on the first line')


args = parser.parse_args()

with open(args.out_fa[0],'w') as out_file:
    out_file.write(">"+ args.name[0]+"\n")
    with open(args.file[0], 'r') as in_file:
        for line in in_file:
            out_file.write(line)
