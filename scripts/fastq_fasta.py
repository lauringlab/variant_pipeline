from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description='This program takes a fastq file and converts it to a fasta file')
parser.add_argument('I', metavar='Input fastq', nargs='+',
                    help='The fastq file that is to be converted')
parser.add_argument('O', metavar='Output fasta', nargs='+',
                    help='The name of the output fasta file')

args = parser.parse_args()
count = SeqIO.convert(args.I[0], "fastq",args.O[0] , "fasta")
print("Converted %i records" % count)
