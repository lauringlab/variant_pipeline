import argparse
import shutil

parser = argparse.ArgumentParser(description='So you\'ve run your analysis and deleted the copied fastq files to save space and now you need to run it again. This script takes in the renaming log file in the form initialfile\ttext\tfinal file and recopies the files.')
parser.add_argument("source",nargs='+',metavar='s',help='The renaming file.')

agrs=args=parser.parse_args()

with open(args.s[0],'r') as file:
    for line in file:
        line=line.split("\t")
        old=line[0]
        new=line[1]
        shutil.copy(old,new)
        print "copying %s to %s" (old,new)
