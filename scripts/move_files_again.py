import argparse
import shutil

parser = argparse.ArgumentParser(description='So you\'ve run your analysis and deleted the copied fastq files to save space and now you need to run it again. This script takes in the renaming log file in the form initialfile\ttext\tfinal file and recopies the files.')
parser.add_argument("s",nargs='+',metavar='s',help='The renaming file.')
parser.add_argument("-r",action="store_true",dest='run', default=False,help = "Boolean switch to actually run the program.")

args=parser.parse_args()

with open(args.s[0],'r') as file:
    for line in file:
        line=line.split("\t")
        old=line[0].strip()
        new=line[2].strip()
        print new
        if args.run==True:
            shutil.copy2(old,new)
        else:
            print "test mode"
        print "copying %s to %s" % (old,new)
