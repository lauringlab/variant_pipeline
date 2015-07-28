#import sys
import os
#import gzip
import argparse
import shutil

parser = argparse.ArgumentParser(description='This program takes hiseq fastq files and renames them as sample.read_direction.#.fastq and keeps a log of the change, it may need to be adjusted if the sequencing core changes their scheme')
parser.add_argument('-s',action='store',dest='s',help='The sorce directory containing the original fastq files')
parser.add_argument('-f',action='store',dest='o',help='The final directory that will hold the renamed fastq files')
parser.add_argument('-k',action='store',dest='key',help='The output csv given by the sequencing core that serves as the key for renaming')
parser.add_argument('-run',action='store_true',dest='test',default=False,help='Boolean switch to run program, without this the program runs in test mode: the log is made but no files are renamed')

args=parser.parse_args()
s=args.s
o=args.f
key=arg.key
test=args.test

if not os.path.exists(o):
    os.makedirs(o)
# input argument is the sample sheet
junk_names = [] # The names given by the sequencing core sampleid_index
new_names = []
# add the bad names from the sample sheet to a list
f =  open(key,"r")
next(f)
for line in f:
    line = line.strip()
    line = line.split(',')
    junk_names.append(line[1]+"_"+line[3])
    new=line[4]#.replace("_","-")
    new_names.append(new)
f.close()

outfile = open(o+'renaming_log.txt','w')

for filename in os.listdir(s):
    name=filename.split("_L")
    bad_name = name[0]
    lane_junk = name[1]
    read_number=lane_junk.split("_R")
    fastq_number=read_number[1][4]
    read_number=read_number[1][0]
    if bad_name in junk_names:
        name_index = junk_names.index(bad_name)
        better_name= new_names[name_index]
        perfect_name= better_name+"."+read_number+"."+fastq_number+".fastq"
        # Write file to new name
        print("COPYING "+ s+filename + " to "+o+perfect_name)
        outfile.write(s+filename + "\t COPIED to \t" + o+perfect_name + "\n")
        if test==True:
            shutil.copy(s+filename,o+perfect_name)

outfile.close()
