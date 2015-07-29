#import sys
import os
import argparse
import shutil

parser = argparse.ArgumentParser(description='This program takes Miseq fastq files and renames them as sample.read_direction.#.fastq and keeps a log of the change')
parser.add_argument('-s',action='store',dest='s',help='The sorce directory containing the original fastq files')
parser.add_argument('-f',action='store',dest='f',help='The final directory that will hold the renamed fastq files')
parser.add_argument('-run',action='store_true',dest='test',default=False,help='Boolean switch to run program, without this the program runs in test mode: the log is made but no files are renamed')
parser.add_argument('-gz',action='store_true',dest='zip',default=False,help='Boolean switch. Activate when working with .gz files')

args=parser.parse_args()
s=args.s
f=args.f
test=args.test
zip=args.zip
# input argument is fastq directory

if not os.path.exists(f):
    os.makedirs(f)

outfile = open(f+'renaming_log.txt','w')
#os.chdir(sys.argv[1])
if test==False:
    print "running in test mode add option -r to run"
for filename in os.listdir(s):
    name=filename.split("_L")
    good_name = name[0].split("_")[0]
    lane_junk = name[1]
    read_number=lane_junk.split("_R")
    fastq_number=read_number[1][4]
    read_number=read_number[1][0]
    if zip==True:
        perfect_name= good_name.replace("-","_")+"."+read_number+"."+fastq_number+".fastq.gz"
    else:
        perfect_name= good_name.replace("-","_")+"."+read_number+"."+fastq_number+".fastq.gz"
    # Write file to new name

    print("COPYING "+ s+filename + " to "+f+perfect_name)
    outfile.write(s+filename + "\t COPIED to \t" + f+perfect_name + "\n")
    if test==True:
        shutil.copy(s+filename,f+perfect_name)


outfile.close()
