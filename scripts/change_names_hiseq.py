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
f=args.f
key=arg.key
test=args.test

if not os.path.exists(f):
    os.makedirs(f)

# add slash to final location if it was forgotten
if f[-1] != '/':
    f=f+'/'
#print f

if s[-1] != '/':
    s=s+'/'
#print s



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

if test==False:
    print "running in test mode add option -run to run"
    
outfile = open(f+'renaming_log.txt','w')

for filename in glob.glob(s + "*.fastq"):
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
            shutil.copy(s+filename,f+perfect_name)
            
for zipfilename in glob.glob(s + "*.fastq.gz"):
    zipfilename= os.path.basename(zipfilename)
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
            shutil.copy(s+filename,f+perfect_name)
            
for zipfile in glob.glob(f + "*.gz"):
    print "unzipping:" + zipfile

    inF = gzip.GzipFile(zipfile, 'rb')
    s = inF.read()
    inF.close()
    #print(os.path.splitext(zipfile)[0])
    outF = file(os.path.splitext(zipfile)[0], 'wb')
    outF.write(s)
    outF.close()
    os.remove(zipfile)


outfile.close()


