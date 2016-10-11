#import sys
import os
import gzip
import argparse
import shutil
import glob
parser = argparse.ArgumentParser(description='This program takes hiseq fastq files and renames them as sample.read_direction.#.fastq and keeps a log of the change, it may need to be adjusted if the sequencing core changes their scheme')
parser.add_argument('-s',action='store',dest='s',help='The sorce directory containing the original fastq files')
parser.add_argument('-f',action='store',dest='f',help='The final directory that will hold the renamed fastq files')
parser.add_argument('-k',action='store',dest='key',help='The output csv given by the sequencing core that serves as the key for renaming')
parser.add_argument('-d',action='store_true',dest='dir',default=False,help='Boolean switch to use of the fastqs are in separate directories in this case the source directory should contain the sample directories which in turn hold the fastq. Some hiseq runs are returned like this.')
parser.add_argument('-run',action='store_true',dest='test',default=False,help='Boolean switch to run program, without this the program runs in test mode: the log is made but no files are renamed')
parser.add_argument('-old',action='store_true',dest='old',default=False,help='Boolean switch to use if the run is was done prior to July 2015 as the sequencing core constantly changes their sample sheets')
parser.add_argument('-l',action='store',dest='l',help='The output file that will keep a log of the name changes')

args=parser.parse_args()
s=args.s
f=args.f
key=args.key
test=args.test
dir=args.dir
old=args.old
log=args.l
#make the output dir if it doesn't exist
if not os.path.exists(f):
    os.makedirs(f)

# add slash to final location if it was forgotten
if f[-1] != '/':
    f=f+'/'
#print f

if s[-1] != '/':
    s=s+'/'
#print s

if dir==True:
    s=s+'*/'

# input argument is the sample sheet
junk_names = [] # The names given by the sequencing core sampleid_index
new_names = [] # the names we want, which are in the same row as the data taht is used to make the bad names. 

# so junk_names[i] should be the junk name for new_names[i] as they both come from the ith column.
# add the bad names from the sample sheet to a list. This uses columns of the sheet to construct the names according to the capricous nature of the sequencing core.Although they may have finally settled down a bit.
names =  open(key,"r")
next(names) # skip the header of the csv
for line in names:
    line = line.strip()
    line = line.split(',')
    if old==True:
        junk_names.append(line[2]+"_"+line[4])
        new=line[5]
    if old==False:
        junk_names.append(line[2])    
        new=line[9]   
    new_names.append(new)
    #print(line[9])
names.close()

if test==False:
    print "running in test mode add option -run to run"
    
outfile = open(log,'w')


# Now to search through the fastq files or fastq.gz and copy them to the new makes
for filename in glob.glob(s + "*.fastq"):
    path=os.path.abspath(filename) # getteh full name and path
    if old==True: # split the fastq name before the extra data added by illumina
        name=filename.split("_L")
    if old==False:
        name=filename.split("_S")
    bad_name = name[0] # the bad name is the first bit
    lane_junk = name[1]
    read_number=lane_junk.split("_R") # info about the read 
    fastq_number=read_number[1][4]
    read_number=read_number[1][0]
    if bad_name in junk_names: # the new name and the bad name have the same index in their respective lists
        name_index = junk_names.index(bad_name)
        better_name= new_names[name_index]
        perfect_name= better_name+"."+read_number+"."+fastq_number+".fastq"
        # Write file to new name
        print("COPYING "+ path + " to "+f+perfect_name)
        outfile.write(path + "\t COPIED to \t" + f+perfect_name + "\n") # log the move in a renaming_log.txt
        if test==True:
            shutil.copy(path,f+perfect_name)
            
#print(junk_names)
#print(new_names)
#This is the same except for zipped files
for zipfilename in glob.glob(s + "*.fastq.gz"):
    path=os.path.abspath(zipfilename)
    zipfilename= os.path.basename(zipfilename)
    if old==True:
        name=zipfilename.split("_L")
    if old==False:
        name=zipfilename.split("_S")
    bad_name = name[0]
    #print(bad_name)
    lane_junk = name[1]
    read_number=lane_junk.split("_R")
    fastq_number=read_number[1][4]
    read_number=read_number[1][0]
    if bad_name in junk_names:
        name_index = junk_names.index(bad_name)
        better_name= new_names[name_index]
        perfect_name= better_name+"."+read_number+"."+fastq_number+".fastq.gz"
        # Write file to new name
        print("COPYING "+ path + " to "+f+perfect_name)
        outfile.write(path + "\t COPIED to \t" + f+perfect_name + "\n")
        if test==True:
            shutil.copy(path,f+perfect_name)
            
#if test==True:
#    for zipfile in glob.glob(f + "*.gz"):
#        print "unzipping:" + zipfile

 #       inF = gzip.GzipFile(zipfile, 'rb')
 #       s = inF.read()
 #       inF.close()
        #print(os.path.splitext(zipfile)[0])
 #       outF = file(os.path.splitext(zipfile)[0], 'wb')
 #       outF.write(s)
 #       outF.close()
 #       os.remove(zipfile)


outfile.close()


