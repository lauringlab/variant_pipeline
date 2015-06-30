import sys
import os
import gzip
# input argument is the sample sheet
junk_names = [] # The names given by the sequencing core sampleid_index
new_names = []
# add the bad names from the sample sheet to a list
f =  open(sys.argv[1],"r") 
next(f)
for line in f:
    line = line.strip()
    line = line.split(',')
    junk_names.append(line[1]+"_"+line[3])
    new=line[4]#.replace("_","-")
    new_names.append(new)
f.close()

outfile = open('renaming_files.txt','w')
os.chdir(sys.argv[2])
for filename in os.listdir(sys.argv[2]):
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
        print("WRITING "+ filename + " to "+perfect_name)
        outfile.write(filename + "\t WRITTEN to \t" + perfect_name + "\n")
        os.rename(filename,perfect_name)
outfile.close()

