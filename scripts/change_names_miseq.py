import sys
import os
# input argument is fastq directory

outfile = open('renaming_files.txt','w')
os.chdir(sys.argv[1])
for filename in os.listdir(sys.argv[1]):
    name=filename.split("_L")
    good_name = name[0].split("_")[0]
    lane_junk = name[1]
    read_number=lane_junk.split("_R")
    fastq_number=read_number[1][4]
    read_number=read_number[1][0]
    perfect_name= good_name.replace("-","_")+"."+read_number+"."+fastq_number+".fastq"
    # Write file to new name
    print("WRITING "+ filename + " to "+perfect_name)
    outfile.write(filename + "\t WRITTEN to \t" + perfect_name + "\n")
    os.rename(filename,perfect_name)
outfile.close()

