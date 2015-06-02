#!/usr/bin/python

# adapted from  http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python/24652788#24652788

# Usage input the bash style reg ex using *s and this will combine all files together only including the header from the first one

import glob
import os
import sys

os.chdir(sys.argv[1])
ending=sys.argv[2]

csv_files=[]
for file in glob.glob("*"+ending):
    csv_files.append(file)
print(csv_files)

fout=open(sys.argv[3],"w")
# first file:
with open(csv_files[0],"r") as first_file:
    print("working with " + csv_files[0])
    header=first_file.readline().strip()+"\tSample\n"
    fout.write(header)
    next(first_file)
    sample=csv_files[0].split('.')[0]
    for line in first_file:
        fout.write(line.strip()+'\t'+sample+'\n')

# now the rest:    
for num in range(1,len(csv_files)):
    print("working with " + csv_files[num])
    with open(csv_files[num],"r") as next_file:
        next_file.next() # skip the header
        sample=csv_files[num].split('.')[0]
        for line in next_file:
            fout.write(line.strip()+'\t'+sample+'\n')
fout.close()
