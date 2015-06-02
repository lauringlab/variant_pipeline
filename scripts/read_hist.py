#!/usr/bin/python

# The goal here is to bin the mapping quality of the reads for each variant called in the sample to save on space and speed in future work

#The csvs are of the form Sample,Mutation,Category, Mapq

import sys
import csv
from collections import Counter



with open(sys.argv[1]) as infile:
    sample=infile.readlines()[1].split(",")[0]

mutations=set()
with open(sys.argv[1]) as infile:
    next(infile)
    for line in infile:
        mutations.add(line.split(",")[1])
mutations=list(mutations)
#print(mutations)
with open(sys.argv[2],'w') as outfile:
    outfile.write("Sample,Mutation,Category,Variable,Qauntity,Count\n")

        
    for mut in mutations:
        #print(mut)
        mapq=[]
        phred=[]
        Read_pos=[]
        category=[]
        with open(sys.argv[1]) as infile:
            next(infile)
            for line in infile:
                mutation=line.split(",")[1].strip()
                Q = line.split(",")[3].strip()
                Pos = line.split(",")[4].strip()
                Red = line.split(",")[3].strip()
                if mutation == mut:
                    mapq.append(Q)
                    Read_pos.append(Pos)
                    phred.append(Red)
                    category=line.split(",")[2].strip()
        mapq=Counter(mapq)
        phred=Counter(phred)
        Read_pos=Counter(Read_pos)

        for qual in mapq:
            outfile.write(sample+","+mut+","+category+","+"MapQ"+","+qual+","+str(mapq[qual])+'\n')
        for qual in Read_pos:
            outfile.write(sample+","+mut+","+category+","+"Read_pos"+","+qual+","+str(Read_pos[qual])+'\n')
        for qual in phred:
            outfile.write(sample+","+mut+","+category+","+"Phred"+","+qual+","+str(phred[qual])+'\n')

            
            
            
            