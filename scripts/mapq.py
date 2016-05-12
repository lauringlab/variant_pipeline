#!/usr/bin/python

# Usage :  This script is used to find quality information regarding a list of variants. It takes as inputs a csv of variants (as outputed by the deepSNV script) and the appropriate bam file.  

# inputs file.csv file.bam output.sum.csv output.reads.csv
#output.sum.csv just adds info to the summary csv that is used as input from deepSNV
#output.read.csv is a file that contains on each read for each variant.  It counts how many reads for a given variant have a given value for a given variable  

import sys
import numpy as np
import pysam
from collections import Counter

csv=sys.argv[1].split(".csv")[0]
csv=csv.split("/")[-1]
bam=sys.argv[2].split(".bam")[0]
bam=bam.split("/")[-1]

#Ensure the bam and csv file match
if csv==bam:
    print("working with csv: "+csv + " and bam: " +bam)
else:
   print( "bam:"+bam +" does not match csv: " +csv)
   sys.exit(1)
   print("We quit so this should not print")





#Put header on the output.read  file
with open(sys.argv[4],"w") as read_csv:
    read_csv.write("Sample,Mutation,Variable,Value,Count\n")


with pysam.AlignmentFile(sys.argv[2], "rb") as bamfile: # open the bam file
    with open(sys.argv[1],'r') as in_var: #open the input summary file
        with open(sys.argv[3],"w") as outfile_csv: # open the output summary file

            header=in_var.readline().strip()+",\"MapQ\"" +",\"Read_pos\""+",\"Phred\""+'\n' # add the mapq column ect. to the header of the input csv
            outfile_csv.write(header) # write the header to the output file


            for line in in_var: # each line is a variant
                mapq=[]# This will hold a list of the mapping qualtties that map to the variant
                phred=[] #This will hold a list of the phred that map to the variant
                Read_pos=[]  # This will hold a list of position relative to the read
                line=line.strip()
                data=line # save the line as it is to append the data to later
                line=line.replace('"','').split(",") # get rid of quotes and split by ,
                chr=line[0] # the chr or segment column
                pos=int(line[1]) #The position relative to the chromosome
                py_pos=pos-1 # deepSNV starts at 1 python starts at 0
                ref=line[2] # the reference column
                var=line[3] # The variant column
                sample= line[16] # The sample column
                mutation= line[17] # mutation name
               # print(mutation)
                for pileupcolumn in bamfile.pileup(chr,py_pos,py_pos+1,truncate=True,stepper="all",max_depth=1E6): # look at the range py_pos to py_pos+1 and all the bases that align to those positions. Think of a column on each position like the samtools output of mpileup
                    if pileupcolumn.pos==py_pos: #If the position is the one we want
                        for pileupread in pileupcolumn.pileups: #Now for each read in that column
                            if  not pileupread.is_del and not pileupread.is_refskip:
                                called_base=pileupread.alignment.query_sequence[pileupread.query_position] # what is the called base at the position
                                called_phred=pileupread.alignment.query_qualities[pileupread.query_position] # what is the phred of that base
                                if called_phred>=30 and called_base==var: # change this if you change the phred cut off in deepSNV. deepSNV only looks a phred>30. and we only want those calles that match the variant.
                                    mapq.append(pileupread.alignment.mapping_quality)# add the mapping quality of the read to the list
                                    phred.append(called_phred)
                                    Read_pos.append(pileupread.query_position)
                mean_map=np.mean(mapq) # get the averages of this quality scores
                mean_phred=np.mean(phred)
                mean_Read_pos=np.mean(Read_pos)
                if mean_map==[]:
                    print( "OOPS didn't find the variant looks like you didn't fix the bug") # Just a check to make sure we are finding every variant we expect to find
                    sys.exit(1)
# add the new data to the summary file
                data=data+","+str(mean_map)+","+str(mean_Read_pos)+","+str(mean_phred)+'\n'
                outfile_csv.write(data) # write to the summary file


#count up the occurances of each value for each quality metric
                mapq=Counter(mapq)
                phred=Counter(phred)
                Read_pos=Counter(Read_pos)
# Write this data to the output.read file
                with open(sys.argv[4],"a") as read_csv:

                    for qual in mapq:
                        read_csv.write(sample+","+mutation+","+"MapQ"+","+str(qual)+","+str(mapq[qual])+'\n')
                    for qual in Read_pos:
                        read_csv.write(sample+","+mutation+","+"Read_pos"+","+str(qual)+","+str(Read_pos[qual])+'\n')
                    for qual in phred:
                        read_csv.write(sample+","+mutation+","+"Phred"+","+str(qual)+","+str(phred[qual])+'\n')
