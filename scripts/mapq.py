#!/usr/bin/python

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






#Put header on true false and exp csv files
with open(sys.argv[4],"w") as read_csv:
    read_csv.write("Sample,Mutation,Variable,Value,Count\n")


with pysam.AlignmentFile(sys.argv[2], "rb") as bamfile:
    with open(sys.argv[1],'r') as in_var:
        with open(sys.argv[3],"w") as outfile_csv:

            header=in_var.readline().strip()+",\"MapQ\"" +",\"Read_pos\""+",\"Phred\""+'\n' # add the mapq column to the header of the input csv
            outfile_csv.write(header) # write the header to the output file


            for line in in_var:
                mapq=[]# This will hold a list of the mapping qualtties that map to the variant
                phred=[] #This will hold a list of the phred that map to the variant
                Read_pos=[]  # This will hold a list of position relative to the read
                line=line.strip()
                data=line # save the line as it is to append the data to later
                line=line.replace('"','').split(",")
                chr=line[0]
                pos=int(line[1])
                py_pos=pos-1
                ref=line[2]
                var=line[3]
                sample= line[16]
                mutation= line[17] # mutation name
               # print(mutation)
                for pileupcolumn in bamfile.pileup(chr,py_pos,py_pos+1,truncate=True,stepper="all",max_depth=1E6):
                    if pileupcolumn.pos==py_pos:
                        for pileupread in pileupcolumn.pileups:
                            if  not pileupread.is_del and not pileupread.is_refskip:
                                called_base=pileupread.alignment.query_sequence[pileupread.query_position]
                                called_phred=pileupread.alignment.query_qualities[pileupread.query_position]
                                if called_phred>30 and called_base==var: # change this if you change the phred cut off in deepSNV
                                    mapq.append(pileupread.alignment.mapping_quality)
                                    phred.append(called_phred)
                                    Read_pos.append(pileupread.query_position)
                mean_map=np.mean(mapq)
                mean_phred=np.mean(phred)
                mean_Read_pos=np.mean(Read_pos)
                if mean_map==[]:
                    print( "OOPS didn't find the variant looks like you didn't fix the bug")
                    sys.exit(1)

                data=data+","+str(mean_map)+","+str(mean_Read_pos)+","+str(mean_phred)+'\n'
                outfile_csv.write(data) # write to the summary file



                mapq=Counter(mapq)
                phred=Counter(phred)
                Read_pos=Counter(Read_pos)
                with open(sys.argv[4],"a") as read_csv:

                    for qual in mapq:
                        read_csv.write(sample+","+mutation+","+"MapQ"+","+str(qual)+","+str(mapq[qual])+'\n')
                    for qual in Read_pos:
                        read_csv.write(sample+","+mutation+","+"Read_pos"+","+str(qual)+","+str(Read_pos[qual])+'\n')
                    for qual in phred:
                        read_csv.write(sample+","+mutation+","+"Phred"+","+str(qual)+","+str(phred[qual])+'\n')
