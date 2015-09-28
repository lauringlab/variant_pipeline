#!/usr/bin/python

import sys
import numpy as np
import pysam
from collections import Counter
import vcf 


input=sys.argv[1].split(".")[0]
input=input.split("/")[-1]
bam=sys.argv[2].split(".")[0]
bam=bam.split("/")[-1]

#Ensure the bam and csv file match
#if input==bam:
#        print("working with vcf: "+input + " and bam: " +bam)
#else:
#   print( "bam:"+bam +" does not match vcf: " +input)
#   sys.exit(1)






#Put header on true false and exp csv files

in_var= vcf.Reader(open(sys.argv[1], 'r'))
variants=list(in_var)
with pysam.AlignmentFile(sys.argv[2], "rb") as bamfile:

#with pysam.AlignmentFile('../data/5_5.removed.bam', "rb") as bamfile:
	for record in variants:
		#print(record)
		mapq=[]# This will hold a list of the mapping qualtties that map to the variant
		phred=[] #This will hold a list of the phred that map to the variant
		Read_pos=[]  # This will hold a list of position relative to the read
		chr=record.CHROM
		pos=int(record.POS)
		py_pos=pos-1
		ref=record.REF
		var=record.ALT[0]
		for pileupcolumn in bamfile.pileup(chr,py_pos,py_pos+1,truncate=True,stepper="all",max_depth=1E6):
			if pileupcolumn.pos==py_pos:
				for pileupread in pileupcolumn.pileups:
					if  not pileupread.is_del and not pileupread.is_refskip:
						called_base=pileupread.alignment.query_sequence[pileupread.query_position]
						called_phred=pileupread.alignment.query_qualities[pileupread.query_position]
						if called_phred>25 and called_base==var: # change this if you change the phred cut off in deepSNV
							mapq.append(pileupread.alignment.mapping_quality)
							phred.append(called_phred)
							Read_pos.append(pileupread.query_position)
		mean_map=np.mean(mapq)
		mean_phred=np.mean(phred)
		mean_Read_pos=np.mean(Read_pos)
		if mean_map==[]:
			print( "OOPS didn't find the variant looks like you didn't fix the bug")
			sys.exit(1)

		record.INFO.update({'MAPQ':mean_map, 'Read_pos':mean_Read_pos, 'phred':mean_phred})
#         print record.INFO
#         print 'done with \n'
#         print record
print "done updating"
iter(variants)
vcf_writer = vcf.Writer(open(sys.argv[3], 'w'), in_var)
csv=open(sys.argv[4],'w')
header="CHROM,POS,REF,ALT,ID"
for entry in record.INFO.keys():
    header=header+","+entry
csv.write(header+"\n")
for record in variants:
    line=record.CHROM+","+str(record.POS)+","+record.REF+","+str(record.ALT[0])+","+input
    for entry in record.INFO:
        line=line+","+str(record.INFO[entry]).replace(",",":")
    csv.write(line+"\n")
    vcf_writer.write_record(record)    
csv.close()
vcf_writer.close()
    
    
    
