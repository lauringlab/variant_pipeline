import sys
import numpy as np
import pysam
from collections import Counter
import vcf
import time


## date and time representation


in_var= vcf.Reader(open(sys.argv[1], 'r'))
out_csv=open(sys.argv[2],'w')

# write meta data
for meta in in_var.metadata:
    out_csv.write('#'+meta+':'+str(in_var.metadata[meta])+'\n')

# write INFO
for info in in_var.infos:
    out_csv.write('#'+info+':'+str(in_var.infos[info])+'\n')

# write filters

for filter in in_var.filters:
    out_csv.write('#'+filter+':'+str(in_var.filters[filter])+'\n')

#write the header
header='CHROM,POS,ID,REF,ALT,QUAL,FILTER,'

first=vcf.Reader(open(sys.argv[1], 'r'))
first=first.next() # first varaint to get the infos that are present
print first
for x in first.INFO:
    header=header+x+","

out_csv.write(header+'\n')

for record in in_var:
    print record
    line=str(record.CHROM)+','+str(record.POS)+','+ str(record.ID)+','+str(record.REF)+','+str(record.ALT[0])+','+str(record.QUAL)+','+str(record.FILTER)+','
    for y in record.INFO:
        line=line+str(record.INFO[y]).replace(",",":")+","
    print line
    out_csv.write(line+'\n')

out_csv.close()
