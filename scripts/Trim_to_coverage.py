#!/usr/bin/env python

import sys

file = open(sys.argv[1],'r')
out = open(sys.argv[2],'w')

out.write('Segment,Position,Coverage\n')

for line in file:

    line = line.split('\t')
    out.write(line[0]+','+line[1]+','+line[3]+'\n')
out.close()
file.close()
