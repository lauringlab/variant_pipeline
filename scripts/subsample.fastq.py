# From http://seqanswers.com/forums/showthread.php?t=12070

import sys, random, itertools
import HTSeq

fraction = float( sys.argv[1] )
in1 = iter( HTSeq.FastqReader( sys.argv[2] ) )
in2 = iter( HTSeq.FastqReader( sys.argv[3] ) )
out1 = open( sys.argv[4], "w" )
out2 = open( sys.argv[5], "w" )

for read1, read2 in itertools.izip( in1, in2 ):
   if random.random() < fraction:
      read1.write_to_fastq_file( out1 )
      read2.write_to_fastq_file( out2 )

out1.close()
out2.close()
