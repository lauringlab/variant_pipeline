#!/bin/bash

EXPECTED_ARGS=2
if [ $# -ne $EXPECTED_ARGS ]; then
    echo "Usage: `basename $0` {input_dir} {output_dir}"
    echo "Example `basename $0` /projectA/input_fastq projectA"
    echo "  input_dir should be absolute path to dir containing fastq pairs (left and right)"
    exit 1
fi

INPUT_DIR=`readlink -f $1`
OUTPUT_DIR=`readlink -f $2`
BIN_DIR=`dirname $0` 
SCRIPT_DIR=`readlink -f $BIN_DIR/../scripts/`
LIB_DIR=`readlink -f $BIN_DIR/../lib`
BPIPE_COMMAND=$LIB_DIR/bpipe-0.9.8.2/bin/bpipe

echo "Processing fastqs from [$INPUT_DIR]."
echo "Results will be saved to [$OUTPUT_DIR]."

(
mkdir -p $OUTPUT_DIR || exit 1
cd $OUTPUT_DIR
cp -a $SCRIPT_DIR/variantPipeline.bpipe* . 

#throttled to 8 processors to be a good neighbor.
#note that running unthrottled can result in errors when bpipe overallocates threads/memory
time $BPIPE_COMMAND run -n 8 variantPipeline.bpipe $INPUT_DIR/*.fastq

) 
