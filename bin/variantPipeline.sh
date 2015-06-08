#!/bin/bash

EXPECTED_ARGS=4
if [ $# -ne $EXPECTED_ARGS ]; then
    echo "Usage: `basename $0` {input_dir} {output_dir} {reference} {plasmid control name}"
    echo "Example `basename $0` /projectA/input_fastq projectA /bowtie2/pWH2000 PR8_con"
    echo "  input_dir should be absolute path to dir containing fastq pairs (left and right)"
    exit 1
fi

INPUT_DIR=`readlink -f $1`
OUTPUT_DIR=`readlink -f $2`
REF=`readlink -f $3`
CONTROL=$4
BIN_DIR=`dirname $0` 
SCRIPT_DIR=`readlink -f $BIN_DIR/../scripts/`
LIB_DIR=`readlink -f $BIN_DIR/../lib`
BPIPE_COMMAND=$LIB_DIR/bpipe-0.9.8.7/bin/bpipe

echo "Processing fastqs from [$INPUT_DIR]."
echo "Results will be saved to [$OUTPUT_DIR]."
echo " Using [$REF] for a reference."

mkdir -p $OUTPUT_DIR || exit 1
cd $OUTPUT_DIR
cp -a $SCRIPT_DIR/variantPipeline.bpipe* . 


# add variables to config reference to config file

sed -i '7iREFERENCE="'$REF'"' variantPipeline.bpipe.config.groovy
sed -i '8iREFERENCE_FA="'$REF'.fa"' variantPipeline.bpipe.config.groovy
sed -i '9iSCRIPTS="'$SCRIPT_DIR'"' variantPipeline.bpipe.config.groovy
sed -i '10iLIBRARY_LOCATION="'$LIB_DIR'"' variantPipeline.bpipe.config.groovy
sed -i '11iCONTROL="'$CONTROL'"' variantPipeline.bpipe.config.groovy
#sed -i '12iCONTROL="'$LIB_DIR'"' variantPipeline.bpipe.config.groovy

#throttled to 8 processors to be a good neighbor.
#note that running unthrottled can result in errors when bpipe overallocates threads/memory
# Run the control first so it is ready to be used as a control in the deepSNV step of the test samples
time $BPIPE_COMMAND run -n 8  variantPipeline.bpipe.groovy $INPUT_DIR/*.fastq
