#!/bin/bash

EXPECTED_ARGS=4
if [ $# -ne $EXPECTED_ARGS ]; then
    echo "Usage: `basename $0` {reference_fasta} {test.bam} {control.bam} {output_dir}"
    echo "Example: `basename $0` projectA/bowtie2/pr8_01.fasta ../test.bam ../control.bam outputA"
    exit 1
fi

REFERENCE_FASTA=$1
TEST_BAM=$2
CONTROL_BAM=$3
OUTPUT_DIR=$4

if [ ! -f $REFERENCE_FASTA ]; then 
    echo "ERROR: Specified reference fasta file [$REFERENCE_FASTA] does not exist. Review params and try again."
    exit 1
elif [ ! -f $TEST_BAM ]; then 
    echo "ERROR: Specified test bam file [$TEST_BAM] does not exist. Review params and try again."
    exit 1
elif [ ! -f $CONTROL_BAM ]; then 
    echo "ERROR: Specified control bam file [$CONTROL_BAM] does not exist. Review params and try again."
    exit 1    
fi



BIN_DIR=`dirname $0`
SCRIPT_DIR=$BIN_DIR/../scripts/ 
LIBRARY_LOCATION=$SCRIPT_DIR/../lib
R_SCRIPT=$SCRIPT_DIR/deepSNV.r
(
mkdir -p $OUTPUT_DIR || exit 1
cd $OUTPUT_DIR
date
set -x
Rscript --vanilla --slave $R_SCRIPT $LIBRARY_LOCATION $REFERENCE_FASTA $TEST_BAM $CONTROL_BAM
)
