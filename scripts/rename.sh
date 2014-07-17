#!/bin/sh


output=`echo $1 | sed 's/\([A-Z].*\)_S.*_R\([0-9]\).*/\1.\2.fastq/'` 


if [ "$1" != "$output"] ; then

    mv $1 $output;

fi

