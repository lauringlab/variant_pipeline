#!/bin/sh


dir=$1

for file in $dir/*.fastq ; do
    output=`echo $file | sed 's/\(.*\)_.*_.*_R\([0-9]\)_.*fastq/\1.\2.fastq/'`
    output=`echo $output| sed 's/-/_/'`
    test "$file" != "$output" && mv "$file" "$output" || echo " No changes to make " 

done








