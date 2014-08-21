#!/bin/sh


dir=$1

for file in $dir/*.fastq ; do
    output=`echo $file | sed 's/\(.*\)_.*_R\([0-9]\).*/\1.\2.fastq/'`
    test "$file" != "$output" && mv "$file" "$output" || echo " No changes to make " 

done








