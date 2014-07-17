#!/bin/sh


start=$1

for dir in $start/*; do
    outdir=`echo $dir | sed 's/-/_/g'`
    test "$dir" != "$outdir" && mv "$dir" "$outdir";dir=$outdir || echo " No changes to make to $dir"
    for file in $dir/*.* ; do
        output=`echo $file | sed 's/-/_/g'`
        test "$file" != "$output" && mv "$file" "$output" || echo " No changes to make in $dir" 

    done
done






