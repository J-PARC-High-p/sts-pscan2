#!/bin/bash


FILES=`cat $1`

for FILE in $FILES ; do
    echo $FILE
    root -q -b show_hist.C\(\"$FILE\"\)
done

