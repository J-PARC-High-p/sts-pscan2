#!/bin/bash


FILES=`cat $1`

for FILE in $FILES ; do
    echo $FILE
    root -l <<EOF
execution("$FILE")
.q
EOF
done

