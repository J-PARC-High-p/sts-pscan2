#!/bin/bash

DIR=${1:-}

start_time=$(date "+%s")

if [ -e $DIR ]; then
echo $DIR
root -b <<EOF
.L execution.C
execution_multi("$DIR")
.q
EOF
root -b <<EOF
.L DirPscanAna_yamada_20240402.C
DirPscanAna("$DIR")
.q
EOF
	
fi

end_time=$(date "+%s")
run_time=$(($end_time - $start_time))
echo "実効時間 : $run_time s"
