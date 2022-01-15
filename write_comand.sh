#!/bin/bash
dir_List=$(ls -l ./ |awk '/^d/ {print $9}')
echo $dir_List
for i in $dir_List
do
    cd $i
    dir=$(ls -l ./ |awk '!/^d/ {print $9}' |grep -i "fasta" )
    #echo $dir
    mkdir /Bigdata/crispr/lmy/results/$i
    for j in $dir
    do
        echo "bash /Bigdata/crispr/lmy/scripts/pipe_lmyMOD.sh /Bigdata/crispr/lmy/CasPDB_contig/2A/"$i $j "/Bigdata/crispr/lmy/results/"$i" 1 > /Bigdata/crispr/lmy/results/$i/$i.$j.out 2>/Bigdata/crispr/lmy/results/$i/$i.$j.out.err "
    done
    cd -
done
