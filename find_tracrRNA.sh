#!/bin/sh

#useage sh find_tracrRNA.sh tracrRNA.fasta query_repeat word_size reward penalty gapoepn gapextend task format
tracrRNA_p=$1	#tracrRNA_potential.fasta
query=$2	#repeat.fasta
word_size=$3	#7
reward=$4	#1
penalty=$5	#-2
gapopen=$6	#1
gapextend=$7	#2
task=$8		#blastn-short
format=$9 #7

makeblastdb -in $tracrRNA_p -dbtype nucl -logfile ${tracrRNA_p/.fasta/.log}
blastn -query $query  -db $tracrRNA_p -word_size $word_size -reward $reward -penalty $penalty -gapopen $gapopen -gapextend $gapextend -task $task -outfmt $format
