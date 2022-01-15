#!/bin/sh


#useage sh find_protospacer.sh db query word_size reward penalty gapoepn gap extend matrix task 
databse=$1	#protospacer_database
query=$2	#multi-spacers fasta file
word_size=$3	#18
reward=$4	#1
penalty=$5	#-1
gapopen=$6	#5
gapextend=$7	#2
#matrix=$8	#PAM-30 query_length < 35
task=$8		#megablast
format=$9

blastn -query $query  -db $databse -word_size $word_size -reward $reward -penalty $penalty -gapopen $gapopen -gapextend $gapextend -task $task -outfmt $format
