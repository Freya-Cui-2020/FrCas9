CRISPRCasFinder="CRISPRCasFinder.pl -so /Bigdata/crispr/app/CRISPRCasFinder-2/sel392v2.so -keep -ccc -cas -ccvr -log -cf CasFinder-2.0.2 -force -html" 
check_step_1="./check_crisprcasfinder.py"
find_protospacer="./find_protospacers.sh"
extract_flanking="./deal_blast.py"
check_array_direction="./check_array_direction.py"
find_tracrRNA="./find_tracrRNA.sh"
confirm_range="./confirm_range.py"
protodb="/Bigdata/crispr/database/CRISPR_pipeline/Database/protospacer/version/without_virus_nodup/withou_virus_v2019622_nodup.fasta"

#genome must be included in fasta_dir
#The program will make a new directory under $result_dir
fasta_dir=$1
genome=$2
result_dir=$3

name=${genome/.fasta/}

seqtk comp ${fasta_dir}/${genome}|awk '{print $1"\t"$2}'> ${fasta_dir}/${name}.bed

$CRISPRCasFinder -in ${fasta_dir}/${genome} -outdir ${result_dir}/${name}_result \
-drpt ${result_dir}/${name}_result/supplementary_files/repeatDirection.tsv \
-rpts ${result_dir}/${name}_result/supplementary_files/Repeat_List.csv \
-dbc ${result_dir}/${name}_result/supplementary_files/CRISPR_crisprdb.csv

python3 $check_step_1 ${result_dir}/${name}_result >${result_dir}/${name}_result/check_step_1.log

check_step_1_result=`grep 'find crisprs' ${result_dir}/${name}_result/check_step_1.log`

if [ "$check_step_1_result" ];
then
    echo "start write files in TSV/"
    spacers_file=`ls ${result_dir}/${name}_result/CRISPRFinderProperties/*_properties/Spacers/Spacers_*[^fasta]`
    cas_report=`ls ${result_dir}/${name}_result/TSV/Cas_REPORT.tsv`
    crispr_report=`ls ${result_dir}/${name}_result//TSV/Crisprs_REPORT.tsv`
    count=0
    best_array=""
    best_array_contig=""
    for i in ${spacers_file[@]}
    do
        if [ "$count" -lt `grep '^>' $i|wc -l|awk '{print $1}'` ];
        then
          best_array=$i
          best_array_contig=`basename -s "_properties" ${i/Spacers*/}`
          count=`grep '^>' $i|wc -l|awk '{print $1}'`
        fi
        echo "$find_protospacer $protodb $i 18 1 -1 5 2 megablast 5 > ${i}_blast_result"
	sh $find_protospacer $protodb $i 18 1 -1 5 2 megablast 5 > ${i}_blast_result
        echo "python $extract_flanking ${i}_blast_result ${protodb} ${i}_Up.fasta ${i}_Down.fasta ${i}_summary.txt ${i}_Medium.fasta"
	echo "$best_array"
	python $extract_flanking ${i}_blast_result ${protodb} ${i}_Up.fasta ${i}_Down.fasta ${i}_summary.txt ${i}_Medium.fasta
        
        #build spacer file
        awk '{print $1}' ${i}_Up.fasta > ${i}_Up_line1.fasta
        awk '{print $1}' ${i}_Medium.fasta > ${i}_Medium_line1.fasta
        paste -d "" ${i}_Up_line1.fasta ${i}_Medium_line1.fasta ${i}_Down.fasta > ${i}_SPACERS.fasta
        
        drid=${i//Spacers/DRs}
        python $check_array_direction $drid ${drid}_repeat >${drid}.log
	echo ${best_array}
    echo -e "${best_array_contig}\t${count}\t${drid}"
    done
    cds=`ls ${result_dir}/${name}_result/Prodigal/prodigal_${best_array_contig}/${best_array_contig}.gff`
    strand_crisprs=`awk -F'\t' '{print $2}' $(echo ${best_array//Spacers/DRs}).log`
    row_num_crisprs=${best_array##*_}
    
    #START generating bedfiles
    echo "start generating bedfiles"
    bed1=`awk -F'\t' -v n=$row_num_crisprs -v c=$best_array_contig -v s=$strand_crisprs 'BEGIN{ac=0}$2~c{ac++;arr[ac]=$2"\t"$6"\t"$7"\t"s"\tarray_"n}END{print arr[n]}' $crispr_report`
    best_array_contig=`echo -e $bed1|awk '{print $1}'`
    bed2=`awk -F'\t' -v c=$best_array_contig '$1~c && $6!=""{print c"\t"$6"\t"$7"\t"$8"\t"$2}' $cas_report`  
    echo -e "${bed1}\n${bed2}"|sort -k2 -n|awk '{print $1"\t"$2-1"\t"$3-1"\t"$4"\t"$5}'> ${result_dir}/${name}_result/${name}.bed
    bedtools complement -i ${result_dir}/${name}_result/${name}.bed -g ${fasta_dir}/${name}.bed>${result_dir}/${name}_result/tmp.bed
    mkdir -p ${result_dir}/${name}_result/tracrRNA
    awk 'BEGIN{while((getline t < ARGV[1]) > 0)last++;close(ARGV[1])}NR==1{if($2<$3-500){print $1"\t"$3-1-500+1"\t"$3-1}else{print $1"\t"$2"\t"$3-1}}NR>1&&NR<last{print $1"\t"$2+1"\t"$3-1}NR==last{if($3>$2+500){print $1"\t"$2+1"\t"$2+1+500-1}else{print $1"\t"$2+1"\t"$3}}' ${result_dir}/${name}_result/tmp.bed>${result_dir}/${name}_result/${name}_tracrRNA_potential.bed

    rm ${result_dir}/${name}_result/tmp.bed
    bedtools getfasta -fi ${fasta_dir}/${genome} -bed ${result_dir}/${name}_result/${name}_tracrRNA_potential.bed >${result_dir}/${name}_result/tracrRNA/${name}_tracrRNA_potential.fasta
    sh $find_tracrRNA ${result_dir}/${name}_result/tracrRNA/${name}_tracrRNA_potential.fasta $(echo ${best_array//Spacers/DRs})_repeat 7 1 -2 1 2 blastn-short 5 >${result_dir}/${name}_result/tracrRNA/${name}_tracrRNA_blast_result
    r1=`awk -F'\t' 'NR==1{print $2}' ${result_dir}/${name}_result/${name}_tracrRNA_potential.bed`
    r2=`awk -F'\t' 'END{print $3}' ${result_dir}/${name}_result/${name}_tracrRNA_potential.bed`
    echo "python $confirm_range $r1 $r2 ${result_dir}/${name}_result/tracrRNA/${name}_tracrRNA_blast_result $cds"
    python $confirm_range $r1 $r2 ${result_dir}/${name}_result/tracrRNA/${name}_tracrRNA_blast_result $cds    
fi
