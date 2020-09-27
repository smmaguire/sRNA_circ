#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 1
#$ -m e

data_path=/mnt/nanopore/Sean_Maguire2/sRNA_nanopore_final/FC2_mirX/basecalled/
mkdir -p /home/smaguire/work/sRNA_circ/spades/aggregated_mirx_data/
agg_path=/home/smaguire/work/sRNA_circ/spades/aggregated_mirx_data/

for j in $( ls -d $data_path"barcode"*); do
    barc=${j#$data_path}
    mkdir -p $agg_path$barc
    file_count=0
    for i in $( ls $j/*.fastq); do
        cat $i >> $agg_path$barc/"agg_file_"$file_count".fastq"
        line_count=`wc -l $agg_path$barc/"agg_file_"$file_count".fastq" | cut -f1 -d' '`
        echo $line_count
        if [ $line_count -gt 100000 ]; then
            ((file_count++))
        fi
    done
done
