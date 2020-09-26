#!/bin/bash

data_path=/mnt/nanopore/Sean_Maguire2/sRNA_nanopore_final/FC1_test_olio_mouse_brain/Final_seq_sRNA_naopore2/basecalled_fastq/
agg_path=/home/smaguire/work/sRNA_circ/spades/data2/aggregated_data/

for j in $( ls -d $data_path"barcode"*); do
    barc=${j#$data_path}
    mkdir -p $agg_path$barc
    file_count=0
    count=0
    for i in $( ls $j/*.fastq); do
        ((count++))
        cat $i >> $agg_path$barc/"agg_file_"$file_count".fastq"
        if [ $count -gt 100 ]; then
            count=0
            ((file_count++))
        fi
    done
done
