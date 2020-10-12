#!/bin/bash
data_path=/home/smaguire/work/sRNA_circ/spades/aggregated_mirx_data/

for j in $( ls -d $data_path"barcode"*); do
    barc=${j#$data_path}
    for i in $( ls $j/*.fastq); do
        new_name=${i#$j/}
        new_name=${new_name%.fastq}
        echo "barcode"
        echo $barc
        echo "name"
        echo $new_name
        qsub -v reads=$i,name=$new_name,barcode=$barc spades.sh
    done
done
