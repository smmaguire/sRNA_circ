#!/bin/bash
data_path=/home/smaguire/work/sRNA_circ/spades/data2/
output_main=/home/smaguire/work/sRNA_circ/spades/output/
for j in $( ls -d $data_path"barcode"*); do
    barc=${j#$data_path}
    mkdir -p $output_main$barc
    output_dir=$output_main$barc
    for i in $( ls -d $j/*"output"); do
        echo "output dir"
        echo $output_dir
        echo "input directory"
        echo $i
        qsub -v input_dir=$i,output_dir_current=$output_dir parse_gb.sh  
    done
done
