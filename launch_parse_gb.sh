#!/bin/bash
data_path=/home/smaguire/work/sRNA_circ/spades/human_brain_output/
output_main=/home/smaguire/work/sRNA_circ/spades/human_brain_output/output/
adapter_path=/home/smaguire/work/sRNA_circ/spades/adapter.fa

for j in $( ls -d $data_path"barcode"*); do
    barc=${j#$data_path}
    mkdir -p $output_main$barc
    output_dir=$output_main$barc
    for i in $( ls -d $j/*"output"); do
        echo "output dir"
        echo $output_dir
        echo "input directory"
        echo $i
        new_name=${i#$j/}
        new_name=${new_name%"_output"}
        echo $new_name
        qsub -v input_dir=$i,output_dir_current=$output_dir,file_name=$new_name,adapter=$adapter_path,bc=$barc parse_gb.sh
    done
done