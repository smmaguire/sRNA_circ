#!/bin/bash
data_path=/home/smaguire/work/sRNA_circ/spades/output/
adapter_path=/home/smaguire/work/sRNA_circ/spades/adapter.fa
for j in $( ls -d $data_path"barcode"*); do
     barc=${j#$data_path}
     qsub -v current_dir=$j,ap=$adapter_path,bc=$barc trim_process_reads.sh
    done
