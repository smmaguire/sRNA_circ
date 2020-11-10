#!/bin/bash
data_path=/home/smaguire/work/sRNA_circ/spades/human_brain_output/output/

for j in $( ls -d $data_path"barcode"*); do
    qsub -v current_dir=$j map_human_genome.sh
done
