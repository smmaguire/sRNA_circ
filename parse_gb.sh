#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 1
#$ -m e

source activate sRNA_circ_spades

Rscript parse_gb.R ${input_dir} ${output_dir_current} ${file_name}

echo "run_now???"