#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 4
#$ -m e

source activate sRNA_circ_spades
cd ${current_dir}
rm *aggregated*
cat *.fasta >> "aggregated_"${bc}"_output.fasta"
cat *.tab >> "aggregated_"${bc}"_output.tab"
flexbar -r "aggregated_"${bc}"_output.fasta" \
-a ${ap} \
--adapter-trim-end ANY \
--min-read-length 10 \
-R "aggregated_"${bc}"trimmed_output.fasta" \
--adapter-error-rate 0.2 \
--adapter-min-overlap 6
Rscript /home/smaguire/work/sRNA_circ/spades/scripts/sRNA_circ/process_trim.R ${current_dir}
