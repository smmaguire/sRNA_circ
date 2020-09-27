#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 1
#$ -m e

source activate sRNA_circ_spades

echo "#################################################"
echo "Starting to parse the genbank"
echo "file is: "${file_name}
echo "#################################################"

#Rscript parse_gb.R ${input_dir} ${output_dir_current} ${file_name}

cd ${output_dir_current}

echo "#################################################"
echo "Running flexbar"
echo "#################################################"

# flexbar -r "temp_fa_"${file_name}".fasta" \
# -a ${adapter} \
# --adapter-trim-end ANY \
# --min-read-length 10 \
# -R ${file_name}"_trimmed_output.fasta" \
# --adapter-error-rate 0.2 \
# --adapter-min-overlap 6

echo "#################################################"
echo "Getting line count"
echo "#################################################"

line_count=`wc -l ${file_name}"_trimmed_output.fasta" | cut -f1 -d' '`
echo $line_count

echo "#################################################"
echo "Getting alignment stats"
echo "#################################################"

Rscript /home/smaguire/work/sRNA_circ/spades/scripts/sRNA_circ/process_trim_v3.R ${file_name}"_trimmed_output.fasta" ${file_name} ${line_count}


