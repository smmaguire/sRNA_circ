#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 1
#$ -m e

source activate bbtools
#source activate sRNA_circ_spades

#echo "#################################################"
#echo "Starting to parse the genbank"
#echo "file is: "${file_name}
#echo "#################################################"

#Rscript parse_gb.R ${input_dir} ${output_dir_current} ${file_name}

cd ${output_dir_current}

#echo "#################################################"
#echo "Running flexbar"
#echo "#################################################"

# flexbar -r "temp_fa_"${file_name}".fasta" \
# -a ${adapter} \
# --adapter-trim-end ANY \
# --min-read-length 10 \
# -R ${file_name}"_trimmed_output.fasta" \
# --adapter-error-rate 0.2 \
# --adapter-min-overlap 6

# echo "#################################################"
# echo "Getting line count"
# echo "#################################################"

# line_count=`wc -l ${file_name}"_trimmed_output.fasta" | cut -f1 -d' '`
# echo $line_count

# echo "#################################################"
# echo "Getting alignment stats"
# echo "#################################################"

# Rscript /home/smaguire/work/sRNA_circ/spades/scripts/sRNA_circ/process_trim_v3.R ${file_name}"_trimmed_output.fasta" ${file_name} ${line_count}

#echo "#################################################"
#echo "Mapping MirXplore"
#echo "#################################################"

#source activate shortstack
#index="/home/smaguire/work/unblock_remakes/data/mirexplore/genome/miRexplore"
##bowtie -n 1 -l 10 -m 100 -k 1 --sam --best --strata $index $downsample_dir/${name}.fastq | samtools view -b - | samtools sort -o $mapped_dir${name}.bam -
#bowtie -v 2 -m 100 -k 1 --sam --best --strata $index -f ${file_name}"_trimmed_output.fasta" | samtools view -b - | samtools sort -o ${file_name}"_mapped.bam" -

#bbmap.sh in=${file_name}"_trimmed_output.fasta" covstats=${file_name}"_covstats.txt" ref=/home/smaguire/work/unblock_remakes/data/mirexplore/mwulf_data/miRxplore.fasta out=${file_name}"_mapped.bam" nodisk local

# echo "#################################################"
# echo "Counting miRNAs"
# echo "#################################################"

# samtools index ${file_name}"_mapped.bam"
# samtools idxstats ${file_name}"_mapped.bam" > ${file_name}"_miRNA_counts.txt"
# samtools flagstat ${file_name}"_mapped.bam" > ${file_name}"_mapping_stats.txt"


echo "#################################################"
echo "Map Human"
echo "#################################################"

bbmap.sh in=${file_name}"_trimmed_output.fasta" \
ref=/home/smaguire/work/human_miRNA/star_pipeline/genome/GRCh38.fasta \
path=/home/smaguire/work/sRNA_circ/spades/genome_ref \
k=10 \
local \
minid=0.9 \
out=${file_name}.sam bamscript=bs.sh; sh bs.sh



