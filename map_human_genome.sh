#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 16
#$ -m e


source activate bbtools

cd ${current_dir}

file=$( ls *"final_trimmed_output.fasta")
bc=${file%"final_trimmed_output.fasta"}

echo "#################################################"
echo "Map Human"
echo "#################################################"

bbduk.sh in=$file out=$bc"_final_trimmed_filtered.fasta" minlen=16 overwrite=T

bbmap.sh in=$bc"_final_trimmed_filtered.fasta" \
ref=/home/smaguire/work/human_miRNA/star_pipeline/genome/GRCh38.fasta \
path=/home/smaguire/work/sRNA_circ/spades/genome_ref \
k=12 \
local \
minid=0.9 \
out=$bc"_mapped.sam" \
overwrite=T

echo "#################################################"
echo "Count miRNAs"
echo "#################################################"

source activate htseq
samtools view -b $bc"_mapped.sam" | samtools sort - > $bc"_mapped_sorted.bam"
samtools index $bc"_mapped_sorted.bam"

htseq-count --nonunique=all -f bam $bc"_mapped_sorted.bam" /home/smaguire/work/unblock_remakes/data/cancer_samples/genome/human_miRbase22.gtf > ${file_name}"_miRNA_counts.txt"
htseq-count --nonunique=all -f bam $bc"_mapped_sorted.bam" /home/smaguire/work/sRNA_circ/spades/genome_ref/GRCh38_tRNA.gtf > ${file_name}"_tRNA_counts.txt"
