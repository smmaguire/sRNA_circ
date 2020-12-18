#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 4
#$ -m e


source activate bbtools

cd ${current_dir}

file=$( ls *"final_trimmed_output.fasta")
bc=${file%"final_trimmed_output.fasta"}

# echo "#################################################"
# echo "Map Human"
# echo "#################################################"

# ##### Trim short ##########
# bbduk.sh in=$file out=$bc"_final_trimmed_filtered.fasta" minlen=16 overwrite=T

# # bbmap.sh in=$bc"_final_trimmed_filtered.fasta" \
# # ref=/home/smaguire/work/human_miRNA/star_pipeline/genome/GRCh38.fasta \
# # path=/home/smaguire/work/sRNA_circ/spades/genome_ref \
# # k=12 \
# # local \
# # minid=0.9 \
# # out=$bc"_mapped.sam" \
# # overwrite=T

# source activate shortstack
# index=/home/smaguire/work/unblock_remakes/data/cancer_samples/genome/igenome/Homo_sapiens/NCBI/GRCh38/Sequence/BowtieIndex/genome
# bowtie -v 3 -l 10 -m 100 -k 1 -f --sam --best --strata $index $bc"_final_trimmed_filtered.fasta" $bc"_mapped.sam"

# echo "#################################################"
# echo "Count miRNAs"
# echo "#################################################"

 source activate htseq
# samtools view -b $bc"_mapped.sam" | samtools sort - > $bc"_mapped_sorted.bam"
# samtools index $bc"_mapped_sorted.bam"

# htseq-count --nonunique=all -f bam -a 0 $bc"_mapped_sorted.bam" /home/smaguire/work/unblock_remakes/data/cancer_samples/genome/human_miRbase22.gtf > $bc"_miRNA_counts.txt"
# htseq-count --nonunique=all -f bam -a 0 $bc"_mapped_sorted.bam" /home/smaguire/work/sRNA_circ/spades/genome_ref/GRCh38_tRNA.gtf > $bc"_tRNA_counts.txt"
htseq-count --nonunique=all -f bam -a 0 $bc"_mapped_sorted.bam" /home/smaguire/work/sRNA_circ/spades/genome_ref/gencode.v33.annotation.gtf > $bc"_all_non_coding_counts.txt"

