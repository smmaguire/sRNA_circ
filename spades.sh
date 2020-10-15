#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 16
#$ -m e

barcode=${barcode}
work_dir="/home/smaguire/work/sRNA_circ/spades/"

mkdir -p $work_dir"human_brain_output/"$barcode
new_data_dir=$work_dir"human_brain_output/"$barcode
mkdir -p $new_data_dir/${name}
mkdir -p $new_data_dir/${name}_output
cd $new_data_dir/${name}

######### Filter Fastq and change to fasta for blast
source activate detect_periods
NanoFilt -l 1000 -q 7 --headcrop 100 ${reads} | reformat.sh in=stdin.fastq out=${name}.fasta overwrite=T --ignorebadquality qin=64

######### Spades
source activate spades

python /home/smaguire/work/sRNA_circ/spades/SPADE/SPADE.py \
-in ${name}.fasta \
-Nw 500 \
-v N \
-n 16 \
-Nk 5 \
-Ns 10 \
-Ng 200 \
-Nm 200 \
-Nu 0.5 \
--delete

######## Move all the output files into a new folder
count=0
find . -type f -name 'repeat.gbk' |
    while IFS= read file_name; do
        mv "$file_name" "$new_data_dir/${name}_output/${file_name##*\/}_${count}"
        ((count=count+1))
    done

####### Delete everything
rm -r *
