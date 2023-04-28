#!/bin/bash
#SBATCH -p RM
#SBATCH -t 2:00:00
#SBATCH --ntasks-per-node=64

reference_genome=$1
keyword=$2

# activate local batch environment
source ~/.bashrc
echo "Bash environment activated!"

# activate conda environment for the pipeline
conda activate annotation
echo "Conda environment activated!"

# change direction to data location
cd $3

# start alignment
minimap2 -a $reference_genome trim_split_"$keyword"_pacbio.fasta > "$keyword"_aligned_pacbio.sam
samtools sort "$keyword"_aligned_pacbio.sam -o "$keyword"_aligned_pacbio.bam

# deactivate conda environment
conda deactivate
