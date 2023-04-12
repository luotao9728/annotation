#!/bin/bash

# Take inputs
echo "**************************************************************"
echo "Please enter the file name with the following prompts:"
read -p "Reference genome (fasta/fna): " ref_genome
read -p "illumina RNA-seq forward (fastq): " ill_for
read -p "illumina RNA-seq reverse (fastq): " ill_rev
read -p "PacBio RNA-seq (fastq): " pacbio
read -p "Reference annotation (gtf): " ref_annotation
echo "**************************************************************"

echo "The file name of reference genome is: $ref_genome"
echo "The file name of illumina RNA-seq (forward) is: $ill_for"
echo "The file name of illumina RNA-seq (reverse) is: $ill_rev"
echo "The file name of PacBio RNA-seq is: $pacbio"
echo "The file name of reference annotation is: $ref_annotation"

#SBATCH -N 2
#SBATCH -p RM
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node=128

# Initiate bash environment
source ~/.bashrc

# Activate conda environment
conda activate annotation

# Start pipeline
echo "Start annotaion pipeline!"

# Start trimming illumina RNA-seq
echo "Start trimming illumina RNA-seq!"
trimmomatic PE -phred33 $ill_for $ill_rev output_forward_paired.fastq output_forward_unpaired.fastq output_reverse_paired.fastq output_reverse_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
echo "Finish trimming illumina RNA-seq!"

# Start correcting PacBio RNA-seq
echo "Start trimming PacBio RNA-seq!"
lordec-correct -2 /ocean/projects/bio230007p/tluo1/data_processing/trimmed_liver_illumina_R2.fastq -k 19 -s 3 -i /ocean/projects/bio230007p/tluo1/data/liver_pacbio/ERR9764407.fastq -o /ocean/projects/bio230007p/tayers/liver_pacbio_lordec_corrected.fasta
echo "Finish trimming PacBio RNA-seq!"

conda deactivate
