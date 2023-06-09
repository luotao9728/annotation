#!/bin/bash
#SBATCH -p RM
#SBATCH -t 10:00:00
#SBATCH --ntasks-per-node=64

# store inputs
reference_genome=$1
forward_short_reads=$2
reverse_short_reads=$3
long_reads=$4
reference_annotation=$5
keyword=$6

# activate local batch environment
source ~/.bashrc
echo "Bash environment activated!"

# activate conda environment for the pipeline
conda activate annotation
echo "Conda environment activated!"

# change direction to data location
cd $7

# sickle
echo "*******************************************************************************************"
echo "Start trimming illumina short RNA-seq!"

sickle pe -f $forward_short_reads -r $reverse_short_reads -t sanger -o trimmed_"$keyword"_illumina_R1.fastq -p trimmed_"$keyword"_illumina_R2.fastq -s temp.fastq
seqtk seq -a trimmed_"$keyword"_illumina_R1.fastq > trimmed_"$keyword"_illumina_R1.fa
seqtk seq -a trimmed_"$keyword"_illumina_R2.fastq > trimmed_"$keyword"_illumina_R2.fa
rm temp.fastq

echo "Finish trimming illumina short RNA-seq!"
echo "*******************************************************************************************"

# LorDEC
## lordec-correct
echo "*******************************************************************************************"
echo "Start correcting PacBio long RNA-seq!"

lordec-correct -2 trimmed_"$keyword"_illumina_R1.fastq trimmed_"$keyword"_illumina_R2.fastq -k 19 -s 3 -i $long_reads -o "$keyword"_pacbio_lordec_corrected.fasta

echo "Finish correcting PacBio long RNA-seq!"
echo "*******************************************************************************************"

## lordec-trim-split
echo "*******************************************************************************************"
echo "Start trim-splitting corrected PacBio long RNA-seq!"

lordec-trim-split -i "$keyword"_pacbio_lordec_corrected.fasta -o trim_split_"$keyword"_pacbio.fasta

echo "Finish trim-splitting corrected PacBio long RNA-seq!"
echo "*******************************************************************************************"

# alignment
echo "*******************************************************************************************"
echo "Start alignment"

minimap2 -a $reference_genome trim_split_"$keyword"_pacbio.fasta > alignments/"$keyword"_aligned_pacbio.sam
cd alignments
hisat2 -f -x "$keyword"_index/"$keyword" -1 ../trimmed_"$keyword"_illumina_R1.fa -2 ../trimmed_"$keyword"_illumina_R2.fa -S "$keyword"_aligned_illumina.sam

echo "Finish alignment!"
echo "*******************************************************************************************"

# SamTools
echo "*******************************************************************************************"
echo "Start converting SAM files to BAM files!"

samtools sort "$keyword"_aligned_illumina.sam -o "$keyword"_aligned_illumina.bam
samtools sort "$keyword"_aligned_pacbio.sam -o "$keyword"_aligned_pacbio.bam

echo "Finish converting SAM files to BAM files!"
echo "*******************************************************************************************"

# Stringtie
echo "*******************************************************************************************"
echo "Start annotaion!"

stringtie --mix -G ../"$reference_annotation" -o ../"$keyword"_annotation.gtf "$keyword"_aligned_illumina.bam "$keyword"_aligned_pacbio.bam

echo "Finish annotation!"
echo "*******************************************************************************************"

# deactivate conda environment
conda deactivate
