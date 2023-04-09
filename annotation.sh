#!/bin/bash

read ref_genome ill_for ill_rev pacbio ref_annotation

echo The file name of reference genome is: $ref_genome
echo The file name of illumina RNA-seq (forward) is: $ill_for
echo The file name of illumina RNA-seq (reverse) is: $ill_rev
echo The file name of PacBio RNA-seq is: $pacbio
echo The file name of reference annotation is: $ref_annotation

echo Start annotaion pipeline!

echo Start trimming illumina RNA-seq!

# trimmomatic PE -phred33 $ill_for $ill_rev output_forward_paired.fastq output_forward_unpaired.fastq output_reverse_paired.fastq output_reverse_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
