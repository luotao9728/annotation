#!/bin/bash

# see files in directory
ls -l

# get input files
echo "*******************************************************************************************"
echo "Please input the reference genome:"
read reference_genome

# check if the files exist, if not echo not exist and exit
if [ ! -f "$reference_genome" ]; then
    echo "$reference_genome does not exist"
    exit 1
fi

echo "Please input the forward short reads file:"
read forward_short_reads
if [ ! -f "$forward_short_reads" ]; then
    echo "$forward_short_reads does not exist"
    exit 1
fi

echo "Please input the reverse short reads file:"
read reverse_short_reads
if [ ! -f "$reverse_short_reads" ]; then
    echo "$reverse_short_reads does not exist"
    exit 1
fi

echo "Please input the long reads file:"
read long_reads
if [ ! -f "$long_reads" ]; then
    echo "$long_reads does not exist"
    exit 1
fi

echo "Please input the reference annotation file:"
read reference_annotation
if [ ! -f "$reference_annotation" ]; then
    echo "$reference_annotation does not exist"
    exit 1
fi

# ask the user to input a keyword for the output file name
echo "Please enter a keyword (e.g. liver or heart) for the output file name:"
read keyword

cwd=$(pwd)

sbatch build_index.sh $reference_genome $keyword $cwd
sbatch pipeline.sh $reference_genome $forward_short_reads $reverse_short_reads $long_reads $reference_annotation $keyword $cwd