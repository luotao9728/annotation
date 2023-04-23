#!/bin/bash
#SBATCH -p RM
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=64

# store inputs
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
mkdir alignments
cd alignments
mkdir "$keyword"_index
cd "$keyword"_index

# build index
echo "*******************************************************************************************"
echo "Start building index from reference genome!"

hisat2-build -p 64 ../../"$reference_genome" $keyword

echo "Start building index from reference genome!"
echo "*******************************************************************************************"

# deactivate conda environment
conda deactivate
