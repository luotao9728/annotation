# Project: Complete a functional annotation of the horse genome
## 03-713 Bioinformatics Data Integration Practicum
## Team-2: Taylor Ayers, Tao Luo, Lilin Huang, Sarah Oladejo


### Prepare the environment
1. Make sure Anaconda3 is installed in your computer
2. Make sure in your current working environment has the following packages:
* Trimmomatic
* LoRDEC
* Magic-BLAST
* SamTool
* StringTie
* Salmon
3. Alternatively, you could directly create a new working conda environment using the following command 
(make sure you have annotation.yml file in your working directory):
> conda env create -n annotation --file annotation.yml

### Requirements for input files
1. Reference genome: fasta/fna
2. illumina RNA-seq (forward/reverse): fastq
3. PacBio RNA-seq: fastq
4. Reference annotation: gtf

### Instruction for the pipeline
1. Make sure you have the environment (with all packages) ready.
2. Make sure you have all the input files ready.
3. Execute the command and follow the prompt:
> ./annotation.sh
4. Be patient. The annotation process may take a long time. Have a great day! :)
