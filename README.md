# Project: Complete a functional annotation of the horse genome
## 03-713 Bioinformatics Data Integration Practicum
## Team-2: Taylor Ayers, Tao Luo, Lilin Huang, Sarah Oladejo

### Prepare the directory
> git clone https://github.com/luotao9728/annotation   

This directory contains files:  
start_pipeline.sh, build_index.sh, pipeline.sh, annotation.yml

### Prepare the environment
1. Make sure Anaconda3 is installed on your computer
2. Make sure in your current working environment has the following packages:
* Trimmomatic &emsp; (Trim illumina short reads)
* LoRDEC &emsp; (Fix long reads by short reads)
* hisat2       &emsp;   (Alignment)
* seqtk      &emsp;     (Convert FASTA and FASTQ format)
* SamTool     &emsp;    (Sort and Convert sam to bam)
* StringTie    &emsp;   (Annotation)
* Salmon       &emsp;   (Quantify expression)
3. Alternatively, you could directly create a new working conda environment using the following command 
(make sure you have annotation.yml file in your working directory):
> conda env create -n annotation --file annotation.yml
> 
> conda activate annotation

### Requirements for input files
1. Reference genome: fasta/fna
2. illumina RNA-seq (forward/reverse): fastq
3. PacBio RNA-seq: fastq
4. Reference annotation: gtf
5. Keyword: name of this annotation
The input files should be in the annotation directory

### Instruction for the pipeline
1. Make sure you have the environment (with all packages) ready.
2. Make sure you have all the input files ready.
3. Execute the command and follow the prompt:
> bash start_pipeline.sh
4. Follow the instructions to enter the file names.
5. Be patient. The annotation process may take a long time. Have a great day! :)
