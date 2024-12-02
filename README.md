# RNA-Seq Analysis
This repository contains scripts related to RNA-Seq data processing, including downloading data, quality control, trimming, and alignment to a reference genome.

## Installation
I used mamba to set up my environments, such as 
`mamba create -n RNASeq -c bioconda fastqc star sra-tools trimmomatic samtools bedtools`

## Steps
The script performs the following steps for each sample:
1. Download Sequencing Data: Uses **fastq-dump** to download the sample data.
2. Quality Control: Runs **FastQC** to assess the quality of the raw reads.
3. Trimming Reads: Utilizes **Trimmomatic** to trim raw reads and remove adapter sequences.
4. Quality Control on Trimmed Reads: Runs **FastQC** again on the trimmed reads.
5. Alignment: Aligns the trimmed reads to the reference genome using **STAR**.
6. Cleanup: Removes the original fastq files to save disk space.

## Notes
You may need to adjust paths in the script depending on your local setup.
