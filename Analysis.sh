##Objectives: To analyse the RNA seq data of 3 replicates of N2_eft-2 and wild type N2_control

## Step 1: Set up working directory, create and activate conda environment
mkdir 1_RNASeq
cd 1_RNASeq
conda create -n RNAseq
conda activate RNAseq

##Step 2: Install sra-tools
conda install bioconda::sra-tools
##Installation failed, downloadeded directly from github, moved it to the working directory and extract
mv /Users/geokf/Downloads/sratoolkit.3.1.1-mac-x86_64.tar.gz .
tar -xvzf sratoolkit.3.1.1-mac-x86_64.tar.gz

##Step 3: Download one sample using FASTQ-DUMP
./sratoolkit.3.1.1-mac-x86_64/bin/fastq-dump --split-files SRR30171833

##Step 4: Download ref genome from igenomes websites
wget https://s3.amazonaws.com/igenomes.illumina.com/Caenorhabditis_elegans/Ensembl/WBcel235/Caenorhabditis_elegans_Ensembl_WBcel235.tar.gz
tar -xvzf Caenorhabditis_elegans_Ensembl_WBcel235.tar.gz

##Step 5: Install fastqc
conda install bioconda::fastqc

# Step 6: run FASTQC
fastqc -t 8 SRR30171833_1.fastq SRR30171833_2.fastq

##Step 7: Install trimmomatic
conda install bioconda::trimmomatic

##Step 8: Trim raw reads using Trimmmomatic
trimmomatic PE -threads 8 SRR30171833_1.fastq SRR30171833_2.fastq Paired_trim_1.fastq Unpar_1.fastq Paired_trim_2.fastq Unpar_2.fastq HEADCROP:10

## Step 9: Install Star aligner
conda install bioconda::star
### fail to download, download from github - ver 2.7.11b
wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz
tar -xzf 2.7.11b.tar.gz
cd STAR-2.7.11be
### fail to execute ver 2.7.11b, download older ver from github - ver 2.7.10b
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR-2.7.10b

##Step 10: Generate star index
./STAR_2.7.10b/MacOSX_x86_64/STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir STAR_index_149 \
--genomeFastaFiles /Users/geokf/1_RNASeq/Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
--sjdbGTFfile /Users/geokf/1_RNASeq/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
--sjdbOverhang 149

##Step 11: Aligned trimmed reads to ref genome
./STAR_2.7.10b/MacOSX_x86_64/STAR --runThreadN 6 \
--genomeDir STAR_index_149 \
--readFilesIn Paired_trim_1.fastq Paired_trim_2.fastq \
--outFileNamePrefix results/test1 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts

##Step 12: rename results dir to SRR30171833
mv results SRR30171833

##Step 13: Clean up unnecessary file
rm *.fastq

##Step 14: write a loop for the remaining 5 files (loop the process of data download, fastqc, trimming, fastqc after trimming, alignment and remove fastqc)
for srr in SRR30171818 SRR30171834 SRR30171833 SRR30171824; do 
    echo "Start to download ${srr}";
    ./sratoolkit.3.1.1-mac-x86_64/bin/fastq-dump --split-files $srr
    echo "Starting FastQC ${srr}";
    fastqc -t 8 ${srr}_1.fastq ${srr}_2.fastq
    echo "Trimming ${srr}";
    trimmomatic PE -threads 8 ${srr}_1.fastq ${srr}_2.fastq ${srr}_Paired_trim_1.fastq ${srr}_Unpar_1.fastq ${srr}_Paired_trim_2.fastq ${srr}_Unpar_2.fastq HEADCROP:10
    echo "Starting FastQC ${srr}";
    fastqc -t 8 ${srr}_Paired_trim_1.fastq ${srr}_Paired_trim_2.fastq
    echo "Starting Alignment ${srr}";
    ./STAR_2.7.10b/MacOSX_x86_64/STAR --runThreadN 6 \
        --genomeDir STAR_index_149 \
        --readFilesIn ${srr}_Paired_trim_1.fastq ${srr}_Paired_trim_2.fastq \
        --outFileNamePrefix results/${srr}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
    echo "Removing ${srr}";
    rm *.fastq 
done

##SRR30171820 RNA-seq of N2_eft-2_1
##SRR30171818 RNA-seq of N2_eft-2_3
##SRR30171834 RNA-seq of N2_control_1
##SRR30171833 RNA-seq of N2_control_2
##SRR30171824 RNA-seq of N2_control_3

    ./STAR_2.7.10b/MacOSX_x86_64/STAR --runThreadN 6 \
        --genomeDir STAR_index_149 \
        --readFilesIn ${srr}_Paired_trim_1.fastq ${srr}_Paired_trim_2.fastq \
        --outFileNamePrefix results/${srr}_ \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts
    echo "Removing ${srr}";
    rm *.fastq 