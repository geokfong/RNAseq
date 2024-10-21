#!/bin/bash

# Part 1: Installation 
mamba install -y bioconda::sra-tools bioconda::fastqc bioconda::trimmomatic bioconda::star || {
    echo "Installation failed for STAR, downloading directly from GitHub..."
    mv /Users/geokf/Downloads/sratoolkit.3.1.1-mac-x86_64.tar.gz .
    tar -xvzf sratoolkit.3.1.1-mac-x86_64.tar.gz
}

# Part 2: Download reference genome
echo "Downloading reference genome..."
wget https://s3.amazonaws.com/igenomes.illumina.com/Caenorhabditis_elegans/Ensembl/WBcel235/Caenorhabditis_elegans_Ensembl_WBcel235.tar.gz
tar -xvzf Caenorhabditis_elegans_Ensembl_WBcel235.tar.gz

# Part 3: Generate STAR index
echo "Generating STAR index..."
STAR_DIR="./STAR-2.7.10b/MacOSX_x86_64"
GENOME_DIR="STAR_index_149"
mkdir -p "$GENOME_DIR"

"$STAR_DIR"/STAR --runThreadN 6 \
    --runMode genomeGenerate \
    --genomeDir "$GENOME_DIR" \
    --genomeFastaFiles ./Caenorhabditis_elegans/Ensembl/WBcel235/Sequence/WholeGenomeFasta/genome.fa \
    --sjdbGTFfile ./Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Genes/genes.gtf \
    --sjdbOverhang 149

# Part 4: Loop through SRR samples
srr_samples=(SRR30171820 SRR30171818 SRR30171824 SRR30171833 SRR30171834)

for srr in "${srr_samples[@]}"; do 
    echo "Processing sample: ${srr}"

    # Step 4.1: Download the sample using fastq-dump
    ./sra-tools/bin/fastq-dump --split-files "$srr"

    # Step 4.2: Run FastQC on raw reads
    echo "Running FastQC on raw reads for ${srr}..."
    fastqc -t 8 "${srr}_1.fastq" "${srr}_2.fastq"

    # Step 4.3: Trim raw reads using Trimmomatic
    echo "Trimming reads for ${srr}..."
    trimmomatic PE -threads 8 "${srr}_1.fastq" "${srr}_2.fastq" \
        "${srr}_Paired_trim_1.fastq" "${srr}_Unpar_1.fastq" \
        "${srr}_Paired_trim_2.fastq" "${srr}_Unpar_2.fastq" HEADCROP:10

    # Step 4.4: Run FastQC on trimmed reads
    echo "Running FastQC on trimmed reads for ${srr}..."
    fastqc -t 8 "${srr}_Paired_trim_1.fastq" "${srr}_Paired_trim_2.fastq"

    # Step 4.5: Align trimmed reads to reference genome
    echo "Aligning reads for ${srr}..."
    "$STAR_DIR"/STAR --runThreadN 6 \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "${srr}_Paired_trim_1.fastq" "${srr}_Paired_trim_2.fastq" \
        --outFileNamePrefix results/"${srr}_" \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts

    # Step 4.6: Clean up fastq files
    echo "Removing fastq files for ${srr}..."
    rm "${srr}"_*.fastq
done

echo "RNA-seq analysis completed!"
