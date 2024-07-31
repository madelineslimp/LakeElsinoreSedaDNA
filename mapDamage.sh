#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <reference_fasta> <fastq_R1> <fastq_R2>"
    exit 1
fi

# Assign input parameters to variables
reference_fasta=$1
fastq_R1=$2
fastq_R2=$3

# Debugging output: print the input files
echo "Reference FASTA: $reference_fasta"
echo "FASTQ R1: $fastq_R1"
echo "FASTQ R2: $fastq_R2"

# Activate the virtual environment (if needed)
#source /scratch1/RachelMeyer/MAPDAMAGE_WORK/mapdamage_venv/bin/activate

# FastQ trimming
fastp -l 30 -q 25 --low_complexity_filter --detect_adapter_for_pe -i "${fastq_R1}" -I "${fastq_R2}" -o "${fastq_R1%.fastq.gz}.trimmed.fastq.gz" -O "${fastq_R2%.fastq.gz}.trimmed.fastq.gz"

# Check if fastp succeeded
if [ $? -ne 0 ]; then
    echo "fastp failed"
    exit 1
fi

# Prepare the reference file (indexing)
ref_base=$(basename "$reference_fasta" .fasta)
bowtie2-build "$reference_fasta" "${ref_base}_index"

# Map to reference
fastq_R1t="${fastq_R1%.fastq.gz}.trimmed.fastq.gz"
fastq_R2t="${fastq_R2%.fastq.gz}.trimmed.fastq.gz"
sam_output="${ref_base}_output.sam"
bam_output="${ref_base}_output.bam"
sorted_bam="${ref_base}_output.sorted.bam"
rg_bam="${ref_base}_output.sorted.rg.bam"

# Check if trimmed files exist
if [ ! -f "$fastq_R1t" ]; then
    echo "Trimmed file $fastq_R1t not found"
    exit 1
fi

if [ ! -f "$fastq_R2t" ]; then
    echo "Trimmed file $fastq_R2t not found"
    exit 1
fi

# Mapping step
bowtie2 -x "${ref_base}_index" -1 "$fastq_R1t" -2 "$fastq_R2t" -S "$sam_output"
if [ $? -ne 0 ]; then
    echo "bowtie2 failed"
    exit 1
fi

# Convert SAM to BAM
samtools view -bS "$sam_output" > "$bam_output"
if [ $? -ne 0 ]; then
    echo "samtools view failed"
    exit 1
fi

# Sort BAM file
samtools sort "$bam_output" -o "$sorted_bam"
if [ $? -ne 0 ]; then
    echo "samtools sort failed"
    exit 1
fi

# Index BAM file
samtools index "$sorted_bam"
if [ $? -ne 0 ]; then
    echo "samtools index failed"
    exit 1
fi

# Flagstat BAM file
samtools flagstat "$sorted_bam"
if [ $? -ne 0 ]; then
    echo "samtools flagstat failed"
    exit 1
fi

# Add read groups using picard from conda environment
picard AddOrReplaceReadGroups I="$sorted_bam" O="$rg_bam" RGID=1 RGLB=library1 RGPL=illumina RGPU=unit1 RGSM="${ref_base}"
if [ $? -ne 0 ]; then
    echo "Picard AddOrReplaceReadGroups failed"
    exit 1
fi

# Run mapDamage
mapDamage -i "$rg_bam" -r "$reference_fasta"
if [ $? -ne 0 ]; then
    echo "mapDamage failed"
    exit 1
fi
