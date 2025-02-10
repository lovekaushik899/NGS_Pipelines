#!/bin/bash

# Load Conda
source /home/user/anaconda3/etc/profile.d/conda.sh

# Set directories
RAW_READS_DIR="./"
TRIMMED_READS_DIR="./trimmed_reads"
QC_RESULTS_DIR="./qc_results"
ALIGNMENT_DIR="./alignment_results"
LOG_FILE="pipeline.log"
READS_FILE="reads.txt"

# Create necessary directories
mkdir -p "$TRIMMED_READS_DIR" "$QC_RESULTS_DIR" "$ALIGNMENT_DIR"

# Logging function
log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log "Starting NGS pipeline..."

# Step 1: Quality Check of Raw Reads
log "Running FastQC on raw reads..."
conda activate base || { log "Failed to activate base environment"; exit 1; }
while read -r SAMPLE; do
    fastqc -t 20 -o "$QC_RESULTS_DIR" "$RAW_READS_DIR/${SAMPLE}_1.fastq" "$RAW_READS_DIR/${SAMPLE}_2.fastq" || log "FastQC failed for $SAMPLE"
done < "$READS_FILE"
multiqc -o "$QC_RESULTS_DIR" "$QC_RESULTS_DIR" || log "MultiQC failed"

# Step 2: Trimming using Trim Galore
log "Running Trim Galore..."
conda activate trim_env_env || { log "Failed to activate trim_env_env"; exit 1; }
while read -r SAMPLE; do
    trim_galore --paired --cores 20 -o "$TRIMMED_READS_DIR" "$RAW_READS_DIR/${SAMPLE}_1.fastq" "$RAW_READS_DIR/${SAMPLE}_2.fastq" || { log "Trim Galore failed for $SAMPLE"; continue; }

    # Check if trimmed files were generated
    if [[ ! -f "$TRIMMED_READS_DIR/${SAMPLE}_1_val_1.fq" || ! -f "$TRIMMED_READS_DIR/${SAMPLE}_2_val_2.fq" ]]; then
        log "Trim Galore output missing for $SAMPLE, skipping..."
        continue
    fi
done < "$READS_FILE"

# Step 3: Quality Check on Trimmed Reads
log "Running FastQC on trimmed reads..."
while read -r SAMPLE; do
    fastqc -t 8 -o "$QC_RESULTS_DIR" "$TRIMMED_READS_DIR/${SAMPLE}_1_val_1.fq" "$TRIMMED_READS_DIR/${SAMPLE}_2_val_2.fq" || log "FastQC on trimmed reads failed for $SAMPLE"
done < "$READS_FILE"
multiqc -o "$QC_RESULTS_DIR" "$QC_RESULTS_DIR" || log "MultiQC on trimmed reads failed"

# Step 4: Aligning Reads using HiSAT2
log "Aligning reads with HiSAT2..."
conda activate hisat2_env || { log "Failed to activate hisat2_env"; exit 1; }
while read -r SAMPLE; do
    FILE_1="$TRIMMED_READS_DIR/${SAMPLE}_1_val_1.fq"
    FILE_2="$TRIMMED_READS_DIR/${SAMPLE}_2_val_2.fq"

    if [[ ! -f "$FILE_1" || ! -f "$FILE_2" ]]; then
        log "Skipping $SAMPLE: one or both trimmed read files are missing."
        continue
    fi

    hisat2 -p 20 -x genome_index -1 "$FILE_1" -2 "$FILE_2" -S "$ALIGNMENT_DIR/${SAMPLE}.sam" || log "HiSAT2 failed for $SAMPLE"
done < "$READS_FILE"

# Step 5: Convert SAM to BAM
log "Converting SAM to BAM..."
while read -r SAMPLE; do
    samtools view -@ 20 -bS "$ALIGNMENT_DIR/${SAMPLE}.sam" > "$ALIGNMENT_DIR/${SAMPLE}.bam" || log "SAM to BAM conversion failed for $SAMPLE"
    rm -f "$ALIGNMENT_DIR/${SAMPLE}.sam"
done < "$READS_FILE"

# Step 6: Sort BAM Files
log "Sorting BAM files..."
while read -r SAMPLE; do
    samtools sort -@ 20 -o "$ALIGNMENT_DIR/${SAMPLE}_sorted.bam" "$ALIGNMENT_DIR/${SAMPLE}.bam" || log "Sorting failed for $SAMPLE"
    rm -f "$ALIGNMENT_DIR/${SAMPLE}.bam"
done < "$READS_FILE"

# Step 7: Index BAM Files
log "Indexing BAM files..."
while read -r SAMPLE; do
    samtools index "$ALIGNMENT_DIR/${SAMPLE}_sorted.bam" || log "Indexing failed for $SAMPLE"
done < "$READS_FILE"

log "Pipeline completed successfully!"

