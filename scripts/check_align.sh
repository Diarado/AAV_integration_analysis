#!/bin/bash

# Define your paths and files
AAV_REF="pssAAV-CB-EGFP ARM.fa"
BAM_DIR="barcode13"
ALIGN_DIR="aligned_bams"
OUTPUT_CSV="mapped_samples.csv"
mkdir -p "$ALIGN_DIR"

# Initialize the CSV file
echo "Sample_Name, Mapped_Reads" > "$OUTPUT_CSV"

# Ensure AAV reference is indexed
if [ ! -f "$AAV_REF.bwt" ]; then
    echo "Indexing AAV reference..."
    bwa index "$AAV_REF"
    if [ $? -ne 0 ]; then
        echo "Error indexing AAV reference. Exiting."
        exit 1
    fi
fi

for BAM_FILE in "$BAM_DIR"/*.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    echo "Processing $SAMPLE_NAME..."

    # Convert BAM to FASTQ
    samtools fastq "$BAM_FILE" > "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
    if [ $? -ne 0 ]; then
        echo "Error converting $SAMPLE_NAME BAM to FASTQ. Skipping."
        continue
    fi

    # Align the FASTQ file to the AAV reference
    bwa mem "$AAV_REF" "$ALIGN_DIR/${SAMPLE_NAME}.fastq" > "$ALIGN_DIR/${SAMPLE_NAME}_aligned.sam"
    if [ $? -ne 0 ]; then
        echo "Error aligning $SAMPLE_NAME. Skipping."
        continue
    fi

    # Convert SAM to BAM
    samtools view -S -b "$ALIGN_DIR/${SAMPLE_NAME}_aligned.sam" > "$ALIGN_DIR/${SAMPLE_NAME}_aligned.bam"
    if [ $? -ne 0 ]; then
        echo "Error converting SAM to BAM for $SAMPLE_NAME. Skipping."
        continue
    fi

    # Sort and index the BAM file
    samtools sort "$ALIGN_DIR/${SAMPLE_NAME}_aligned.bam" -o "$ALIGN_DIR/${SAMPLE_NAME}_sorted.bam"
    if [ $? -ne 0 ]; then
        echo "Error sorting BAM for $SAMPLE_NAME. Skipping."
        continue
    fi

    samtools index "$ALIGN_DIR/${SAMPLE_NAME}_sorted.bam"

    # Check alignment statistics
    samtools flagstat "$ALIGN_DIR/${SAMPLE_NAME}_sorted.bam" > "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt"
    if [ $? -ne 0 ]; then
        echo "Error generating flagstat for $SAMPLE_NAME. Skipping."
        continue
    fi

    # Get the number of mapped reads from the flagstat output
    MAPPED=$(grep "mapped (" "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" | head -1 | cut -d ' ' -f 1)

    # Ensure we correctly check if MAPPED is greater than 0 and save to CSV
    if [[ "$MAPPED" -gt 0 ]]; then
        echo "$SAMPLE_NAME has mapped reads. Proceeding with analysis."
        echo "$SAMPLE_NAME, $MAPPED" >> "$OUTPUT_CSV"
    else
        echo "$SAMPLE_NAME has no mapped reads. Skipping."
    fi

done
