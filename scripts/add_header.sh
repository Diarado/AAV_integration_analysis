#!/bin/bash

# Define paths
BAM_DIR="barcode13"
HEADER_BAM_DIR="barcode13_header"
HG38_REF="hg38.fa"
AAV_REF="pssAAV-CB-EGFP ARM.fa"

# Create output directory if not exists
mkdir -p "$HEADER_BAM_DIR"

# Ensure both references are indexed
if [ ! -f "${HG38_REF}.fai" ]; then
    echo "Indexing human genome reference (hg38)..."
    samtools faidx "$HG38_REF"
fi

if [ ! -f "${AAV_REF}.fai" ]; then
    echo "Indexing AAV reference..."
    samtools faidx "$AAV_REF"
fi

# Generate combined header from hg38 and AAV references
echo "Generating combined header from hg38 and AAV..."
samtools view -H "$HG38_REF" > combined_header.sam
samtools view -H "$AAV_REF" >> combined_header.sam

# Process each BAM file in the barcode13 folder
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extract the sample name
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    echo "Processing $SAMPLE_NAME..."
    
    # Reheader the BAM file and save to the new directory
    HEADER_BAM_FILE="$HEADER_BAM_DIR/${SAMPLE_NAME}_with_header.bam"
    echo "Adding header to BAM file and saving to: $HEADER_BAM_FILE"
    samtools reheader combined_header.sam "$BAM_FILE" > "$HEADER_BAM_FILE"
    
    # Check if the reheadering succeeded
    if [ $? -ne 0 ]; then
        echo "Error reheadering $SAMPLE_NAME. Skipping."
        continue
    fi
    
    # Index the reheadered BAM file
    echo "Indexing $HEADER_BAM_FILE..."
    samtools index "$HEADER_BAM_FILE"
    
    # Check if indexing succeeded
    if [ $? -ne 0 ]; then
        echo "Error indexing $HEADER_BAM_FILE. Skipping."
        continue
    fi
    
    echo "Finished processing $SAMPLE_NAME."
done

# Clean up
rm combined_header.sam

echo "All BAM files processed successfully!"

