#!/bin/bash

# Define paths
barcode_folder="barcode13"
aligned_folder="aligned_bams"
ref_genome="combined_ref.fa"

# Create aligned_bams folder if it doesn't exist
mkdir -p $aligned_folder

# Loop through all BAM files in the barcode13 folder
for bam_file in $barcode_folder/*.bam; do
    # Extract the base name of the file (without path and extension)
    base_name=$(basename $bam_file .bam)
    
    # Align BAM file to the combined reference genome
    bwa mem $ref_genome $bam_file > $aligned_folder/${base_name}_with_header.sam

    # Convert SAM to BAM, sort, and index
    samtools view -Sb $aligned_folder/${base_name}_with_header.sam | samtools sort -o $aligned_folder/${base_name}_aligned.bam
    samtools index $aligned_folder/${base_name}_aligned.bam

    echo "Processed $bam_file"
done

echo "All BAM files have been aligned, sorted, and indexed."

