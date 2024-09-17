#!/bin/bash

# Directory paths
ZIP_DIR="/mnt/d/downloads"  # Directory containing the zip files
EXTRACT_DIR="/mnt/d/downloads"  # Temporary extraction directory
ALIGNED_DIR="/mnt/d/Jiahe/IU/AAV/HeLa_project/fastq_files"  # Destination directory

# Step 1: Unzip all zip files
for file in "$ZIP_DIR"/*.zip; do
    unzip "$file" -d "$EXTRACT_DIR"
done

# Step 2: Extract all tar.gz (or other compressed files) in extracted directory
for archive in "$EXTRACT_DIR"/*.tar.gz; do
    tar -xvf "$archive" -C "$EXTRACT_DIR"
done

# Step 3: Move everything to aligned_bams directory
mv "$EXTRACT_DIR"/* "$ALIGNED_DIR"

echo "All files unzipped, extracted, and moved to aligned_bams!"

