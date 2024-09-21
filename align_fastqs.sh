#!/bin/bash

# Define paths
fastq_folder="fastq_files"      
aligned_folder="aligned_bams"
ref_genome="combined_ref.fa"
threads=4                       # Set the number of threads to use

# Create aligned_bams folder if it doesn't exist
mkdir -p $aligned_folder

# Loop through all FASTQ files in the folder
for fastq_file in $fastq_folder/*.fastq; do
    # Extract the base name of the file (without path and extension)
    base_name=$(basename $fastq_file .fastq)

    # Check if the BAM file already exists
    if [ -f $aligned_folder/${base_name}_with_header.bam ]; then
        echo "Skipping $fastq_file, BAM file already exists"
        continue
    fi

    # Check if FASTQ file is corrupt by checking the first and last lines
    if ! head -n 1 "$fastq_file" >/dev/null || ! tail -n 1 "$fastq_file" >/dev/null; then
        echo "Skipping $fastq_file, it appears to be corrupt"
        continue
    fi
    
    # Align FASTQ file to the combined reference genome using specified threads
    bwa mem -t $threads $ref_genome $fastq_file > $aligned_folder/${base_name}_with_header.sam
    if [ $? -ne 0 ]; then
        echo "Error aligning $fastq_file"
        continue
    fi

    # Convert SAM to BAM, sort, and index using specified threads
    samtools view -@ $threads -Sb $aligned_folder/${base_name}_with_header.sam | \
    samtools sort -@ $threads -o $aligned_folder/${base_name}_with_header.bam
    if [ $? -ne 0 ]; then
        echo "Error processing BAM for $fastq_file"
        continue
    fi
    samtools index $aligned_folder/${base_name}_with_header.bam
    if [ $? -ne 0 ]; then
        echo "Error indexing BAM for $fastq_file"
        continue
    fi

    # Optional: Remove the SAM file to save space
    rm $aligned_folder/${base_name}_with_header.sam

    echo "Processed $fastq_file"
done

echo "All FASTQ files have been aligned, sorted, and indexed."

