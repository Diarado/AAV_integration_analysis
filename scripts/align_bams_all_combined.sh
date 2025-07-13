#!/bin/bash

# BAM to Combined Reference Alignment Script (Hybrid Approach)
# 
# This script processes BAM files containing sequencing reads and re-aligns them 
# to a combined reference genome (human hg38 + AAV vector). This combines the 
# best parts of the previous scripts for AAV integration analysis.
#
# Workflow:
# 1. Creates a combined reference (hg38 + AAV) if it doesn't exist
# 2. Indexes the combined reference for alignment
# 3. Converts each BAM file to FASTQ format (reusing Script 3 approach)
# 4. Re-aligns reads to the combined reference using BWA-MEM
# 5. Converts back to sorted, indexed BAM format
# 6. Outputs files suitable for downstream AAV integration site analysis
#
# Input: BAM files in 'barcode03/' directory
# Output: Re-aligned BAM files in 'aligned_bams_barcode03/' directory with '_with_header.bam' suffix

# Define paths
BAM_DIR="../data/barcode03"                    # Input BAM files directory
ALIGN_DIR="../data/aligned_bams_barcode03"               # Output directory
HG38_REF="../data/hg38/hg38.fa"                    # Human genome reference
AAV_REF="../data/barcode03_ref/37825-AAV2-CAG-EGFP.gb"  # AAV reference
COMBINED_REF="../data/barcode03_ref/combined_ref.fa"         # Combined reference
OUTPUT_CSV="../data/alignment_summary.csv"     # Summary statistics
threads=4                              # Number of threads to use

# Create output directory
mkdir -p "$ALIGN_DIR"

# Initialize summary CSV
echo "Sample_Name,Total_Reads,Mapped_Reads,Mapping_Rate" > "$OUTPUT_CSV"

# Check if combined reference exists, if not create it
if [ ! -f "$COMBINED_REF" ]; then
    echo "Combined reference not found. Creating combined reference..."
    
    # Check if hg38 reference exists
    if [ ! -f "$HG38_REF" ]; then
        echo "Error: hg38.fa not found at expected path: $HG38_REF"
        echo "Current working directory: $(pwd)"
        echo "Please check the path and ensure the file exists."
        exit 1
    fi
    
    # Check if AAV reference exists
    if [ ! -f "$AAV_REF" ]; then
        echo "Error: AAV reference not found at $AAV_REF"
        echo "Current working directory: $(pwd)"
        echo "Please check the path and ensure the file exists."
        exit 1
    fi
    
    # Create directory for combined reference if it doesn't exist
    mkdir -p "$(dirname "$COMBINED_REF")"
    
    # Combine the references
    echo "Combining hg38 and AAV references..."
    cat "$HG38_REF" "$AAV_REF" > "$COMBINED_REF"
    
    if [ $? -ne 0 ]; then
        echo "Error creating combined reference. Exiting."
        exit 1
    fi
    
    echo "Combined reference created successfully at: $COMBINED_REF"
fi

# Ensure combined reference is indexed
if [ ! -f "${COMBINED_REF}.bwt" ]; then
    echo "Indexing combined reference for BWA..."
    bwa index "$COMBINED_REF"
    if [ $? -ne 0 ]; then
        echo "Error indexing combined reference. Exiting."
        exit 1
    fi
fi

if [ ! -f "${COMBINED_REF}.fai" ]; then
    echo "Indexing combined reference for samtools..."
    samtools faidx "$COMBINED_REF"
fi

# Check if BAM directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "Error: BAM directory not found at $BAM_DIR"
    echo "Current working directory: $(pwd)"
    echo "Please check the path and ensure the directory exists."
    exit 1
fi

# Process each BAM file (adapted from Script 3 approach)
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Check if any BAM files exist
    if [ ! -f "$BAM_FILE" ]; then
        echo "No BAM files found in $BAM_DIR"
        exit 1
    fi
    
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    # Check if output already exists
    if [ -f "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam" ]; then
        echo "Skipping $SAMPLE_NAME, output BAM file already exists"
        continue
    fi
    
    echo "Processing $SAMPLE_NAME..."

    # Convert BAM to FASTQ (from Script 3)
    echo "Converting BAM to FASTQ..."
    samtools fastq "$BAM_FILE" > "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
    if [ $? -ne 0 ]; then
        echo "Error converting $SAMPLE_NAME BAM to FASTQ. Skipping."
        continue
    fi

    # Check if FASTQ is empty
    if [ ! -s "$ALIGN_DIR/${SAMPLE_NAME}.fastq" ]; then
        echo "Generated FASTQ file is empty for $SAMPLE_NAME. Skipping."
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
        continue
    fi

    # Align to combined reference (modified from Script 2)
    echo "Aligning to combined reference..."
    bwa mem -t $threads "$COMBINED_REF" "$ALIGN_DIR/${SAMPLE_NAME}.fastq" > "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
    if [ $? -ne 0 ]; then
        echo "Error aligning $SAMPLE_NAME. Skipping."
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
        continue
    fi

    # Convert SAM to BAM, sort and index (from Script 2)
    echo "Converting to BAM and sorting..."
    samtools view -@ $threads -Sb "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam" | \
    samtools sort -@ $threads -o "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam"
    if [ $? -ne 0 ]; then
        echo "Error processing BAM for $SAMPLE_NAME. Skipping."
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq" "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
        continue
    fi

    # Index the BAM file
    echo "Indexing BAM file..."
    samtools index "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam"
    if [ $? -ne 0 ]; then
        echo "Error indexing BAM for $SAMPLE_NAME. Skipping."
        continue
    fi

    # Generate alignment statistics (from Script 3)
    echo "Generating alignment statistics..."
    samtools flagstat "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam" > "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt"
    if [ $? -eq 0 ]; then
        # Extract statistics
        TOTAL_READS=$(grep "in total" "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" | cut -d ' ' -f 1)
        MAPPED_READS=$(grep "mapped (" "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" | head -1 | cut -d ' ' -f 1)
        
        # Calculate mapping rate
        if [[ "$TOTAL_READS" -gt 0 ]]; then
            MAPPING_RATE=$(echo "scale=2; $MAPPED_READS * 100 / $TOTAL_READS" | bc -l 2>/dev/null || echo "0")
        else
            MAPPING_RATE="0"
        fi
        
        # Add to summary CSV
        echo "$SAMPLE_NAME,$TOTAL_READS,$MAPPED_READS,${MAPPING_RATE}%" >> "$OUTPUT_CSV"
        
        echo "$SAMPLE_NAME: $MAPPED_READS/$TOTAL_READS mapped (${MAPPING_RATE}%)"
    fi

    # Clean up intermediate files
    rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq" "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
    
    echo "Successfully processed $SAMPLE_NAME"
done

echo "All BAM files processed successfully!"
echo "Summary statistics saved to: $OUTPUT_CSV"

# Check if combined reference exists, if not create it
if [ ! -f "$COMBINED_REF" ]; then
    echo "Combined reference not found. Creating combined reference..."
    
    # Check if hg38 reference exists
    if [ ! -f "$HG38_REF" ]; then
        echo "Error: hg38.fa not found. Please ensure it exists in the current directory."
        exit 1
    fi
    
    # Check if AAV reference exists
    if [ ! -f "$AAV_REF" ]; then
        echo "Error: AAV reference not found at $AAV_REF"
        exit 1
    fi
    
    # Combine the references
    echo "Combining hg38 and AAV references..."
    cat "$HG38_REF" "$AAV_REF" > "$COMBINED_REF"
    
    if [ $? -ne 0 ]; then
        echo "Error creating combined reference. Exiting."
        exit 1
    fi
    
    echo "Combined reference created successfully."
fi

# Ensure combined reference is indexed
if [ ! -f "${COMBINED_REF}.bwt" ]; then
    echo "Indexing combined reference for BWA..."
    bwa index "$COMBINED_REF"
    if [ $? -ne 0 ]; then
        echo "Error indexing combined reference. Exiting."
        exit 1
    fi
fi

if [ ! -f "${COMBINED_REF}.fai" ]; then
    echo "Indexing combined reference for samtools..."
    samtools faidx "$COMBINED_REF"
fi

# Process each BAM file (adapted from Script 3 approach)
for BAM_FILE in "$BAM_DIR"/*.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    # Check if output already exists
    if [ -f "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam" ]; then
        echo "Skipping $SAMPLE_NAME, output BAM file already exists"
        continue
    fi
    
    echo "Processing $SAMPLE_NAME..."

    # Convert BAM to FASTQ (from Script 3)
    echo "Converting BAM to FASTQ..."
    samtools fastq "$BAM_FILE" > "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
    if [ $? -ne 0 ]; then
        echo "Error converting $SAMPLE_NAME BAM to FASTQ. Skipping."
        continue
    fi

    # Check if FASTQ is empty
    if [ ! -s "$ALIGN_DIR/${SAMPLE_NAME}.fastq" ]; then
        echo "Generated FASTQ file is empty for $SAMPLE_NAME. Skipping."
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
        continue
    fi

    # Align to combined reference (modified from Script 2)
    echo "Aligning to combined reference..."
    bwa mem -t $threads "$COMBINED_REF" "$ALIGN_DIR/${SAMPLE_NAME}.fastq" > "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
    if [ $? -ne 0 ]; then
        echo "Error aligning $SAMPLE_NAME. Skipping."
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
        continue
    fi

    # Convert SAM to BAM, sort and index (from Script 2)
    echo "Converting to BAM and sorting..."
    samtools view -@ $threads -Sb "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam" | \
    samtools sort -@ $threads -o "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam"
    if [ $? -ne 0 ]; then
        echo "Error processing BAM for $SAMPLE_NAME. Skipping."
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq" "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
        continue
    fi

    # Index the BAM file
    echo "Indexing BAM file..."
    samtools index "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam"
    if [ $? -ne 0 ]; then
        echo "Error indexing BAM for $SAMPLE_NAME. Skipping."
        continue
    fi

    # Generate alignment statistics (from Script 3)
    echo "Generating alignment statistics..."
    samtools flagstat "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam" > "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt"
    if [ $? -eq 0 ]; then
        # Extract statistics
        TOTAL_READS=$(grep "in total" "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" | cut -d ' ' -f 1)
        MAPPED_READS=$(grep "mapped (" "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" | head -1 | cut -d ' ' -f 1)
        
        # Calculate mapping rate
        if [[ "$TOTAL_READS" -gt 0 ]]; then
            MAPPING_RATE=$(echo "scale=2; $MAPPED_READS * 100 / $TOTAL_READS" | bc -l 2>/dev/null || echo "0")
        else
            MAPPING_RATE="0"
        fi
        
        # Add to summary CSV
        echo "$SAMPLE_NAME,$TOTAL_READS,$MAPPED_READS,${MAPPING_RATE}%" >> "$OUTPUT_CSV"
        
        echo "$SAMPLE_NAME: $MAPPED_READS/$TOTAL_READS mapped (${MAPPING_RATE}%)"
    fi

    # Clean up intermediate files
    rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq" "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
    
    echo "Successfully processed $SAMPLE_NAME"
done

echo "All BAM files processed successfully!"
echo "Summary statistics saved to: $OUTPUT_CSV"