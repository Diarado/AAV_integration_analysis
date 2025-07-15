#!/bin/bash

# BAM to Combined Reference Alignment Script (Fixed Version)
# 
# This script processes BAM files containing sequencing reads and re-aligns them 
# to a combined reference genome (human hg38 + AAV vector).

# Define paths
BAM_DIR="../data/barcode03"                    # Input BAM files directory
ALIGN_DIR="../data/aligned_bams_barcode03"     # Output directory
HG38_REF="../data/hg38/hg38.fa"               # Human genome reference
AAV_REF_GB="../data/barcode03_ref/37825-AAV2-CAG-EGFP.gb"  # AAV reference (GenBank)
AAV_REF_FA="../data/barcode03_ref/AAV2-CAG-EGFP.fa"        # AAV reference (FASTA)
COMBINED_REF="../data/barcode03_ref/combined_ref.fa"       # Combined reference
OUTPUT_CSV="../data/alignment_summary.csv"     # Summary statistics
threads=4                                       # Number of threads to use

# Create output directory
mkdir -p "$ALIGN_DIR"
mkdir -p "$(dirname "$COMBINED_REF")"

# Initialize summary CSV
echo "Sample_Name,Total_Reads,Mapped_Reads,Mapping_Rate" > "$OUTPUT_CSV"

# Function to convert GenBank to FASTA
convert_genbank_to_fasta() {
    local gb_file=$1
    local fa_file=$2
    
    echo "Converting GenBank format to FASTA..."
    
    # Check if we have seqret (EMBOSS) available
    if command -v seqret &> /dev/null; then
        seqret -sequence "$gb_file" -outseq "$fa_file" -osformat fasta
    # Check if we have any2fasta available
    elif command -v any2fasta &> /dev/null; then
        any2fasta "$gb_file" > "$fa_file"
    # Use Python with Biopython if available
    elif command -v python3 &> /dev/null; then
        python3 -c "
from Bio import SeqIO
import sys
try:
    records = SeqIO.parse('$gb_file', 'genbank')
    SeqIO.write(records, '$fa_file', 'fasta')
    print('Successfully converted GenBank to FASTA')
except ImportError:
    print('Error: Biopython not installed. Please install with: pip install biopython')
    sys.exit(1)
except Exception as e:
    print(f'Error converting file: {e}')
    sys.exit(1)
"
    else
        echo "Error: No suitable tool found to convert GenBank to FASTA"
        echo "Please install one of: EMBOSS (seqret), any2fasta, or Biopython"
        return 1
    fi
    
    return $?
}

# Check and convert AAV reference if needed
if [ ! -f "$AAV_REF_FA" ]; then
    if [ -f "$AAV_REF_GB" ]; then
        echo "AAV reference found in GenBank format. Converting to FASTA..."
        convert_genbank_to_fasta "$AAV_REF_GB" "$AAV_REF_FA"
        if [ $? -ne 0 ] || [ ! -f "$AAV_REF_FA" ]; then
            echo "Error: Failed to convert GenBank to FASTA format"
            echo "Please manually convert $AAV_REF_GB to FASTA format"
            exit 1
        fi
    else
        echo "Error: AAV reference not found at $AAV_REF_GB"
        exit 1
    fi
fi

# Function to validate FASTA format
validate_fasta() {
    local fasta_file=$1
    local file_desc=$2
    
    echo "Validating $file_desc FASTA format..."
    
    # Check if file exists and is not empty
    if [ ! -s "$fasta_file" ]; then
        echo "Error: $file_desc file is empty or doesn't exist"
        return 1
    fi
    
    # Check FASTA format (should start with >)
    if ! head -n 1 "$fasta_file" | grep -q "^>"; then
        echo "Error: $file_desc doesn't appear to be in FASTA format"
        return 1
    fi
    
    # Check for common issues
    if grep -q -E "^[^>ACGTUNRYSWKMBDHV\-\n\r]" "$fasta_file"; then
        echo "Warning: $file_desc contains unexpected characters"
    fi
    
    echo "$file_desc format validation passed"
    return 0
}

# Check if combined reference exists before creating
if [ -f "$COMBINED_REF" ]; then
    echo "Combined reference already exists at $COMBINED_REF"
    
    # Validate the existing combined reference
    validate_fasta "$COMBINED_REF" "existing combined reference"
    if [ $? -ne 0 ]; then
        echo "Existing combined reference failed validation. Recreating..."
        rm -f "$COMBINED_REF" "${COMBINED_REF}."*
    fi
fi

# Check if combined reference needs to be created
if [ ! -f "$COMBINED_REF" ]; then
    echo "Creating combined reference..."
    
    # Validate hg38 reference
    if [ ! -f "$HG38_REF" ]; then
        echo "Error: hg38.fa not found at expected path: $HG38_REF"
        exit 1
    fi
    validate_fasta "$HG38_REF" "hg38 reference"
    if [ $? -ne 0 ]; then
        exit 1
    fi
    
    # Validate AAV reference
    validate_fasta "$AAV_REF_FA" "AAV reference"
    if [ $? -ne 0 ]; then
        exit 1
    fi
    
    # Combine the references
    echo "Combining hg38 and AAV references..."
    cat "$HG38_REF" > "$COMBINED_REF"
    echo "" >> "$COMBINED_REF"  # Ensure newline between files
    cat "$AAV_REF_FA" >> "$COMBINED_REF"
    
    if [ $? -ne 0 ]; then
        echo "Error creating combined reference"
        exit 1
    fi
    
    # Validate combined reference
    validate_fasta "$COMBINED_REF" "combined reference"
    if [ $? -ne 0 ]; then
        echo "Error: Combined reference validation failed"
        rm -f "$COMBINED_REF"
        exit 1
    fi
    
    echo "Combined reference created successfully"
fi

# Index combined reference for BWA
if [ ! -f "${COMBINED_REF}.bwt" ]; then
    echo "Indexing combined reference for BWA..."
    bwa index "$COMBINED_REF" 2>&1 | tee bwa_index.log
    if [ $? -ne 0 ] || [ ! -f "${COMBINED_REF}.bwt" ]; then
        echo "Error: BWA indexing failed. Check bwa_index.log for details"
        exit 1
    fi
    echo "BWA indexing completed successfully"
fi

# Index combined reference for samtools
if [ ! -f "${COMBINED_REF}.fai" ]; then
    echo "Indexing combined reference for samtools..."
    samtools faidx "$COMBINED_REF" 2>&1 | tee samtools_index.log
    if [ $? -ne 0 ] || [ ! -f "${COMBINED_REF}.fai" ]; then
        echo "Error: Samtools indexing failed. Check samtools_index.log for details"
        exit 1
    fi
    echo "Samtools indexing completed successfully"
fi

# Check if BAM directory exists
if [ ! -d "$BAM_DIR" ]; then
    echo "Error: BAM directory not found at $BAM_DIR"
    exit 1
fi

# Count BAM files
BAM_COUNT=$(find "$BAM_DIR" -name "*.bam" -type f | wc -l)
if [ "$BAM_COUNT" -eq 0 ]; then
    echo "Error: No BAM files found in $BAM_DIR"
    exit 1
fi
echo "Found $BAM_COUNT BAM files to process"

# Process each BAM file
processed=0
failed=0
skipped=0

for BAM_FILE in "$BAM_DIR"/*.bam; do
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam)
    
    # Check if output already exists
    if [ -f "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam" ]; then
        echo "Skipping $SAMPLE_NAME, output BAM file already exists"
        
        # If stats exist, add them to the summary
        if [ -f "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" ]; then
            TOTAL_READS=$(grep "in total" "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" | cut -d ' ' -f 1)
            MAPPED_READS=$(grep "mapped (" "$ALIGN_DIR/${SAMPLE_NAME}_flagstat.txt" | head -1 | cut -d ' ' -f 1)
            
            if [[ "$TOTAL_READS" -gt 0 ]]; then
                MAPPING_RATE=$(echo "scale=2; $MAPPED_READS * 100 / $TOTAL_READS" | bc -l 2>/dev/null || echo "0")
            else
                MAPPING_RATE="0"
            fi
            
            echo "$SAMPLE_NAME,$TOTAL_READS,$MAPPED_READS,${MAPPING_RATE}%" >> "$OUTPUT_CSV"
        fi
        
        ((skipped++))
        continue
    fi
    
    echo ""
    echo "Processing $SAMPLE_NAME..."
    echo "----------------------------------------"

    # Convert BAM to FASTQ
    echo "Converting BAM to FASTQ..."
    samtools fastq "$BAM_FILE" > "$ALIGN_DIR/${SAMPLE_NAME}.fastq" 2>"$ALIGN_DIR/${SAMPLE_NAME}_fastq.err"
    if [ $? -ne 0 ]; then
        echo "Error converting $SAMPLE_NAME BAM to FASTQ. Check ${SAMPLE_NAME}_fastq.err"
        ((failed++))
        continue
    fi

    # Check if FASTQ is empty
    if [ ! -s "$ALIGN_DIR/${SAMPLE_NAME}.fastq" ]; then
        echo "Generated FASTQ file is empty for $SAMPLE_NAME. Skipping."
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq"
        ((failed++))
        continue
    fi

    # Get read count from FASTQ
    READ_COUNT=$(grep -c "^@" "$ALIGN_DIR/${SAMPLE_NAME}.fastq" || echo "0")
    echo "Found $READ_COUNT reads in FASTQ"

    # Align to combined reference
    echo "Aligning to combined reference..."
    bwa mem -t $threads "$COMBINED_REF" "$ALIGN_DIR/${SAMPLE_NAME}.fastq" \
        > "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam" 2>"$ALIGN_DIR/${SAMPLE_NAME}_bwa.err"
    if [ $? -ne 0 ]; then
        echo "Error aligning $SAMPLE_NAME. Check ${SAMPLE_NAME}_bwa.err"
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq" "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
        ((failed++))
        continue
    fi

    # Convert SAM to BAM, sort and index
    echo "Converting to BAM and sorting..."
    samtools view -@ $threads -Sb "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam" 2>"$ALIGN_DIR/${SAMPLE_NAME}_view.err" | \
    samtools sort -@ $threads -o "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam" 2>"$ALIGN_DIR/${SAMPLE_NAME}_sort.err"
    if [ $? -ne 0 ]; then
        echo "Error processing BAM for $SAMPLE_NAME. Check view/sort error files"
        rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq" "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
        ((failed++))
        continue
    fi

    # Index the BAM file
    echo "Indexing BAM file..."
    samtools index "$ALIGN_DIR/${SAMPLE_NAME}_with_header.bam"
    if [ $? -ne 0 ]; then
        echo "Error indexing BAM for $SAMPLE_NAME"
        ((failed++))
        continue
    fi

    # Generate alignment statistics
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
        
        echo "Alignment stats: $MAPPED_READS/$TOTAL_READS mapped (${MAPPING_RATE}%)"
    fi

    # Clean up intermediate files
    rm -f "$ALIGN_DIR/${SAMPLE_NAME}.fastq" "$ALIGN_DIR/${SAMPLE_NAME}_with_header.sam"
    rm -f "$ALIGN_DIR/${SAMPLE_NAME}_"*.err  # Clean up error files if successful
    
    echo "Successfully processed $SAMPLE_NAME"
    ((processed++))
done

echo ""
echo "========================================="
echo "Processing complete!"
echo "Successfully processed: $processed samples"
echo "Failed: $failed samples"
echo "Skipped (already aligned): $skipped samples"
echo "Summary statistics saved to: $OUTPUT_CSV"
echo "========================================="

# Additional diagnostics if all samples failed
if [ "$processed" -eq 0 ] && [ "$failed" -gt 0 ] && [ "$skipped" -eq 0 ]; then
    echo ""
    echo "WARNING: All samples failed to process!"
    echo "Common issues to check:"
    echo "1. Verify the combined reference was created correctly"
    echo "2. Check BWA index files exist: ${COMBINED_REF}.bwt, ${COMBINED_REF}.pac, etc."
    echo "3. Check samtools index exists: ${COMBINED_REF}.fai"
    echo "4. Review error files in $ALIGN_DIR for specific error messages"
fi