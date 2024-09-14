# Load required libraries
print("Loading required libraries...")
library(ShortRead)
library(GenomicAlignments)
library(Rsamtools)
library(tidyverse)
library(yaml)

# Set constants
ALIGN_DIR <- "/mnt/d/jiahe/IU/AAV/HeLa_project/aligned_bams_test"
AAV_REF_NAME <- "pssAAV-CB-EGFP"  
OUTPUT_CSV <- "/mnt/d/jiahe/IU/AAV/HeLa_project/output/aav_reads_locations.csv"

# Retrieve aligned BAM files
print("Retrieving BAM files...")
bam_files <- list.files(ALIGN_DIR, pattern = "_with_header\\.bam$", full.names = TRUE)
print(paste("Total number of BAM files found:", length(bam_files)))

# Function to find reads aligned to AAV and save their locations
findAAVReads <- function(aligned_bam) {
  sample_name <- basename(aligned_bam)
  sample_name <- sub("_with_header\\.bam$", "", sample_name)
  print(paste("Processing sample:", sample_name))
  
  if (!file.exists(aligned_bam)) {
    stop(paste("BAM file not found:", aligned_bam))
  }
  
  # Load alignments
  print("Loading alignments...")
  param <- ScanBamParam(what = c("qname", "rname", "pos", "cigar"))
  reads <- readGAlignments(aligned_bam, param = param)
  
  if (length(reads) == 0) {
    print(paste("No alignments found in BAM file:", aligned_bam))
    return(NULL)
  }
  
  # Print all reference names to verify
  ref_names <- unique(seqnames(reads))
  print("Reference names found in BAM file:")
  print(ref_names)
  
  # Find reads aligned to AAV (verify the correct reference name)
  aav_reads <- reads[seqnames(reads) == AAV_REF_NAME]
  
  if (length(aav_reads) == 0) {
    print(paste("No AAV reads found in sample:", sample_name))
    return(NULL)
  }
  
  # Create a data frame with relevant information
  read_data <- data.frame(
    Sample = sample_name,
    Read_Name = mcols(aav_reads)$qname,
    AAV_Start = start(aav_reads),
    AAV_End = end(aav_reads)
  )
  
  print(paste("Found", nrow(read_data), "AAV reads in sample:", sample_name))
  return(read_data)
}

# Initialize a list to collect all AAV reads
all_aav_reads <- list()

# Loop over BAM files and collect AAV reads
for (aligned_bam in bam_files) {
  print("------------------------------------------------------------")
  aav_reads <- findAAVReads(aligned_bam)
  if (!is.null(aav_reads)) {
    all_aav_reads[[length(all_aav_reads) + 1]] <- aav_reads
  }
  print("------------------------------------------------------------")
}

# Combine all results into a single data frame
if (length(all_aav_reads) > 0) {
  all_aav_reads_df <- bind_rows(all_aav_reads)
  print(paste("Total AAV reads found across all samples:", nrow(all_aav_reads_df)))
  
  # Save to CSV
  print(paste("Saving results to:", OUTPUT_CSV))
  write_csv(all_aav_reads_df, OUTPUT_CSV)
} else {
  print("No AAV reads found in any sample.")
}

print("AAV read analysis completed.")

