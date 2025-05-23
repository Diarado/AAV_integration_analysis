# Load required libraries
print("Loading required libraries...")
library(GenomicAlignments)
library(Rsamtools)
library(tidyverse)

# Set constants
ALIGN_DIR <- "D:/jiahe/IU/AAV/HeLa_project/aligned_bams"
AAV_REF_NAME <- "pssAAV-CB-EGFP"  
OUTPUT_CSV <- "D:/jiahe/IU/AAV/HeLa_project/output/aav_reads_locations_with_readlen.csv"

# Retrieve aligned BAM files
print("Retrieving BAM files...")
bam_files <- list.files(ALIGN_DIR, pattern = "_with_header\\.bam$", full.names = TRUE)
print(paste("Total number of BAM files found:", length(bam_files)))

# Function to compute alignment length on the reference
compute_alignment_length <- function(cigar_string) {
  cigarWidthAlongReferenceSpace(cigar_string)
}

# Function to find reads aligned to AAV and save their locations
findAAVReads <- function(aligned_bam) {
  sample_name <- basename(aligned_bam)
  sample_name <- sub("_with_header\\.bam$", "", sample_name)
  print(paste("Processing sample:", sample_name))
  
  # Load alignments
  print("Loading alignments...")
  param <- ScanBamParam(what = c("qname", "rname", "pos", "cigar"), tag = "SA")
  reads <- tryCatch({
    readGAlignments(aligned_bam, param = param)
  }, error = function(e) {
    print(paste(sample_name, "is null, skipping"))
    return(NULL)
  })
  
  # Check if reads are NULL or empty
  if (is.null(reads) || length(reads) == 0) {
    print(paste("No alignments found in BAM file:", aligned_bam, "Skipping sample..."))
    return(NULL)
  }
  
  # Find reads aligned to AAV
  aav_reads <- reads[seqnames(reads) == AAV_REF_NAME]
  
  if (length(aav_reads) == 0) {
    print(paste("No AAV reads found in sample:", sample_name))
    return(NULL)
  }
  
  # Create a data frame with relevant information
  aav_data <- data.frame(
    Sample = rep(sample_name, length(aav_reads)),
    Read_Name = mcols(aav_reads)$qname,
    AAV_Start = start(aav_reads),
    AAV_End = end(aav_reads),
    Host_Chromosome = rep(NA, length(aav_reads)),  # Initialize with NA
    Host_Start = rep(NA, length(aav_reads)),       # Initialize with NA
    Host_Read_Start = rep(NA, length(aav_reads)),  # New column for host read start
    Host_Read_End = rep(NA, length(aav_reads)),    # New column for host read end
    stringsAsFactors = FALSE
  )
  
  # Extract SA tag details (supplementary alignments) and filter out AAV vector matches
  sa_tags <- mcols(aav_reads)$SA
  
  if (!is.null(sa_tags)) {
    sa_list <- lapply(sa_tags, function(tag) {
      if (!is.null(tag) && is.character(tag)) {
        return(strsplit(tag, ";"))
      } else {
        return(NULL)
      }
    })
    
    for (i in seq_along(sa_list)) {
      current_sa <- sa_list[[i]]
      
      # Ensure the list is not NULL and has entries
      if (!is.null(current_sa) && length(current_sa) > 0) {
        entries <- unlist(current_sa)
        for (entry in entries) {
          if (length(entry) > 0 && !is.na(entry) && entry != "") {
            sa_fields <- strsplit(entry, ",")[[1]]
            if (length(sa_fields) >= 4 && !grepl(AAV_REF_NAME, sa_fields[1])) {  # Ensure we have enough fields and not AAV_REF_NAME
              host_start <- as.numeric(sa_fields[2])
              strand <- sa_fields[3]
              cigar_string <- sa_fields[4]
              host_alignment_length <- compute_alignment_length(cigar_string)
              host_end <- host_start + host_alignment_length - 1
              aav_data$Host_Chromosome[i] <- sa_fields[1]
              aav_data$Host_Start[i] <- host_start  # Leave Host_Start as it is
              aav_data$Host_Read_Start[i] <- host_start  # New column value
              aav_data$Host_Read_End[i] <- host_end      # New column value
              break  # We found a host alignment, no need to look further
            }
          }
        }
      }
    }
  }
  
  return(aav_data)
}

# Process each BAM file and collect results
all_aav_data <- lapply(bam_files, findAAVReads)
all_aav_data_df <- bind_rows(all_aav_data)

# Write results to CSV
if (nrow(all_aav_data_df) > 0) {
  write_csv(all_aav_data_df, OUTPUT_CSV)
  print(paste("Results saved to:", OUTPUT_CSV))
} else {
  print("No AAV reads found in any sample.")
}

print("AAV read analysis completed.")
