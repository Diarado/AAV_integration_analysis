# Load necessary libraries
library(Rsamtools)
library(Biostrings)
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)

# Set working directory
setwd("D:/Jiahe/IU/AAV/HeLa_project")

# Read the updated gRNAs.xlsx file with 'id' column
gRNAs <- read_excel("gRNAs.xlsx")

# Generate sequences to search for
gRNAs <- gRNAs %>%
  mutate(
    first17 = substr(sequence, 1, 17),
    last3 = substr(sequence, nchar(sequence) - 2, nchar(sequence)),
    last3_PAM = paste0(last3, PAM)
  )

# Create a mapping from sequences to gRNA ids
sequence_id_map <- gRNAs %>%
  select(id, first17, last3_PAM) %>%
  pivot_longer(cols = c(first17, last3_PAM), names_to = "type", values_to = "sequence") %>%
  group_by(sequence) %>%
  summarise(gRNA_id = paste(unique(id), collapse = ",")) %>%
  ungroup()

# Convert the mapping to a named vector for efficient lookup
sequence_to_id <- setNames(sequence_id_map$gRNA_id, sequence_id_map$sequence)

# Function to process a BAM file
process_bam_file <- function(bam_file, sequence_to_id, min_alignment_quality = 20, min_read_length = 40) {
  sample_name <- basename(bam_file)
  sample_name <- sub("_with_header.bam$", "", sample_name)
  
  # Open the BAM file with additional parameters for quality
  bam <- BamFile(bam_file)
  
  # Include mapq (mapping quality) in the parameters
  param <- ScanBamParam(
    what = c("qname", "seq", "rname", "pos", "mapq", "cigar"),
    mapqFilter = min_alignment_quality
  )
  
  # Read BAM data
  bam_data <- scanBam(bam, param = param)[[1]]
  
  # Create a data frame from the BAM data with quality information
  bam_df <- data.frame(
    Read_Name = bam_data$qname,
    Sequence = as.character(bam_data$seq),
    Ref_Name = as.character(bam_data$rname),
    Position = bam_data$pos,
    MapQ = bam_data$mapq,
    CIGAR = bam_data$cigar,
    stringsAsFactors = FALSE
  )
  
  # Remove reads with NA sequences or below minimum length
  bam_df <- bam_df[!is.na(bam_df$Sequence) & 
                     nchar(bam_df$Sequence) >= min_read_length, ]
  
  # Initialize a list to store results
  results_list <- list()
  
  # Get unique read names
  read_names <- unique(bam_df$Read_Name)
  
  # Progress indicator
  total_reads <- length(read_names)
  cat(sprintf("Processing %d reads for %s\n", total_reads, sample_name))
  progress_step <- max(1, floor(total_reads/10))
  
  # Loop through each read
  for (i in seq_along(read_names)) {
    # Progress update
    if (i %% progress_step == 0) {
      cat(sprintf("Processed %d%% of reads...\n", floor(i/total_reads * 100)))
    }
    
    read_name <- read_names[i]
    read_records <- bam_df[bam_df$Read_Name == read_name, ]
    
    # Skip if less than 2 alignments (not a junction)
    if (nrow(read_records) < 2) {
      next
    }
    
    read_seq <- read_records$Sequence[1]
    read_len <- nchar(read_seq)
    
    # Skip if read is too short
    if (read_len < min_read_length) {
      next
    }
    
    # Get the first and last 20 bp of the read
    first20 <- substr(read_seq, 1, min(20, read_len))
    last20 <- substr(read_seq, max(1, read_len - 19), read_len)
    
    # Strict matching using exact matches
    matched_sequences <- sequence_id_map$sequence[
      sapply(sequence_id_map$sequence, function(seq) {
        # Check for exact matches at the start or end
        (nchar(seq) <= nchar(first20) && startsWith(first20, seq)) || 
          (nchar(seq) <= nchar(last20) && endsWith(last20, seq))
      })
    ]
    
    if (length(matched_sequences) > 0) {
      # Validate junction: check for both AAV and host genome alignments
      has_aav <- any(grepl("AAV", read_records$Ref_Name, ignore.case = TRUE))
      has_host <- any(grepl("^chr", read_records$Ref_Name, ignore.case = TRUE))
      
      if (!(has_aav && has_host)) {
        next  # Skip if not a true junction
      }
      
      # Get unique gRNA ids for matched sequences
      matched_gRNA_ids <- unique(unlist(strsplit(sequence_to_id[matched_sequences], ",")))
      gRNA_ids_str <- paste(matched_gRNA_ids, collapse = ",")
      
      # Initialize alignment variables
      AAV_Start <- NA
      AAV_End <- NA
      Host_Chromosome <- NA
      Host_Start <- NA
      AAV_MapQ <- NA
      Host_MapQ <- NA
      
      # Process alignments with quality checks
      for (j in seq_len(nrow(read_records))) {
        ref_name <- read_records$Ref_Name[j]
        position <- read_records$Position[j]
        mapq <- read_records$MapQ[j]
        
        if (grepl("AAV", ref_name, ignore.case = TRUE)) {
          # Update AAV alignment if better quality
          if (is.na(AAV_MapQ) || mapq > AAV_MapQ) {
            AAV_Start <- position
            AAV_End <- position  # Will be updated with CIGAR if needed
            AAV_MapQ <- mapq
          }
        } else if (grepl("^chr", ref_name, ignore.case = TRUE)) {
          # Update host genome alignment if better quality
          if (is.na(Host_MapQ) || mapq > Host_MapQ) {
            Host_Chromosome <- ref_name
            Host_Start <- position
            Host_MapQ <- mapq
          }
        }
      }
      
      # Only store results if both alignments meet quality threshold
      if (!is.na(AAV_MapQ) && !is.na(Host_MapQ) && 
          AAV_MapQ >= min_alignment_quality && 
          Host_MapQ >= min_alignment_quality) {
        
        # Store the result
        result <- data.frame(
          Sample = sample_name,
          Read_Name = read_name,
          AAV_Start = ifelse(is.na(AAV_Start), "NA", AAV_Start),
          AAV_End = ifelse(is.na(AAV_End), "NA", AAV_End),
          AAV_MapQ = AAV_MapQ,
          Host_Chromosome = ifelse(is.na(Host_Chromosome), "NA", Host_Chromosome),
          Host_Start = ifelse(is.na(Host_Start), "NA", Host_Start),
          Host_MapQ = Host_MapQ,
          gRNA_id = gRNA_ids_str,
          Read_Length = read_len,
          stringsAsFactors = FALSE
        )
        
        results_list[[length(results_list) + 1]] <- result
      }
    }
  }
  
  # Combine results into a data frame
  if (length(results_list) > 0) {
    results_df <- do.call(rbind, results_list)
    
    # Add additional summary statistics
    results_df <- results_df %>%
      group_by(gRNA_id) %>%
      mutate(
        Support_Reads = n(),
        Unique_Junction_Sites = n_distinct(paste(Host_Chromosome, Host_Start, AAV_Start))
      ) %>%
      ungroup()
    
  } else {
    results_df <- NULL
  }
  
  cat(sprintf("Completed processing %s\n", sample_name))
  if (!is.null(results_df)) {
    cat(sprintf("Found %d junction reads\n", nrow(results_df)))
  }
  
  return(results_df)
}

# Get list of BAM files
bam_files <- list.files(path = "aligned_bams_test", pattern = "_with_header.bam$", full.names = TRUE)

# Initialize a list to collect all results
all_results_list <- list()

# Process each BAM file
for (bam_file in bam_files) {
  cat("Processing", bam_file, "\n")
  bam_results <- process_bam_file(bam_file, sequence_to_id)
  if (!is.null(bam_results)) {
    all_results_list[[length(all_results_list) + 1]] <- bam_results
  }
}

# Combine all results into a single data frame
if (length(all_results_list) > 0) {
  all_results <- do.call(rbind, all_results_list)
} else {
  all_results <- data.frame(
    Sample = character(),
    Read_Name = character(),
    AAV_Start = character(),
    AAV_End = character(),
    Host_Chromosome = character(),
    Host_Start = character(),
    gRNA_id = character(),
    stringsAsFactors = FALSE
  )
}

# Write the results to a CSV file
write.csv(all_results, file = "cut_sites_with_gRNA_id.csv", row.names = FALSE)

cat("Analysis complete. Results saved to cut_sites_with_gRNA_id.csv\n")









# 
# 
# 
# # statistical analysis
# # 4. Define the get_mode Function
# get_mode <- function(v) {
#   uniqv <- unique(v)
#   tab <- tabulate(match(v, uniqv))
#   mode_values <- uniqv[tab == max(tab)]
#   return(mode_values)
# }
# 
# # 5. Calculate Basic Statistics
# max_value <- max(difference, na.rm = TRUE)
# min_value <- min(difference, na.rm = TRUE)
# median_value <- median(difference, na.rm = TRUE)
# mean_value <- mean(difference, na.rm = TRUE)
# sd_value <- sd(difference, na.rm = TRUE)
# quantiles <- quantile(difference, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
# 
# # 6. Calculate Mode
# mode_value <- get_mode(difference)
# 
# # 7. Display the Results
# cat("Statistics for AAV_End - AAV_Start:\n")
# cat("Maximum:", max_value, "\n")
# cat("Minimum:", min_value, "\n")
# cat("Median:", median_value, "\n")
# cat("Mean:", mean_value, "\n")
# cat("Standard Deviation:", sd_value, "\n")
# cat("Quantiles (25%, 50%, 75%):\n")
# print(quantiles)
# cat("Mode:", paste(mode_value, collapse = ", "), "\n")
# 
# # Filter for barcode13 (where '13' appears in the specified position)
# barcode13_cnt <- data |>
#   filter(!is.na(Host_Chromosome) & substr(Sample, 22, 23) == "13")
# 
# # Filter for barcode14 (where '14' appears in the specified position)
# barcode14_cnt <- data |>
#   filter(!is.na(Host_Chromosome) & substr(Sample, 22, 23) == "14")
# 
# 
# write_csv(barcode13_cnt, "D:/Jiahe/IU/AAV/HeLa_project/output/aav_13_integrations.csv")
# write_csv(barcode14_cnt, "D:/Jiahe/IU/AAV/HeLa_project/output/aav_14_integrations.csv")
# find gRNA seq