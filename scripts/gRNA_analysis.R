library(Rsamtools)
library(Biostrings)
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(stringi)

# Set working directory
setwd("D:/Jiahe/IU/AAV/HeLa_project")

gRNAs <- read_excel("gRNAs.xlsx")
barcodes <- read_excel("Nanopore Native barcodes list.xlsx")

# Generate sequences to search for and add cut site information
gRNAs <- gRNAs %>%
  mutate(
    first17 = substr(sequence, 1, 17),
    last3 = substr(sequence, nchar(sequence) - 2, nchar(sequence)),
    last3_PAM = paste0(last3, PAM),
    cut_site_position = nchar(sequence) - 3
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

# Function to count mismatches between two sequences
count_mismatches <- function(seq1, seq2) {
  # Ensure sequences are same length
  if (nchar(seq1) != nchar(seq2)) return(Inf)
  # Count positions where characters differ
  sum(strsplit(seq1, '')[[1]] != strsplit(seq2, '')[[1]])
}

# Function to check for barcode matches with tolerance
find_barcode_match <- function(sequence, barcode_list, max_mismatches = 2) {
  for (barcode in barcode_list) {
    if (nchar(sequence) >= nchar(barcode)) {
      mismatches <- count_mismatches(
        substr(sequence, 1, nchar(barcode)),
        barcode
      )
      if (mismatches <= max_mismatches) return(TRUE)
    }
  }
  return(FALSE)
}

trim_adapters <- function(read_seq, barcode_seqs_forward, barcode_seqs_reverse) {
  # Convert input sequence to DNAString object
  if (is.character(read_seq)) {
    read_seq <- DNAString(read_seq)
  }
  
  # Convert barcodes to DNAStringSet objects
  if (is.character(barcode_seqs_forward)) {
    barcode_seqs_forward <- DNAStringSet(barcode_seqs_forward)
  }
  # if (is.character(barcode_seqs_reverse)) {
  #   barcode_seqs_reverse <- DNAStringSet(barcode_seqs_reverse)
  # }
  
  # Helper function to trim AT sequences
  trim_AT_sequence <- function(seq, from_start = TRUE) {
    if (length(seq) == 0) return(seq)
    
    seq_chars <- strsplit(as.character(seq), "")[[1]]
    
    if (from_start) {
      # Find first position that's not A or T
      non_at_pos <- which(!seq_chars %in% c("A", "T"))[1]
      
      if (!is.na(non_at_pos)) {
        return(subseq(seq, start = non_at_pos))
      }
    } else {
      # Find last position that's not A or T
      non_at_pos <- rev(which(!seq_chars %in% c("A", "T")))[1]
      
      if (!is.na(non_at_pos)) {
        return(subseq(seq, end = non_at_pos))
      }
    }
    return(seq)
  }
  
  # Process head of sequence
  original_length <- length(read_seq)
  head_processed <- FALSE
  
  # Iterate through barcodes
  for (i in 1:length(barcode_seqs_forward)) {
    barcode <- barcode_seqs_forward[[i]]
    matches <- matchPattern(barcode, read_seq, max.mismatch = 1)
    
    if (length(matches) > 0) {
      first_match_end <- end(matches)[1]
      if (start(matches)[1] == 1) {  # Found at start
        remaining_seq <- subseq(read_seq, start = first_match_end + 1)
        read_seq <- trim_AT_sequence(remaining_seq, from_start = TRUE)
        head_processed <- TRUE
        break
      }
    }
  }
  
  # Only process tail if head processing didn't consume entire sequence
  if (!head_processed || length(read_seq) > 0) {
    for (i in 1:length(barcode_seqs_forward)) {
      barcode <- barcode_seqs_forward[[i]]
      matches <- matchPattern(barcode, read_seq, max.mismatch = 1)
      
      if (length(matches) > 0) {
        last_match_idx <- length(matches)
        if (end(matches)[last_match_idx] == length(read_seq)) {  # Found at end
          remaining_seq <- subseq(read_seq, end = start(matches)[last_match_idx] - 1)
          read_seq <- trim_AT_sequence(remaining_seq, from_start = FALSE)
          break
        }
      }
    }
  }
  
  # Convert back to character if needed
  return(as.character(read_seq))
}

# Function to find sequence matches with tolerance
find_sequence_matches <- function(target_seq, reference_seq, max_mismatches = 2) {
  target_len <- nchar(target_seq)
  ref_len <- nchar(reference_seq)
  
  if (target_len < ref_len) return(FALSE)
  
  # Check start of sequence
  start_seq <- substr(target_seq, 1, ref_len)
  if (count_mismatches(start_seq, reference_seq) <= max_mismatches) return(TRUE)
  
  # Check end of sequence
  end_seq <- substr(target_seq, target_len - ref_len + 1, target_len)
  if (count_mismatches(end_seq, reference_seq) <= max_mismatches) return(TRUE)
  
  return(FALSE)
}

# Function to process a BAM file
process_bam_file <- function(bam_file, sequence_to_id, gRNAs_info = gRNAs, min_alignment_quality = 20, min_read_length = 40) {
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
    
    # Trim adapter sequences
    read_seq <- trim_adapters(read_seq, 
                              barcodes$`Forward sequence`, 
                              barcodes$`Reverse sequence`)
    
    read_len <- nchar(read_seq)
    if (read_len < min_read_length) next
    
    # Get the first and last 30 bp of the read
    first30 <- substr(read_seq, 1, min(30, read_len))
    last30 <- substr(read_seq, max(1, read_len - 29), read_len)
    
    # Find matches with tolerance
    matched_sequences <- sequence_id_map$sequence[
      sapply(sequence_id_map$sequence, function(seq) {
        find_sequence_matches(first30, seq) || 
          find_sequence_matches(last30, seq)
      })
    ]
    
    if (length(matched_sequences) > 0) {
      # Validate junction: check for both AAV and host genome alignments
      has_aav <- any(grepl("AAV", read_records$Ref_Name, ignore.case = TRUE))
      has_host <- any(grepl("^chr", read_records$Ref_Name, ignore.case = TRUE))
      
      if (!(has_aav && has_host)) {
        next  # Skip if not a true junction
      }
      
      # Get matched gRNA information
      matched_gRNA_ids <- unique(unlist(strsplit(sequence_to_id[matched_sequences], ",")))
      
      # For each matched gRNA, process the cut site information
      for (gRNA_id in matched_gRNA_ids) {
        gRNA_info <- gRNAs_info[gRNAs_info$id == gRNA_id, ]
        cut_site_pos <- gRNA_info$cut_site_position
        
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
            if (is.na(AAV_MapQ) || mapq > AAV_MapQ) {
              AAV_Start <- position
              AAV_End <- position
              AAV_MapQ <- mapq
            }
          } else if (grepl("^chr", ref_name, ignore.case = TRUE)) {
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
          
          # Calculate distance from cut site
          distance_from_cut <- abs(as.numeric(AAV_Start) - as.numeric(cut_site_pos))
          
          result <- data.frame(
            Sample = sample_name,
            Read_Name = read_name,
            gRNA_id = gRNA_id,
            gRNA_Sequence = gRNA_info$sequence,
            Cut_Site_Position = cut_site_pos,
            AAV_Start = AAV_Start,
            AAV_End = AAV_End,
            AAV_MapQ = AAV_MapQ,
            Host_Chromosome = Host_Chromosome,
            Host_Start = Host_Start,
            Host_MapQ = Host_MapQ,
            Distance_From_Cut = distance_from_cut,
            Read_Length = read_len,
            stringsAsFactors = FALSE
          )
          
          results_list[[length(results_list) + 1]] <- result
        }
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
        Unique_Junction_Sites = n_distinct(paste(Host_Chromosome, Host_Start, AAV_Start)),
        Mean_Distance_From_Cut = mean(Distance_From_Cut),
        Reads_Within_10bp_Of_Cut = sum(Distance_From_Cut <= 10)
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
  
  # Create a summary data frame
  summary_results <- all_results %>%
    group_by(gRNA_id, gRNA_Sequence) %>%
    summarize(
      Total_Reads = n(),
      Unique_Integration_Sites = n_distinct(paste(Host_Chromosome, Host_Start)),
      Mean_Distance_From_Cut = mean(Distance_From_Cut),
      Median_Distance_From_Cut = median(Distance_From_Cut),
      Reads_Within_10bp = sum(Distance_From_Cut <= 10),
      Mean_AAV_MapQ = mean(AAV_MapQ),
      Mean_Host_MapQ = mean(Host_MapQ),
      .groups = 'drop'
    ) %>%
    arrange(desc(Total_Reads))
  
  # Write both detailed and summary results
  write.csv(all_results, file = "detailed_cut_sites_with_gRNA.csv", row.names = FALSE)
  write.csv(summary_results, file = "cut_sites_summary_by_gRNA.csv", row.names = FALSE)
} else {
  all_results <- data.frame(
    Sample = character(),
    Read_Name = character(),
    gRNA_id = character(),
    gRNA_Sequence = character(),
    Cut_Site_Position = numeric(),
    AAV_Start = character(),
    AAV_End = character(),
    Host_Chromosome = character(),
    Host_Start = character(),
    Distance_From_Cut = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(all_results, file = "detailed_cut_sites_with_gRNA.csv", row.names = FALSE)
}

cat("Analysis complete. Results saved to detailed_cut_sites_with_gRNA.csv and 
    cut_sites_summary_by_gRNA.csv\n")