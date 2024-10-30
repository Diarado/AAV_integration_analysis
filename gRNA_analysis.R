# Load necessary libraries
library(Rsamtools)
library(Biostrings)
library(dplyr)
library(stringr)

# Set working directory if needed
setwd("D:/Jiahe/IU/AAV/HeLa_project")

# Read the gRNAs.csv file
gRNAs <- read.csv("gRNAs.csv", sep="\t", stringsAsFactors = FALSE)

# Generate sequences to search for
gRNAs$first17 <- substr(gRNAs$sequence, 1, 17)
gRNAs$last3 <- substr(gRNAs$sequence, nchar(gRNAs$sequence) - 2, nchar(gRNAs$sequence))
gRNAs$last3_PAM <- paste0(gRNAs$last3, gRNAs$PAM)

# Create a vector of sequences to search for
search_sequences <- unique(c(gRNAs$first17, gRNAs$last3_PAM))

# Function to process a BAM file
process_bam_file <- function(bam_file, search_sequences) {
  sample_name <- basename(bam_file)
  sample_name <- sub("_with_header.bam$", "", sample_name)
  
  # Open the BAM file
  bam <- BamFile(bam_file)
  
  # Define parameters to extract necessary information
  param <- ScanBamParam(what = c("qname", "seq", "rname", "pos"))
  bam_data <- scanBam(bam, param = param)[[1]]
  
  # Create a data frame from the BAM data
  bam_df <- data.frame(
    Read_Name = bam_data$qname,
    Sequence = as.character(bam_data$seq),
    Ref_Name = as.character(bam_data$rname),
    Position = bam_data$pos,
    stringsAsFactors = FALSE
  )
  
  # Remove reads with NA sequences (unmapped reads)
  bam_df <- bam_df[!is.na(bam_df$Sequence), ]
  
  # Initialize a list to store results
  results_list <- list()
  
  # Get unique read names
  read_names <- unique(bam_df$Read_Name)
  
  # Loop through each read
  for (read_name in read_names) {
    read_records <- bam_df[bam_df$Read_Name == read_name, ]
    read_seq <- read_records$Sequence[1]  # Assuming sequence is the same across alignments
    read_len <- nchar(read_seq)
    
    # Get the first and last 20 bp of the read
    first20 <- substr(read_seq, 1, min(20, read_len))
    last20 <- substr(read_seq, max(1, read_len - 19), read_len)
    
    # Check for matches in the ends of the read
    match_found <- FALSE
    for (search_seq in search_sequences) {
      if (str_detect(first20, fixed(search_seq)) || str_detect(last20, fixed(search_seq))) {
        match_found <- TRUE
        break
      }
    }
    
    if (match_found) {
      # Initialize variables
      AAV_Start <- NA
      AAV_End <- NA
      Host_Chromosome <- NA
      Host_Start <- NA
      
      # Check for alignments to AAV and host
      for (i in seq_len(nrow(read_records))) {
        ref_name <- read_records$Ref_Name[i]
        position <- read_records$Position[i]
        
        if (grepl("AAV", ref_name, ignore.case = TRUE)) {
          # Alignment to AAV
          if (is.na(AAV_Start) || position < AAV_Start) {
            AAV_Start <- position
          }
          if (is.na(AAV_End) || position > AAV_End) {
            AAV_End <- position
          }
        } else if (grepl("^chr", ref_name, ignore.case = TRUE)) {
          # Alignment to host genome
          Host_Chromosome <- ref_name
          if (is.na(Host_Start) || position < Host_Start) {
            Host_Start <- position
          }
        }
      }
      
      # Store the result
      result <- data.frame(
        Sample = sample_name,
        Read_Name = read_name,
        AAV_Start = AAV_Start,
        AAV_End = AAV_End,
        Host_Chromosome = Host_Chromosome,
        Host_Start = Host_Start,
        stringsAsFactors = FALSE
      )
      results_list[[length(results_list) + 1]] <- result
    }
  }
  
  # Combine results into a data frame
  if (length(results_list) > 0) {
    results_df <- do.call(rbind, results_list)
  } else {
    results_df <- NULL
  }
  
  return(results_df)
}

# Get list of BAM files
bam_files <- list.files(path = "aligned_bams", pattern = "_with_header.bam$", full.names = TRUE)

# Initialize a list to collect all results
all_results_list <- list()

# Process each BAM file
for (bam_file in bam_files) {
  cat("Processing", bam_file, "\n")
  bam_results <- process_bam_file(bam_file, search_sequences)
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
    AAV_Start = numeric(),
    AAV_End = numeric(),
    Host_Chromosome = character(),
    Host_Start = numeric(),
    stringsAsFactors = FALSE
  )
}

# Write the results to a CSV file
write.csv(all_results, file = "cut_sites.csv", row.names = FALSE)

cat("Analysis complete. Results saved to cut_sites.csv\n")










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