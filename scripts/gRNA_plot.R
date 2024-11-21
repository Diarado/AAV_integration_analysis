# Load required libraries
library(dplyr)
library(ggplot2)
library(Biostrings)
library(readr)

# File paths - modify these to match your file locations
AAV_REF_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/pssAAV-CB-EGFP_ARM.fa"
AAV_READS_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/output/aav_reads_locations.csv"
CUT_SITES_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/output/detailed_cut_sites_with_gRNA.csv"

# Function to read and process the data
load_and_process_data <- function(aav_reads_path, cut_sites_path, aav_ref_path) {
  # Read AAV reference genome and get length
  aav_seq <- readDNAStringSet(aav_ref_path)
  aav_length <- width(aav_seq)[1]
  
  # Read AAV reads data and ensure proper data types
  aav_reads <- read.csv(aav_reads_path, stringsAsFactors = FALSE) %>%
    filter(!is.na(Host_Chromosome) & !is.na(Host_Start)) %>%
    select(Sample, Read_Name, AAV_Start, AAV_End, Host_Chromosome, Host_Start) %>%
    mutate(
      Sample = as.character(Sample),
      Read_Name = as.character(Read_Name),
      AAV_Start = as.numeric(AAV_Start),
      AAV_End = as.numeric(AAV_End),
      Host_Chromosome = as.character(Host_Chromosome),
      Host_Start = as.numeric(Host_Start)
    )
  
  # Read cut sites data and ensure proper data types
  cut_sites <- read.csv(cut_sites_path, stringsAsFactors = FALSE) %>%
    select(Sample, Read_Name, gRNA_id, Cut_Site_Position, AAV_Start, AAV_End) %>%
    mutate(
      Sample = as.character(Sample),
      Read_Name = as.character(Read_Name),
      gRNA_id = as.numeric(gRNA_id),
      Cut_Site_Position = as.numeric(Cut_Site_Position),
      AAV_Start = as.numeric(AAV_Start),
      AAV_End = as.numeric(AAV_End)
    )
  
  # Print data structure for debugging
  cat("AAV Reads structure:\n")
  str(aav_reads)
  cat("\nCut Sites structure:\n")
  str(cut_sites)
  
  # Verify data overlap
  overlapping_reads <- intersect(unique(aav_reads$Read_Name), unique(cut_sites$Read_Name))
  cat("\nNumber of overlapping reads:", length(overlapping_reads), "\n")
  
  if(length(overlapping_reads) == 0) {
    stop("No overlapping reads found between datasets")
  }
  
  return(list(
    aav_reads = aav_reads,
    cut_sites = cut_sites,
    aav_length = aav_length
  ))
}

# Function to merge and prepare data for plotting
prepare_plot_data <- function(aav_reads, cut_sites) {
  # Merge the data
  merged_data <- aav_reads %>%
    inner_join(cut_sites, by = c("Read_Name", "Sample"), suffix = c(".aav", ".cut"))
  
  # Convert to data frame and handle duplicates
  merged_data <- as.data.frame(merged_data) %>%
    distinct(Read_Name, .keep_all = TRUE) %>%
    mutate(
      Absolute_Cut_Position = AAV_Start.aav + Cut_Site_Position - 1,
      Read_Index = row_number()  # Add index for vertical positioning
    )
  
  return(merged_data)
}

# Function to create the visualization
create_aav_alignment_plot <- function(merged_data, aav_length) {
  ggplot() +
    # Base genome line
    geom_segment(
      aes(x = 0, xend = aav_length, y = 0, yend = 0),
      color = "black",
      size = 1
    ) +
    
    # Plot reads as gray segments
    geom_segment(
      data = merged_data,
      aes(
        x = AAV_Start.aav,
        xend = AAV_End.aav,
        y = Read_Index,
        yend = Read_Index
      ),
      color = "gray50",
      size = 0.5
    ) +
    
    # Plot cut sites as red marks
    geom_segment(
      data = merged_data,
      aes(
        x = Absolute_Cut_Position,
        xend = Absolute_Cut_Position,
        y = Read_Index - 0.2,
        yend = Read_Index + 0.2
      ),
      color = "red",
      size = 0.5
    ) +
    
    # Customize the plot appearance
    labs(
      title = "AAV Reads Alignment with Cut Sites",
      x = "AAV Genome Position",
      y = ""
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 10)
    ) +
    scale_x_continuous(
      limits = c(0, aav_length),
      breaks = seq(0, aav_length, by = 1000)
    )
}

# Function to generate summary statistics
generate_summary_stats <- function(merged_data) {
  summary_stats <- data.frame(
    total_reads = length(unique(merged_data$Read_Name)),
    unique_cut_sites = length(unique(merged_data$Cut_Site_Position)),
    mean_read_length = mean(merged_data$AAV_End.aav - merged_data$AAV_Start.aav),
    median_read_length = median(merged_data$AAV_End.aav - merged_data$AAV_Start.aav),
    unique_grna_ids = length(unique(merged_data$gRNA_id))
  )
  
  return(summary_stats)
}

# Main execution
main <- function() {
  tryCatch({
    # Load and process data
    cat("Loading and processing data...\n")
    data <- load_and_process_data(AAV_READS_PATH, CUT_SITES_PATH, AAV_REF_PATH)
    
    # Prepare data for plotting
    cat("Preparing plot data...\n")
    plot_data <- prepare_plot_data(data$aav_reads, data$cut_sites)
    
    # Generate summary statistics
    cat("\nGenerating summary statistics...\n")
    stats <- generate_summary_stats(plot_data)
    print(stats)
    
    # Create visualization
    cat("\nCreating visualization...\n")
    plot <- create_aav_alignment_plot(plot_data, data$aav_length)
    
    # Save the plot
    cat("\nSaving plot...\n")
    ggsave(
      "aav_alignment_plot.pdf",
      plot,
      width = 12,
      height = 8,
      units = "in",
      dpi = 300
    )
    
    # Display the plot
    print(plot)
    
    cat("\nAnalysis complete!\n")
    
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    print(str(e))
  })
}

# Run the analysis
main()