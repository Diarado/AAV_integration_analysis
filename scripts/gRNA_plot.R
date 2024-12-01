library(dplyr)
library(ggplot2)
library(Biostrings)
library(readr)
library(RColorBrewer)  

AAV_REF_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/pssAAV-CB-EGFP_ARM.fa"
AAV_READS_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/output/aav_reads_locations.csv" 
CUT_SITES_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/output/detailed_cut_sites_with_gRNA.csv"

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
    distinct(Read_Name, .keep_all = TRUE) %>%
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

prepare_plot_data <- function(aav_reads, cut_sites) {
  # Merge the data
  merged_data <- aav_reads %>%
    inner_join(cut_sites, by = c("Read_Name", "Sample"), suffix = c(".aav", ".cut"))
  
  # Get total number of unique reads
  total_unique_reads <- length(unique(merged_data$Read_Name))
  
  # Calculate position-specific read counts and percentages
  position_data <- merged_data %>%
    # Ensure we're counting unique reads
    distinct(Read_Name, AAV_Start.aav, AAV_End.aav, Cut_Site_Position, gRNA_id) %>%
    # Calculate absolute cut position
    mutate(Absolute_Cut_Position = AAV_Start.aav + Cut_Site_Position - 1) %>%
    # Group by position and gRNA
    group_by(Absolute_Cut_Position, gRNA_id) %>%
    # Count reads at each position
    summarise(
      read_count = n(),
      percentage = (n() / total_unique_reads) * 100,
      .groups = 'drop'
    ) %>%
    # Convert gRNA_id to factor for proper color mapping
    mutate(gRNA_id = factor(gRNA_id))
  
  return(position_data)
}

create_aav_alignment_plot <- function(position_data, aav_length) {
  # Get number of unique gRNA IDs
  n_grnas <- length(unique(position_data$gRNA_id))
  
  # Create color palette
  if(n_grnas <= 8) {
    colors <- brewer.pal(max(n_grnas, 3), "Set2")
  } else {
    colors <- c(
      brewer.pal(8, "Set2"),
      brewer.pal(8, "Set1"),
      brewer.pal(8, "Dark2")
    )[1:n_grnas]
  }
  
  # Create the plot
  ggplot(position_data, aes(x = Absolute_Cut_Position, y = percentage)) +
    # Add bars
    geom_bar(stat = "identity", position = "stack", width = 50) +  # Adjust width as needed
    
    # Customize the appearance
    scale_fill_manual(
      values = colors,
      name = "gRNA ID"
    ) +
    
    labs(
      title = "AAV Cut Site Distribution",
      subtitle = paste("Percentage of reads at each position", sep=""),
      x = "AAV Genome Position",
      y = "On-Target reads (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(
      limits = c(0, aav_length),
      breaks = seq(0, aav_length, by = 1000)
    ) +
    scale_y_continuous(
      limits = c(0, NA),  # Start at 0, auto-scale upper limit
      expand = expansion(mult = c(0, 0.05))  # Remove padding at bottom, add small padding at top
    )
}

generate_summary_stats <- function(position_data) {
  # Calculate gRNA-specific statistics
  grna_stats <- position_data %>%
    group_by(gRNA_id) %>%
    summarise(
      total_reads = sum(read_count),
      mean_position = weighted.mean(Absolute_Cut_Position, read_count),
      peak_position = Absolute_Cut_Position[which.max(percentage)],
      peak_percentage = max(percentage)
    )
  
  # Calculate overall statistics
  overall_stats <- data.frame(
    total_positions = length(unique(position_data$Absolute_Cut_Position)),
    total_reads = sum(position_data$read_count),
    max_percentage = max(position_data$percentage),
    unique_grna_ids = length(unique(position_data$gRNA_id))
  )
  
  return(list(
    overall_stats = overall_stats,
    grna_stats = grna_stats
  ))
}

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
    cat("\nOverall Statistics:\n")
    print(stats$overall_stats)
    cat("\ngRNA-specific Statistics:\n")
    print(stats$grna_stats)
    
    # Create visualization
    cat("\nCreating visualization...\n")
    plot <- create_aav_alignment_plot(plot_data, data$aav_length)
    
    # Save the plot
    cat("\nSaving plot...\n")
    ggsave(
      "aav_alignment_histogram.pdf",
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
unqiue_reads <- unique(all_results$Read_Name)
