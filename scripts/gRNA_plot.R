library(dplyr)
library(ggplot2)
library(Biostrings)
library(readr)
library(RColorBrewer)  

AAV_REF_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/pssAAV-CB-EGFP_ARM.fa"
AAV_READS_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/output/aav_reads_locations.csv" 

# CUT_SITES_PATH_unfiltered <- "D:/Jiahe/IU/AAV/HeLa_project/output/detailed_cut_sites_with_gRNA.csv"
# 
# all_reads <- read.csv(AAV_READS_PATH)
# all_reads <- all_reads |>
#   filter(!is.na(Host_Chromosome))
# 
# df <- read.csv(CUT_SITES_PATH_unfiltered)
# df <- df |>
#   filter(Read_Name %in% all_reads$Read_Name)
# 
# write.csv(df, "D:/Jiahe/IU/AAV/HeLa_project/output/detailed_cut_sites_with_gRNA_filtered.csv")

CUT_SITES_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/output/detailed_cut_sites_with_gRNA_filtered.csv"

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

create_read_alignment_plot <- function(all_reads, filtered_reads, aav_length) {
  # Get unique Read_Names from filtered_reads
  unique_filtered_reads <- filtered_reads %>%
    distinct(Read_Name)
  
  # Filter reads that are in the filtered dataset, keeping only unique reads
  plot_reads <- all_reads %>%
    filter(Read_Name %in% unique_filtered_reads$Read_Name) %>%
    # Add a y-position for each read (stacking them vertically)
    group_by(Read_Name) %>%
    mutate(
      y_position = cur_group_id() * -1  # Negative to plot below the reference line
    )
  
  # Create the reference genome data
  reference_genome <- data.frame(
    x = c(1, aav_length),
    y = c(0, 0)
  )
  
  # Create the plot
  ggplot() +
    # Add the reference genome line
    geom_line(
      data = reference_genome,
      aes(x = x, y = y),
      color = "black",
      size = 1
    ) +
    # Add read alignment lines
    geom_segment(
      data = plot_reads,
      aes(
        x = AAV_Start,
        xend = AAV_End,
        y = y_position,
        yend = y_position
      ),
      color = "lightblue",
      size = 0.5,
      alpha = 0.6
    ) +
    # Customize the appearance
    labs(
      title = "AAV Read Alignments",
      subtitle = paste("Total unique reads:", nrow(unique_filtered_reads)),
      x = "AAV Genome Position",
      y = "Read Number"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_x_continuous(
      limits = c(0, aav_length),
      breaks = seq(0, aav_length, by = 1000)
    ) +
    # Hide y-axis labels since they're just read indices
    scale_y_continuous(labels = NULL)
}

# Update the main function to include the new visualization
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
    
    # Create histogram visualization
    cat("\nCreating histogram visualization...\n")
    histogram_plot <- create_aav_alignment_plot(plot_data, data$aav_length)
    
    # Create read alignment visualization
    cat("\nCreating read alignment visualization...\n")
    alignment_plot <- create_read_alignment_plot(all_reads, df, data$aav_length)
    
    # Save the plots
    cat("\nSaving plots...\n")
    ggsave(
      "aav_histogram.pdf",
      histogram_plot,
      width = 12,
      height = 8,
      units = "in",
      dpi = 300
    )
    
    ggsave(
      "aav_read_alignment.pdf",
      alignment_plot,
      width = 12,
      height = 8,
      units = "in",
      dpi = 300
    )
    
    # Display the plots
    print(histogram_plot)
    print(alignment_plot)
    
    cat("\nAnalysis complete!\n")
    
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    print(str(e))
  })
}

main()