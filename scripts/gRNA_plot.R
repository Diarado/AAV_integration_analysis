library(dplyr)
library(ggplot2)
library(Biostrings)
library(readr)
library(RColorBrewer)  

AAV_REF_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/pssAAV-CB-EGFP_ARM.fa"
AAV_READS_PATH <- "D:/Jiahe/IU/AAV/HeLa_project/output/aav_reads_locations.csv" 

CUT_SITES_PATH_unfiltered <- "D:/Jiahe/IU/AAV/HeLa_project/output/detailed_cut_sites_with_gRNA.csv"

all_reads <- read.csv(AAV_READS_PATH)
all_reads <- all_reads |>
  filter(!is.na(Host_Chromosome))

# this is filtered_reads
df <- read.csv(CUT_SITES_PATH_unfiltered)
df <- df |>
  filter(Read_Name %in% all_reads$Read_Name)

write.csv(df, "D:/Jiahe/IU/AAV/HeLa_project/output/detailed_cut_sites_with_gRNA_filtered.csv")

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
  # Read and extract the AAV sequence
  aav_seq <- as.character(readDNAStringSet(AAV_REF_PATH)[[1]])
  # Extract the relevant region (1281-1453) and convert to lowercase
  itr_seq <- tolower(substr(aav_seq, 1281, 1453))
  
  # Create sequence labels (10 bases per label)
  seq_labels <- strsplit(itr_seq, "")[[1]]
  seq_df <- data.frame(
    x = 1281:1453,
    label = seq_labels,
    stringsAsFactors = FALSE
  )
  
  # Get unique Read_Names from filtered_reads
  unique_filtered_reads <- filtered_reads %>%
    distinct(Read_Name)
  
  # Filter reads that are in the filtered dataset, keeping only unique reads
  # print(unique_filtered_reads)
  # print(all_reads)
  plot_reads <- all_reads %>%
    filter(AAV_Start >= 1281 & AAV_End <= 1453) %>%
    filter(Read_Name %in% unique_filtered_reads$Read_Name) %>%
    arrange(Host_Start) %>%    
    mutate(
      y_position = row_number()  
    )
  print(plot_reads)
  # Calculate y-range for the plot
  max_y <- max(plot_reads$y_position)
  
  # Define ITR feature coordinates and colors - ensure each feature has a unique color
  itr_features <- list(
    list(start = 1281, end = 1453, name = "5' ITR", color = "#FFE4B5", layer = 1),
    list(start = 1285.5, end = 1301.5, name = "RBE 5' End", color = "#FFD700", layer = 2),
    list(start = 1310.5, end = 1317.5, name = "C'", color = "#DDA0DD", layer = 2),
    list(start = 1322.5, end = 1329.5, name = "C", color = "#9370DB", layer = 2),
    list(start = 1332.5, end = 1339.5, name = "B", color = "#90EE90", layer = 2),
    list(start = 1339.5, end = 1344.5, name = "RBE'", color = "#FFFFE0", layer = 2),
    list(start = 1344.5, end = 1351.5, name = "B'", color = "#98FB98", layer = 2),
    list(start = 1360.5, end = 1376.5, name = "RBE 3' End", color = "#DAA520", layer = 2),
    list(start = 1376.5, end = 1393.5, name = "A", color = "#ADD8E6", layer = 2),
    list(start = 1391.5, end = 1393.5, name = "Trs", color = "#FF6B6B", layer = 3),
    list(start = 1393.5, end = 1413.5, name = "D'", color = "orange", layer = 2)
  )
  
  # Create feature rectangles data frame
  feature_rects <- do.call(rbind, lapply(itr_features, function(f) {
    data.frame(
      xmin = f$start,
      xmax = f$end,
      ymin = -0.15 * max_y + (f$layer - 1) * -0.05 * max_y,
      ymax = -0.10 * max_y + (f$layer - 1) * -0.05 * max_y,
      name = f$name,
      color = f$color,
      stringsAsFactors = FALSE
    )
  }))
  
  # Create color vector that matches the features exactly
  feature_colors <- setNames(
    sapply(itr_features, function(x) x$color),
    sapply(itr_features, function(x) x$name)
  )
  
  # Create the plot
  p <- ggplot() +
    # Add read alignment lines
    geom_segment(
      data = plot_reads,
      aes(
        x = AAV_Start,
        xend = AAV_End,
        y = y_position,
        yend = y_position
      ),
      color = "blue",
      size = 1,
      alpha = 0.6
    ) +
    scale_y_discrete() + 
    # Add sequence labels
    geom_text(
      data = seq_df,
      aes(x = x, y = -0.05 * max_y, label = label),
      size = 2.5,
      family = "mono"
    ) +
    # Add feature rectangles
    geom_rect(
      data = feature_rects,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = name),
      alpha = 0.7
    ) +
    # Add feature labels
    geom_text(
      data = feature_rects,
      aes(x = (xmin + xmax)/2, y = ymin - 0.02 * max_y, label = name),
      size = 2.5,
      angle = 45,
      hjust = 0
    ) +
    # Customize the appearance with exact color matching
    scale_fill_manual(values = feature_colors) +
    labs(
      title = "5' ITR gRNA-cut AAV Read Alignments",
      x = "AAV Genome Position",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    ) +
    scale_x_continuous(
      limits = c(1281, 1453),
      breaks = seq(1280, 1460, by = 20)
    ) +
    scale_y_continuous(
      limits = c(-0.35 * max_y, max_y),
      labels = NULL
    )
  
  return(p)
}

# Update the main function to avoid opening new devices
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
    print(histogram_plot)
    
    # Create read alignment visualization
    cat("\nCreating read alignment visualization...\n")
    alignment_plot <- create_read_alignment_plot(all_reads, df, data$aav_length)
    print(alignment_plot)
    
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
    
    cat("\nAnalysis complete!\n")
    
  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    print(str(e))
  })
}

main()