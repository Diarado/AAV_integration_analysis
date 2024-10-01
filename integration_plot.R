library(ggplot2)
library(scales)

# Read chromosome information
chromosome_info <- read.table("D:/Jiahe/IU/AAV/HeLa_project/hg38.fa.fai", header = FALSE, stringsAsFactors = FALSE)
colnames(chromosome_info) <- c("chrom", "length", "offset", "line_bases", "line_bytes")

# Read integration data
data <- read.csv("D:/Jiahe/IU/AAV/HeLa_project/output/aav_reads_locations.csv", stringsAsFactors = FALSE)

# Filter out entries with NA in Host_Chromosome or Host_Start
integration_data <- integration_data <- data |>
  filter(
    !is.na(Host_Chromosome),
    !is.na(Host_Start),
    (AAV_Start <= 2000 & 1269 <= AAV_End) &
    (AAV_Start < 1413 & 1413 < AAV_End)
      
  )

# Define the main chromosomes
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y")) # Include all variations and alts

# Function to extract the main chromosome name
extract_main_chrom <- function(chrom_name) {
  match <- regexpr("^chr[0-9XY]+", chrom_name)
  if (match[1] != -1) {
    substr(chrom_name, match[1], attr(match, "match.length"))
  } else {
    NA
  }
}

# Add 'chrom_main' to chromosome_info and integration_data
chromosome_info$chrom_main <- sapply(chromosome_info$chrom, extract_main_chrom)
integration_data$chrom_main <- sapply(integration_data$Host_Chromosome, extract_main_chrom)

# Filter chromosomes to include only those with chrom_main in standard chromosomes
chromosome_info <- subset(chromosome_info, chrom_main %in% standard_chromosomes)
integration_data <- subset(integration_data, chrom_main %in% standard_chromosomes)

# Order chromosomes
chromosome_info$chrom_main <- factor(chromosome_info$chrom_main, levels = standard_chromosomes)
integration_data$chrom_main <- factor(integration_data$chrom_main, levels = standard_chromosomes)

# Assign numeric y positions to chromosomes
chromosome_info$y_pos <- as.numeric(chromosome_info$chrom_main)
integration_data$y_pos <- as.numeric(integration_data$chrom_main)

# Prepare data for plotting
# Aggregate lengths by chrom_main (take maximum length)
chromosomes_df <- aggregate(length ~ chrom_main + y_pos, data = chromosome_info, max)
chromosomes_df$start <- 0
chromosomes_df$end <- chromosomes_df$length

# Create the plot
ggplot() +
  # Plot chromosomes as horizontal bars
  geom_segment(data = chromosomes_df,
               aes(x = start, xend = end, y = y_pos, yend = y_pos),
               size = 5, color = "grey") +
  # Plot integration sites as vertical lines
  geom_segment(data = integration_data,
               aes(x = Host_Start, xend = Host_Start, y = y_pos - 0.2, yend = y_pos + 0.2),
               color = "red", size = 0.5) +
  # Adjust scales and labels
  scale_x_continuous(
    labels = scales::comma_format(scale = 1/1000), 
    name = "Position (kbp)"                          
  ) +
  scale_y_continuous(breaks = chromosomes_df$y_pos, labels = chromosomes_df$chrom_main) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),           
    panel.background = element_blank(),     
    axis.ticks.y = element_blank(),         
    axis.text.y = element_text(size = 10)   
  ) +
  ylab("Chromosome") +
  ggtitle("AAV Integration Sites in Human Chromosomes (Those Includes Both ITR and CB)")
