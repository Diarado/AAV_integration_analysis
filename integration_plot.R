# Load necessary libraries
library(ggplot2)
library(scales)

# Read chromosome information
chromosome_info <- read.table("D:/Jiahe/IU/AAV/HeLa_project/hg38.fa.fai", header = FALSE, stringsAsFactors = FALSE)
colnames(chromosome_info) <- c("chrom", "length", "offset", "line_bases", "line_bytes")

# Read integration data
data <- read.csv("D:/Jiahe/IU/AAV/HeLa_project/output/aav_reads_locations_test.csv", stringsAsFactors = FALSE)

# Filter out entries with NA in Host_Chromosome or Host_Start
integration_data <- subset(data, !is.na(Host_Chromosome) & !is.na(Host_Start))

# Keep only standard chromosomes (chr1 to chr22, chrX, chrY)
standard_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
chromosome_info <- subset(chromosome_info, chrom %in% standard_chromosomes)
integration_data <- subset(integration_data, Host_Chromosome %in% standard_chromosomes)

# Order chromosomes
chromosome_info$chrom <- factor(chromosome_info$chrom, levels = standard_chromosomes)
integration_data$Host_Chromosome <- factor(integration_data$Host_Chromosome, levels = standard_chromosomes)

# Assign numeric y positions to chromosomes
chromosome_info$y_pos <- as.numeric(chromosome_info$chrom)
integration_data$y_pos <- as.numeric(integration_data$Host_Chromosome)

# Prepare data for plotting
chromosomes_df <- data.frame(
  chrom = chromosome_info$chrom,
  start = 0,
  end = chromosome_info$length,
  y_pos = chromosome_info$y_pos
)

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
  scale_x_continuous(labels = comma) +
  scale_y_continuous(breaks = chromosomes_df$y_pos, labels = chromosomes_df$chrom) +
  theme_minimal() +
  xlab("Position (bp)") +
  ylab("Chromosome") +
  ggtitle("AAV Integration Sites on Human Chromosomes")
