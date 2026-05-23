library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(tidyr)

setwd("E:/AAV/Hela_Project_cleaned/AAV_integration_analysis")
data <- read.csv("CRISPR_paper/AAV_distribution_11_24.csv", stringsAsFactors = FALSE) |>
  drop_na(Host_Chromosome, Host_Start)
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Distance threshold (in base pairs)
# 300bp is quite small - typical promoter regions are 1-2kb upstream
# Consider 2000-5000bp for regulatory regions, or 10000bp for broader context
MAX_DISTANCE <- 5000

# Find nearby gene for each integration site
find_nearby_gene <- function(chr, pos) {
  # Create integration site as genomic range
  site <- GRanges(chr, IRanges(pos, pos))
  
  # Check if site overlaps any gene
  overlaps <- findOverlaps(site, genes)
  if (length(overlaps) > 0) {
    gene_ids <- genes[subjectHits(overlaps)]$gene_id
    symbols <- mapIds(org.Hs.eg.db, gene_ids, "SYMBOL", "ENTREZID", multiVals = "first")
    return(paste(na.omit(symbols), collapse = ";"))
  }
  
  # Find nearest downstream gene (higher position on same chromosome)
  downstream <- genes[seqnames(genes) == chr & start(genes) > pos]
  if (length(downstream) == 0) return(NA)
  
  # Get distance to nearest gene
  distances <- start(downstream) - pos
  min_distance <- min(distances)
  
  # Check if within threshold
  if (min_distance > MAX_DISTANCE) return(NA)
  
  nearest_gene <- downstream[which.min(distances)]
  symbol <- mapIds(org.Hs.eg.db, nearest_gene$gene_id, "SYMBOL", "ENTREZID")
  return(as.character(symbol))
}

# Build output dataframe
output <- data.frame(
  chromosome = data$Host_Chromosome,
  position = data$Host_Start,
  nearby_genes = mapply(find_nearby_gene, data$Host_Chromosome, data$Host_Start)
)

# Save results
write.csv(output, "CRISPR_paper/AAV_integration_nearby_genes.csv", row.names = FALSE)
cat("Complete! Analyzed", nrow(output), "integration sites\n")
cat("Using distance threshold:", MAX_DISTANCE, "bp\n")
cat("Sites with nearby genes:", sum(!is.na(output$nearby_genes)), "\n")