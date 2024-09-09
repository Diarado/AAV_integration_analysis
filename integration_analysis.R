# Required Libraries
library(ShortRead)
library(parallel)
library(GenomicRanges)
library(tidyverse)
library(yaml)

options(stringsAsFactors = FALSE)

configFile <- "/mnt/d/Jiahe/IU/AAV/HeLa_project/config.yaml"
config <- read_yaml(configFile)

# Create Output Directory
dir.create(config$outputDir, showWarnings = FALSE)

# Start log file
logFile <- file.path(config$outputDir, 'logs', 'log')
dir.create(file.path(config$outputDir, 'logs'), showWarnings = FALSE)  # Ensure logs dir exists before writing log
write(date(), file = logFile)

# Read in the sample configuration (sampleConfigFile contains the input BAM files)
samples <- read_delim(config$sampleConfigFile, delim = ',', col_names = TRUE, col_types = cols())

# Set up parallel processing
cluster <- makeCluster(config$alignment.CPUs)
clusterExport(cluster, c('config', 'samples', 'logFile'))  # Export logFile to parallel environment

# Align the AAV-Aligned Reads to the Human Genome
alignReads <- function(sampleRow) {
  # Load necessary libraries
  library(ShortRead)
  library(tidyverse)
  library(GenomicRanges)
  
  # Define the BAM File and Reference Genome based on the sampleRow
  bamFile <- file.path("barcode13", paste0("PAW40361_pass_barcode13_04266714_eb72dc1a_", sampleRow$sample, ".bam"))
  refGenome <- sampleRow$refGenomeBLATdb
  
  # Debug: Print the sample being processed
  print(paste("Processing sample:", sampleRow$sample))
  print(paste("BAM File:", bamFile))
  print(paste("Reference Genome:", refGenome))
  
  # Define the output file paths
  alignedSamFile <- file.path(config$outputDir, paste0(sampleRow$sample, "_aligned.sam"))
  sortedBamFile <- file.path(config$outputDir, paste0(sampleRow$sample, "_sorted.bam"))
  
  # Debug: Print the output file paths
  print(paste("Aligned SAM File:", alignedSamFile))
  print(paste("Sorted BAM File:", sortedBamFile))
  
  # BWA Command (if needed for realignment, otherwise skip this)
  bwaCmd <- paste("bwa mem", refGenome, bamFile, ">", alignedSamFile)
  print(paste("Running BWA Command:", bwaCmd))
  system(bwaCmd)
  
  # SAM to BAM and Sort Command
  samtoolsViewSortCmd <- paste("samtools view -Sb", alignedSamFile, "| samtools sort -o", sortedBamFile)
  print(paste("Running SAM to BAM and Sort Command:", samtoolsViewSortCmd))
  system(samtoolsViewSortCmd)
  
  # Log completion
  write(paste("Completed alignment for sample:", sampleRow$sample), file = logFile, append = TRUE)
  
  # Debug: Print completion
  print(paste("Completed alignment for sample:", sampleRow$sample))
}

# Run alignment for each sample sequentially for debugging
# For parallel execution, uncomment the following line:
# invisible(parLapply(cluster, split(samples, 1:nrow(samples)), alignReads))

# Sequential debugging for one sample
alignReads(samples[1, ])

# Stop the cluster
stopCluster(cluster)

# Process aligned reads to identify integration sites
identifyIntegrationSites <- function(sampleRow) {
  # Load the aligned BAM file
  bamFile <- file.path("barcode13", 
                       paste0("PAW40361_pass_barcode13_04266714_eb72dc1a_", 
                              sampleRow$sample, ".bam"))
  print(paste0("PAW40361_pass_barcode13_04266714_eb72dc1a_", 
               sampleRow$sample, ".bam"))
  # Use GenomicRanges to identify integration sites
  reads <- readGAlignments(bamFile)
  integrationSites <- reduce(granges(reads))
  
  # Save integration sites
  save(integrationSites, file = file.path(config$outputDir, paste0(sampleRow$sample, "_integration_sites.RData")))
  
  # Log completion
  write(paste("Completed integration site identification for sample:", sampleRow$sample), file = logFile, append = TRUE)
}

# Run integration site identification for each sample sequentially for debugging
identifyIntegrationSites(samples[1, ])

write("Integration site analysis completed.", file = logFile, append = TRUE)

