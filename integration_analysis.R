# Required Libraries
library(ShortRead)
library(parallel)
library(GenomicRanges)
library(tidyverse)
library(yaml)
library(GenomicAlignments)

options(stringsAsFactors = FALSE)

# Load the config file
configFile <- "/mnt/d/Jiahe/IU/AAV/HeLa_project/config.yaml"
config <- read_yaml(configFile)

# Debug: Print the config file path and its content
print("Config file loaded successfully:")
print(config)

# Create Output Directory
print("Creating output directory...")
dir.create(config$outputDir, showWarnings = FALSE)

# Start log file
logFile <- file.path(config$outputDir, 'logs', 'log')
dir.create(file.path(config$outputDir, 'logs'), showWarnings = FALSE)  # Ensure logs dir exists before writing log
write(date(), file = logFile)

# Read in the sample configuration (sampleConfigFile contains the input BAM files)
print("Loading sample configuration from:")
print(config$sampleConfigFile)
samples <- read_delim(config$sampleConfigFile, delim = ',', col_names = TRUE, col_types = cols())

# Debug: Show the structure of the loaded sample configuration
print("Sample configuration loaded:")
print(samples)

# Load gRNA sequences (gRNAs.csv)
print("Loading gRNA sequences from gRNAs.csv...")
gRNAs <- read_delim("/mnt/d/Jiahe/IU/AAV/HeLa_project/gRNAs.csv", delim = ',', col_names = TRUE)
print("gRNA sequences loaded:")
print(gRNAs)

# Load AAV reference sequence (pssAAV-CB-EGFP ARM.fa)
print("Loading AAV reference sequence from pssAAV-CB-EGFP ARM.fa...")
AAV_sequence <- readDNAStringSet("/mnt/d/Jiahe/IU/AAV/HeLa_project/pssAAV-CB-EGFP ARM.fa")
print("AAV reference sequence loaded.")

# Load barcode information (barcodes.csv)
print("Loading barcode information from barcodes.csv...")
barcodes <- read_delim("/mnt/d/Jiahe/IU/AAV/HeLa_project/barcodes.csv", delim = ',', col_names = TRUE)
print("Barcode information loaded.")
print(barcodes)

# Set up parallel processing (DISABLED for debugging)
#print("Setting up parallel processing...")
#cluster <- makeCluster(config$alignment.CPUs)
#clusterExport(cluster, c('config', 'samples', 'gRNAs', 'AAV_sequence', 'barcodes', 'logFile'))

# Function to process aligned reads and identify integration sites
identifyIntegrationSites <- function(sampleRow) {
  print(paste("Processing sample row:", sampleRow$sample))  # Debug: Print sample details
  
  bamFile <- file.path("barcode13", paste0("PAW40361_pass_barcode13_04266714_eb72dc1a_", sampleRow$sample, ".bam"))
  print(paste("Processing BAM file:", bamFile))
  
  if (!file.exists(bamFile)) {
    stop(paste("BAM file not found:", bamFile))
  }
  
  print("Loading aligned reads from BAM file...")
  tryCatch({
    reads <- readGAlignments(bamFile)
    
    if (length(reads) == 0) {
      stop("No alignments found in BAM file")
    }
    
    print(paste("Number of reads loaded:", length(reads)))
    
    print("Identifying integration sites...")
    integrationSites <- GenomicRanges::reduce(granges(reads))
    
    print(paste("Number of integration sites identified:", length(integrationSites)))
    rDataFile <- file.path(config$outputDir, paste0(sampleRow$sample, "_integration_sites.RData"))
    print(paste("Saving integration sites to:", rDataFile))
    save(integrationSites, file = rDataFile)
    
    if (file.exists(rDataFile)) {
      print(paste("RData file saved successfully:", rDataFile))
    } else {
      stop("Failed to save the RData file")
    }
    
    write(paste("Completed integration site identification for sample:", sampleRow$sample), file = logFile, append = TRUE)
    print(paste("Integration site identification completed for sample:", sampleRow$sample))
    
  }, error = function(e) {
    print(paste("Error during integration site identification:", e$message))
  })
}

# Running in a loop (sequential for debugging)
print("Starting integration site analysis for each sample...")
for (i in 1:nrow(samples)) {
  print(paste("Running integration for sample:", samples$sample[i]))  # Debug: Notify function execution
  identifyIntegrationSites(samples[i, ])
}

# Stop cluster (DISABLED for now)
# print("Stopping the cluster...")
# stopCluster(cluster)

# Write final completion message to log file
write("Integration site analysis completed.", file = logFile, append = TRUE)
print("Integration site analysis completed.")
