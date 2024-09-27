# check whether barcode 13 or 14 has more integrations
# Read in the data
data <- read.csv("D:/Jiahe/IU/AAV/HeLa_project/output/aav_reads_locations.csv", stringsAsFactors = FALSE)

# Filter for barcode13 (where '13' appears in the specified position)
barcode13_cnt <- data |>
  filter(!is.na(Host_Chromosome) & substr(Sample, 22, 23) == "13")

# Filter for barcode14 (where '14' appears in the specified position)
barcode14_cnt <- data |>
  filter(!is.na(Host_Chromosome) & substr(Sample, 22, 23) == "14")



# find gRNA seq