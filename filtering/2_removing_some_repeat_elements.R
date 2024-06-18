library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("mutations_homopolymers_removed.csv")

ltr <- read.table("data/ltr.txt")
simple <- read.table("data/simple_repeat.txt")
low <- read.table("data/low_complexity.txt")
sat <- read.table("data/satellite.txt")

repeatmask <- bind_rows(ltr, simple,low,sat)

rm(low,ltr,sat,simple)

names(repeatmask) <- c("CHROM", "START", "END", "NAME", "TYPE")

repeatmask <- repeatmask %>% filter(CHROM %in% dat$CHROM)

# Convert data frames to data.tables
setDT(dat)
setDT(repeatmask)

# Create an empty list to store the filtered results
filtered_results <- list()

# Loop through each unique chromosome in dat
for (chromosome in unique(dat$CHROM)) {
  # Subset dat and repeatmask for the current chromosome
  dat_subset <- dat[CHROM == chromosome]
  repeatmask_subset <- repeatmask[CHROM == chromosome]
  
  # Perform the filtering
  filtered_results[[chromosome]] <- dat_subset[!repeatmask_subset, on = .(CHROM, POS >= START, POS <= END)]
}

# Combine the filtered results into a single data.table
filtered_dat <- rbindlist(filtered_results)

filtered_dat

write.csv(filtered_dat, "mutations_repeats_homopolymers_removed.csv", row.names = F)
