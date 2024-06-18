library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("removing_regions_overlapping_founder_SVs.csv")

olf <- read.csv("data/MGI_olf_mask.csv") %>%
  select(-X) %>%
  mutate(Chr = paste0("chr", Chr))

names(olf) <- c("CHROM", "START", "END")

olf <- olf %>% filter(CHROM %in% dat$CHROM)

# Convert data frames to data.tables
setDT(dat)
setDT(olf)

# Create an empty list to store the filtered results
filtered_results <- list()

# Loop through each unique chromosome in dat
for (chromosome in unique(dat$CHROM)) {
  # Subset dat and olf for the current chromosome
  dat_subset <- dat[CHROM == chromosome]
  olf_subset <- olf[CHROM == chromosome]
  
  # Perform the filtering
  filtered_results[[chromosome]] <- dat_subset[!olf_subset, on = .(CHROM, POS >= START, POS <= END)]
}

# Combine the filtered results into a single data.table
filtered_dat <- rbindlist(filtered_results)

filtered_dat
write.csv(filtered_dat, "mutations_overlap_olf_removed.csv", row.names = F)
