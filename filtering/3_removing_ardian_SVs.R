library(tidyverse)
library(ggpubr)
library(data.table)
setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("mutations_repeats_homopolymers_removed.csv")
ardian <-bind_rows(read.table("data/combined_founder_and_cc_insdel.vcf", header = F), 
                   read.table("data/combined_founder_and_cc_inv.vcf", header = F))

ardian_mult <- filter(ardian,
                     grepl(";", V3)) %>% 
  separate(V3, sep=";", into=c("var1", "var2", "var3", "var4", "var5")) %>% 
  select(starts_with("var")) %>%
  pivot_longer(1:last_col()) %>% 
  select(-name) %>% 
  separate(value, sep="-", into=c("CHROM_ID", "START", "TYPE", "LEN"), convert = T) %>% 
  filter(!is.na(LEN))
  

ardian <- ardian %>% 
  select(CHROM=V1, ID=V3) %>% 
  filter(!grepl(";", ID)) %>% 
  separate(ID, sep="-", into=c("CHROM_ID", "START", "TYPE", "LEN"), convert = T) %>%
  mutate( END = START + LEN) %>% 
  filter(!CHROM %in% c("chrX")) %>% 
  bind_rows(ardian_mult)

ardian_buffer <- ardian %>% 
  mutate(START = START - 10,
         END = END + 10) %>% distinct()

# Convert data frames to data.tables
setDT(dat)
setDT(ardian_buffer)

# Create an empty list to store the filtered results
filtered_results <- list()

# Loop through each unique chromosome in dat
for (chromosome in unique(dat$CHROM)) {
  # Subset dat and ardian_buffer for the current chromosome
  dat_subset <- dat[CHROM == chromosome]
  ardian_subset <- ardian_buffer[CHROM == chromosome]
  
  # Perform the filtering
  filtered_results[[chromosome]] <- dat_subset[!ardian_subset, on = .(CHROM, POS >= START, POS <= END)]
}

# Combine the filtered results into a single data.table
filtered_dat <- rbindlist(filtered_results)

write.csv(filtered_dat, "removing_regions_overlapping_founder_SVs.csv", row.names = F)

nrow(dat) - nrow(filtered_dat)
