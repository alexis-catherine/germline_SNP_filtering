library(tidyverse)  
library(ggpubr)
library(data.table)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("mutations_hard_filtering.csv")
vcf <- read.table("data/truvari_merge.vcf", header = T) %>% 
  filter(CHROM %in% paste0("chr", 1:19)) %>%
  pivot_longer(10:last_col(), names_to="STRAIN") %>% 
  filter(grepl("^1/1|^0/1", value)) %>% 
  mutate(TYPE = gsub(".*SVTYPE=", "", INFO) %>% 
           gsub(";.*", "", .),
         LEN = gsub(".*SVLEN=", "", INFO) %>% 
           gsub(";.*", "", .),
         STRAIN=gsub("_.*","",STRAIN))

filtered_by_strain <- NULL

for(i in unique(dat$Strain)){
  strain_vcf <- filter(vcf, STRAIN == i)
  strain_dat <- filter(dat, Strain == i)
  

  strain_vcf <- strain_vcf %>% 
    mutate(LEN=ifelse(grepl("BND", LEN), 0, LEN)) %>% 
    select(CHROM,
           START=POS,
           TYPE,
           LEN) %>% 
    mutate(LEN = as.numeric(LEN) %>% abs()) %>%
    mutate(END = START + LEN)%>% 
    mutate(START = START - 500,
           END = END + 500) 
  
  strain_vcf_buffer <- strain_vcf %>%  distinct()
  
  # Convert data frames to data.tables
  setDT(strain_dat)
  setDT(strain_vcf_buffer)
  
  # Create an empty list to store the filtered results
  filtered_results <- list()
  
  # Loop through each unique chromosome in dat
  for (chromosome in unique(strain_dat$CHROM)) {
    # Subset dat and ardian_buffer for the current chromosome
    strain_dat_subset <- strain_dat[CHROM == chromosome]
    strain_vcf_subset <- strain_vcf_buffer[CHROM == chromosome]
    
    # Perform the filtering
    filtered_results[[chromosome]] <- strain_dat_subset[!strain_vcf_subset, on = .(CHROM, POS >= START, POS <= END)]
  }
  
  # Combine the filtered results into a single data.table
  filtered_strain_dat <- rbindlist(filtered_results)
  filtered_by_strain <- rbind(filtered_by_strain, filtered_strain_dat)
}



write.csv(filtered_by_strain, "removing_regions_overlapping_lumpy_sv.csv", row.names = F)
