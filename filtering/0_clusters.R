library(tidyverse)
library(ggpubr)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <-read.csv("gatk_graph_deep.csv") %>% 
  mutate(context_id=paste(CHROM, POS, sep=":")) %>% 
  dplyr::select(-X)

test_10 <- dat %>%
  arrange(CHROM, POS) %>%  # Arrange the data by chromosome and position
  group_by(Strain, CHROM) %>%  # Group by strain and chromosome
  mutate(prev_pos = lag(POS),  # Create a column for previous position
         diff = ifelse(is.na(prev_pos), NA, POS - prev_pos),  # Calculate difference
         keep = ifelse(is.na(diff) | diff >= 10, TRUE, FALSE)) %>%  # Determine whether to keep the SNP
  filter(keep) %>%  # Filter out the SNPs to keep
  select(-prev_pos, -diff, -keep)

test <- dat %>%
  arrange(CHROM, POS) %>%  # Arrange the data by chromosome and position
  group_by(Strain, CHROM) %>%  # Group by strain and chromosome
  mutate(prev_pos = lag(POS),  # Create a column for previous position
         diff = ifelse(is.na(prev_pos), NA, POS - prev_pos),  # Calculate difference
         keep = ifelse(is.na(diff) | diff >= 100, TRUE, FALSE)) %>%  # Determine whether to keep the SNP
  filter(keep) %>%  # Filter out the SNPs to keep
  select(-prev_pos, -diff, -keep)

test <- dat %>%
  arrange(CHROM, POS) %>%  # Arrange the data by chromosome and position
  group_by(Strain, CHROM) %>%  # Group by strain and chromosome
  mutate(prev_pos = lag(POS),  # Create a column for previous position
         diff = ifelse(is.na(prev_pos), NA, POS - prev_pos),  # Calculate difference
         keep = ifelse(is.na(diff) | diff >= 150, TRUE, FALSE)) %>%  # Determine whether to keep the SNP
  filter(keep) %>%  # Filter out the SNPs to keep
  select(-prev_pos, -diff, -keep)


write.csv(test, "clusters_removed_100.csv")
write.csv(test, "clusters_removed_150.csv")
write.csv(test_10, "clusters_removed_10.csv")
