library(tidyverse)
library(ggpubr)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

context <- read.csv("data/mutation_context.csv")
dat <-read.csv("clusters_removed_10.csv") %>% 
  mutate(context_id=paste(CHROM, POS, sep=":")) %>% 
  dplyr::select(-X)

dat <- merge(dat, context, by.x="context_id", by.y="ID")

dat <- dat %>% 
  mutate(putative_illumina_homopolymer_strong=grepl("^A{10}\\.|^T{10}\\.|T{10}$|A{10}$", context_21),
         putative_illumina_homopolymer=grepl("^AAAAA\\.|^TTTTT\\.|TTTTT$|AAAAA$", context_11),
         putative_illumina_homopolymer_weak = grepl("^AA\\.|^TT\\.|TT$|AA$", context_5))


dat <- dat %>% filter(! CHROM %in% c("chrY", "chrX"))

dat <- dat %>%
  filter(!putative_illumina_homopolymer)

write.csv(dat, "mutations_homopolymers_removed.csv", row.names = F)