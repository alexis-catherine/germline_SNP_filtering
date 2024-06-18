library(tidyverse)
library(ggpubr)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.table("data/MGI_Features.txt", sep="\t", head=T) %>% 
  filter(Build=="GRCm39",
         Chr %in% c(1:19,"X","Y"),
         Start != "syntenic",
         !grepl("cM", Start)) %>% 
  select(Chr, Start, End) %>% 
  mutate(Start = as.integer(Start)-100) %>% 
  mutate(End = as.integer(End)+100) %>% 
  distinct()

write.csv(dat, "data/MGI_olf_mask.csv")
