library(tidyverse)
library(ggpubr)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("removing_regions_overlapping_lumpy_sv.csv")

####Filtering problematic megabase regions####
bins <- seq(0, max(dat$POS) + 1000000, by = 1000000)

filtered_summary <- dat %>% 
  group_by(CHROM) %>% 
  mutate(bin = Hmisc::cut2(as.integer(POS), bins, right = FALSE)) %>%
  group_by(CHROM, bin) %>%
  summarize(count = n(),
            mean_qual_dv = min(QUAL_dv)) %>% 
  arrange(-count)

ggplot(filtered_summary, aes(x=count, y=mean_qual_dv)) + 
  geom_point()
 
dat <- dat %>%
  mutate(bin = Hmisc::cut2(as.integer(POS), bins, right = FALSE)) %>%
  merge(filtered_summary, by = c("bin", "CHROM"), all.x=T) %>% 
  filter(count <= quantile(filtered_summary$count, probs=0.95)[1])%>% 
  select(-bin, -count)

####Filtering problematic 6.3 kilobase regions (SINE)####
bins <- seq(0, max(dat$POS) + 6300, by = 6300)

filtered_summary <- dat %>% 
  select(CHROM, POS) %>%
  group_by(CHROM) %>% 
  mutate(bin = Hmisc::cut2(as.integer(POS), bins, right = FALSE)) %>%
  group_by(CHROM, bin) %>%
  summarize(count = n()) %>% 
  arrange(-count)

ggplot(filtered_summary, aes(x=count, 
                             fill=ifelse(count>quantile(filtered_summary$count, probs=0.95)[1], 
                                         CHROM,NA))) + 
  geom_histogram() + scale_y_log10() +
  xlab("Number of DNM / 6.3 kb") +
  labs (fill="CHROM")

dat <- dat %>%
  mutate(bin = Hmisc::cut2(as.integer(POS), bins, right = FALSE)) %>%
  merge(filtered_summary, by = c("bin", "CHROM"), all.x=T) %>% 
  filter(count <= quantile(filtered_summary$count, probs=0.95)[1]) %>% 
  select(-bin, -count)

####Filtering problematic 200 bp (TE)####
bins <- seq(0, max(dat$POS) + 200, by = 200)

filtered_summary <- dat %>% 
  select(CHROM, POS) %>%
  group_by(CHROM) %>% 
  mutate(bin = Hmisc::cut2(as.integer(POS), bins, right = FALSE)) %>%
  group_by(CHROM, bin) %>%
  summarize(count = n()) %>% 
  arrange(-count)

dat <- dat %>%
  mutate(bin = Hmisc::cut2(as.integer(POS), bins, right = FALSE)) %>%
  merge(filtered_summary, by = c("bin", "CHROM"), all.x=T) %>% 
  filter(count <= quantile(filtered_summary$count, probs=0.95)[1])%>% 
  select(-bin, -count)

write.csv(dat, "removing_problematic_mB.csv")
