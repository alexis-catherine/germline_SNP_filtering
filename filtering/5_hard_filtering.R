library(tidyverse)
library(ggpubr)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("mutations_overlap_olf_removed.csv")

high_depth_gk <- quantile(dat$DP_gk, probs = 0.99)
high_depth_dv <- quantile(dat$DP_dv, probs = 0.99)
high_depth_gt <- quantile(dat$DP_gt, probs = 0.99)

test <- dat %>%
  mutate(AA_score =gsub(";.*", "", INFO_gt) %>% 
           gsub("AAScore=", "", .) %>% as.numeric()) %>% 
  filter(AA_score > 0.5) %>% 
  filter(QUAL_gt >= 255) %>%
  separate(AD_dv, sep=",", convert = T, into=c("AD_ref_dv", "AD_alt_dv"))%>%
  separate(AD_gk, sep=",", convert = T, into=c("AD_ref_gk", "AD_alt_gk"))%>%
  separate(AD_gt, sep=",", convert = T, into=c("AD_ref_gt", "AD_alt_gt"))%>%
  filter(!grepl("RankSum", INFO_gk)) %>% ## Removes variants where GATK had reads supporting the ALT in other samples %>%
  separate(INFO_gk, sep=";", into=c("AC_gk", "AF_gk", "AN_gk", "DP_info_gk", "ExcessHet_gk",
                                    "FS_gk", "Inbreeding_Coeff_gk", "MLEAC_gk", "MLEAF_gk", "MQ_gk",
                                    "QD_gk", "SOR_gk")) %>%
  mutate_at(vars(ends_with("gk")), ~gsub(".*=", "", .)) %>%
  select(-AC_gk, -AF_gk) %>%
  filter(as.numeric(QD_gk) > 2.0,
         as.numeric(QUAL_gk) >= 30,
         as.numeric(SOR_gk) < 2.5,
         as.numeric(FS_gk) <= 60,
         as.numeric(MQ_gk) >= 40,
         as.numeric(AN_gk) ==140) %>%
  filter(AD_ref_dv == 0 & AD_ref_gk == 0 & AD_ref_gt == 0)  %>%
  filter(AD_alt_dv >= 10 & AD_alt_gk >= 10 & AD_alt_gt >= 10) %>%
  filter(GQ_dv >= 30 & GQ_gt >= 30 & GQ_gk >=30) %>% 
  filter(QUAL_dv > 40) %>% 
  filter(DP_dv <= high_depth_dv, 
         DP_gk <= high_depth_gk,
         DP_gt <= high_depth_gt)

write.csv(test, "mutations_hard_filtering.csv")

