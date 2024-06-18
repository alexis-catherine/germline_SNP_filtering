library(tidyverse)


setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

names <- c("CHROM",	
           "POS",
           "ID",
           "REF",
           "ALT",
           "QUAL",
           "FILTER",
           "INFO",	
           "FORMAT",
           "CC001",
           "CC002",
           "CC003",
           "CC004",
           "CC005",
           "CC006",
           "CC007",
           "CC008",
           "CC009",
           "CC010",
           "CC011",
           "CC012",
           "CC013",
           "CC015",
           "CC016",
           "CC017",
           "CC018",
           "CC020",
           "CC021",
           "CC022",
           "CC023",
           "CC024",
           "CC025",
           "CC027",
           "CC028",
           "CC029",
           "CC030",
           "CC031",
           "CC032",
           "CC033",
           "CC034",
           "CC035",
           "CC036",
           "CC037",
           "CC038",
           "CC039",
           "CC040",
           "CC041",
           "CC042",
           "CC043",
           "CC044",
           "CC045",
           "CC046",
           "CC047",
           "CC050",
           "CC051",
           "CC052",
           "CC053",
           "CC055",
           "CC056",
           "CC057",
           "CC058",
           "CC059",
           "CC060",
           "CC061",
           "CC062",
           "CC063",
           "CC065",
           "CC068",
           "CC071",
           "CC072",
           "CC073",
           "CC074",
           "CC075",
           "CC078",
           "CC079",
           "CC080",
           "CC081",
           "CC082",
           "CC083")
graphtyper <- read.table("data/0000_filtered.doubleton.vcf.recode.vcf", header=F)
deepvar <- read.table("data/0001_filtered.doubleton.vcf.recode.vcf", header=F)
gatk <-read.table("data/0002_filtered.doubleton.vcf.recode.vcf", header=F)

names(deepvar) <- names
names(gatk) <- names
names(graphtyper) <- names 

deepvar <- deepvar %>% 
  pivot_longer(starts_with("CC"), names_to="Strain") %>% 
  filter(!grepl("^0\\/0|^\\./\\.|^0/|^\\./|^1/2", value)) %>%
  separate(value, sep=":",
           into=c("GT_dv", "DP_dv", "AD_dv", "GQ_dv", "PL_dv", "RNC_dv"), 
           convert=T)%>%
  filter(RNC_dv == "..") %>% 
  filter(QUAL > 25)
gc()

paste(deepvar$CHROM, deepvar$POS) %>% unique() %>% length()

gatk <- gatk %>% 
  pivot_longer(starts_with("CC"), names_to="Strain") %>% 
  filter(!grepl("^0\\/0|^\\./\\.|^0/|^0\\||^\\./|^0\\|0|^\\.|^1/2|^1/3", value)) %>%
  separate(value, sep=":", into=c("GT_gk", "AD_gk", "DP_gk", "GQ_gk", "PGT_gk", "PID_gk",
                                  "PL_gk", "PS_gk"), convert=T)
gc()

graphtyper <- graphtyper %>% 
  pivot_longer(starts_with("CC"), names_to="Strain") %>% 
  filter(!grepl("^0\\/0|^\\./\\.|^0/|^\\./|^2/3|^1/2|^1/3", value)) %>%
  separate(value, sep=":", into=c("GT_gt", "AD_gt", "MD_gt", "DP_gt", "GQ_gt", "PL_gt"), convert=T)
gc()


deepvar_uniq <- deepvar %>% group_by(ID) %>% tally() %>% arrange(n) %>% 
  filter(n==1)
deepvar <- deepvar %>% filter(ID %in% deepvar_uniq$ID)

gatk_uniq <- gatk %>% group_by(ID=paste(CHROM, POS)) %>% tally() %>% arrange(n) %>% 
  filter(n==1)
gatk <- gatk %>% filter(paste(CHROM, POS) %in% gatk_uniq$ID)

graphtyper_uniq <- graphtyper %>%
  group_by(ID) %>% tally() %>% 
  arrange(n) %>% 
  filter(n==1)
graphtyper <- graphtyper %>% filter(ID %in% graphtyper_uniq$ID)

together <- bind_rows(select(deepvar, CHROM, POS, Strain, ALT),
                    select(graphtyper, CHROM, POS, Strain, ALT),
                    select(gatk, CHROM, POS, Strain, ALT)) %>% 
  group_by(CHROM,POS,Strain, ALT) %>% 
  summarise(n=n()) %>% 
  arrange(-n) %>% 
  filter(n==3) %>% 
  filter(nchar(ALT) == 1) %>% 
  mutate(ID = paste(CHROM, POS, Strain, ALT))

graphtyper <- graphtyper %>% 
  filter(paste(CHROM, POS, Strain, ALT) %in% together$ID)
gatk <- gatk %>% 
  filter(paste(CHROM, POS, Strain, ALT) %in% together$ID)
deepvar <- deepvar %>% 
  filter(paste(CHROM, POS, Strain, ALT) %in% together$ID)

rm(deepvar_uniq, gatk_uniq, graphtyper_uniq, together, names)
gc()

deepvar <- deepvar %>% 
  select(-ID, 
         -FILTER,
         -FORMAT,
         -GT_dv) %>% 
  select(Strain,
         CHROM, 
         POS,
         REF,
         ALT,
         INFO_dv=INFO,
         QUAL_dv=QUAL,
         everything())

gatk <- gatk %>% 
  select(-ID, 
         -FILTER,
         -FORMAT,
         -GT_gk) %>% 
  select(Strain,
         CHROM, 
         POS,
         REF,
         ALT,
         INFO_gk=INFO,
         QUAL_gk=QUAL,
         everything())

graphtyper <- graphtyper %>% 
  select(-ID, 
         -FILTER,
         -FORMAT,
         -GT_gt) %>% 
  select(Strain,
         CHROM, 
         POS,
         REF,
         ALT,
         INFO_gt=INFO,
         QUAL_gt=QUAL,
         everything())

full <- merge(deepvar,
              gatk, 
              by=c("Strain", "CHROM", "POS", "REF", "ALT")) %>% 
  merge(graphtyper, 
        by=c("Strain", "CHROM", "POS", "REF", "ALT"))

full <- full %>% select(-RNC_dv)

rm(deepvar, gatk,graphtyper)
gc()


strains <- read.csv("data/cc_genome_stats.csv")
strains$Strain <- gsub("\\/.*", "", strains$Strain)

full <- merge(strains, full, by="Strain") 

write.csv(full, "gatk_graph_deep.csv")
