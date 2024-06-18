library(tidyverse)


setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("gatk_graph_deep.csv", header=T) %>% 
  select(-X)

dat_full <- dat %>% 
  mutate(INFO=INFO_dv,
         GT="1/1",
         FILTER="PASS") %>% 
  mutate(FORMAT="GT") %>% 
  select(CHROM, POS, ID, REF,ALT,QUAL=QUAL_dv,FILTER,INFO,FORMAT,CCFAKE=GT) %>% 
  write.table(file=paste0("data/cc_fake.vcf"),
              col.names = F, row.names = F,
              quote=F, sep="\t")

for(i in unique(dat$Strain)) {
dat %>% 
    filter(Strain==i) %>% 
    mutate(INFO=INFO_dv,
           GT="1/1",
           FILTER="PASS") %>% 
    mutate(FORMAT="GT") %>% 
    select(CHROM, POS, ID, REF,ALT,QUAL=QUAL_dv,FILTER,INFO,FORMAT,Strain, GT) %>%
    pivot_wider(names_from = Strain, values_from=GT) %>% 
    write.table(file=paste0("data/sample_vcf/",i, ".vcf"),
                col.names = F, row.names = F,
                quote=F, sep="\t")
}
