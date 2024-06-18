library(tidyverse)
library(ggpubr)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("removing_problematic_mB.csv")

dat %>%
  mutate(type = paste0(REF, ">", ALT)) %>% 
  mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                          type %in% c("T>C", "A>G") ~ "T>C|A>G",
                          type %in% c("T>A", "A>T") ~ "T>A|A>T",
                          type %in% c("G>T", "C>A") ~ "C>A|G>T",
                          type %in% c("G>C", "C>G") ~ "C>G|G>C",
                          type %in% c("G>A", "C>T") ~ "C>T|G>A")) %>% 
  ggplot(aes(x=trinucleotide)) + 
  facet_wrap(~type) +
  geom_bar() +
  theme_classic() +
  rotate_x_text()

dat <- dat %>% group_by(Strain, Generation) %>%
  mutate(type = paste0(REF, ">", ALT)) %>% 
  mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                          type %in% c("T>C", "A>G") ~ "T>C|A>G",
                          type %in% c("T>A", "A>T") ~ "T>A|A>T",
                          type %in% c("G>T", "C>A") ~ "C>A|G>T",
                          type %in% c("G>C", "C>G") ~ "C>G|G>C",
                          type %in% c("G>A", "C>T") ~ "C>T|G>A")) %>% 
  group_by(Strain, Generation, type) %>%
  summarise(n=n()) %>% 
  mutate(total_n = sum(n),
         total=sum(n)/Generation,
         n=n/Generation) %>% 
  pivot_wider(names_from=type, values_from=n, values_fill = 0)

######################################
#calculate the expected number of mutations in sequenced cc lines

mutation_rate = 0.5e-09 #assume equal mutation rate in males and females
n_gens_inbreeding = 37
length = 2900000000 #size of the region to simulate
nReps = 100
genome_size = 2900000000

##########################################

mutation_simulations <- function(mutation_rate,
                                 n_gens_inbreeding,
                                 nReps,
                                 genome_size,
                                 lengths){
  #determine the fraction of de novo SNPs that should be het versus homozygous
  
  fraction_het = numeric(nReps)
  fraction_hom = numeric(nReps)
  
  number_het = numeric(nReps)
  number_hom = numeric(nReps)
  
  for (rep in 1:nReps) {
    
    #start with 4 G2:F1 haplotypes, each carrying de novo mutations inherited from either mom or dad
    r = rpois(2, 2*length*mutation_rate) # 2 parents/gen * 0.5 chance of mutation transfer at any generation * 2 generations
    f.prev.hap1 = runif(r[1], min = 1, max = length)
    f.prev.hap2 = runif(r[2], min = 1, max = length)
    
    r = rpois(2, 2*length*mutation_rate)
    m.prev.hap1 = runif(r[1], min = 1, max = length)
    m.prev.hap2 = runif(r[2], min = 1, max = length)
    
    for (g in 1:n_gens_inbreeding) {
      #randomly sample 1 male and 1 female chromosome to simulate female
      r = runif(1)
      if (r < 0.5) {f.curr.hap1 = f.prev.hap1}
      if (r >= 0.5) {f.curr.hap1 = f.prev.hap2}
      
      r = runif(1)
      if (r < 0.5) {f.curr.hap2 = m.prev.hap1}
      if (r >= 0.5) {f.curr.hap2 = m.prev.hap2}
      
      #mutate the haplotypes
      r = rpois(2, length*mutation_rate)
      pos = runif(r[1], min = 1, max = length)
      f.curr.hap1 = c(f.curr.hap1, pos)
      pos = runif(r[2], min = 1, max = length)
      f.curr.hap2 = c(f.curr.hap2, pos)
      
      #randomly sample 1 male and 1 female chromosome to simulate male
      r = runif(1)
      if (r < 0.5) {m.curr.hap1 = f.prev.hap1}
      if (r >= 0.5) {m.curr.hap1 = f.prev.hap2}
      
      r = runif(1)
      if (r < 0.5) {m.curr.hap2 = m.prev.hap1}
      if (r >= 0.5) {m.curr.hap2 = m.prev.hap2}
      
      #mutate the haplotypes
      r = rpois(2, length*mutation_rate)
      pos = runif(r[1], min = 1, max = length)
      m.curr.hap1 = c(m.curr.hap1, pos)
      pos = runif(r[2], min = 1, max = length)
      m.curr.hap2 = c(m.curr.hap2, pos)
      
      #re-assign current haplotypes as previous haplotypes
      m.prev.hap1 = m.curr.hap1
      m.prev.hap2 = m.curr.hap2
      f.prev.hap1 = f.curr.hap1
      f.prev.hap2 = f.curr.hap2
    }
    
    #calculate the fraction of het sites in the male
    n_sites = length(m.prev.hap1)+length(m.prev.hap2)
    n_het = length(intersect(m.prev.hap1, m.prev.hap2))
    number_het[rep] = n_het
    number_hom[rep] = n_sites-n_het
    fraction_het[rep] = n_het/n_sites
    fraction_hom[rep] = 1-(n_het/n_sites)
    
  }
  return(number_hom)
}

for(l in 1:nrow(dat)){
  out <- mutation_simulations(mutation_rate,
                              n_gens_inbreeding=dat$Generation[l],
                              nReps=1000,
                              genome_size,
                              lengths) %>% data.frame()
  p <- ggplot(out, aes(x=.)) + 
    geom_histogram() + 
    geom_vline(xintercept=dat$n[l])
  ggsave(plot= p, 
         filename=paste0("simulations/",
                         dat$Strain[l],
                         ".png"),
         height=3, width=4)
  dat$zscore[l] <- (dat$total_n[l]-mean(out$.))/sd(out$.)
  }

ggplot(dat, aes(x=zscore)) +
  geom_histogram() + 
  theme_classic()

ggplot(dat, aes(x=Generation, y=zscore)) +
  geom_point() +
  geom_smooth(method="lm") + 
  theme_classic() + 
  stat_cor()

ggplot(dat, aes(x=mut, y=zscore)) +
  geom_point() +
  geom_smooth(method="lm") + 
  theme_classic() + 
  stat_cor()

dat
write.csv(dat, "new_mutation_Rates_05_15_2024.csv", row.names = F)
