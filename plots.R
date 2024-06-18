library(tidyverse)
library(ggpubr)
library(ggridges)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")

dat <- read.csv("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/mutations_hard_filtering.csv")

dat %>% group_by(Strain, Generation) %>% 
  summarise(n=n()) %>% 
  mutate(mut=n/Generation)

dat %>% 
  group_by(Strain,Founders,Generation, Coverage) %>% 
  summarise(n=n()) %>% 
  mutate(mut=n/Generation) %>%
  ggplot(aes(x=reorder(Strain,mut), y=mut)) +
  geom_bar(stat="identity") +
  theme_classic() +
#  geom_hline(yintercept = 2) +
  xlab("Strain") +
#  geom_hline(yintercept = 5) +
  ylab("Per Generation Mutation Rate") +
  rotate_x_text()
ggsave("strain_mut_bar.png", height=5, width=8)

dat %>% 
  group_by(Strain,Founders,Generation, Coverage) %>% 
  summarise(n=n()) %>% 
  mutate(mut=n/Generation) %>%
  ggplot(aes(x=Founders, y=mut)) + 
  geom_point() + 
  theme_classic() + 
  geom_smooth(method="lm") +
  stat_cor()

ggplot(dat %>%
       mutate(type = paste0(REF, ">", ALT)) %>% 
         mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                                 type %in% c("T>C", "A>G") ~ "T>C|A>G",
                                 type %in% c("T>A", "A>T") ~ "T>A|A>T",
                                 type %in% c("G>T", "C>A") ~ "C>A|G>T",
                                 type %in% c("G>C", "C>G") ~ "C>G|G>C",
                                 type %in% c("G>A", "C>T") ~ "C>T|G>A")) , aes(x=trinucleotide)) +
  geom_bar() +
  facet_wrap(~type) + theme_classic() +rotate_x_text()



dat %>%
  group_by(Strain,Founders,Generation, Coverage, CHROM) %>% 
  summarise(n=n()) %>% 
  mutate(mut=n/Generation) %>%
  ggplot(aes(y=reorder(Strain,mut), x=mut)) +
  geom_point(size=1) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_vline(xintercept = 5/19)+
  geom_vline(xintercept = 2/19)+
  ylab("Strain") +
  xlab("Per Generation Mutation Rate") + xlim(0,10)
ggsave("strain_mut_box.png", height=8, width=5)

dat %>%
  group_by(Strain,Founders,Generation, Coverage, CHROM) %>% 
  
  summarise(n=n()) %>% 
  mutate(mut=n/Generation) %>%
  ggplot(aes(x=reorder(CHROM,mut), y=mut)) +
  geom_point() +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  geom_hline(yintercept = 5/19)+
  geom_hline(yintercept = 2/19)+
  ylab("Strain") +
  xlab("Per Generation Mutation Rate")+
  scale_y_log10()
#ggsave("strain_mut_box.png", height=8, width=5)

dat %>% 
  mutate(type = paste0(REF, ">", ALT)) %>% 
  mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                          type %in% c("T>C", "A>G") ~ "T>C|A>G",
                          type %in% c("T>A", "A>T") ~ "T>A|A>T",
                          type %in% c("G>T", "C>A") ~ "C>A|G>T",
                          type %in% c("G>C", "C>G") ~ "C>G|G>C",
                          type %in% c("G>A", "C>T") ~ "C>T|G>A")) %>% 
  group_by(Strain,Founders,Generation, CHROM) %>% 
  summarise(n=n()) %>% 
  mutate(mut=n/Generation) %>%
  summarise(mut=median(mut)) %>% 
  ggplot(aes(x=reorder(Strain,mut), y=mut*22)) +
  geom_bar(stat="identity") +
  theme_classic() +
  geom_hline(yintercept = 5)+
  geom_hline(yintercept = 2)+
  xlab("Strain") +
  ylab("Median Per Generation Mutation Rate") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("median_for_poster.png", height=5, width=10)

dat %>% 
  mutate(type = paste0(REF, ">", ALT)) %>% 
  mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                          type %in% c("T>C", "A>G") ~ "T>C|A>G",
                          type %in% c("T>A", "A>T") ~ "T>A|A>T",
                          type %in% c("G>T", "C>A") ~ "C>A|G>T",
                          type %in% c("G>C", "C>G") ~ "C>G|G>C",
                          type %in% c("G>A", "C>T") ~ "C>T|G>A")) %>% 
  group_by(Strain,Founders,Generation, type, CHROM) %>% 
  summarise(n=n()) %>% 
  summarise(n=mean(n)) %>% 
  mutate(n=n/Generation) %>% 
  pivot_wider(names_from=type, values_from=n) -> mut_medians
write.csv(mut_medians,"type_mutation_rate_Apr_not_final.csv", row.names = F)

dat %>% 
  mutate(type = paste0(REF, ">", ALT)) %>% 
  mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                          type %in% c("T>C", "A>G") ~ "T>C|A>G",
                          type %in% c("T>A", "A>T") ~ "T>A|A>T",
                          type %in% c("G>T", "C>A") ~ "C>A|G>T",
                          type %in% c("G>C", "C>G") ~ "C>G|G>C",
                          type %in% c("G>A", "C>T") ~ "C>T|G>A")) %>% 
  group_by(Strain,Founders,Generation, Coverage, CHROM, type) %>% 
  summarise(n=n()) %>% 
  mutate(mut=n/Generation) %>%
  summarise(mut=median(mut)) %>% 
  write.csv("mutation_rate_Apr_not_final.csv", row.names = F)



dat %>%
  mutate(type = paste0(REF, ">", ALT)) %>% 
  mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                          type %in% c("T>C", "A>G") ~ "T>C|A>G",
                          type %in% c("T>A", "A>T") ~ "T>A|A>T",
                          type %in% c("G>T", "C>A") ~ "C>A|G>T",
                          type %in% c("G>C", "C>G") ~ "C>G|G>C",
                          type %in% c("G>A", "C>T") ~ "C>T|G>A",
  )) %>% 
  group_by(Strain, type) %>% 
  summarise(n=n()) %>% group_by(Strain) %>% 
  mutate(fract=n/sum(n)) %>%  
  ggplot(aes(x=type,y=fract)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha=0.5) +
  theme_classic() + rotate_x_text()
ggsave("type_box.png", height=4,width=3)

dat %>%
  mutate(type = paste0(REF, ">", ALT)) %>% 
  mutate(type = case_when(type %in% c("T>G", "A>C") ~ "T>G|A>C",
                          type %in% c("T>C", "A>G") ~ "T>C|A>G",
                          type %in% c("T>A", "A>T") ~ "T>A|A>T",
                          type %in% c("G>T", "C>A") ~ "C>A|G>T",
                          type %in% c("G>C", "C>G") ~ "C>G|G>C",
                          type %in% c("G>A", "C>T") ~ "C>T|G>A",
  )) %>% 
  group_by(Strain, type) %>% 
  summarise(n=n()) %>%
  group_by(Strain) %>% 
  mutate(fract=n/sum(n)) %>%  
  ggplot(aes(x=Strain,y=fract,fill=type))+
  geom_bar(stat="identity") +
  theme_classic() +
  ylab("Fraction of DNMs") +
  xlab("Strain") +
  scale_fill_manual(values=unname(qtl2::CCcolors))  +
  scale_x_discrete(guide = guide_axis(angle = 90)) 
ggsave("strain_fract.png", height=5, width=8)


dat %>% filter(Strain == "CC038", 
               CHROM == "chr2") %>% 
  ggplot(aes(x=POS)) +
  geom_histogram() 

dat %>% filter(Strain == "CC038") %>% 
  ggplot(aes(x=DP_gt, y=CHROM, fill=ifelse(CHROM=="chr2", "chr2", NA))) +
  geom_density_ridges_gradient() +
  theme_classic() +
  theme(legend.position = "none")

dat %>% filter(Strain == "CC038") %>% 
  ggplot(aes(x=QUAL_gt, y=CHROM, fill=ifelse(CHROM=="chr2", "chr2", NA))) +
  geom_density_ridges_gradient() +
  theme_classic() +
  theme(legend.position = "none")

dat %>% filter(Strain == "CC038") %>% 
  ggplot(aes(x=GQ_gt, y=CHROM, fill=ifelse(CHROM=="chr2", "chr2", NA))) +
  geom_density_ridges_gradient() +
  theme_classic() +
  theme(legend.position = "none")

dat %>% filter(Strain == "CC038") %>% 
  filter(CHROM == "chr2") %>% 
  ggplot(aes(x=POS, y=DP_gt)) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point()

dat %>%
  filter(CHROM == "chr2") %>% 
  ggplot(aes(x=POS, y=DP_gt, color=Strain)) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point()
dat %>%
  filter(CHROM == "chr2") %>% 
  ggplot(aes(x=POS, y=QUAL_gt, color=Strain)) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_point()

dat %>%
  filter(CHROM == "chr2") %>% 
  filter(Strain == "CC038") %>%
  ggplot(aes(x=POS, y=GQ_gt, color=Strain)) +
  theme_classic() +
  theme(legend.position = "none") +
  geom_jitter(size=0.05) +
  geom_hline(data=dat %>%
               filter(Strain == "CC038") %>%
               filter(CHROM != "chr2") %>% 
               group_by(Strain) %>% 
               summarise(me=mean(GQ_gt)),
             aes(yintercept=me)) +
  geom_hline(data=dat %>%
               filter(Strain == "CC038") %>%
               filter(CHROM != "chr2") %>% 
               group_by(Strain) %>% 
               summarise(me=mean(GQ_gt),
                         sd=sd(GQ_gt)),
             aes(yintercept=me+sd)) +
  geom_hline(data=dat %>%
               filter(Strain == "CC038") %>%
               filter(CHROM != "chr2") %>% 
               group_by(Strain) %>% 
               summarise(me=mean(GQ_gt),
                         sd=sd(GQ_gt)),
             aes(yintercept=me-sd)) +
  facet_wrap(~Strain, scales = "free_y")

dat %>%
  filter(CHROM == "chr2") %>% 
  group_by(Strain) %>% 
  summarise(me=mean(GQ_gt)) %>% 
  arrange(me)

ggplot(dat, aes(x=DP_dv)) + geom_histogram()

