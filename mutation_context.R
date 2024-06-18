library(tidyverse)
library(SomaticSignatures)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/")
ref="~/Downloads/mm39.fa"
genome <- FaFile("~/Downloads/mm39.fa")


vcfRange <- readVcfAsVRanges("data/cc_fake.vcf",ref)


sca_motifs_tri = mutationContext(vcfRange,
                                      genome,k=3, check=T)
sca_motifs_5 = mutationContext(vcfRange,
                                 genome,k=5, check=T)

sca_motifs_11 = mutationContext(vcfRange,
                               genome,k=11, check=T)
sca_motifs_21 = mutationContext(vcfRange,genome,k=21, check=T)

context <- data.frame(ID=as.character(sca_motifs_tri),
                      trinucleotide=sca_motifs_tri$context %>% as.character(),
                      context_5 = sca_motifs_5$context %>% as.character(),
                      context_11 = sca_motifs_11$context %>% as.character(),
                      context_21 = sca_motifs_21$context %>% as.character())
context
write.csv(context, "data/mutation_context.csv", row.names=F)
