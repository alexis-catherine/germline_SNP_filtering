library(tidyverse)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SomaticSignatures")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(SomaticSignatures)
library(Rsamtools)
library(VariantAnnotation)
library(deconstructSigs)
library(ggfortify)
library(ggrepel)

setwd("/Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/Apr_2024/data/sample_vcf/")

files <- list.files(pattern="*.vcf")
ref="~/Downloads/mm39.fa"
genome <- FaFile("~/Downloads/mm39.fa")

for( i in 1:70){
  vcfRange <- readVcfAsVRanges(files[i],ref)

  
  sca_motifs = tryCatch(mutationContext(vcfRange,
                                        genome),
                        error=function(e) e)
  
  sca_mm = motifMatrix(sca_motifs,
                       group = "sampleNames",
                       normalize = TRUE) #+ 0.0000000001
  plotMutationSpectrum(sca_motifs, "sampleNames") + 
    theme_classic() + 
    scale_fill_viridis_d()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0("spectrum/", gsub(".vcf", "", files[i]),"ivf_spectrum.png"),
         height=3*vcfRange@sampleNames %>% unique() %>% length(),
         width=10,limitsize = F)
  
  write.csv(sca_mm, paste0("spectrum/", gsub(".vcf", "", files[i]), "IVF_spectrum.csv"))
}

setwd("spectrum/")
files<- list.files(pattern="*spectrum.csv")
theData_list<-lapply(files, read.csv)
dat <-bind_cols(theData_list)
titles <- dat$X...1
dat <- dat %>% 
  dplyr::select(-starts_with("X."))
rownames(dat) <- titles

sigs_mm <- as.matrix(dat)
n_sigs=5

sigs_pca <- identifySignatures(sigs_mm,
                               n_sigs,
                               pcaDecomposition)
#sigs_nmf = identifySignatures(sigs_mm, n_sigs, nmfDecomposition)

sig_input <- t(sigs_mm)%>% as.data.frame()
names(sig_input) <- names(signatures.nature2013)

weights <- data.frame()
weights_nature2013 <- data.frame()

for(j in 1:70){
  plot <- whichSignatures(tumor.ref=sig_input,
                          sample.id=rownames(sig_input)[j],
                          signatures.ref = signatures.cosmic)
  weights <- bind_rows(weights, plot$weights)
  
  plot <- whichSignatures(tumor.ref=sig_input,
                            sample.id=rownames(sig_input)[j],
                            signatures.ref = signatures.nature2013)
  weights_nature2013 <- bind_rows(weights_nature2013, plot$weights)
}
gc()

write.csv(weights_nature2013, "../../signatures_nature2013_across_all_samples.csv")
write.csv(weights, "../../signatures_across_all_samples.csv")

write.csv(sigs_mm, "../../trinucleotide_context.csv")

setwd("..")
files <- list.files(pattern="*.vcf")
for( i in 1:70){
  vcfRange <- readVcfAsVRanges(files[i],ref)
  
  
  sca_motifs = tryCatch(mutationContext(vcfRange,
                                        genome, k=5, check=T),
                        error=function(e) e)
  
  sca_mm = motifMatrix(sca_motifs,
                       group = "sampleNames",
                       normalize = TRUE) #+ 0.0000000001
  plotMutationSpectrum(sca_motifs, "sampleNames") + 
    theme_classic() + 
    scale_fill_viridis_d()+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(paste0("spectrum/", gsub(".vcf", "", files[i]),"ivf_spectrum_5.png"),
         height=3*vcfRange@sampleNames %>% unique() %>% length(),
         width=10,limitsize = F)
  
  write.csv(sca_mm, paste0("spectrum/", gsub(".vcf", "", files[i]), "IVF_spectrum_5.csv"))
}

