#qtlMapping_sexRatio.R

library(qtl2)
library(data.table)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
#setwd("C:/Users/alexi/Box/Dumont_Lab/Research/CC_reproduction/")
setwd("//Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_genome_projects/")
colors = brewer.pal(n = 8, name="Dark2")

load("//Users/garrea/Library/CloudStorage/Box-Box/Dumont_Lab (Alexis C. Garretson)/Research/CC_reproduction/reproductive_success_and_decline/Reproductive_Decline/orig/CC_mapping.RData")

setwd("Apr_2024/")
strains <- srd$geno$`1` %>% row.names()
pheno <- read.csv("new_mutation_Rates_05_15_2024.csv") %>% select(-Generation)
merge(pheno, data.frame(strains = strains,
                      short_strains = gsub("\\/.*", "",
                                           strains)),
      by.x="Strain", by.y="short_strains") %>%
  select(Strain = strains,
         everything(),
         -Strain) %>%
  write.csv("combined_pheno.csv", row.names = F)

srd$pheno = read_pheno("combined_pheno.csv") 

srd$pheno

mainDir=getwd()

scan_out <- scan1(pr,
                  srd$pheno, 
                  kinship_loco,
                  Xcovar=Xcovar, 
                  cores = 2)  
out_perm <- scan1perm(pr, srd$pheno,
                      kinship_loco,
                      Xcovar=Xcovar, n_perm=100,
                      perm_Xsp=TRUE,
                      chr_lengths=chr_lengths(srd$gmap),
                      cores=2)
saveRDS(out_perm, "out_perm.rds")

thr = summary(out_perm)
out_ymx <- maxlod(scan_out) # overall maximum LOD score

phenotype = "mut"
png(paste0(phenotype, ".png"))
for(i in 1:ncol(srd$pheno)){
par(mar=c(5.1, 4.1, 1.1, 1.1))
  plot(x = scan_out,
       map = srd$gmap, 
       lodcolumn=i,
       main = paste(colnames(srd$pheno)[i]),
  ylim=c(0,thr$A[i]))#, round(thr$A[[i]],3)))
abline(h = thr$A[i], col = "red", lwd = 2)
}

plot(x = scan_out,
     map = srd$gmap, 
     lodcolumn=7,
     main = paste(colnames(srd$pheno)[7]),
     ylim=c(0,(thr$A[7])+5))#, round(thr$A[[i]],3)))
abline(h = thr$A[i], col = "red", lwd = 2)

abline(h = summary(out_perm,0.05)$A[1], col = "maroon", lwd = 2)
abline(h = summary(out_perm,0.1)$A[1], col = "maroon", lwd = 2)
abline(h = summary(out_perm,0.15)$A, col = "darkred", lwd = 2)
dev.off()

thr = summary(out_perm)
thr_0.1 = summary(out_perm,0.1)
thr_0.5 = summary(out_perm,0.5)
out_ymx <- maxlod(scan_out) # overall maximum LOD score


peak_Mbp <- find_peaks(scan1_output = scan_out, map = srd$pmap, peakdrop=1.8, prob=0.95) %>% 
  filter(lod > 6) %>% 
  filter(lod == max(lod))
peak_Mbp$ci_hi <- ifelse(peak_Mbp$ci_lo==peak_Mbp$ci_hi, peak_Mbp$ci_hi+1, peak_Mbp$ci_hi)

coef_cX <- scan1coef(apr[,peak_Mbp$chr], srd$pheno[,1], kinship_loco[[peak_Mbp$chr]])
# plot_coefCC(coef_cX,
#             list(`5`=srd$gmap$`5`[srd$gmap$`5`<42.944&srd$gmap$`5`>38.427]), 
#             scan1_output=scan_out, bgcolor="gray95",
#             legend="topright")
peak_gmap <- find_peaks(scan1_output = scan_out, map = srd$pmap, peakdrop=1.8, prob=0.95) %>% 
  filter(lod > 1) %>% 
  filter(lod == max(lod))
peak_gmap$ci_hi <- ifelse(peak_gmap$ci_lo==peak_gmap$ci_hi, peak_gmap$ci_hi+3, peak_gmap$ci_hi)

region <- list(srd$gmap[[as.character(peak_gmap$chr)]][srd$gmap[[as.character(peak_gmap$chr)]]
                                                       < peak_gmap$ci_hi & srd$gmap[[as.character(peak_gmap$chr)]]>peak_gmap$ci_lo])

names(region) <- as.character(peak_gmap$chr)
png(paste0(2, "_short_coeff.png"))
par(mar=c(8.1, 8.1, 1.1, 1.1))
plot_coefCC(coef_cX,
            region, 
            scan1_output=scan_out, 
            bgcolor="gray95",
            legend="bottomleft")
dev.off()

png(paste0(i, "_coeff.png"))
par(mar=c(8.1, 8.1, 1.1, 1.1))
plot_coefCC(coef_cX,
            srd$gmap, 
            scan1_output=scan_out, 
            bgcolor="gray95",
            legend="topright")
dev.off()

query_variants <- create_variant_query_func("../../cc_variants.sqlite")

find_peaks(scan1_output = scan_out, map = srd$gmap, peakdrop=1.8, prob=0.95) %>% filter(lod > 6)

variants <- query_variants(peak_Mbp$chr, peak_Mbp$ci_lo, peak_Mbp$ci_hi)
out_snps <- scan1snps(pr, 
                      srd$pmap,
                      srd$pheno, kinship_loco[[peak_Mbp$chr]],
                      query_func=query_variants,
                      chr=peak_Mbp$chr, 
                      start=peak_Mbp$ci_lo,
                      end=peak_Mbp$ci_hi, 
                      keep_all_snps=TRUE)

query_genes <- create_gene_query_func("../../../Research/mouse_genes_mgi.sqlite")

genes <- query_genes(peak_Mbp$chr[1], 
                     peak_Mbp$ci_lo[1], 
                     peak_Mbp$ci_hi[1]) %>% 
  bind_rows(query_genes(peak_Mbp$chr[2], 
                        peak_Mbp$ci_lo[2], 
                        peak_Mbp$ci_hi[2])) %>% 
  bind_rows(query_genes(peak_Mbp$chr[3], 
                        peak_Mbp$ci_lo[3], 
                        peak_Mbp$ci_hi[3])) %>% 
  bind_rows(query_genes(peak_Mbp$chr[4], 
                        peak_Mbp$ci_lo[4], 
                        peak_Mbp$ci_hi[4])) %>%
  bind_rows(query_genes(peak_Mbp$chr[5], 
                        peak_Mbp$ci_lo[5], 
                        peak_Mbp$ci_hi[5])) 
png(paste0(i, "_snps.png"))
par(mar=c(8.1, 8.1, 1.1, 1.1))
plot(out_snps$lod, 
     out_snps$snpinfo,
     drop_hilit=0.05, 
     genes=genes, color="black",
     show_all_snps=T)
dev.off()

write.csv(genes, "genes.csv", row.names = F)
write.csv(out_snps$snpinfo, "snps.csv", row.names = F)
write.csv(out_snps$lod, "snps_lod.csv")
setwd("..")

gc()

print(i)
print(phenotype)


apr <- genoprob_to_alleleprob(pr)

coef_cX <- scan1coef(apr[,"2"], srd$pheno[,"Litter_Size_Dam_Age_individual"], kinship_loco[["2"]])
plot_coefCC(coef_cX,
            srd$gmap["2"], scan1_output=scan_out, bgcolor="gray95",
            legend="bottomleft")

query_variants <- create_variant_query_func("D:/Classes_Fall_2020/churchill/Data/sql_mouse_data/cc_variants.sqlite")

find_peaks(scan1_output = scan_out, map = srd$gmap, peakdrop=1.8, prob=0.95) %>% filter(lod > 6)

peak_Mbp <- max(scan_out, srd$pmap)
variants <- query_variants(peak_Mbp$chr, peak_Mbp$pos-1, peak_Mbp$pos+1)
out_snps <- scan1snps(pr, 
                      srd$pmap,
                      srd$pheno, kinship_loco[["2"]],
                      query_func=query_variants,
                      chr=peak_Mbp$chr, 
                      start=peak_Mbp$pos-1,
                      end=peak_Mbp$pos+5, 
                      keep_all_snps=TRUE)


query_genes <- create_gene_query_func("D:/Classes_Fall_2020/churchill/Data/sql_mouse_data/mouse_genes_mgi.sqlite")
genes <- query_genes(peak_Mbp$chr,
                     peak_Mbp$pos-1, 
                     peak_Mbp$pos+5)

plot(out_snps$lod, 
     out_snps$snpinfo,
     drop_hilit=0.5, 
     genes=genes, color="black")

out_gwas <- scan1snps(pr,
                      srd$pmap,
                      srd$pheno,
                      kinship = kinship_loco,
                      query_func = query_variants,
                      cores=1)

plot(out_gwas$lod,
     out_gwas$snpinfo, 
     altcol="green4", 
     gap=0)

write.csv(genes, "genes_decline_litter_size_dam_Age.csv", row.names = F)
write.csv(out_snps$snpinfo, "snps_decline_litter_size_dam_age.csv", row.names = F)
write.csv(out_snps$lod, "snps_lod_decline_litter_size_dam_age.csv")


ggplot(dat, aes(x=Litter_Size_Dam_Age_individual)) +
  geom_histogram(bins =10) +
  theme_classic() +
  ylab("Number of Strains") +
  xlab("Litter Size ~ Dam Age (Slope)")
ggsave("litter_size_dam_age_hist.png", height=5, width=5)


correlations <- read.csv("output/mutation_rates_07062021.csv") %>% merge(dat, by.x="Strain", by.y="strain") %>%
  select(-Strain) %>% 
  select(-starts_with("age"))
  
library(ggcorrplot)
library(ggpubr)
ggcorrplot::ggcorrplot(cor(correlations)[1:7,-c(1:7)], 
                       )

ggplot(correlations %>% 
       filter(mutation_rate < 15), aes(x=mutation_rate, y=Litter_Size_Dam_Age_individual)) +
  geom_point() +
  stat_cor() +
  geom_smooth(method="lm", se=F) +
  ylab("Litter Size ~ Dam Age") + 
  xlab("Mutation Rate")
