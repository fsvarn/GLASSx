library(tidyverse)
library(ssgsea.GBM.classification)
library(GSVA)
library(RColorBrewer)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/downstream/myeloid/glass_transcriptional_subtype_comparisons.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/downstream/myeloid/tcga_transcriptional_subtype_comparisons.txt"
myinf3 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/downstream/myeloid/klemm_transcriptional_subtype_comparisons.txt"

dat1 <- read.delim(myinf1)
dat2 <- read.delim(myinf2)
dat3 <- read.delim(myinf3)

subdat1 <- dat1[which(dat1[,"subtype"] == "Mesenchymal"),]
subdat2 <- dat2[which(dat2[,"subtype"] == "Mesenchymal"),]
subdat3 <- dat3[which(dat3[,"subtype"] == "Mesenchymal"),]

subdat1[which(subdat1[,"median_effect"] < 0), "q.value"] <-  subdat1[which(subdat1[,"median_effect"] < 0), "q.value"] * -1
subdat2[which(subdat2[,"median_effect"] < 0), "q.value"] <-  subdat2[which(subdat2[,"median_effect"] < 0), "q.value"] * -1
subdat3[which(subdat3[,"median_effect"] < 0), "q.value"] <-  subdat3[which(subdat3[,"median_effect"] < 0), "q.value"] * -1

examine <- data.frame(subdat1[,"q.value"], subdat2[,"q.value"], subdat3[,"q.value"])
colnames(examine) <- c("glass","tcga","klemm")

glass_pos <- rownames(examine[which(examine[,"glass"] > 0 & abs(examine[,"glass"]) < 0.1),])
glass_neg <- rownames(examine[which(examine[,"glass"] < 0 & abs(examine[,"glass"]) < 0.1),])
glass_sig <- c(glass_pos, glass_neg)

tcga_pos <- rownames(examine[which(examine[,"tcga"] > 0 & abs(examine[,"tcga"]) < 0.1),])
tcga_neg <- rownames(examine[which(examine[,"tcga"] < 0 & abs(examine[,"tcga"]) < 0.1),])
tcga_sig <- c(tcga_pos, tcga_neg)

klemm_pos <- rownames(examine[which(examine[,"klemm"] > 0 & abs(examine[,"klemm"]) < 0.1),])
klemm_neg <- rownames(examine[which(examine[,"klemm"] < 0 & abs(examine[,"klemm"]) < 0.1),])
klemm_sig <- c(klemm_pos, klemm_neg)

csx <- c(intersect(glass_pos, tcga_pos),  intersect(glass_neg, tcga_neg))
valpo <- intersect(intersect(glass_pos, klemm_pos), intersect(tcga_pos, klemm_pos))
valne <- intersect(intersect(glass_neg, klemm_neg), intersect(tcga_neg, klemm_neg))


# Mesenchymal
# > valpo
# [1] "15" "20" "28"
# > valne
# [1] "26" "35"