###################################################
# Final necrosis vs CIBERSORT plots for TCGA
# Author: Frederick Varn
# Date: 2022.01.19
# Figures S2D- TCGA
##################################################

library(odbc)
library(DBI)
library(tidyverse)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

myinf1 <- "data/dataset/tcga/TCGA_ivygap_fractions.txt"
myinf2 <- "data/dataset/tcga/tcga_extent_of_nec_sample.txt"
myinf3 <- "data/dataset/tcga/tcga_extent_of_angio.txt"
myinf4 <- "data/dataset/tcga/tcga_extent_of_nec_case.txt"

dat <- read.delim(myinf1)
nec <- read.delim(myinf2)
ang <- read.delim(myinf3)
nec2 <- read.delim(myinf4)

dat[,"Mixture"] <- gsub("\\.","-",dat$Mixture)
dat[,"Sample"] <- sapply(strsplit(dat$Mixture,"-"),function(x)paste(x[1:4],collapse="-"))
dat[,"Case"] <- sapply(strsplit(dat$Mixture,"-"),function(x)paste(x[1:3],collapse="-"))

nec_comp <- dat %>%
inner_join(nec, by = "Sample")

cor(nec_comp$Markup.Percentage...., nec_comp$CTpan)
cor(nec_comp$Markup.Percentage...., nec_comp$CTmvp)
cor(nec_comp$Markup.Percentage...., nec_comp$CT)
cor(nec_comp$Markup.Percentage...., nec_comp$LE)


ang_comp <- dat %>%
inner_join(ang, by = "Case")

cor(ang_comp$Markup.Percentage...., ang_comp$CTpan)
cor(ang_comp$Markup.Percentage...., ang_comp$CTmvp)
cor(ang_comp$Markup.Percentage...., ang_comp$CT)
cor(ang_comp$Markup.Percentage...., ang_comp$LE)


nec_comp2 <- dat %>%
  inner_join(nec2, by = "Case")

cor(nec_comp2$Markup.Percentage...., nec_comp2$CTpan)
cor(nec_comp2$Markup.Percentage...., nec_comp2$CTmvp)
cor(nec_comp2$Markup.Percentage...., nec_comp2$CT)
cor(nec_comp2$Markup.Percentage...., nec_comp2$LE)

nec_comp[,"shape"] <- "txt"
pdf("figures/pathology/necrosis_cor_tcga.pdf", height = 1, width = 1)
ggplot(data = nec_comp, aes(x = Markup.Percentage...., y = CTpan*100)) +
geom_point(aes(shape = shape),size = 1) +
geom_smooth(method="lm",se = FALSE) +
#labs(x = "Necrosis % (Pathology)", y = "CTpan % (CIBERSORT)") +
theme_classic() +
theme(axis.text = element_text(size=7),
axis.title= element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position = "none")
dev.off()


pdf("figures/pathology/angiogenesis_cor_tcga.pdf", height = 2, width = 2)
ggplot(data = ang_comp, aes(x = Markup.Percentage...., y = CTmvp*100)) +
geom_point() +
geom_smooth(method="lm",se = FALSE) +
labs(x = "Angiogenesis % (Pathology)", y = "MVP % (CIBERSORT)") +
theme_classic() +
theme(axis.text = element_text(size=7),
axis.title= element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"))
dev.off()