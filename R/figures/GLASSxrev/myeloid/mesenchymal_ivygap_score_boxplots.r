###################################################
# Create boxplot examining how the myeloid mesenchymal signature is distributed in Ivy GAP features
# Author: Frederick Varn
# Date: 2021.12.15
# Figure 5E
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

#######################################################

rm(list=ls())

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/ivygap_scgp/CIBERSORTxHiRes_ivygap_myeloid_Window48.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/ivygap_scgp/CIBERSORTxHiRes_ivygap_differentiated_tumor_Window48.txt"
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/columns-samples.csv"

# Read in myeloid cell information. Get the TCGA score
geps <- read.delim(myinf1, row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("X","",colnames(geps))

# Check how TCGA signature changes:
tcga_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v4.txt")
tcga_sig <- rownames(tcga_dat %>% filter(sig & eff > log2(1.5)))

tcga_sig <- intersect(rownames(geps), tcga_sig)
mes_score <- apply(geps[tcga_sig,], 2, mean)

# Examine how these variables are regionally distributed
samps <- read.csv(myinf3,stringsAsFactor=FALSE)

#Simplify analysis for now: Reference histology
samp_struct <- samps$structure_name
names(samp_struct) <- samps$rna_well_id

tumor_id <- samps$tumor_id
names(tumor_id) <- samps$rna_well_id

plot_dat <- data.frame(tumor_id, samp_struct, mes_score)
colnames(plot_dat) <- c("tumor_id", "structure", "mes_score")

plot_dat <- plot_dat[grep("reference histology", plot_dat[,"structure"]),]
plot_dat[,"structure"] <- gsub(" sampled by reference histology","",plot_dat[,"structure"])

plot_dat <- plot_dat %>%
mutate(structure = recode(structure, "Cellular Tumor" = "CT", "Leading Edge" = "LE", "Infiltrating Tumor" = "IT",
	"Pseudopalisading cells around necrosis" = "PAN", "Microvascular proliferation" = "MVP")) %>%
mutate(structure = as_factor(structure)) %>%
mutate(structure = fct_relevel(structure, 
"LE", "IT", "CT", "PAN", "MVP")) 

plot_dat %>%
group_by(structure) %>%
summarise(median(mes_score))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_ivygap_mesenchymal_sig_by_region.pdf",width=1.6,height=1.5)
ggplot(plot_dat, aes(structure, mes_score, fill = structure)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
labs(y = "Mesenchymal signature") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7, hjust = 0.5),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=c("LE" = "#009999", "IT"= "#ad4597","CT" = "#01b050", "PAN"="#02ffcc", "MVP"="#c00000")) 
dev.off()
