
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

# Read in myeloid cell information. Get the OSM ligand expression and TCGA score
geps <- read.delim(myinf1, row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("X","",colnames(geps))

osm <- as.numeric(geps["OSM",])
names(osm) <- colnames(geps)

# Check how TCGA signature changes:
tcga_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v3.txt")
tcga_sig <- rownames(tcga_dat %>% filter(q.val < 0.05 & eff > log2(1.1)))

tcga_sig <- intersect(rownames(geps), tcga_sig)
mes_score <- apply(geps[tcga_sig,], 2, mean)

# Get the OSMR receptor information from differentiated tumors
geps <- read.delim(myinf2, row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("X","",colnames(geps))

osmr <- as.numeric(geps["OSMR",])
names(osmr) <- colnames(geps)

# Examine how these variables are regionally distributed
samps <- read.csv(myinf3,stringsAsFactor=FALSE)

#Simplify analysis for now: Reference histology
samp_struct <- samps$structure_name
names(samp_struct) <- samps$rna_well_id

tumor_id <- samps$tumor_id
names(tumor_id) <- samps$rna_well_id

plot_dat <- data.frame(tumor_id, samp_struct, mes_score, osm, osmr)
colnames(plot_dat) <- c("tumor_id", "structure", "mes_score", "osm", "osmr")

plot_dat <- plot_dat[grep("reference histology", plot_dat[,"structure"]),]
plot_dat[,"structure"] <- gsub(" sampled by reference histology","",plot_dat[,"structure"])

plot_dat <- plot_dat %>%
mutate(structure = recode(structure, "Cellular Tumor" = "CT", "Leading Edge" = "LE", "Infiltrating Tumor" = "IT",
	"Pseudopalisading cells around necrosis" = "PAN", "Microvascular proliferation" = "MVP")) %>%
mutate(structure = as_factor(structure)) %>%
mutate(structure = fct_relevel(structure, 
"LE", "IT", "CT", "PAN", "MVP")) 
# mutate(structure = recode(structure, "Cellular Tumor" = "Cellular tumor", "Leading Edge" = "Leading edge", "Infiltrating Tumor" = "Infiltrating tumor")) %>%
# mutate(structure = as_factor(structure)) %>%
# mutate(structure = fct_relevel(structure, 
# "Leading edge", "Infiltrating tumor", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation"))

plot_dat %>%
group_by(structure) %>%
summarise(median(mes_score))

plot_dat %>%
group_by(structure) %>%
summarise(median(osm))

plot_dat %>%
group_by(structure) %>%
summarise(median(osmr))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_ivygap_mesenchymal_sig_by_region.pdf",width=1.93,height=1.5)
ggplot(plot_dat, aes(structure, mes_score, fill = structure)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
labs(y = "Mesenchymal signature") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7, angle = 45, hjust = 1),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=c("LE" = "#009999", "IT"= "#ad4597","CT" = "#01b050", "PAN"="#02ffcc", "MVP"="#c00000")) 
dev.off()

plot_dat2 <- plot_dat %>% 
pivot_longer(c(osm, osmr), names_to = "gene", values_to = "expr") %>%
mutate(gene = recode(gene, "osm" = "OSM", "osmr" = "OSMR"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_ivygap_osm_osmr_by_region.pdf",width=3.5,height=4)
ggplot(plot_dat2, aes(structure, expr, fill = structure)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
labs(y = "Inferred expression") +
facet_grid(gene~.) +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7, angle = 45, hjust = 1),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") +
scale_fill_manual(values=c("Leading edge" = "#009999", "Infiltrating tumor"= "#ad4597","Cellular tumor" = "#01b050", "Pseudopalisading cells around necrosis"="#02ffcc", "Microvascular proliferation"="#c00000"))  +
coord_cartesian(ylim = c(0,2))
dev.off()


wilcox.test(plot_dat %>% filter(structure == "Cellular tumor") %>% .$osmr, 
			plot_dat %>% filter(structure == "Microvascular proliferation") %>% .$osm)