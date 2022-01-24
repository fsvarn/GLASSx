##################################################
# Compare recurrent stem-like neoplastic signature in core vs rim single neoplastic cells
# Author: Frederick Varn
# Date: 2022.01.10
# Figure 4D
##################################################

library(tidyverse)
library(Seurat)

#######################################################
rm(list=ls())
set.seed(11)

# Read in core/rim data

mydat <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/garofano_2021/GSE117891_all_6148.umi.count.matrix.tsv", row.names = 1, stringsAsFactor = FALSE)
myseurat <- CreateSeuratObject(counts = mydat)
myseurat <- NormalizeData(myseurat, normalization.method = "LogNormalize", scale.factor = 10000)
tpm <- myseurat[["RNA"]]@data
tpm <- data.frame(tpm)
info <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/garofano_2021/GSE117891_info.txt", stringsAsFactor=FALSE)
tumor.cells <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/garofano_2021/tumor_cells.txt",stringsAsFactor=FALSE)
tpm <- tpm[,tumor.cells[,1]]
cell_source <- sapply(strsplit(colnames(tpm), "_"), function(x)x[1])

# Read in signatures

myDir1 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/analysis/"
myinf1 <- dir(myDir1)
myinf1 <- myinf1[grep("GLASS_idhwt_", myinf1)]
myinf1 <- paste(myDir1, myinf1, sep = "")

diff_sig <- read.delim(myinf1[1])
stem_sig <- read.delim(myinf1[3])

diff_sig <- rownames(diff_sig %>% filter(q.val < 0.05, eff > 0))
stem_sig <- rownames(stem_sig %>% filter(q.val < 0.05, eff > 0))
inter_sig <- intersect(diff_sig, stem_sig)

stem_sig_nona <- intersect(stem_sig, rownames(tpm))
stem_tpm <- tpm[stem_sig_nona,]
stem_tpm <- apply(stem_tpm, 2, mean)

diff_sig_nona <- intersect(diff_sig, rownames(tpm))
diff_tpm <- tpm[diff_sig_nona,]
diff_tpm <- apply(diff_tpm, 2, mean)

inter_sig_nona <- intersect(inter_sig, rownames(tpm))
inter_tpm <- tpm[inter_sig_nona,]
inter_tpm <- apply(inter_tpm, 2, mean)

#############################################

# Compare the stem-like sig differences
#---------------------------------------

info[,"tumor_id"] <- substring(info$Biopsy,1, nchar(info$Biopsy)-2)
tumor_id <- unique(info$tumor_id)
p.val <- eff <- n1 <- n2 <- rep(0, length(tumor_id))
for(i in 1:length(tumor_id))
{
	sub_info <- info[grep(paste(tumor_id[i],"P",sep=""), info$Biopsy),]
	core <- sub_info %>% filter(Tumor.location == "Core") %>% .$Biopsy
	rim <- sub_info %>% filter(Tumor.location == "Rim") %>% .$Biopsy

	g1 <- stem_tpm[which(cell_source %in% core)]
	g2 <- stem_tpm[which(cell_source %in% rim)]
	
	p.val[i] <- wilcox.test(g1,g2)$p.value
	eff[i] <- log2(mean(g2)/mean(g1))
	n1[i] <- length(g1)
	n2[i] <- length(g2)
}
res <- data.frame(tumor_id, p.val, eff, n1, n2)
sig_res <- res %>% filter(p.val < 0.05)
sig_res

# Across all cells
core <- info %>% filter(Tumor.location == "Core") %>% .$Biopsy
rim <- info %>% filter(Tumor.location == "Rim") %>% .$Biopsy

g1 <- stem_tpm[which(cell_source %in% core)]
g2 <- stem_tpm[which(cell_source %in% rim)]

full_p.val <- wilcox.test(g1,g2)$p.value
full_eff <- log2(mean(g2)/mean(g1))
data.frame(full_p.val, full_eff)

# Plot intersection data
cell_name <- c(names(g1), names(g2))
score <- c(g1, g2)
region <- c(rep("Core", length(g1)), rep("Rim", length(g2)))
plot_res <- data.frame(cell_name, score, region)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/garofano_core_vs_rim_stem_neuronal.pdf",width=1.75, height=2.2)
ggplot(plot_res, aes(region, score)) + 
geom_violin(aes(fill=region))  +
geom_boxplot(outlier.size = 0.1, width = 0.5)  +
scale_fill_manual(values = c("Rim" = "#009999", "Core" = "#01b050")) +
labs(y = "Signature score") +
theme_bw() +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x =  element_text(size=7),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size = 7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
coord_cartesian(ylim = c(0,0.22))
dev.off()

# Repeat with diff_sig
#---------------------------------------

info[,"tumor_id"] <- substring(info$Biopsy,1, nchar(info$Biopsy)-2)
tumor_id <- unique(info$tumor_id)
p.val <- eff <- n1 <- n2 <- rep(0, length(tumor_id))
for(i in 1:length(tumor_id))
{
	sub_info <- info[grep(paste(tumor_id[i],"P",sep=""), info$Biopsy),]
	core <- sub_info %>% filter(Tumor.location == "Core") %>% .$Biopsy
	rim <- sub_info %>% filter(Tumor.location == "Rim") %>% .$Biopsy

	g1 <- diff_tpm[which(cell_source %in% core)]
	g2 <- diff_tpm[which(cell_source %in% rim)]
	
	p.val[i] <- wilcox.test(g1,g2)$p.value
	eff[i] <- log2(mean(g2)/mean(g1))
	n1[i] <- length(g1)
	n2[i] <- length(g2)
}
res <- data.frame(tumor_id, p.val, eff, n1, n2)
sig_res <- res %>% filter(p.val < 0.05)
sig_res

# Across all cells
core <- info %>% filter(Tumor.location == "Core") %>% .$Biopsy
rim <- info %>% filter(Tumor.location == "Rim") %>% .$Biopsy

g1 <- diff_tpm[which(cell_source %in% core)]
g2 <- diff_tpm[which(cell_source %in% rim)]

full_p.val <- wilcox.test(g1,g2)$p.value
full_eff <- log2(mean(g2)/mean(g1))
data.frame(full_p.val, full_eff)

# Repeat with intersection sig
#---------------------------------------

info[,"tumor_id"] <- substring(info$Biopsy,1, nchar(info$Biopsy)-2)
tumor_id <- unique(info$tumor_id)
p.val <- eff <- n1 <- n2 <- rep(0, length(tumor_id))
for(i in 1:length(tumor_id))
{
	sub_info <- info[grep(paste(tumor_id[i],"P",sep=""), info$Biopsy),]
	core <- sub_info %>% filter(Tumor.location == "Core") %>% .$Biopsy
	rim <- sub_info %>% filter(Tumor.location == "Rim") %>% .$Biopsy

	g1 <- inter_tpm[which(cell_source %in% core)]
	g2 <- inter_tpm[which(cell_source %in% rim)]
	
	p.val[i] <- wilcox.test(g1,g2)$p.value
	eff[i] <- log2(mean(g2)/mean(g1))
	n1[i] <- length(g1)
	n2[i] <- length(g2)
}
res <- data.frame(tumor_id, p.val, eff, n1, n2)
sig_res <- res %>% filter(p.val < 0.05)
sig_res

# Across all cells
core <- info %>% filter(Tumor.location == "Core") %>% .$Biopsy
rim <- info %>% filter(Tumor.location == "Rim") %>% .$Biopsy

g1 <- inter_tpm[which(cell_source %in% core)]
g2 <- inter_tpm[which(cell_source %in% rim)]

full_p.val <- wilcox.test(g1,g2)$p.value
full_eff <- log2(mean(g2)/mean(g1))
data.frame(full_p.val, full_eff)

# Plot intersection data
cell_name <- c(names(g1), names(g2))
score <- c(g1, g2)
region <- c(rep("Core", length(g1)), rep("Rim", length(g2)))
plot_res <- data.frame(cell_name, score, cell_source)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/garofano_core_vs_rim_neuronal.pdf",width=1.75, height=2.56)
ggplot(plot_res, aes(region, score, fill = cell_source)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") 
dev.off()