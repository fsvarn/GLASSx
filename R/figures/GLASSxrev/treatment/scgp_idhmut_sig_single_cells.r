###################################################
# Compare IDH-mutant cell state-specific gene expression in IDH-mutant single-cells
# Author: Frederick Varn
# Date: 2022.01.05
# Figure S4I
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)

#######################################################
rm(list=ls())
# Location of the SCGP single cell data
myinf1 <- "/projects/verhaak-lab/scgp/data/10x/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds"

load(myinf1)

# Remove all cells with no expression
#sums <- apply(log2cpm, 2, sum)	
#range(sums)
# [1]  2453.769 45315.686
# No cells with 0 expression

# Remove QC genes:
qc_genes <- c("ENSGGENES","ENSGUMI","ENSGMITO", "ENSGSEQSAT","ENSGSAMP") 
log2cpm <- log2cpm[-which(rownames(log2cpm) %in% qc_genes),]
featuredata <- featuredata[-which(rownames(featuredata) %in% qc_genes),]

# Annotate clusters using previous definitions
clust_annot = tsne.data %>%
	rownames_to_column('cell') %>%
	mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
						  `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
						  `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) %>%
	column_to_rownames('cell')


		
# Get sample names
sample_id <- sapply(strsplit(rownames(clust_annot), "-"), "[[", 3)
sample_id <- recode(sample_id, "0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006",
					"5" = "SM011", "6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018")
clust_annot[,"sample_id"] <- sample_id
names(sample_id) <- rownames(clust_annot)


# Convert ensembl ID to gene symbol
featuredata <- featuredata[rownames(log2cpm),]

gene <- featuredata[,"Associated.Gene.Name"]
names(gene) <- rownames(featuredata)

rownames(log2cpm) <- gene

#  Read in post-treatment stem cell signature
myDir2 <- "data/res/CIBERSORTx/analysis/"
myinf2 <- dir(myDir2)
myinf2 <- myinf2[grep("idhmut", myinf2)]
myinf2 <- myinf2[grep("_postreatment_result", myinf2)]
mytag <- myinf2[grep("stemcell|differentiated",myinf2)]		# These are the signatures we have some confidence in
myinf2 <- paste(myDir2, mytag, sep = "/")
mytag <- gsub("_postreatment_result.txt","",mytag)
mytag <- gsub("GLASS_idhmut_","",mytag)

p.val <- eff <- n1 <- n2 <- rep(0, length(myinf2))
plotList <- upList <- dnList <- list()

# Treatment (don't use this)
# for(i in 1:length(myinf2))
# {
# 
# 	sigtable <- read.delim(myinf2[i])	
# 	mysig <- rownames(sigtable %>% filter(sig, eff < 0))
# 
# 	mysig <- intersect(mysig, rownames(log2cpm))
# 
# 	g1 <- clust_annot %>%
# 			filter(sample_id %in% c("SM001", "SM002", "SM008"), cell_type == mytag[i]) %>%
# 			rownames(.)
# 	g2 <- clust_annot %>%
# 			filter(sample_id %in% c("SM004", "SM015","UC917"), cell_type == mytag[i]) %>%
# 			rownames(.)
# 
# 	vars <- apply(log2cpm[mysig,c(g1,g2)], 1, var)
# 	if(sum(vars==0) > 0){
# 		mysig <- mysig[-which(vars==0)]}
# 
# 	g1_mean <- apply(log2cpm[mysig,g1], 2, mean)
# 	g2_mean <- apply(log2cpm[mysig,g2], 2, mean)
# 	
# 	score <- c(g1_mean, g2_mean)
# 	status <- c(rep("Recurrent",length(g1_mean)), rep("Initial", length(g2_mean)))
# 	cell_state <- rep(mytag[i], length(score))
# 	plotList[[i]] <- data.frame(score, status, cell_state)
# 
# 	p.val[i] <- wilcox.test(g1_mean, g2_mean)$p.value
# 	eff [i] <- median(g1_mean) - median(g2_mean)
# 
# }
# 
# res <- data.frame(mytag, p.val, eff)

# Grade (use this)
for(i in 1:length(myinf2))
{

	sigtable <- read.delim(myinf2[i])	
	mysig <- rownames(sigtable %>% filter(sig, eff > 0))

	mysig <- intersect(mysig, rownames(log2cpm))

	g1 <- clust_annot %>%
			filter(sample_id %in% c("SM001", "UC917", "SM002", "SM008"), cell_type == mytag[i]) %>%
			rownames(.)
	g2 <- clust_annot %>%
			filter(sample_id %in% c("SM004", "SM015"), cell_type == mytag[i]) %>%
			rownames(.)

	vars <- apply(log2cpm[mysig,c(g1,g2)], 1, var)
	if(sum(vars==0) > 0){
		mysig <- mysig[-which(vars==0)]}

	g1_mean <- apply(log2cpm[mysig,g1], 2, mean)
	g2_mean <- apply(log2cpm[mysig,g2], 2, mean)
	
	score <- c(g1_mean, g2_mean)
	status <- c(rep("Grade III",length(g1_mean)), rep("Grade II", length(g2_mean)))
	cell_state <- rep(mytag[i], length(score))
	upList[[i]] <- data.frame(score, status, cell_state)

	p.val[i] <- wilcox.test(g1_mean, g2_mean)$p.value
	eff [i] <- median(g1_mean) - median(g2_mean)
	
	n1[i] <- length(g1)
	n2[i] <- length(g2)

}

up_res <- data.frame(mytag, p.val, eff)

# We don't focus on the down-regulated signature in the revision
# for(i in 1:length(myinf2))
# {
# 
# 	sigtable <- read.delim(myinf2[i])	
# 	mysig <- rownames(sigtable %>% filter(sig, eff < 0))
# 
# 	mysig <- intersect(mysig, rownames(log2cpm))
# 
# 	g1 <- clust_annot %>%
# 			filter(sample_id %in% c("SM001", "UC917", "SM002", "SM008"), cell_type == mytag[i]) %>%
# 			rownames(.)
# 	g2 <- clust_annot %>%
# 			filter(sample_id %in% c("SM004", "SM015"), cell_type == mytag[i]) %>%
# 			rownames(.)
# 
# 	vars <- apply(log2cpm[mysig,c(g1,g2)], 1, var)
# 	if(sum(vars==0) > 0){
# 		mysig <- mysig[-which(vars==0)]}
# 
# 	g1_mean <- apply(log2cpm[mysig,g1], 2, mean)
# 	g2_mean <- apply(log2cpm[mysig,g2], 2, mean)
# 	
# 	score <- c(g1_mean, g2_mean)
# 	status <- c(rep("Grade III",length(g1_mean)), rep("Grade II", length(g2_mean)))
# 	cell_state <- rep(mytag[i], length(score))
# 	dnList[[i]] <- data.frame(score, status, cell_state)
# 
# 	p.val[i] <- wilcox.test(g1_mean, g2_mean)$p.value
# 	eff [i] <- median(g1_mean) - median(g2_mean)
# 
# }
# 
# dn_res <- data.frame(mytag, p.val, eff)


upTab <- do.call(rbind, upList)

plot_res <- upTab
plot_res <- plot_res %>%
			filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor")) %>%
			mutate(cell_state = recode(cell_state, "differentiated_tumor" = "Diff.-like", "stemcell_tumor" = "Stem-like")) %>%
			mutate(status = fct_relevel(status, "Grade II","Grade III"))
			
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_idhmut_grade_boxplot_v2.pdf", width=1.65,height=1.8)
ggplot(data = plot_res %>% filter(cell_state == "Diff.-like" | cell_state == "Stem-like"), aes(x = cell_state, y = score, fill = status)) +
geom_boxplot(outlier.size=1,colour="black") +
scale_fill_manual(values=c("#27408B", "#CD4F39")) +
labs(y = "Signature score") +
theme_classic() +
theme(plot.title = element_text(size= 7, hjust = 0.5),
	axis.text.x = element_text(size=7, hjust = 0.5),
	axis.text.y = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim=c(0,0.52))
dev.off()


