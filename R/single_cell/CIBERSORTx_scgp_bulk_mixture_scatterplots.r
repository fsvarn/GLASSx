###################################################
# Validate the CIBERSORTx  bulk RNA SCGP proportions against the SCGP single cell proportions 
# Updated: 2020.06.12
# Author: Frederick Varn
##################################################

library(tidyverse)
library(readxl)
library(RColorBrewer)
library(gridExtra)
library(scales)

#myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_bulk_mix_cibersortx_prop_06092020.txt"
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/scgp_bulk/CIBERSORTxGEP_scgp_Fractions-Adjusted.txt"
myinf2 <- "/projects/verhaak-lab/scgp/data/10x/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds"
myinf3 <- "/projects/verhaak-lab/scgp/data/clinical/scgp-subject-metadata.xlsx"

prop <-read.delim(myinf1,row.names=1)

# Convert rownames format to match the single cell format
rownames(prop) <- sapply(strsplit(rownames(prop),"\\."),function(x)paste(x[2],x[3],sep=""))

load(myinf2)
sample_metadata <- read_excel(myinf3)

# Annotate clusters using previous definitions
clust_annot = tsne.data %>%
	rownames_to_column('cell') %>%
	mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                              `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                              `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) %>%
	column_to_rownames('cell')

# Rename sample_id using the information in the sample_metadata file
sample_id <- sapply(strsplit(rownames(clust_annot), "-"), "[[", 3)
sample_id <- recode("0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006", "5" = "SM011",
					"6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018", sample_id)
clust_annot[,"sample_id"] <- sample_id

# Calculate the true proportions in the single cell data:
uni_cell_type <- unique(clust_annot[,"cell_type"])
uni_samp <- unique(sample_id)
truth <- matrix(0, nrow = length(uni_samp), ncol = length(uni_cell_type))
rownames(truth) <- uni_samp
colnames(truth) <- uni_cell_type

for(i in 1:length(uni_samp))
{
	sub_annot <- clust_annot[which(clust_annot[,"sample_id"] == uni_samp[i]),]
	cell_counts <- table(sub_annot[,"cell_type"])
	cell_props <- cell_counts/nrow(sub_annot)
	truth[i,names(cell_props)] <- cell_props
}

# Get rid of unnecessary columns and rows that did not have both types of data
comxx <- intersect(colnames(prop), colnames(truth))
prop <- prop[,comxx]
truth <- truth[,comxx]

comxx <- intersect(rownames(prop), rownames(truth))
prop <- prop[comxx,]
truth <- truth[comxx,]

# Get correlation coefficients
mycor <- rep(0, ncol(prop))
names(mycor) <- colnames(prop)
for(i in 1:ncol(prop))
{
	mycor[i] <- round(cor(prop[,i], truth[,i], method="p"),2)
}


# Recalculate the proportions based on tumor and immune sorting
tumor_prop <- prop[,grep("_tumor",colnames(prop))]
tumor_truth <- truth[,grep("_tumor",colnames(truth))]

tumor_prop <- t(apply(tumor_prop, 1, function(x)(x/sum(x))))
tumor_truth <- t(apply(tumor_truth, 1, function(x)(x/sum(x))))

tumor_cor <- rep(0, ncol(tumor_prop))
names(tumor_cor) <- colnames(tumor_prop)
for(i in 1:ncol(tumor_prop))
{
	tumor_cor[i] <- round(cor(tumor_prop[,i], tumor_truth[,i], method="p"),2)
}

# Correlate malignant cells in long format
tfull_prop <- tumor_prop %>%  
data.frame() %>%  
rownames_to_column("sample") %>%
pivot_longer(!sample, names_to = "cell_state", values_to = "csx_fraction")

tfull_truth <- tumor_truth %>%
data.frame() %>%  
rownames_to_column("sample") %>%
pivot_longer(!sample, names_to = "cell_state", values_to = "true_fraction")

cor(tfull_truth[,"true_fraction"], tfull_prop[,"csx_fraction"])

# Plot the full value

tfull_plot <- tfull_prop %>%
			 inner_join(tfull_truth) %>%
			 mutate(cell_state = recode(cell_state, 
					"differentiated_tumor" = "Diff.-like", "prolif_stemcell_tumor" = "Prolif. stem-like",
					"stemcell_tumor" = "Stem-like")) %>%
			 mutate(cell_state = as_factor(cell_state)) %>%
			 mutate(cell_state = fct_relevel(cell_state, "Prolif. stem-like","Stem-like","Diff.-like")) 						 
						 
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scrna_csx_proportions_full_tumor.pdf",width=1.9,height=1.5)  
ggplot(tfull_plot, aes(x = true_fraction * 100, y = csx_fraction * 100)) +
geom_point(aes(colour=cell_state)) +
geom_smooth(method = lm, se = FALSE) + 
scale_colour_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
theme_classic() +
labs(x = "scRNAseq (%)", y = "CIBERSORTx (%)") +
theme(axis.text = element_text(size=7),
axis.title = element_text(size=7), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.background = element_blank(),
legend.position="none") 
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scrna_csx_proportions_tumor.pdf",width=3,height=1.3)  
ggplot(tfull_plot, aes(x = true_fraction * 100, y = csx_fraction * 100)) +
geom_point() +
geom_smooth(method = lm, se = FALSE, aes(colour = cell_state)) + 
facet_wrap(.~cell_state, scales = "free") + 
scale_colour_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
theme_classic() +
labs(x = "scRNAseq (%)", y = "CIBERSORTx (%)") +
theme(axis.text = element_text(size=7),
axis.title = element_text(size=7), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") 
dev.off()




# Recalculate the proportions based on tumor and immune sorting
cd45 <- c("myeloid", "t_cell", "granulocyte", "dendritic_cell", "b_cell")
normal_prop <- prop[,cd45]
normal_truth <- truth[,cd45]

normal_prop <- t(apply(normal_prop, 1, function(x)(x/sum(x))))
normal_truth <- t(apply(normal_truth, 1, function(x)(x/sum(x))))

normal_cor <- rep(0, ncol(normal_prop))
names(normal_cor) <- colnames(normal_prop)
for(i in 1:ncol(normal_prop))
{
	normal_cor[i] <- round(cor(normal_prop[,i], normal_truth[,i], method="s"),2)
}

# Correlate everything in long format
full_prop <- prop %>%  
rownames_to_column("sample") %>%
pivot_longer(!sample, names_to = "cell_state", values_to = "csx_fraction")


full_truth <- truth %>%
data.frame() %>%  
rownames_to_column("sample") %>%
pivot_longer(!sample, names_to = "cell_state", values_to = "true_fraction")

cor(full_truth[,"true_fraction"], full_prop[,"csx_fraction"])

# Plot the full value

full_plot <- full_prop %>%
			 inner_join(full_truth) %>%
			 mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
					"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
					"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
					"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
					"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
					"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
			 mutate(cell_state = as_factor(cell_state)) %>%
			 mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
											"Oligodendrocyte", 
											"Endothelial", "Pericyte",
											"Fibroblast", 
											"Diff.-like", "Stem-like", "Prolif. stem-like")) 						 
						 
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scrna_csx_proportions_full.pdf",width=1.9,height=1.5)  
ggplot(full_plot, aes(x = true_fraction * 100, y = csx_fraction * 100)) +
geom_point(aes(colour=cell_state)) +
geom_smooth(method = lm, se = FALSE) + 
scale_colour_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
theme_classic() +
labs(x = "scRNAseq (%)", y = "CIBERSORTx (%)") +
theme(axis.text = element_text(size=7),
axis.title = element_text(size=7), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.background = element_blank(),
legend.position="none") 
dev.off()




#colors <- brewer.pal(12,"Paired")
colors <- c("stemcell_tumor" = "#fb6a4a", "myeloid" = "#08519c", "differentiated_tumor" = "#fcbba1", "prolif_stemcell_tumor" = "#a50f15",
			"oligodendrocyte" = "#2ca25f", "t_cell" = "#6baed6", "granulocyte" = "#bdd7e7", "pericyte" = "#fee391", "endothelial" = "#ffffd4", 
			"dendritic_cell" = "#3182bd", "fibroblast" = "#feb24c", "b_cell" = "#eff3ff")

# Create scatterplots
se <- list()
for(i in 1:ncol(prop))
{
	cibersortx_prop <- prop[,i] * 100
	true_prop <- truth[,i] * 100
	plot_cor <- data.frame(cibersortx_prop, true_prop)
	
	mycoef <- format(mycor[i], nsmall = 2)
	names(mycoef) <- NULL
	myr <- deparse((bquote(italic("R") ~" = " ~ .(mycoef))))
	annotation_text <- data.frame(true_prop = 0.75 * max(true_prop), 
	cibersortx_prop = 0.02 * max(cibersortx_prop),
	myr)
	
	cell <- colnames(prop)[i]
	cell <- gsub("_"," ",cell)
	substr(cell,1,1) <- toupper(substr(cell,1,2))
	
	se[[i]] <- ggplot(plot_cor, aes(x=true_prop,y=cibersortx_prop)) +
	geom_point(colour="black",size=1) +
	geom_smooth(colour=colors[i],method='lm',se=FALSE) +
	geom_text(data=annotation_text,label=myr, size=2.5, parse=TRUE) +
	labs(title = cell, x = "True fraction (%)", y = "CIBERSORTx fraction (%)") +
	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.title=element_text(size=7,hjust=0.5),
	axis.title=element_text(size=7),
	axis.text.x=element_text(size=7),
	axis.text.y=element_text(size=7),
	legend.position="none") 
}

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/10x_scgp_cibersortx_bulk_mix_validation.pdf",width=6.5,height=5.2)
grid.arrange(se[[1]],se[[2]],se[[3]],se[[4]],se[[5]],se[[6]],se[[7]],se[[8]],se[[9]],se[[10]],se[[11]],se[[12]],nrow=3)
dev.off()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/10x_scgp_cibersortx_bulk_mix_tumor_validation.pdf",width=4.2,height=1.5)
grid.arrange(se[[3]],se[[1]],se[[4]],nrow=1)
dev.off()


# Look at subtype

myinf4 <- "/projects/verhaak-lab/SCGP-analysis/results/kallisto/kallisto/final/p_result_gene_tpm_matrix_all_samples.gct.txt"
subtype <- read.delim(myinf4)

rownames(subtype) <- sapply(strsplit(rownames(subtype),"\\."),function(x)paste(x[2],x[3],sep=""))
subtype <- subtype[rownames(prop),]

subtype_cor <- cor(subtype[,1:3], prop, method = "p")