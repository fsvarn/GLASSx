###################################################
# Digital cytometry vs flow cytometry signature heatmap
# Author: Frederick Varn
# Date: 2022.01.05
# Figures S4B, S4C
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(viridis)

#######################################################
rm(list=ls())

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze_validation_rev/CIBERSORTxGEP_GLASS_SM_GEPs_HeatMap.txt"

mat <- read.delim(myinf1)

colnames(mat)[1] <- "gene"
colnames(mat)[2:8] <- paste("cibersortx_", colnames(mat)[2:8], sep = "")
colnames(mat)[9:15] <- paste("klemm_", colnames(mat)[9:15], sep = "")
colnames(mat) <- gsub("\\.1","",colnames(mat))

sub_mat <- mat[,c(1:5,9,10)]
colnames(sub_mat)[which(colnames(sub_mat) == "klemm_stemcell_tumor")] <- "klemm_tumor"

# Get order for genes
ord1 <- as.character(sub_mat$gene)
ord1 <- ord1[order(sub_mat[,"cibersortx_differentiated_tumor"], decreasing = TRUE)]
ord1 <- ord1[which(!is.na(sub_mat[ord1,"cibersortx_differentiated_tumor"]))]

ord2 <- as.character(sub_mat$gene)
ord2 <- ord2[order(sub_mat[,"cibersortx_stemcell_tumor"], decreasing = TRUE)]
ord2 <- ord2[which(!is.na(sub_mat[ord2,"cibersortx_stemcell_tumor"]))]
ord2 <- ord2[-which(ord2 %in% ord1)]

ord <- c(ord1, ord2)

ord3 <- as.character(sub_mat$gene)
ord3 <- ord3[order(sub_mat[,"cibersortx_prolif_stemcell_tumor"], decreasing = TRUE)]
ord3 <- ord3[which(!is.na(sub_mat[ord3,"cibersortx_prolif_stemcell_tumor"]))]
ord3 <- ord3[-which(ord3 %in% ord)]

ord <- c(ord, ord3)

ord4 <- as.character(sub_mat$gene)
ord4 <- ord4[order(sub_mat[,"cibersortx_myeloid"], decreasing = TRUE)]
ord4 <- ord4[which(!is.na(sub_mat[ord4,"cibersortx_myeloid"]))]
ord4 <- ord4[-which(ord4 %in% ord)]

ord <- c(ord, ord4)

remain <- as.character(sub_mat$gene[which(!(sub_mat$gene %in% ord))])
ord <- c(ord, remain)

plot_dat <- sub_mat %>%
			pivot_longer(-gene, names_to = "gep", values_to = "expression") %>%
			#mutate(expression = replace_na(expression, replace = 0)) %>%
			mutate(gene = fct_relevel(gene, rev(ord))) %>%
			mutate(dataset = sapply(strsplit(as.character(gep), "_"),function(x)x[1])) %>%
			mutate(dataset = recode(dataset, "cibersortx" = "CIBERSORTx", "klemm" = "Klemm et al")) %>%
			mutate(gep = recode(gep, "cibersortx_differentiated_tumor" = "Diff.-like", "cibersortx_stemcell_tumor" = "Stem-like",
										  "cibersortx_prolif_stemcell_tumor" = "Prolif. stem-like","cibersortx_myeloid" = "Myeloid",
										  "klemm_tumor" = "CD45neg", klemm_myeloid = "Myeloid")) %>%
			mutate(gep = fct_relevel(gep, c("Diff.-like", "Stem-like", 
											"Prolif. stem-like", "CD45neg", "Myeloid")))

thr <- quantile(plot_dat$expression, .99, na.rm=TRUE) 

plot_dat <- plot_dat %>% mutate(expression = replace_na(expression, replace = 0))
plot_dat[which(plot_dat$expression > thr),"expression"] <- thr

# Exclude genes that were not imputed in any of the malignant or myeloid cell states
sub_mat2 <- sub_mat[,2:7]
sub_mat2 <- data.matrix(sub_mat2)
sub_mat2[which(is.na(sub_mat2))] <- 0
non_tumor_myeloid <- apply(sub_mat2[,1:4],1,sum)
non_tumor_myeloid <- names(non_tumor_myeloid[which(non_tumor_myeloid > 0)])

plot_dat <- plot_dat %>% filter(gene %in% non_tumor_myeloid)			

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/csx_klemm_comparison.pdf",width=1.5,height=2)
ggplot(data = plot_dat, aes(x = gep, y = gene)) +
geom_tile(aes(fill=expression)) +
scale_fill_viridis() + 
facet_grid(.~dataset, scale = "free_x", space = "free_x") +
theme_void() +
theme(
axis.text.x = element_text(size=7, angle = 45, hjust = 1, vjust =1),
axis.text.y = element_blank(),
axis.title= element_blank(),
axis.ticks = element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
strip.text = element_text(size=7),
strip.background = element_blank(),
legend.position = "none")
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/csx_klemm_comparison_legend.pdf",width=1.5,height=2)
ggplot(data = plot_dat, aes(x = gep, y = gene)) +
geom_tile(aes(fill=expression)) +
scale_fill_viridis() + 
facet_grid(.~dataset, scale = "free_x", space = "free_x") +
theme_void() +
theme(
axis.text.x = element_text(size=7, angle = 45, hjust = 1, vjust =1),
axis.text.y = element_blank(),
axis.title= element_blank(),
axis.ticks = element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
strip.text = element_text(size=7),
strip.background = element_blank(),
legend.position = "right")
dev.off()

cor_matrix <- cor(sub_mat2,method="p")
cor_matrix <- cor_matrix[grep("cibersortx_",rownames(cor_matrix)),]
cor_matrix <- cor_matrix[,grep("klemm_",colnames(cor_matrix))]

plot_cor <- cor_matrix %>%
			data.frame() %>%
			rownames_to_column("cibersortx") %>%
			pivot_longer(-cibersortx, names_to = "klemm", values_to = "cor") %>%
			mutate(cibersortx = recode(cibersortx, "cibersortx_differentiated_tumor" = "Diff.-like", "cibersortx_stemcell_tumor" = "Stem-like",
									  "cibersortx_prolif_stemcell_tumor" = "Prolif. stem-like","cibersortx_myeloid" = "Myeloid")) %>%
			mutate(klemm = recode(klemm, "klemm_tumor" = "CD45neg", klemm_myeloid = "Myeloid")) %>%
			mutate(cibersortx = fct_relevel(cibersortx, "Myeloid", "Prolif. stem-like","Stem-like","Diff.-like")) %>%
			mutate(klemm = fct_relevel(klemm, "CD45neg","Myeloid"))
									  	

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/csx_klemm_cor.pdf",width=1.1,height=1)
ggplot(data = plot_cor, aes(x = klemm, y = cibersortx)) +
geom_tile(aes(fill=cor)) +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0,  limits = c(-1, 1),
space = "Lab", na.value = "grey50", guide = "colourbar",
aesthetics = "fill") +
theme_void() +
theme(
axis.text.x = element_text(size=7, angle = 45, hjust = 1, vjust =1),
axis.text.y = element_text(size=7, hjust = 1),
axis.title= element_blank(),
axis.ticks = element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
strip.text = element_text(size=7),
strip.background = element_blank(),
legend.position = "none")
dev.off()


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/csx_klemm_cor_legend.pdf",width=1.1,height=1)
ggplot(data = plot_cor, aes(x = klemm, y = cibersortx)) +
geom_tile(aes(fill=cor)) +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0,  limits = c(-1, 1),
space = "Lab", na.value = "grey50", guide = "colourbar",
aesthetics = "fill") +
theme_void() +
theme(
axis.text.x = element_text(size=7, angle = 45, hjust = 1, vjust =1),
axis.text.y = element_text(size=7, hjust = 1),
axis.title= element_blank(),
axis.ticks = element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
strip.text = element_text(size=7),
strip.background = element_blank(),
legend.position = "right")
dev.off()