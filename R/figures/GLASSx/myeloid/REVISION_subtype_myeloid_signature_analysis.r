###################################################
# Create scatterplot correlating the mesenchymal myeloid score with bulk mesenchymal tumor score
# Updated: 2021.01.14
# Author: Frederick Varn
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
	filter(dbCluster == 2) %>%
	column_to_rownames('cell')
	
log2cpm_myeloid <- log2cpm[,rownames(clust_annot)]
	
# 	mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
#                               `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
#                               `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) %>%
# 	column_to_rownames('cell')

# clust_annot file and log2cpm files are in the same order so no need to match the order between them
# cell_type <- clust_annot[,"cell_type"]

# Assign each cell a subtype
# log2cpm_annot <- log2cpm
# colnames(log2cpm_annot) <- cell_type

# Get sample names
sample_id <- sapply(strsplit(rownames(clust_annot), "-"), "[[", 3)
sample_id <- recode(sample_id, "0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006",
					"5" = "SM011", "6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018")
clust_annot[,"sample_id"] <- sample_id
names(sample_id) <- rownames(clust_annot)

# Convert ensembl ID to gene symbol
featuredata <- featuredata[rownames(log2cpm_myeloid),]

gene <- featuredata[,"Associated.Gene.Name"]
names(gene) <- rownames(featuredata)

rownames(log2cpm_myeloid) <- gene

tcga_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_nonmes_myeloid_ts_result.txt", stringsAsFactor=FALSE)

sigs <- unique(tcga_dat[,"analysis"])
cor_coef <- cor_coef.s <- p.value <- rep(0, length(sigs))
all_sigs <- list()
names(cor_coef) <- sigs
for(i in 1:length(sigs))
{
	
	cat("\r",i)
	sub_dat <- tcga_dat %>% filter(analysis == sigs[i])
	sub_sig <- sub_dat %>% filter(q.val < 0.05 & eff > log2(1.1)) %>% .$gene

	sc_sig_score <- apply(log2cpm_myeloid[sub_sig,],2,mean)
	sample_id <- clust_annot[names(sc_sig_score),"sample_id"]
	plot_res <- data.frame(names(sc_sig_score), sc_sig_score, sample_id)
	colnames(plot_res) <- c("cell_id","sc_sig_score","sample_id")


	# Read in bulk subtype results
	con <- DBI::dbConnect(odbc::odbc(), "scgp")

	bulk_res <- dbReadTable(con, Id(schema = "analysis", table="transcriptional_subtype"))
	
	tsub <- ifelse(sigs[i] == "Non-mesenchymal", "Mesenchymal", sigs[i])
	bulk_sig <- bulk_res %>% 
				filter(signature_name == tsub) %>%
				mutate(sample_id = sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[2],x[3],sep=""))) #%>%
				#filter(sample_id %in% c("SM006","SM012","SM017","SM018","SM011"))

	plot_res <- plot_res %>% 
	group_by(sample_id) %>% 
	summarise(mean = mean(sc_sig_score)) %>% 
	ungroup() %>%
	inner_join(bulk_sig, by = "sample_id") %>%
	as.data.frame()

	cor_coef[i] <- cor(plot_res$mean, plot_res$enrichment_score)
	cor_coef.s[i] <- cor(plot_res$mean, plot_res$enrichment_score,method="s")
	p.value[i] <- cor.test(plot_res$mean, plot_res$enrichment_score)$p.value
	
	all_sigs[[i]] <- sub_sig
}
res <- data.frame(cor_coef, cor_coef.s, p.value)
names(all_sigs) <- sigs

# Compare to mesenchymal

mes_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v3.txt", stringsAsFactor=FALSE)
mes_sig <- rownames(mes_dat %>% filter(q.val < 0.05 & eff > log2(1.1)))

sum(mes_sig %in% all_sigs[[1]])	#0
sum(mes_sig %in% all_sigs[[2]])	#0
sum(mes_sig %in% all_sigs[[3]])	#0


# Correlate at single-cell level
non_sig <- tcga_dat %>% filter(analysis == sigs[3], q.val < 0.05 & eff > log2(1.1)) %>% .$gene
mes_sig <- rownames(mes_dat %>% filter(q.val < 0.05 & eff > log2(1.1)))


non_score <- apply(log2cpm_myeloid[non_sig,],2,mean)
mes_score <- apply(log2cpm_myeloid[mes_sig,],2,mean)

cor.test(non_score, mes_score, method="p")
plot_res <- data.frame(non_score, mes_score)

# Compare proneural vs mesenchymal and classical vs mesenchymal and examine overlap. Others (mesenchymal vs proneural/mesenchymal vs classical?)


# Look at ligands

lig_rec <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/ramilowski_2015/receptor_ligand_pairs_fantom5_ramilowski_ncomms_2015.txt",stringsAsFactor=FALSE)
lig_rec <- lig_rec %>%
filter(lig_rec[,"Pair.Evidence"] == "literature supported")
lig_rec <- lig_rec[,c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]

png("/projects/verhaak-lab/GLASS-III/figures/analysis/myeloid_subtype_scatterplot.png")#,width=1.9,height=1.5)  
ggplot(plot_res, aes(x = non_score, y = mes_score)) +
geom_point() +
geom_smooth(method = lm, se = FALSE) + 
theme_classic() +
labs(x = "Mean myeloid cell mes. score", y = "Bulk mes. score") +
theme(axis.text = element_text(size=7),
axis.title = element_text(size=7), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.background = element_blank(),
legend.position="none") #+
#coord_cartesian(xlim=c(-35, 35))
dev.off()
