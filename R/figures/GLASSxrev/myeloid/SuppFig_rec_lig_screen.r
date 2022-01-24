###################################################
# Identify ligand-receptor candidate pairs associated with mesenchymal transition
# Author: Frederick Varn
# Date: 2022.01.08
# Table S9
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name, ts1.enrichment_score AS score_a, ts2.enrichment_score AS score_b, ts2.enrichment_score - ts1.enrichment_score AS change
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE cs.idh_codel_subtype = 'IDHwt' AND ts1.signature_name ='Mesenchymal' AND ts2.signature_name = 'Mesenchymal'
ORDER BY 1"

dat <- dbGetQuery(con,q)

# Read in receptor-ligand interactions

lig_rec <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/ramilowski_2015/receptor_ligand_pairs_fantom5_ramilowski_ncomms_2015.txt",stringsAsFactor=FALSE)
lig_rec <- lig_rec %>%
filter(lig_rec[,"Pair.Evidence"] == "literature supported")
lig_rec <- lig_rec[,c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]

# Myeloid analysis
geps <- read.delim(myinf1[2], row.names=1)
#geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))


lig_my <- unique(intersect(lig_rec[,"Ligand.ApprovedSymbol"], rownames(geps)))
rec_my <- unique(intersect(lig_rec[,"Receptor.ApprovedSymbol"], rownames(geps)))
genes_my <- c(lig_my, rec_my)

my_cor <- my_p <- rep(length(genes_my))
for(i in 1:length(genes_my))
{
	g1 <- as.numeric(geps[genes_my[i],dat[,"tumor_barcode_a"]])
	g2 <- as.numeric(geps[genes_my[i],dat[,"tumor_barcode_b"]])
	tpm_change <- g2 - g1

	my_cor[i] <- cor(tpm_change, dat$change, method="p")
	my_p[i] <-  cor.test(tpm_change, dat$change, method="p")$p.value

}
names(my_cor) <- names(my_p) <- genes_my

# Differentiated-like tumor analysis
geps <- read.delim(myinf1[1], row.names=1)
#geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))


lig_diff <- unique(intersect(lig_rec[,"Ligand.ApprovedSymbol"], rownames(geps)))
rec_diff <- unique(intersect(lig_rec[,"Receptor.ApprovedSymbol"], rownames(geps)))
genes_diff <- c(lig_diff, rec_diff)

diff_cor <- diff_p <- rep(length(genes_diff))
for(i in 1:length(genes_diff))
{
	g1 <- as.numeric(geps[genes_diff[i],dat[,"tumor_barcode_a"]])
	g2 <- as.numeric(geps[genes_diff[i],dat[,"tumor_barcode_b"]])
	tpm_change <- g2 - g1

	diff_cor[i] <- cor(tpm_change, dat$change, method="p")
	diff_p[i] <-  cor.test(tpm_change, dat$change, method="p")$p.value

}
names(diff_cor) <- names(diff_p) <- genes_diff

ligand <- lig_rec[,"Ligand.ApprovedSymbol"]
receptor <- lig_rec[,"Receptor.ApprovedSymbol"]

my_lig_cor <- my_cor[ligand]
my_lig_p <- my_p[ligand]
diff_rec_cor <- diff_cor[receptor]
diff_rec_p <- diff_p[receptor]

my_rec_cor <- my_cor[receptor]
my_rec_p <- my_p[receptor]
diff_lig_cor <- diff_cor[ligand]
diff_lig_p <- diff_p[ligand]

ligand_source <- c(rep("Myeloid", length(my_lig_cor)), rep("Diff", length(my_rec_cor)))
receptor_source <- c(rep("Diff", length(my_lig_cor)), rep("Myeloid", length(my_rec_cor)))

res1 <- data.frame(my_lig_cor, my_lig_p, diff_rec_cor, diff_rec_p)
res2 <- data.frame(my_rec_cor, my_rec_p, diff_lig_cor, diff_lig_p)
colnames(res1) <- colnames(res2) <- c("ligand_cor", "ligand_p", "receptor_cor", "receptor_p")

res <- rbind(res1, res2)
res <- data.frame(ligand, receptor, res, ligand_source, receptor_source)
res$ligand_q <- p.adjust(res$ligand_p, "BH")
res$receptor_q <- p.adjust(res$receptor_p, "BH")

res %>% filter(ligand_cor > 0.4, receptor_cor > 0.4, ligand_q < 0.05, receptor_q < 0.05)

# Plot the longitudinal result
plot_res <- res %>%
			filter(!is.na(receptor_cor), !is.na(ligand_cor)) %>%
			mutate(group = paste(ligand_source, receptor_source, sep = "_")) %>%
			mutate(group = recode(group, "Diff_Myeloid" = "Ligand: Diff.-like, Receptor: Myeloid", "Myeloid_Diff" = "Ligand: Myeloid, Receptor: Diff.-like")) %>%
			mutate(sig = as.character(ligand_q < 0.05 & receptor_q < 0.05 & ligand_cor > 0 & receptor_cor > 0))
			
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/longitudinal_mes_rec_lig_screen.pdf",width=1.6,height=2.95)
ggplot(plot_res, aes(ligand_cor, receptor_cor, colour = sig)) + 
geom_point()  +
geom_hline(yintercept=0, colour = "gray") +
geom_vline(xintercept=0, colour = "gray") +
scale_colour_manual(values = c("TRUE" = "#CD4F39", "FALSE" = "black")) +
labs(x="Myeloid Pearson coefficient", y = "Diff.-like Pearson coefficient") +
facet_wrap(group~., nrow=2) +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") 
dev.off()


# Single cell validation


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
annot = tsne.data %>%
	rownames_to_column('cell') %>%
	mutate(cell_type = recode(dbCluster, `1` = "Diff.-like",  `2` = "Myeloid", `3` = "Stem-like",
							  `4` = "Oligodendrocyte", `5` = "Prolif. stem-like", `6` = "granulocyte", `7` = "endothelial",
							  `8` = "T cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) %>%
	column_to_rownames('cell')
	

# Get sample names
sample_id <- sapply(strsplit(rownames(annot), "-"), "[[", 3)
sample_id <- recode(sample_id, "0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006",
					"5" = "SM011", "6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018")
annot[,"sample_id"] <- sample_id
names(sample_id) <- rownames(annot)

myeloid_cells <- annot %>% filter(cell_type == "Myeloid")
tumor_cells <- annot %>% filter(cell_type == "Diff.-like")

# Convert ensembl ID to gene symbol
featuredata <- featuredata[rownames(log2cpm),]

gene <- featuredata[,"Associated.Gene.Name"]
names(gene) <- rownames(featuredata)

rownames(log2cpm) <- gene

# Read in bulk subtype results
con <- DBI::dbConnect(odbc::odbc(), "scgp")

bulk_res <- dbReadTable(con, Id(schema = "analysis", table="transcriptional_subtype"))
bulk_mes <- bulk_res %>% 
			filter(signature_name == "Mesenchymal") %>%
			mutate(sample_id = sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[2],x[3],sep=""))) %>%
			filter(sample_id %in% c("SM006","SM012","SM017","SM018","SM011"))

# Test receptor-ligand associations

# Filter for mesenchymal myeloid signature only (old way)
# sub_lig <- lig_plot %>% filter(rec_eff > log2(1.1), rec_qval < 0.1 &
# 			lig_cell_state %in% c("Diff.-like", "Stem-like", "Prolif. stem-like")) # Tumor only
# sub_rec <- rec_plot %>% filter(lig_eff > log2(1.1), lig_qval < 0.1 &
# 			rec_cell_state %in% c("Diff.-like", "Stem-like", "Prolif. stem-like")) # Tumor only

hits <- res %>% filter(ligand_cor > 0, receptor_cor > 0, ligand_q < 0.05, receptor_q < 0.05)
hits <- hits[which(hits$ligand %in% rownames(log2cpm) & hits$receptor %in% rownames(log2cpm)),]
ligand_cor <- ligand_p <- receptor_cor <- receptor_p <- rep(NA, nrow(hits))
for(i in 1:nrow(hits))
{
	cat("\r", i)
	
	myeloid_gene <- ifelse(hits[i,"ligand_source"]=="Myeloid", as.character(hits[i,"ligand"]), as.character(hits[i,"receptor"]))
	tumor_gene <- ifelse(hits[i,"ligand_source"]=="Diff", as.character(hits[i,"ligand"]), as.character(hits[i,"receptor"]))
	
	myeloid_cpm <- t(log2cpm[myeloid_gene,rownames(myeloid_cells)]) %>% data.frame()
	myeloid_cpm[,"sample_id"] <- myeloid_cells[,"sample_id"]
	colnames(myeloid_cpm) <- c("cpm","sample_id")
	
	myeloid_cor <- myeloid_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor(mean_cpm, enrichment_score, method="p")) %>%
	as.numeric()
	
	myeloid_p <- myeloid_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor.test(mean_cpm, enrichment_score, method="p")$p.value) %>%
	as.numeric()
	
	#tumor_cells <- annot %>% filter(cell_type == rec_lig[i,"tumor_class"])

	tumor_cpm <- t(log2cpm[tumor_gene,rownames(tumor_cells)]) %>% data.frame()
	tumor_cpm[,"sample_id"] <- tumor_cells[,"sample_id"]
	colnames(tumor_cpm) <- c("cpm","sample_id")
	
	
	tumor_cor <- tumor_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor(mean_cpm, enrichment_score, method="p")) %>%
	as.numeric()
	
	tumor_p <- tumor_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor.test(mean_cpm, enrichment_score, method="p")$p.value) %>%
	as.numeric()
	
	if(hits[i,"ligand_source"]=="Myeloid")
	{
		ligand_cor[i] <- myeloid_cor
		ligand_p[i] <- myeloid_p
		receptor_cor[i] <- tumor_cor
		receptor_p[i] <- tumor_p
	}
	if(hits[i,"receptor_source"]=="Myeloid")
	{
		ligand_cor[i] <- tumor_cor
		ligand_p[i] <- tumor_p
		receptor_cor[i] <- myeloid_cor
		receptor_p[i] <- myeloid_p
	}
}
v_res <- data.frame(hits[,c(1,2,7,8)], ligand_cor, ligand_p, receptor_cor, receptor_p)
test <- apply(data.frame(ligand_cor, receptor_cor), 1, mean)

v_res <- v_res %>% filter(!is.na(myeloid_cor), !is.na(tumor_cor))
v_res[order(test,decreasing=TRUE),]

plot_res2 <- v_res %>%
			filter(!is.na(receptor_cor), !is.na(ligand_cor)) %>%
			mutate(group = paste(ligand_source, receptor_source, sep = "_")) %>%
			mutate(group = recode(group, "Diff_Myeloid" = "Ligand: Diff.-like, Receptor: Myeloid", "Myeloid_Diff" = "Ligand: Myeloid, Receptor: Diff.-like")) %>%
			mutate(highlight = as.character(ligand_cor > 0.7 & receptor_cor > 0.7))
			
#Get numbers in each ligand/receptor expression group
plot_res2 %>% filter(receptor_cor > 0 & ligand_cor > 0) %>% .$group %>% table()

# Create supplementary table of candidates that passed both criteria
pass <- plot_res2 %>%
		filter(receptor_cor > 0 & ligand_cor > 0)  %>%
		dplyr::select(ligand, receptor, ligand_source, receptor_source, 
		ligand_cor, receptor_cor)
pass_join <- pass %>%
		inner_join(plot_res, by = c("ligand", "receptor", "ligand_source", "receptor_source"), suffix=c("_sc","_long")) %>%
		dplyr::select(-c(group,sig)) %>%
		dplyr::select(ligand, receptor, ligand_source, receptor_source, 
			   ligand_cor_long, ligand_p, ligand_q, receptor_cor_long, receptor_p, receptor_q, 
			   ligand_cor_sc, receptor_cor_sc)
pj_ord <- apply(data.frame(pass_join$ligand_cor_sc, pass_join$receptor_cor_sc), 1, mean)
pass_join <- pass_join[order(pass_join$ligand_source, pj_ord,decreasing=TRUE),]

write.table(pass_join, "data/res/CIBERSORTx/analysis/ligand_receptor_screen.txt",sep="\t",quote=FALSE,row.names=FALSE)
