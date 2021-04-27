library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)
library(topGO)

#######################################################
rm(list=ls())

# CIBERSORTx cell-specific profiles
#------------------------------
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated|t_cell|oligodendrocyte",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_TCGA_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

# TCGA subtyping information
#------------------------------

# Read clinical data
clin_path <- "/projects/verhaak-lab/USERS/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt"

info <- read.delim(clin_path)
mycase <- info %>% 
		filter(IDH.codel.subtype == "IDHwt") %>%
		dplyr::select(Case) %>%
		.$Case %>%
		as.character()

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/p_result_TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct.txt")
rownames(subtype_ssgsea) <- gsub("\\.","-",rownames(subtype_ssgsea))

aliquot_barcode <- rep(rownames(subtype_ssgsea),3)
signature_name <- c(rep("Proneural",nrow(subtype_ssgsea)),rep("Classical",nrow(subtype_ssgsea)),rep("Mesenchymal",nrow(subtype_ssgsea)))
enrichment_score <- c(subtype_ssgsea[,"Proneural"],subtype_ssgsea[,"Classical"],subtype_ssgsea[,"Mesenchymal"])
p_value <- c(subtype_ssgsea[,"Proneural_pval"],subtype_ssgsea[,"Classical_pval"],subtype_ssgsea[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(aliquot_barcode,signature_name,enrichment_score,p_value)
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"aliquot_barcode"]),]

transcriptional_subtype[,"signif"] <- transcriptional_subtype[,"p_value"] < 0.05

sig_sub <- transcriptional_subtype %>%
		   group_by(aliquot_barcode, signif) %>%
		   summarise(signature_name = paste(signature_name, collapse=",")) %>%
		   data.frame()
# Pull out significant subtypes first
sig_sub1 <- sig_sub %>%
		    filter(signif)
signif_ali <- sig_sub1[,"aliquot_barcode"]

sig_sub2 <- sig_sub %>%
			filter(!(aliquot_barcode %in% signif_ali) & !signif)
sig_sub <- rbind(sig_sub1, sig_sub2)
  
sig_sub[which(sig_sub[,"signature_name"] == "Proneural,Classical,Mesenchymal"),"signature_name"] <- "Mixed"

sig_sub[,"case_barcode"] <- sapply(strsplit(as.character(sig_sub[,"aliquot_barcode"]),"-"),function(x)paste(x[1:3],collapse="-"))
sig_sub <- sig_sub %>%
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case")) %>%
		   inner_join((transcriptional_subtype %>% filter(signature_name=="Mesenchymal"))[,c("aliquot_barcode","enrichment_score")], by = "aliquot_barcode")

# Read in TCGA mesenchymal myeloid sig
#------------------------------
tcga_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v3.txt")
mes_genes <- rownames(tcga_dat %>% filter(sig, eff > log2(1)))

##################################################
# Identify ligand-receptor pairs
##################################################

mygeps <- read.delim(myinf1[2],row.names=1)

lig_rec <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/ramilowski_2015/receptor_ligand_pairs_fantom5_ramilowski_ncomms_2015.txt",stringsAsFactor=FALSE)
lig_rec <- lig_rec %>%
filter(lig_rec[,"Pair.Evidence"] == "literature supported")
lig_rec <- lig_rec[,c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]

lig_rec <- lig_rec[which(lig_rec[,"Ligand.ApprovedSymbol"] %in% rownames(mygeps) & 
						 lig_rec[,"Receptor.ApprovedSymbol"] %in% rownames(mygeps)),]

# Select the mesenchymal myeloid-specific ligands and receptors
cand_lig_pairs <- lig_rec[which(lig_rec[,"Ligand.ApprovedSymbol"] %in% mes_genes),]
cand_rec_pairs <- lig_rec[which(lig_rec[,"Receptor.ApprovedSymbol"] %in% mes_genes),]

# Identify where the complementary receptors are expressed
log2fc <- function(expr, mes_status){
	log2(median(expr[which(mes_status == "Mes")])/median(expr[which(mes_status == "Non-mes")]))}

nonmyelinf <- myinf1[-2]
rec_list <- lig_list <- list()
for(i in 1:length(nonmyelinf))
{
	geps <- read.delim(nonmyelinf[i],row.names=1)
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	geps <- log10(geps+1)
	
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]

	# IDHwt only
 	barcodes <- sig_sub %>% filter(IDH.codel.subtype=="IDHwt") %>% .$aliquot_barcode
 	geps <- geps[,barcodes]
	
	# Compare receptor expression between Mes and non-mes IDHwt samples
	recs <- intersect(cand_lig_pairs[,"Receptor.ApprovedSymbol"], rownames(geps))
	rec_geps <- geps[recs,]
	
	rec_res <- rec_geps %>% 
	rownames_to_column("receptor") %>%
	pivot_longer(!receptor, names_to = "aliquot_barcode", values_to = "expr") %>%
	inner_join(sig_sub, by = "aliquot_barcode") %>%
	filter(IDH.codel.subtype == "IDHwt") %>%
	filter(signature_name == "Mesenchymal" | !grepl("Mesenchymal", signature_name)) %>%
	mutate(mes_status = ifelse(grepl("Mesenchymal",signature_name), "Mes", "Non-mes")) %>%
	mutate(mes_status = as_factor(mes_status)) %>%
	mutate(mes_status = fct_relevel(mes_status, "Non-mes", "Mes")) %>%
	group_by(receptor) %>%
	summarise(p.val = wilcox.test(expr ~ mes_status)$p.value, eff = log2fc(expr, mes_status)) %>%
	as.data.frame()

	rec_res[,"q.val"] <- p.adjust(rec_res$p.val,"BH")

	# Compare ligand expression between Mes and non-mes IDHwt samples	
	ligs <- intersect(cand_rec_pairs[,"Ligand.ApprovedSymbol"], rownames(geps))
	lig_geps <- geps[ligs,]
	
	lig_res <- lig_geps %>% 
	rownames_to_column("ligand") %>%
	pivot_longer(!ligand, names_to = "aliquot_barcode", values_to = "expr") %>%
	inner_join(sig_sub, by = "aliquot_barcode") %>%
	filter(IDH.codel.subtype == "IDHwt") %>%
	filter(signature_name == "Mesenchymal" | !grepl("Mesenchymal", signature_name)) %>%
	mutate(mes_status = ifelse(grepl("Mesenchymal",signature_name), "Mes", "Non-mes")) %>%
	mutate(mes_status = as_factor(mes_status)) %>%
	mutate(mes_status = fct_relevel(mes_status, "Non-mes", "Mes")) %>%
	group_by(ligand) %>%
	summarise(p.val = wilcox.test(expr ~ mes_status)$p.value, eff = log2fc(expr, mes_status)) %>%
	as.data.frame()
	
	lig_res[,"q.val"] <- p.adjust(lig_res$p.val,"BH")
	
	mytag <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga//CIBERSORTxHiRes_TCGA_","",nonmyelinf[i])
	mytag <- gsub("_Window48.txt","",mytag)
	rec_res[,"cell_state"] <- mytag
	lig_res[,"cell_state"] <- mytag
	
	rec_list[[i]] <- rec_res
	lig_list[[i]] <- lig_res
}

rec_plot <- do.call(rbind,rec_list)
lig_plot <- do.call(rbind,lig_list)

# Fix the myeloid signature table
tcga_dat <- tcga_dat %>%
rownames_to_column("gene_symbol")

# Integrate with the myeloid ligand results
rec_plot <- rec_plot %>%
inner_join(cand_lig_pairs, by = c("receptor" = "Receptor.ApprovedSymbol")) %>%
inner_join(tcga_dat, by = c("Ligand.ApprovedSymbol" = "gene_symbol")) %>%
mutate(pair = paste(receptor, Ligand.ApprovedSymbol, sep = "_")) %>%
dplyr::select(pair = pair, receptor = receptor, ligand = Ligand.ApprovedSymbol, rec_pval = p.val.x, rec_qval = q.val.x, rec_eff = eff.x, lig_pval = p.val.y, lig_qval = q.val.y, lig_eff = eff.y,  rec_cell_state = cell_state) %>%
add_column(lig_cell_state = "Myeloid") %>%
mutate(rec_cell_state = recode(rec_cell_state, 
		"differentiated_tumor" = "Diff.-like", "oligodendrocyte" = "Oligodendrocyte",
		"prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
filter(rec_qval < 0.05 & lig_qval < 0.05) %>%
arrange(desc(rec_cell_state))


#mutate(sig = ifelse(rec_pval < 0.05 & lig_pval < 0.05, rec_cell_state, "nonsig"))
		
rec_plot[which(rec_plot$rec_eff > log2(1.1) & rec_plot$lig_eff > log2(1.1)),]

lig_plot <- lig_plot %>%
inner_join(cand_rec_pairs, by = c("ligand" = "Ligand.ApprovedSymbol")) %>%
inner_join(tcga_dat, by = c("Receptor.ApprovedSymbol" = "gene_symbol")) %>%
mutate(pair = paste(Receptor.ApprovedSymbol, ligand, sep = "_")) %>%
dplyr::select(pair = pair, receptor = Receptor.ApprovedSymbol, ligand = ligand, rec_pval = p.val.y, rec_qval = q.val.y, rec_eff = eff.y, lig_pval = p.val.x, lig_qval = q.val.y, lig_eff = eff.x, lig_cell_state = cell_state) %>%
add_column(rec_cell_state = "Myeloid")  %>%
mutate(lig_cell_state = recode(lig_cell_state, 
		"differentiated_tumor" = "Diff.-like", "oligodendrocyte" = "Oligodendrocyte",
		"prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
#mutate(lig_cell_state = fct_relevel(lig_cell_state, "T cell", "Oligodendrocyte", "Prolif. stem-like", "Stem-like", "Diff.-like")) %>%
filter(rec_qval < 0.05 & lig_qval < 0.05) %>%
arrange(desc(lig_cell_state))

#mutate(sig = ifelse(rec_pval < 0.05 & lig_pval < 0.05, lig_cell_state, "nonsig"))


lig_plot[which(lig_plot$rec_eff > log2(1.1) & lig_plot$lig_eff > log2(1.1)),]


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_myel_lig_environ_rec.pdf",width=1.75,height=1.4)
ggplot(rec_plot %>% filter(rec_eff >= 0 & lig_eff >= 0), aes(lig_eff, rec_eff, colour = rec_cell_state)) + 
geom_point(alpha = 0.7)  +
scale_colour_manual(values=c("T cell" = "#6baed6", "Oligodendrocyte" = "#2ca25f", 
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
labs(x="Myeloid ligand log2(FC)", y = "Non-myeloid receptor log2(FC)") +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
geom_vline(xintercept=log2(1.1),linetype=2,size = 0.5) +
geom_hline(yintercept=log2(1.1),linetype=2,size = 0.5) + 
coord_cartesian(ylim = c(0, 1), xlim=c(0, 1))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_myel_rec_environ_lig.pdf",width=1.75,height=1.4)
ggplot(lig_plot %>% filter(rec_eff >= 0 & lig_eff >= 0), aes(rec_eff, lig_eff, colour = lig_cell_state)) + 
geom_point(alpha = 0.7)  +
scale_colour_manual(values=c("T cell" = "#6baed6", "Oligodendrocyte" = "#2ca25f", 
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
labs(x="Myeloid receptor log2(FC)", y = "Non-myeloid ligand log2(FC)") +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
geom_vline(xintercept=log2(1.1),linetype=2, size = 0.5) +
geom_hline(yintercept=log2(1.1),linetype=2, size = 0.5) +
coord_cartesian(ylim = c(0, 1), xlim=c(0, 1))
dev.off()

##################################################
# Read in SCGP 
##################################################

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
			mutate(sample_id = sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[2],x[3],sep=""))) #%>%
			#filter(sample_id %in% c("SM006","SM012","SM017","SM018","SM011"))

# Test receptor-ligand associations

# Filter for mesenchymal myeloid signature only (old way)
# sub_lig <- lig_plot %>% filter(rec_eff > log2(1.1), rec_qval < 0.1 &
# 			lig_cell_state %in% c("Diff.-like", "Stem-like", "Prolif. stem-like")) # Tumor only
# sub_rec <- rec_plot %>% filter(lig_eff > log2(1.1), lig_qval < 0.1 &
# 			rec_cell_state %in% c("Diff.-like", "Stem-like", "Prolif. stem-like")) # Tumor only

# Both components significant
sub_lig <- lig_plot %>% filter(rec_eff > log2(1.1), rec_qval < 0.1, lig_eff > log2(1.1), lig_qval < 0.1,
			lig_cell_state %in% c("Diff.-like", "Stem-like", "Prolif. stem-like")) # Tumor only
sub_rec <- rec_plot %>% filter(lig_eff > log2(1.1), lig_qval < 0.1, rec_eff > log2(1.1), rec_qval < 0.1,
			rec_cell_state %in% c("Diff.-like", "Stem-like", "Prolif. stem-like")) # Tumor only

myeloid_genes <- c(sub_lig$receptor, sub_rec$ligand) #c("CD44", "GALR2","ITGB1", "LTBR", "OSMR", "CCR10", "IL7", "TSLP")
tumor_genes <- c(sub_lig$ligand, sub_rec$receptor) #c("COL14A1", "GAL", "ICAM4", "TNFSF14", "OSM", "CCL2", "IL7R", "IL7R")
tumor_class <- c(sub_lig$lig_cell_state, sub_rec$rec_cell_state)#c("differentiated_tumor", "differentiated_tumor", "differentiated_tumor", "differentiated_tumor", "differentiated_tumor", "stemcell_tumor", "differentiated_tumor", "differentiated_tumor")
rec_lig <- data.frame(myeloid_genes, tumor_genes, tumor_class, stringsAsFactors=FALSE)

myeloid_cells <- annot %>% filter(cell_type == "Myeloid")
myeloid_cor <- myeloid_p <- tumor_cor <- tumor_p <- rep(NA, nrow(rec_lig))
for(i in 1:nrow(rec_lig))
{
	cat("\r", i)
	myeloid_cpm <- t(log2cpm[rec_lig[i,"myeloid_genes"],rownames(myeloid_cells)]) %>% data.frame()
	myeloid_cpm[,"sample_id"] <- myeloid_cells[,"sample_id"]
	colnames(myeloid_cpm) <- c("cpm","sample_id")
	
	myeloid_cor[i] <- myeloid_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor(mean_cpm, enrichment_score, method="p")) %>%
	as.numeric()
	
	myeloid_p[i] <- myeloid_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor.test(mean_cpm, enrichment_score, method="p")$p.value) %>%
	as.numeric()
	
	tumor_cells <- annot %>% filter(cell_type == rec_lig[i,"tumor_class"])

	tumor_cpm <- t(log2cpm[rec_lig[i,"tumor_genes"],rownames(tumor_cells)]) %>% data.frame()
	tumor_cpm[,"sample_id"] <- tumor_cells[,"sample_id"]
	colnames(tumor_cpm) <- c("cpm","sample_id")
	
	
	tumor_cor[i] <- tumor_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor(mean_cpm, enrichment_score, method="p")) %>%
	as.numeric()
	
	tumor_p[i] <- tumor_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(cor.test(mean_cpm, enrichment_score, method="p")$p.value) %>%
	as.numeric()
}
res <- data.frame(rec_lig, myeloid_cor, myeloid_p, tumor_cor, tumor_p)
test <- apply(data.frame(myeloid_cor, tumor_cor), 1, mean)
res[order(test,decreasing=TRUE),]

res <- res %>% filter(!is.na(myeloid_cor), !is.na(tumor_cor))

# Plot correlation coefficients

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_rec_lig_cor.pdf",width=3,height=3)
ggplot(res, aes(myeloid_cor, tumor_cor, colour = tumor_class)) + 
geom_point()  +
scale_colour_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
labs(x="Myeloid receptor log2(FC)", y = "Non-myeloid ligand log2(FC)") +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") 
dev.off()



