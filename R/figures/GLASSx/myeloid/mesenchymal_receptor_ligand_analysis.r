
library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)

#######################################################
rm(list=ls())
set.seed(11)
##################################################
# Step 1: Subtype each TCGA sample
##################################################

gct_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct"
clin_path <- "/projects/verhaak-lab/USERS/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt"

info <- read.delim(clin_path)
mycase <- info %>% 
		filter(IDH.codel.subtype == "IDHwt") %>%
		dplyr::select(Case) %>%
		.$Case %>%
		as.character()

# Run Qianghu's transcriptional classifier (only needed once)
#runSsGSEAwithPermutation(gct_path,100)

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
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case"))

##################################################
# Step 2: UMAP of each CIBERSORTx GEP
##################################################

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated|t_cell|oligodendrocyte",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_TCGA_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

diff_list <- embedding <- gep_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	colnames(geps) <- paste(mytag[i], colnames(geps), sep="__")
	geps <- log10(geps+1)
	gep_list[[i]] <- t(geps)

	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]


	embedding[[i]] <- umap(t(geps))$layout


	g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Mesenchymal"),"aliquot_barcode"]]
	g2 <- geps[,sig_sub[which(!grepl("Mesenchymal",sig_sub[,"signature_name"])),"aliquot_barcode"]]

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(g1[j,])
		group2 <- as.numeric(g2[j,])
	
		p.val[j] <- wilcox.test(group1,group2)$p.value
		eff[j] <- log2(median(group2)/median(group1))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	idhwt_res <- data.frame(p.val, q.val, eff)
	idhwt_res <- idhwt_res[order(p.val),]

	idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
	idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1
		
	mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
	
	diff_list[[i]] <- idhwt_res
}

full_embedding <- do.call(rbind, embedding)

plot_res <- full_embedding
aliquot_barcode <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[2])
cell_state <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[1])
plot_res <- data.frame(aliquot_barcode, plot_res, cell_state)
colnames(plot_res) <- c("aliquot_barcode","UMAP1","UMAP2", "cell_state")
plot_res <- plot_res %>%
	inner_join(sig_sub, by = "aliquot_barcode")
rownames(plot_res) <- paste(plot_res[,"cell_state"], plot_res[,"aliquot_barcode"], sep = "__")

# Save the myeloid signature to apply to GLASS
mes_myel_gene_sig <- diff_list[[2]]
write.table(mes_myel_gene_sig, "data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result.txt",sep="\t",quote=FALSE,row.names=TRUE)


# Create a boxplot showing the signature's difference across transcriptional subtypes
mes_geps <- gep_list[[2]]
mes_genes <- rownames(mes_myel_gene_sig[which(mes_myel_gene_sig[,"sig"] & mes_myel_gene_sig[,"eff"] < 0),])
mes_score <- apply(mes_geps[,mes_genes], 1, mean)
names(mes_score) <- sapply(strsplit(names(mes_score), "__"), function(x)x[2])
mes_score <- mes_score[sig_sub[,"aliquot_barcode"]]
mes_score <- data.frame(sig_sub, mes_score)
mes_score <- mes_score %>% filter(IDH.codel.subtype == "IDHwt")

colors <- c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#3F264B")
names(colors) <- levels(factor(mes_score$signature_name))

mes_score <- mes_score %>%
mutate(signature_name = fct_relevel(signature_name, "Proneural","Classical","Mixed","Proneural,Mesenchymal","Classical,Mesenchymal","Mesenchymal"))
reord_colors <- colors[levels(mes_score$signature_name)]

wilcox.test(mes_score %>% filter(signature_name == "Mesenchymal") %>% .$mes_score, mes_score %>% filter(signature_name != "Mesenchymal") %>% .$mes_score, )

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_mesenchymal_sig_by_subtype.pdf",width=3.5,height=2.5)
ggplot(mes_score, aes(signature_name, mes_score, fill = signature_name)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
labs(title = "Mesenchymal signature") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=reord_colors)
dev.off()

# UMAP plot
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_tcga_facet_idhwt_by_batch_umap_test.pdf",width=7,height=3.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(UMAP1, UMAP2, color = signature_name)) + 
geom_point()  +
theme_bw() +
facet_wrap(~ cell_state, scales = 'free', nrow = 2) +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()

##################################################
# Step 3: Ligand-receptor pairs in CIBERSORTx
##################################################


full_geps <- do.call(rbind, gep_list)

lig_rec <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/ramilowski_2015/receptor_ligand_pairs_fantom5_ramilowski_ncomms_2015.txt",stringsAsFactor=FALSE)
lig_rec <- lig_rec %>%
filter(lig_rec[,"Pair.Evidence"] == "literature supported")
lig_rec <- lig_rec[,c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]

lig_rec <- lig_rec[which(lig_rec[,"Ligand.ApprovedSymbol"] %in% colnames(full_geps) & 
						 lig_rec[,"Receptor.ApprovedSymbol"] %in% colnames(full_geps)),]

#mypair <- full_geps[rownames(plot_res),c("OSM","OSMR")]

dl_pval <- dl_eff <- dr_pval <- dr_eff <- ml_pval <- ml_eff <- mr_pval <- mr_eff <- rep(NA, nrow(lig_rec))
for(i in 1:nrow(lig_rec))
{
	cat("\r", i)
	mypair <- as.character(lig_rec[i,])
	names(mypair) <- c("ligand","receptor")
	pair_expr <- full_geps[rownames(plot_res),mypair]
	diff1 <- plot_res %>%
			 filter(signature_name == "Mesenchymal" & cell_state == "differentiated_tumor") %>%
			 rownames()
	diff2 <- plot_res %>%
			 filter(signature_name != "Mesenchymal" & cell_state == "differentiated_tumor") %>%
			 rownames()
	myel1 <- plot_res %>%
			 filter(signature_name == "Mesenchymal" & cell_state == "myeloid") %>%
			 rownames()
	myel2 <- plot_res %>%
			 filter(signature_name != "Mesenchymal" & cell_state == "myeloid") %>%
			 rownames()
	
	# Get p-value/fc for whether ligand is differentially expressed in mes tumor cells relative to non-mes tumor cells
	if(sum(is.na(full_geps[c(diff1,diff2),mypair["ligand"]])) < length(c(diff1,diff2)))
	{
		dl_pval[i] <- wilcox.test(full_geps[diff1,mypair["ligand"]], full_geps[diff2, mypair["ligand"]])$p.value
		dl_eff[i] <- log2(mean(full_geps[diff1,mypair["ligand"]])/mean(full_geps[diff2, mypair["ligand"]]))
		#dl_eff[i] <- median(full_geps[diff1,mypair["ligand"]]) - median(full_geps[diff2, mypair["ligand"]])
	}
	# Get p-value/fc for whether receptor is differentially expressed in mes tumor cells relative to non-mes tumor cells
	if(sum(is.na(full_geps[c(diff1,diff2),mypair["receptor"]])) < length(c(diff1,diff2)))
	{
		dr_pval[i] <- wilcox.test(full_geps[diff1,mypair["receptor"]], full_geps[diff2, mypair["receptor"]])$p.value
		dr_eff[i] <- log2(mean(full_geps[diff1,mypair["receptor"]])/mean(full_geps[diff2, mypair["receptor"]]))		
		#dr_eff[i] <- median(full_geps[diff1,mypair["receptor"]]) - median(full_geps[diff2, mypair["receptor"]])
	}
	# Get p-value/fc for whether ligand is differentially expressed in mes myeloid cells relative to non-mes myeloid cells
	if(sum(is.na(full_geps[c(myel1,myel2),mypair["ligand"]])) < length(c(myel1,diff2)))
	{
		ml_pval[i] <- wilcox.test(full_geps[myel1,mypair["ligand"]], full_geps[myel2, mypair["ligand"]])$p.value
		ml_eff[i] <- log2(mean(full_geps[myel1,mypair["ligand"]])/mean(full_geps[myel2, mypair["ligand"]]))
		#ml_eff[i] <- median(full_geps[myel1,mypair["ligand"]]) - median(full_geps[myel2, mypair["ligand"]])
	}
	# Get p-value/fc for whether receptor is differentially expressed in mes myeloid cells relative to non-mes myeloid cells
	if(sum(is.na(full_geps[c(myel1,myel2),mypair["receptor"]])) < length(c(myel1,myel2)))
	{
		mr_pval[i] <- wilcox.test(full_geps[myel1,mypair["receptor"]], full_geps[myel2, mypair["receptor"]])$p.value
		mr_eff[i] <- log2(mean(full_geps[myel1,mypair["receptor"]])/mean(full_geps[myel2, mypair["receptor"]]))
		#mr_eff[i] <- median(full_geps[myel1,mypair["receptor"]]) - median(full_geps[myel2, mypair["receptor"]])
	}
}

lig_rec_pair <- paste(lig_rec[,"Ligand.ApprovedSymbol"], lig_rec[,"Receptor.ApprovedSymbol"],sep="_")
dl_qval <- p.adjust(dl_pval, "BH")
dr_qval <- p.adjust(dr_pval, "BH")
ml_qval <- p.adjust(ml_pval, "BH")
mr_qval <- p.adjust(mr_pval, "BH")

#Set a threshold that both pairs be over the 75th percentile:
dlmrres <- data.frame(lig_rec_pair, dl_qval, mr_qval, dl_eff, mr_eff)
dlmrres <- dlmrres %>% filter(!is.na(dl_qval) & !is.na(mr_qval))

dlmr_thr <- which(dlmrres[,"dl_eff"] > quantile(dlmrres[,"dl_eff"],.75) & dlmrres[,"mr_eff"] > quantile( dlmrres[,"mr_eff"], .75))
dlmrres[dlmr_thr,]


drmlres <- data.frame(lig_rec_pair, dr_qval, ml_qval, dr_eff, ml_eff)
drmlres <- drmlres %>% filter(!is.na(dr_qval) & !is.na(ml_qval))

drml_thr <- which(drmlres[,"dr_eff"] > quantile(drmlres[,"dr_eff"],.75) & drmlres[,"ml_eff"] > quantile( drmlres[,"ml_eff"], .75))
drmlres[drml_thr,]

# Identify eligible pairs
# res <- res[which((dl_eff > 0 & mr_eff > 0) | (dr_eff > 0 & ml_eff > 0)),]
# sig_res <- res %>%
# filter(dl_qval < 0.01 & mr_qval < 0.01 | dr_qval < 0.01 & ml_qval < 0.01)
# 
# Tumor signalling to myeloid:
# tmres <- res %>%
# 		 filter(dl_qval < 0.01 & mr_qval < 0.01)
# mtres <- res %>%
# 		 filter(dr_qval < 0.01 & ml_qval < 0.01)
# 
# 
# logres <-  -log10(res[,2:5])
# logres <- data.frame(res[,1], logres)
# colnames(logres)[1] <- "lig_rec_pair"


# Plot
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_myel_lig_diff_rec.pdf",width=2.5,height=2)
ggplot(drmlres %>% filter(dr_eff > 0 & ml_eff > 0), aes(dr_eff, ml_eff)) + 
geom_point()  +
theme_bw() +
geom_vline(xintercept=quantile(drmlres[,"dr_eff"],.75)) +
geom_hline(yintercept=quantile(drmlres[,"ml_eff"],.75)) 
dev.off()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_myel_rec_diff_lig.pdf",width=2.5,height=2)
ggplot(dlmrres %>% filter(dl_eff > 0 & mr_eff > 0), aes(dl_eff, mr_eff)) + 
geom_point()  +
theme_bw() +
geom_vline(xintercept=quantile(dlmrres[,"dl_eff"],.75)) +
geom_hline(yintercept=quantile(dlmrres[,"mr_eff"],.75)) 
dev.off()

# 
# 
# 
# Plot
# pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_myel_lig_diff_rec.pdf",width=2.5,height=1.5)
# ggplot(logres, aes(dr_qval, ml_qval)) + 
# geom_point()  +
# theme_bw()
# dev.off()
# 
# logres %>% filter(dr_qval > 15 & ml_qval > 15)
# 
# pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_myel_rec_diff_lig.pdf",width=2.5,height=1.5)
# ggplot(logres, aes(dl_qval, mr_qval)) + 
# geom_point()  +
# theme_bw()
# dev.off()
# 
# logres %>% filter(dl_qval > 9.5 & mr_qval > 10)