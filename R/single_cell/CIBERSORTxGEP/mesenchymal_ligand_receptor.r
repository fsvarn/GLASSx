
library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)
library(gridExtra)

#######################################################
rm(list=ls())
set.seed(11)
##################################################
# Step 1: Subtype each TCGA sample
##################################################

gct_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct"
clin_path <- "/projects/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt"

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

sig_sub <- transcriptional_subtype %>% 
		   filter(p_value < 0.05) 

case_barcode <- sapply(strsplit(as.character(sig_sub[,"aliquot_barcode"]), "-"),function(x)paste(x[1:3],collapse="-"))

sig_sub <- sig_sub[which(case_barcode %in% mycase),] %>%
		   group_by(aliquot_barcode) %>%
		   summarise(signature_name = paste(signature_name, collapse=",")) %>%
		   data.frame()
sig_sub <- sig_sub[-grep(",",sig_sub[,"signature_name"]),]

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
	g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]!="Mesenchymal"),"aliquot_barcode"]]

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

# UMAP plot
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_tcga_facet_idhwt_by_batch_umap.pdf",width=7,height=3.5)
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

mypair <- full_geps[rownames(plot_res),c("OSM","OSMR")]
lig_rec <- data.frame(plot_res, mypair)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_tcga_facet_osm.pdf",width=3.2,height=1.5)
p1 <- ggplot(lig_rec %>% 
	   filter(!grepl("-02A-|-02B-",lig_rec[,"aliquot_barcode"])) %>% 
	   filter(cell_state == "myeloid"), 
	   aes(UMAP1, UMAP2, color = signature_name)) + 
geom_point()  +
theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = "none") +
#facet_wrap(~ cell_state, scales = 'free', nrow = 2) +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))

p2 <- ggplot(lig_rec %>% 
	   filter(!grepl("-02A-|-02B-",lig_rec[,"aliquot_barcode"])) %>% 
	   filter(cell_state == "myeloid"), 
	   aes(UMAP1, UMAP2, color = OSM)) + 
geom_point()  +
theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = "none") +
#facet_wrap(~ cell_state, scales = 'free', nrow = 2) +
scale_color_gradient(low = "#D3D3D3", high = "#a50000", na.value = NA)

grid.arrange(p1,p2,nrow=1)
dev.off()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_tcga_facet_osmr.pdf",width=3.2,height=1.5)
p1 <- ggplot(lig_rec %>% 
	   filter(!grepl("-02A-|-02B-",lig_rec[,"aliquot_barcode"])) %>% 
	   filter(cell_state == "differentiated_tumor"), 
	   aes(UMAP1, UMAP2, color = signature_name)) + 
geom_point()  +
theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = "none") +
#facet_wrap(~ cell_state, scales = 'free', nrow = 2) +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))

p2 <- ggplot(lig_rec %>% 
	   filter(!grepl("-02A-|-02B-",lig_rec[,"aliquot_barcode"])) %>% 
	   filter(cell_state == "differentiated_tumor"), 
	   aes(UMAP1, UMAP2, color = OSMR)) + 
geom_point()  +
theme_bw() +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	  legend.position = "none") +
#facet_wrap(~ cell_state, scales = 'free', nrow = 2) +
scale_color_gradient(low = "#D3D3D3", high = "#a50000", na.value = NA)

grid.arrange(p1,p2,nrow=1)
dev.off()
