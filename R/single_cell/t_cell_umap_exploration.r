
library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)

#######################################################
rm(list=ls())

set.seed(11) #Set seed so that code is consistent
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

##################################################
# Step 2: UMAP T cells, which do not associate with subtype
##################################################

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("t_cell",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_TCGA_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")


geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
colnames(geps) <- paste(mytag, colnames(geps), sep="__")
geps <- log10(geps+1)

rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]

embedding <- umap(t(geps))$layout


plot_res <- embedding
aliquot_barcode <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[2])
cell_state <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[1])
plot_res <- data.frame(aliquot_barcode, plot_res, cell_state)
colnames(plot_res) <- c("aliquot_barcode","UMAP1","UMAP2", "cell_state")
#plot_res <- plot_res %>%
#	inner_join(sig_sub, by = "aliquot_barcode")
rownames(plot_res) <- paste(plot_res[,"cell_state"], plot_res[,"aliquot_barcode"], sep = "__")

# UMAP plot
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_tcga_tcell_idhwt_umap.pdf",width=2.5,height=1.75)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(UMAP1, UMAP2)) + 
geom_point()  +
theme_bw()
#scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()

##################################################
# Step 2: Examine T cell clustering
##################################################

tcell_umap <- plot_res %>% 
		 filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])) %>%
		 filter(cell_state == "t_cell")

tcell_gep <- t(geps)
tcell_gep <- tcell_gep[rownames(tcell_umap),]

umap1_cor <- cor(tcell_gep, tcell_umap[,"UMAP1"], method="s")
umap1_cor <- umap1_cor[order(umap1_cor),]

umap2_cor <- cor(tcell_gep, tcell_umap[,"UMAP2"], method="s")
umap2_cor <- umap2_cor[order(umap2_cor),]

# Compare two groups

g1 <- rownames(tcell_umap[which(tcell_umap[,"UMAP2"] > 0),])
g2 <- rownames(tcell_umap[which(tcell_umap[,"UMAP2"] < 0),])

umap_dif <- apply(tcell_gep, 2, function(x) wilcox.test(x[g1], x[g2])$p.value)
umap_dif <- umap_dif[order(umap_dif)]
umap_dif <- p.adjust(umap_dif, "BH")
mygenes <- names(umap_dif[which(umap_dif < 0.1)])