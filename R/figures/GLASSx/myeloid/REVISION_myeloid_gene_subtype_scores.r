library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)
library(topGO)

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

#transcriptional_subtype[,"signif"] <- transcriptional_subtype[,"p_value"] < 0.05

sig_sub <- transcriptional_subtype
		   
sig_sub[,"case_barcode"] <- sapply(strsplit(as.character(sig_sub[,"aliquot_barcode"]),"-"),function(x)paste(x[1:3],collapse="-"))
sig_sub <- sig_sub %>%
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case")) 

##################################################
# Step 2: Correlate myeloid cells with subtype scores
##################################################

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_TCGA_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
colnames(geps) <- paste(mytag, colnames(geps), sep="__")
colnames(geps) <- gsub("myeloid__","",colnames(geps))
geps <- log10(geps+1)

rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
#geps <- geps[which(vars > 0),]		We want these genes to determine pan-myeloid genes


sigs <- c("Proneural", "Classical", "Mesenchymal")
cor_coef <- p_val <- subtype <- gene <- c()
for(i in 1:length(sigs))
{
	sig_score <- sig_sub %>% filter(signature_name == sigs[i], IDH.codel.subtype == "IDHwt")
	sub_geps <- geps[,sig_score[,"aliquot_barcode"]]

	cor_coef <- c(cor_coef, apply(sub_geps, 1, function(x)cor(x,sig_score$enrichment_score,method="s")))
	p_val <- c(p_val, apply(sub_geps, 1, function(x)cor.test(x,sig_score$enrichment_score,method="s")$p.value))
	q_val <- p.adjust(p_val, "BH")
	gene <- c(gene, rownames(sub_geps))
	subtype <- c(subtype, rep(sigs[i], nrow(sub_geps)))
}

cor_res <- data.frame(gene, cor_coef, p_val, q_val, subtype)

proneural_sig <- cor_res %>% filter(cor_coef > 0, q_val < 0.05, subtype == "Proneural")
classical_sig <- cor_res %>% filter(cor_coef > 0, q_val < 0.05, subtype == "Classical")
mesenchymal_sig <- cor_res %>% filter(cor_coef > 0, q_val < 0.05, subtype == "Mesenchymal")

# Define a basal signature
basal_sig <- cor_res %>% 
			 group_by(gene) %>%
			 summarise(insig = sum(q_val > 0.05 | cor_coef < 0 | is.na(q_val))) %>%
			 filter(insig == 3)


##################################################
# Step 3: Test signatures in single-cell data
##################################################


myinf1 <- "/projects/verhaak-lab/scgp/data/10x/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds"

load(myinf1)

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

# Read in bulk subtype results
con <- DBI::dbConnect(odbc::odbc(), "scgp")

bulk_res <- dbReadTable(con, Id(schema = "analysis", table="transcriptional_subtype"))
bulk_mes <- bulk_res %>% 
			#filter(signature_name == "Mesenchymal") %>%
			mutate(sample_id = sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[2],x[3],sep=""))) #%>%
			#filter(sample_id %in% c("SM006","SM012","SM017","SM018","SM011"))
			
sigs <- c("Proneural", "Classical", "Mesenchymal")
valid_cor <- rep(0, length(sigs))
for(i in 1:length(sigs))
{
	mysig <- cor_res %>% filter(cor_coef > 0.5, q_val < 0.05, subtype == sigs[i]) %>% .$gene %>% as.character()

	sc_mes_score <- apply(log2cpm_myeloid[mysig,],2,mean)
	sample_id <- clust_annot[names(sc_mes_score),"sample_id"]
	plot_res <- data.frame(names(sc_mes_score), sc_mes_score, sample_id)
	colnames(plot_res) <- c("cell_id","sc_mes_score","sample_id")			

	tmp_bulk <- bulk_mes %>% filter(signature_name == sigs[i])
	
	plot_res <- plot_res %>% 
	group_by(sample_id) %>% 
	summarise(mean = mean(sc_mes_score)) %>% 
	ungroup() %>%
	inner_join(tmp_bulk, by = "sample_id") %>%
	as.data.frame()
	
	valid_cor[i] <- cor(plot_res$mean, plot_res$enrichment_score)
}
valid_cor

basal_cor <- rep(0, length(sigs))
for(i in 1:length(sigs))
{
	cat("\r",i)
	mysig <- as.character(basal_sig$gene)

	sc_mes_score <- apply(log2cpm_myeloid[mysig,],2,mean)
	sample_id <- clust_annot[names(sc_mes_score),"sample_id"]
	plot_res <- data.frame(names(sc_mes_score), sc_mes_score, sample_id)
	colnames(plot_res) <- c("cell_id","sc_mes_score","sample_id")			

	tmp_bulk <- bulk_mes %>% filter(signature_name == sigs[i])
	
	plot_res <- plot_res %>% 
	group_by(sample_id) %>% 
	summarise(mean = mean(sc_mes_score)) %>% 
	ungroup() %>%
	inner_join(tmp_bulk, by = "sample_id") %>%
	as.data.frame()
	
	basal_cor[i] <- cor(plot_res$mean, plot_res$enrichment_score)
}
