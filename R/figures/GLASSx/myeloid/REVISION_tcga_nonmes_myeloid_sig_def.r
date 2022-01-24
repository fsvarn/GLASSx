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

##################################################
# Step 2: UMAP of each CIBERSORTx GEP
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
geps <- log10(geps+1)

rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]


# Define the proneural signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Proneural" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(!grepl("Proneural",sig_sub[,"signature_name"]) & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(j in 1:nrow(geps))
{
	group1 <- as.numeric(g1[j,])
	group2 <- as.numeric(g2[j,])

	p.val[j] <- wilcox.test(group1,group2)$p.value
	eff[j] <- log2(median(group1)/median(group2))
	#eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
pro_res <- data.frame(p.val, q.val, eff)
pro_res <- pro_res[order(p.val),]

pro_res[,"logp"] <- -log10(pro_res[,"p.val"])
pro_res[,"sig"] <- pro_res[,"q.val"] < 0.1
	
# Define the classical signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Classical" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(!grepl("Classical",sig_sub[,"signature_name"]) & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(j in 1:nrow(geps))
{
	group1 <- as.numeric(g1[j,])
	group2 <- as.numeric(g2[j,])

	p.val[j] <- wilcox.test(group1,group2)$p.value
	eff[j] <- log2(median(group1)/median(group2))
	#eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
class_res <- data.frame(p.val, q.val, eff)
class_res <- class_res[order(p.val),]

class_res[,"logp"] <- -log10(class_res[,"p.val"])
class_res[,"sig"] <- class_res[,"q.val"] < 0.1

# Define the non-mesenchymal signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Mesenchymal" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(!grepl("Mesenchymal",sig_sub[,"signature_name"]) & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(j in 1:nrow(geps))
{
	group1 <- as.numeric(g1[j,])
	group2 <- as.numeric(g2[j,])

	p.val[j] <- wilcox.test(group1,group2)$p.value
	eff[j] <- log2(median(group2)/median(group1))			# Reverse the numerator/denominator here
	#eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
idhwt_res <- data.frame(p.val, q.val, eff)
idhwt_res <- idhwt_res[order(p.val),]

idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1

pro_res[,"analysis"] <- "Proneural"
class_res[,"analysis"] <- "Classical"
idhwt_res[,"analysis"] <- "Non-mesenchymal"

pro_res <- pro_res %>% rownames_to_column("gene")
class_res <- class_res %>% rownames_to_column("gene")
idhwt_res <- idhwt_res %>% rownames_to_column("gene")

full_res <- rbind(pro_res, class_res, idhwt_res)

# Save the myeloid results to apply to GLASS
write.table(full_res, "data/res/CIBERSORTx/analysis/TCGA_nonmes_myeloid_ts_result.txt",sep="\t",quote=FALSE,row.names=TRUE) 
#v2 is with the directionality so that mes is the "change" and non-mes is the reference
#v3 is the signature developed from IDHwt only
