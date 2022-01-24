library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)
library(topGO)
library(UpSetR)

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


# Define the proneural vs classical signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Proneural" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Classical" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

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
pro_class <- data.frame(p.val, q.val, eff)
pro_class <- pro_class[order(p.val),]

pro_class[,"logp"] <- -log10(pro_class[,"p.val"])
pro_class[,"sig"] <- pro_class[,"q.val"] < 0.1
sig1 <- pro_class %>% filter(sig, eff > log2(1.1))

# Define the proneural vs mesenchymal signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Proneural" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Mesenchymal" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

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
pro_mes <- data.frame(p.val, q.val, eff)
pro_mes <- pro_mes[order(p.val),]

pro_mes[,"logp"] <- -log10(pro_mes[,"p.val"])
pro_mes[,"sig"] <- pro_mes[,"q.val"] < 0.1
sig2 <- pro_mes %>% filter(sig, eff > log2(1.1))


# Define the classical vs proneural signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Classical" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Proneural" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

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
mes_class <- data.frame(p.val, q.val, eff)
class_pro <- class_pro[order(p.val),]

class_pro[,"logp"] <- -log10(class_pro[,"p.val"])
class_pro[,"sig"] <- class_pro[,"q.val"] < 0.1
sig3 <- class_pro %>% filter(sig, eff > log2(1.1))

# Define the classical vs mesenchymal signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Classical" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Mesenchymal" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

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
class_mes <- data.frame(p.val, q.val, eff)
class_mes <- class_mes[order(p.val),]

class_mes[,"logp"] <- -log10(class_mes[,"p.val"])
class_mes[,"sig"] <- class_mes[,"q.val"] < 0.1
sig4 <- class_mes %>% filter(sig, eff > log2(1.1))

# Define the mesenchymal vs proneural signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Mesenchymal" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Proneural" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

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
mes_pro <- data.frame(p.val, q.val, eff)
mes_pro <- mes_pro[order(p.val),]

mes_pro[,"logp"] <- -log10(mes_pro[,"p.val"])
mes_pro[,"sig"] <- mes_pro[,"q.val"] < 0.1
sig5 <- mes_pro %>% filter(sig, eff > log2(1.1))

# Define the mesenchymal vs classical signature
g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Mesenchymal" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Classical" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]

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
mes_class <- data.frame(p.val, q.val, eff)
mes_class <- mes_class[order(p.val),]

mes_class[,"logp"] <- -log10(mes_class[,"p.val"])
mes_class[,"sig"] <- mes_class[,"q.val"] < 0.1
sig6 <- mes_class %>% filter(sig, eff > log2(1.1))

# Make UpSet plot input

# Option 1
mes_base <- list("Proneural" = rownames(sig2), "Classical" = rownames(sig4))
class_base <- list("Proneural" = rownames(sig1), "Mesenchymal" = rownames(sig6))
pro_base <- list("Classical" = rownames(sig3), "Mesenchymal" = rownames(sig5))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/proneural_base_upset.pdf")  
upset(fromList(pro_base), keep.order = TRUE, sets.bar.color = NA, mainbar.y.max=250, empty.intersections="on")
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/classical_base_upset.pdf")  
upset(fromList(class_base), keep.order = TRUE, sets.bar.color = NA, mainbar.y.max=250, empty.intersections="on")
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/mesenchymal_base_upset.pdf")  
upset(fromList(mes_base), keep.order = TRUE, sets.bar.color = NA, mainbar.y.max=250, empty.intersections="on")
dev.off()

# Option 2
pro_base <- list("Classical" = rownames(sig1), "Mesenchymal" = rownames(sig2))
class_base <- list("Proneural" = rownames(sig3), "Mesenchymal" = rownames(sig4))
mes_base <- list("Proneural" = rownames(sig5), "Classical" = rownames(sig6))

dims = 3
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/proneural_comparison_upset.pdf",width=dims, height=dims)  
upset(fromList(pro_base), keep.order = TRUE, sets.bar.color = NA, mainbar.y.max=200, empty.intersections="on")
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/classical_comparison_upset.pdf",width=dims, height=dims)   
upset(fromList(class_base), keep.order = TRUE, sets.bar.color = NA, mainbar.y.max=200, empty.intersections="on")
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/mesenchymal_comparison_upset.pdf",width=dims, height=dims)  
upset(fromList(mes_base), keep.order = TRUE, sets.bar.color = NA, mainbar.y.max=200, empty.intersections="on")
dev.off()

# Save the myeloid results to apply to GLASS
write.table(full_res, "data/res/CIBERSORTx/analysis/TCGA_nonmes_myeloid_ts_result.txt",sep="\t",quote=FALSE,row.names=TRUE) 
#v2 is with the directionality so that mes is the "change" and non-mes is the reference
#v3 is the signature developed from IDHwt only
