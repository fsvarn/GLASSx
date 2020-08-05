###################################################
# Compare macrophage activation profiles in TCGA between transcriptional subtypes
# Updated: 2020.07.22
# Author: Frederick Varn
##################################################
library(tidyverse)
library(ssgsea.GBM.classification)
library(GSVA)
library(RColorBrewer)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.txt"

expr <- read.delim(myinf1,stringsAsFactors=FALSE)

unigene <- table(expr[,"GeneSymbol"])
expr1 <- expr[which(unigene==1),]

# Only one gene symbol, can throw out as its mostly low compared to its duplicate
expr2 <- expr[which(unigene>1),]

rownames(expr1) <- expr1[,1]
expr <- expr1[,-1]
expr <- data.matrix(expr)


##################################################
# Helper functions
##################################################

# GCT creator function (for transcriptional classifier)
#-------------------------------------------

gct <- function(data, gct_out)
{
	Description <- rep(NA,nrow(data))
	data2 <- cbind(rownames(data),Description,data[,1:ncol(data)])
	colnames(data2)[1] <- "NAME"
	write.table(data2, gct_out, sep="\t", quote=FALSE, row.names=FALSE)

	conIn <- file(gct_out, "r")
	rawfile = readLines(conIn)
	close(conIn)

	mytext <- c("#1.2", paste(nrow(data2),"\t",(ncol(data)-1),sep=""),rawfile)
	conOut = file(gct_out, "w")
	writeLines(mytext, conOut)
	close(conOut)
}
	
##################################################
# Step 1: Subtype each TCGA using bulk RNAseq data
##################################################

# Save glioma-specific expr file in GCT format for classifier

tcga_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct"
gct(expr, tcga_path)
	
# Run Qianghu's transcriptional classifier (only needed once)
#runSsGSEAwithPermutation(tcga_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/p_result_TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct.txt")

aliquot_barcode <- rep(rownames(subtype_ssgsea),3)
signature_name <- c(rep("Proneural",nrow(subtype_ssgsea)),rep("Classical",nrow(subtype_ssgsea)),rep("Mesenchymal",nrow(subtype_ssgsea)))
enrichment_score <- c(subtype_ssgsea[,"Proneural"],subtype_ssgsea[,"Classical"],subtype_ssgsea[,"Mesenchymal"])
p_value <- c(subtype_ssgsea[,"Proneural_pval"],subtype_ssgsea[,"Classical_pval"],subtype_ssgsea[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(aliquot_barcode,signature_name,enrichment_score,p_value)
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"aliquot_barcode"]),]

#Calculate simplicity scores and assign subtypes
aliquots <- unique(as.character(transcriptional_subtype[,"aliquot_barcode"]))

simplicity_score <- subtype_class <- rep(NA,length(aliquots))
for(i in 1:length(aliquots))
{
	sub_dat <- transcriptional_subtype[which(transcriptional_subtype[,"aliquot_barcode"] == aliquots[i]),]
	sub_dat[,"p_rank"] <- rank(sub_dat[,"p_value",],ties.method="min")
	
	subtype_class[i] <- paste(sub_dat[which(sub_dat[,"p_value"] == min(sub_dat[,"p_value"])),"signature_name"],collapse=",")
	
	r0 <- sub_dat[which(sub_dat[,"p_rank"] ==1),"p_value"][1]
	ri <- sub_dat[which(sub_dat[,"p_rank"] > 1),"p_value"]
	ri <- ri[order(ri)]
	
	adds <- sum(ri - r0)
	
	d <- abs(outer(ri,ri,"-"))
	diag(d) <- NA
	d[lower.tri(d)] <- NA
	adns <- sum(d,na.rm=TRUE)
	
	rn1 <- sub_dat[which(sub_dat[,"p_rank"] == max(sub_dat[,"p_rank"])),"p_value"][1]
	n1 <- 2	#Number of unique subtypes - 1
	simplicity_score[i] <- (adds - adns) * (rn1 - r0)/n1
}

# Store full results in a data frame
transcriptional_subtype <- data.frame(aliquots, subtype_class, simplicity_score,stringsAsFactors=FALSE)
colnames(transcriptional_subtype) <- c("aliquot_barcode","transcriptional_subtype","simplicity_score")

transcriptional_subtype[,"case_barcode"] <- sapply(strsplit(transcriptional_subtype[,"aliquot_barcode"],"\\."),function(x)paste(x[1:3],collapse="-"))

##################################################
# Step 2: Assess macrophage activation activity on CIBERSORTx macrophage profiles
##################################################

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/CIBERSORTxHiRes_TCGA_myeloid_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
#colnames(geps) <- gsub("\\.","-",colnames(geps))
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]

myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"
mac_module_df <- read.delim(myinf2,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- paste("module_",1:49,sep="")

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)

# Test how many genes are expressed per module
fraction <- rep(0, length(mac_modules))
for(i in 1:length(mac_modules))
{
	mygenes <- mac_modules[[i]]
	denom <- length(mygenes)
	mygenes <- intersect(mygenes,rownames(geps))
	subgep <- geps[mygenes,]
	rem <- apply(subgep, 1, function(x)sum(x!=1))
	mygenes <- mygenes[which(rem == ncol(subgep))]
	fraction[i]<- length(mygenes)/denom
}
names(fraction) <- names(mac_modules)
fraction <- fraction[order(fraction)]


res <- t(gsva(data.matrix(geps), mac_modules, method="ssgsea",parallel.sz=1))

##################################################
# Step 3: Compare subtypes
##################################################

info <- read.delim("/projects/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt",stringsAsFactor=FALSE)
mycase <- info[which(info[,"IDH.codel.subtype"] == "IDHwt"),"Case"]

idhwt <- transcriptional_subtype[which(transcriptional_subtype[,"case_barcode"] %in% mycase),]

subtype <- unique(idhwt[,"transcriptional_subtype"])
subtype <- subtype[-grep(",",subtype)]

eff <- p.val <- q.val <- matrix(0, nrow = ncol(res), ncol = length(subtype))
rownames(eff) <- rownames(p.val) <- rownames(q.val) <- colnames(res)
colnames(eff) <- colnames(p.val) <- colnames(q.val) <- subtype
for(i in 1:length(subtype))
{

	mymes <- idhwt[grep(subtype[i], idhwt[,"transcriptional_subtype"]),"aliquot_barcode"]
	mynonmes <- idhwt[-grep(subtype[i], idhwt[,"transcriptional_subtype"]),"aliquot_barcode"]

	eff[,i] <- apply(res, 2, function(x)median(x[mymes]) - median(x[mynonmes]))
	mypvalue <- apply(res, 2, function(x)wilcox.test(x[mymes],x[mynonmes])$p.value)
	p.val[,i] <- mypvalue
	q.val[,i] <- p.adjust(mypvalue,"BH")
}

# Write results to table
module <- rep(rownames(p.val), 3)
subtype <- rep(subtype, each = nrow(p.val))
signature_genes <- rep(fraction, 3)
p.value <- c(p.val[,1], p.val[,2], p.val[,3])
q.value <- c(q.val[,1], q.val[,2], q.val[,3])
median_effect <- c(eff[,1], eff[,2], eff[,3])
full_result <- data.frame(module, subtype, signature_genes, p.value, q.value, median_effect)

write.table(full_result, "/projects/verhaak-lab/GLASS-III/results/cibersortx/downstream/myeloid/tcga_transcriptional_subtype_comparisons.txt",sep="\t",quote=FALSE,row.names=FALSE)
