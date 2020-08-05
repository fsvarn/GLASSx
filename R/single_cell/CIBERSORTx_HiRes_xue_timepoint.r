###################################################
# Compare how macrophage activation changes pre- and post-treatment in samples exhibiting subtype switches
# Updated: 2020.07.17
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(ggdendro)
library(grid)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS/CIBERSORTxHiRes_GLASS_myeloid_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
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
fraction <- functional_length <- rep(0, length(mac_modules))
for(i in 1:length(mac_modules))
{
	mygenes <- mac_modules[[i]]
	denom <- length(mygenes)
	mygenes <- intersect(mygenes,rownames(geps))
	subgep <- geps[mygenes,]
	rem <- apply(subgep, 1, function(x)sum(x==1))
	mygenes <- mygenes[-which(rem == ncol(subgep))]
	fraction[i]<- length(mygenes)/denom
	functional_length[i] <- length(mygenes)
}
names(fraction) <- names(mac_modules)

module_info <- data.frame(fraction, functional_length)
module_info <- module_info[order(module_info[,"fraction"]),]



res <- t(gsva(data.matrix(geps), mac_modules, method="ssgsea",parallel.sz=1))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT ps.*, ts1.signature_name AS signature_name_a, ts2.signature_name AS signature_name_b
FROM analysis.rna_silver_set ps
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
"
dat <- dbGetQuery(con,q)

# Compare timepoints
time_eff <- apply(res, 2, function(x) median(x[dat[,"tumor_barcode_b"]] - x[dat[,"tumor_barcode_a"]]))
time_pvalue <- apply(res, 2, function(x)wilcox.test(x[dat[,"tumor_barcode_a"]], x[dat[,"tumor_barcode_b"]], paired=TRUE)$p.value)
time_qvalue <- p.adjust(time_pvalue, "BH")

subtype <- unique(dat[,"signature_name_a"])

eff <- p.val <- q.val <- matrix(0, nrow = ncol(res), ncol = length(subtype))
rownames(eff) <- rownames(p.val) <- rownames(q.val) <- colnames(res)
colnames(eff) <- colnames(p.val) <- colnames(q.val) <- subtype
for(i in 1:length(subtype))
{
	sub_dat <- dat[which(dat[,"signature_name_a"] == subtype[i] & dat[,"signature_name_b"] != subtype[i] | dat[,"signature_name_a"] != subtype[i] & dat[,"signature_name_b"] == subtype[i]),]
	ali1 <- sub_dat[which(sub_dat[,"signature_name_a"] == subtype[i]),"tumor_barcode_a"]
	ali2 <- sub_dat[which(sub_dat[,"signature_name_b"] != subtype[i]),"tumor_barcode_b"]
	ali3 <- sub_dat[which(sub_dat[,"signature_name_b"] == subtype[i]),"tumor_barcode_b"]
	ali4 <- sub_dat[which(sub_dat[,"signature_name_a"] != subtype[i]),"tumor_barcode_a"]
	
	ali_in <- c(ali1, ali3)
	ali_out <- c(ali2, ali4)
	
	eff[,i] <- apply(res, 2, function(x) median(x[ali_in] - x[ali_out]))
	mypval <- apply(res, 2, function(x)wilcox.test(x[ali_in], x[ali_out], paired=TRUE)$p.value)
	p.val[,i] <- mypval
	q.val[,i] <- p.adjust(mypval, "BH")
}
q.val[which(q.val[,3] < 0.05),]
eff[which(q.val[,3] < 0.05),]

# Write results to table
module <- rep(rownames(p.val), 3)
subtype <- rep(subtype, each = nrow(p.val))
signature_genes <- rep(fraction, 3)
p.value <- c(p.val[,1], p.val[,2], p.val[,3])
q.value <- c(q.val[,1], q.val[,2], q.val[,3])
median_effect <- c(eff[,1], eff[,2], eff[,3])
full_result <- data.frame(module, subtype, signature_genes, p.value, q.value, median_effect)

write.table(full_result, "/projects/verhaak-lab/GLASS-III/results/cibersortx/downstream/myeloid/glass_transcriptional_subtype_comparisons.txt",sep="\t",quote=FALSE,row.names=FALSE)

# IDHmut analysis (not enough switching)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT ps.*, ts1.signature_name AS signature_name_a, ts2.signature_name AS signature_name_b
FROM analysis.rna_silver_set ps
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE cs.idh_codel_subtype LIKE 'IDHmut%'
"
dat <- dbGetQuery(con,q)

# Compare timepoints
time_eff <- apply(res, 2, function(x) median(x[dat[,"tumor_barcode_b"]] - x[dat[,"tumor_barcode_a"]]))
time_pvalue <- apply(res, 2, function(x)wilcox.test(x[dat[,"tumor_barcode_a"]], x[dat[,"tumor_barcode_b"]], paired=TRUE)$p.value)


subtype <- unique(dat[,"signature_name_a"])

eff <- p.val <- q.val <- matrix(0, nrow = ncol(res), ncol = length(subtype))
rownames(eff) <- rownames(p.val) <- rownames(q.val) <- colnames(res)
colnames(eff) <- colnames(p.val) <- colnames(q.val) <- subtype
for(i in 1:length(subtype))
{
	sub_dat <- dat[which(dat[,"signature_name_a"] == subtype[i] & dat[,"signature_name_b"] != subtype[i] | dat[,"signature_name_a"] != subtype[i] & dat[,"signature_name_b"] == subtype[i]),]
	ali1 <- sub_dat[which(sub_dat[,"signature_name_a"] == subtype[i]),"tumor_barcode_a"]
	ali2 <- sub_dat[which(sub_dat[,"signature_name_b"] != subtype[i]),"tumor_barcode_b"]
	ali3 <- sub_dat[which(sub_dat[,"signature_name_b"] == subtype[i]),"tumor_barcode_b"]
	ali4 <- sub_dat[which(sub_dat[,"signature_name_a"] != subtype[i]),"tumor_barcode_a"]
	
	ali_in <- c(ali1, ali3)
	ali_out <- c(ali2, ali4)
	
	eff[,i] <- apply(res, 2, function(x) median(x[ali_in] - x[ali_out]))
	mypval <- apply(res, 2, function(x)wilcox.test(x[ali_in], x[ali_out], paired=TRUE)$p.value)
	p.val[,i] <- mypval
	q.val[,i] <- p.adjust(mypval, "BH")
}
q.val[which(q.val[,3] < 0.05),]
eff[which(q.val[,3] < 0.05),]


