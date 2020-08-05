###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/CIBERSORTxHiRes_GLASS_myeloid_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cs.idh_codel_subtype
FROM analysis.platinum_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'"
dat <- dbGetQuery(con,q)

g1 <- geps[,dat[,"rna_barcode_a"]]
g2 <- geps[,dat[,"rna_barcode_b"]]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(i in 1:nrow(geps))
{
	group1 <- as.numeric(g1[i,])
	group2 <- as.numeric(g2[i,])
	
	p.val[i] <- wilcox.test(group1,group2)$p.value
	eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
idhwt_res <- data.frame(p.val, q.val, eff)
idhwt_res <- idhwt_res[order(p.val),]


q <- "SELECT ps.rna_barcode_a, signature_name
FROM analysis.platinum_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
WHERE cs.idh_codel_subtype = 'IDHwt'"
dat <- dbGetQuery(con,q)

sub_geps <- geps[,dat[,"rna_barcode_a"]]
g1 <- sub_geps[,dat[which(dat[,"signature_name"]=="Mesenchymal"),"rna_barcode_a"]]
g2 <- sub_geps[,dat[which(dat[,"signature_name"]!="Mesenchymal"),"rna_barcode_a"]]


p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(i in 1:nrow(geps))
{
	group1 <- as.numeric(g1[i,])
	group2 <- as.numeric(g2[i,])
	
	p.val[i] <- wilcox.test(group1,group2)$p.value
	eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
mes_res <- data.frame(p.val, q.val, eff)
mes_res <- mes_res[order(p.val),]
