###################################################
# Test different ways of subsetting the gene expression matrix
# Updated: 2020.07.16
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/CIBERSORTxHiRes_GLASS_myeloid_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps_nona <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps_nona_vars <- geps[which(vars > 0),]

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

res <- gsva(data.matrix(geps), mac_modules, method="ssgsea",parallel.sz=1)
res_nona <- gsva(data.matrix(geps_nona), mac_modules, method="ssgsea",parallel.sz=1)
res_nona_vars <- gsva(data.matrix(geps_nona_vars), mac_modules, method="ssgsea",parallel.sz=1)

# Add back missing modules to res_nona_vars
module_27 <- rep(0, ncol(res_nona_vars))
module_49 <- rep(0, ncol(res_nona_vars))
res_nona_vars <- rbind(res_nona_vars, module_27, module_49)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.rna_barcode_a AS aliquot_barcode, signature_name, 'Initial' AS timepoint,  cs.idh_codel_subtype, al.aliquot_batch
FROM analysis.platinum_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_a
WHERE cs.idh_codel_subtype = 'IDHwt'"
dat <- dbGetQuery(con,q)

mes <- dat[which(dat[,"signature_name"] == "Mesenchymal"), "aliquot_barcode"]
nonmes <- dat[which(dat[,"signature_name"] != "Mesenchymal"), "aliquot_barcode"]

#p.val and p.val_nona are virtually identical. Eliminating genes that are not variable removes too much information
p.val <- p.val_nona <- p.val_nona_vars <- p.val_nona_eff <- rep(0, nrow(res))
for(i in 1:length(mac_modules))
{
	mymod <- names(mac_modules)[i]
	p.val[i] <- wilcox.test(res[mymod,mes], res[mymod,nonmes])$p.value
	p.val_nona[i] <- wilcox.test(res_nona[mymod,mes], res_nona[mymod,nonmes])$p.value
	p.val_nona_eff[i] <- median(res_nona[mymod,mes]) - median(res_nona[mymod,nonmes])
	p.val_nona_vars[i] <- wilcox.test(res_nona_vars[mymod,mes], res_nona_vars[mymod,nonmes])$p.value
}
q.val_nona <- p.adjust(p.val_nona, "BH")
full_res <- data.frame(names(mac_modules), p.val_nona, q.val_nona, p.val_nona_eff)
colnames(full_res) <- c("module","p.val","q.val","eff")
full_res <- full_res[order(full_res[,"p.val"]),]