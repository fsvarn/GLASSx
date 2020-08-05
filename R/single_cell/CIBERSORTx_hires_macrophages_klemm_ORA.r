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
bg <- rownames(geps)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.rna_barcode_a AS aliquot_barcode, signature_name, 'Initial' AS timepoint,  cs.idh_codel_subtype, al.aliquot_batch
FROM analysis.platinum_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_a
WHERE cs.idh_codel_subtype = 'IDHwt'"
dat <- dbGetQuery(con,q)


sub_geps <- geps[,dat[,"aliquot_barcode"]]
p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(i in 1:nrow(geps))
{
	mygene <- as.numeric(sub_geps[i,dat[,"aliquot_barcode"]])
	mes <- as.numeric(dat[,"signature_name"] == "Mesenchymal")
	testdat <- cbind(dat, mygene)
	
	p.val[i] <- summary(lm(mygene ~ mes + aliquot_batch, data = testdat))$coefficients[38]
	eff[i] <- summary(lm(mygene ~ mes + aliquot_batch, data = testdat))$coefficients[2]
}
q.val <- p.adjust(p.val,"BH")
mes_res <- data.frame(p.val, q.val, eff)
mes_res <- mes_res[order(p.val),]

sig_genes <- rownames(mes_res)[which(mes_res[,"q.val"] < 0.05)]

# Examine how they associate with the Xue macrophage modules

module_inf <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"

mac_module_df <- read.delim(module_inf,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- gsub("X","module",colnames(mac_module_df))

ora.p.value <- rep(0,ncol(mac_module_df))
names(ora.p.value) <- colnames(mac_module_df)
for(i in 1:ncol(mac_module_df))
{
	mac_mod <- mac_module_df[,i]
	mac_mod <- intersect(mac_mod, bg)
	
	deg_in_mod_in <- length(intersect(sig_genes, mac_mod))
	deg_in_mod_out <- length(sig_genes) -  deg_in_mod_in
	deg_out_mod_in <- length(mac_mod) - deg_in_mod_in
	deg_out_mod_out <- length(bg) - deg_in_mod_in - deg_in_mod_out - deg_out_mod_in
	
	ct <- matrix(c(deg_in_mod_out, deg_in_mod_in, deg_out_mod_out, deg_out_mod_in), nrow=2, ncol=2)
	ora.p.value[i] <- fisher.test(ct)$p.value
}
ora.p.value <- ora.p.value[order(ora.p.value)]
ora.q.value <- p.adjust(ora.p.value,method="BH")

# Validate using ssGSEA