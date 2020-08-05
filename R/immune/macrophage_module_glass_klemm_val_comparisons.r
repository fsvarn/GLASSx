
library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(ggdendro)
library(grid)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/CIBERSORTxHiRes_GLASS_myeloid_Window48.txt"

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

res <- gsva(data.matrix(geps), mac_modules, method="ssgsea",parallel.sz=1)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
WITH long_pairs AS
(
	SELECT ps.rna_barcode_a AS aliquot_barcode, signature_name, 'Initial' AS timepoint,  cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
	JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_a
	UNION
	SELECT ps.rna_barcode_b AS aliquot_barcode, signature_name, 'Recurrent' AS timepoint, cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_b
	JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_b
)
SELECT * 
FROM long_pairs
ORDER BY aliquot_batch
"
dat <- dbGetQuery(con,q)

module <- rownames(res)
res <- data.frame(module, res,stringsAsFactors=FALSE)
colnames(res) <- gsub("\\.","-",colnames(res))
plot_res <- res %>%
 			pivot_longer(-module, names_to = "aliquot_barcode", values_to = "es") %>%
 			inner_join(dat, "aliquot_barcode") %>%
 			as.data.frame()

validation_profiles <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/macrophage_subtype_res.txt",stringsAsFactor=FALSE)

test <- validation_profiles %>%
	pivot_wider(id_cols="module", names_from = "aliquot_barcode", values_from = "es") %>%
	as.data.frame()
rownames(test) <- test[,1]
test <- test[,2:ncol(test)]
	
target_res <- data.matrix(res[,-1])
cor.res <- cor(target_res, test, method="s")

# Microglia vs macropahges
g1 <- apply(cor.res[,grep("_mdm",colnames(cor.res))],2,mean)
g2 <- apply(cor.res[,grep("_mg",colnames(cor.res))],2,mean)

wilcox.test(g1,g2)

# Try looking at correlation coefficients by IDH subtype and by transcriptional subtype