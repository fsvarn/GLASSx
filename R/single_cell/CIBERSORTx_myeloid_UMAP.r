###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(umap)

# DO THIS ALL BY BATCH

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

# run UMAP algorithm
embedding <- umap(t(geps))

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")


q <- "SELECT ps.rna_barcode_a AS aliquot_barcode, signature_name, 'Initial' AS timepoint,  cs.idh_codel_subtype, al.aliquot_batch
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
"
dat <- dbGetQuery(con,q)

plot_res <- data.frame(embedding$layout )
aliquot_barcode <- rownames(plot_res)
aliquot_barcode <- gsub("\\.","-",aliquot_barcode)
rownames(plot_res) <- NULL
plot_res <- cbind(aliquot_barcode, plot_res)
colnames(plot_res) <- c("aliquot_barcode","UMAP1","UMAP2")
	plot_res %>%
	inner_join(dat, by = "aliquot_barcode") %>%
	filter(idh_codel_subtype == 'IDHwt')

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_SCGP_hires_myeloid_idhwt_by_batch_umap.pdf",width=3,height=3)
ggplot(plot_res, aes(UMAP1, UMAP2, color = aliquot_batch, shape = timepoint)) + 
geom_point()  +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()
	
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_SCGP_hires_myeloid_idhwt_sm_umap.pdf",width=3,height=3)
plot_res %>%
inner_join(dat, by = "aliquot_barcode") %>%
filter(idh_codel_subtype == 'IDHwt', aliquot_batch == 'GLSS-SM-RNA') %>%
ggplot(aes(UMAP1, UMAP2, color = signature_name, shape = timepoint)) + 
geom_point()  +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_SCGP_hires_myeloid_idhwt_cu_umap.pdf",width=3,height=3)
plot_res %>%
inner_join(dat, by = "aliquot_barcode") %>%
filter(idh_codel_subtype == 'IDHwt', aliquot_batch == 'GLSS-CU-RNA') %>%
ggplot(aes(UMAP1, UMAP2, color = signature_name, shape = timepoint)) + 
geom_point()  +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_SCGP_hires_myeloid_idhwt_tcga_umap.pdf",width=3,height=3)
plot_res %>%
inner_join(dat, by = "aliquot_barcode") %>%
filter(idh_codel_subtype == 'IDHwt', aliquot_batch == 'TCGA-GB-RNA') %>%
ggplot(aes(UMAP1, UMAP2, color = signature_name, shape = timepoint)) + 
geom_point()  +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()