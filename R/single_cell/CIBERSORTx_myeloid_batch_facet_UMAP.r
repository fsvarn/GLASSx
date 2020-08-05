###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(umap)


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

q <- "
WITH long_pairs AS
(
	SELECT ps.rna_barcode_a AS aliquot_barcode, signature_name, 'Initial' AS timepoint,  cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
	JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_a
	WHERE cs.idh_codel_subtype = 'IDHwt'
	UNION
	SELECT ps.rna_barcode_b AS aliquot_barcode, signature_name, 'Recurrent' AS timepoint, cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_b
	JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_b
	WHERE cs.idh_codel_subtype = 'IDHwt'
)
SELECT * 
FROM long_pairs
ORDER BY aliquot_batch
"
dat <- dbGetQuery(con,q)

batches <- unique(dat[,"aliquot_batch"])

batch_list <- list()
for(i in 1:length(batches))
{
	sub_dat <- dat %>%
		filter(aliquot_batch == batches[i])
	sub_geps <- geps[,sub_dat[,"aliquot_barcode"]]
	
	neighbors <- max(c(3, 0.2 * ncol(sub_geps)))
	
	custom.settings = umap.defaults
	custom.settings$n_neighbors = neighbors
	
	# run UMAP algorithm
	embedding <- umap(t(sub_geps), config = custom.settings)
	
	plot_res <- data.frame(embedding$layout )
	aliquot_barcode <- rownames(plot_res)
	aliquot_barcode <- gsub("\\.","-",aliquot_barcode)
	rownames(plot_res) <- NULL
	plot_res <- cbind(aliquot_barcode, plot_res)
	colnames(plot_res) <- c("aliquot_barcode","UMAP1","UMAP2")
	plot_res <-	plot_res %>%
		inner_join(sub_dat, by = "aliquot_barcode") 
	
	batch_list[[i]] <- plot_res	
}

batch_plot <- do.call(rbind,batch_list)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_SCGP_hires_myeloid_batch_facet_umap.pdf",width=7,height=7)
ggplot(batch_plot, aes(UMAP1, UMAP2, color = signature_name, shape = timepoint)) + 
geom_point()  +
facet_wrap(vars(aliquot_batch), nrow = 4) +
scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()
