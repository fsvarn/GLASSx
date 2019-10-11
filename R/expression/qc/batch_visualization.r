library(DBI)
library(odbc)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

rm(list=ls())

#Load and prepare gene expression matrix
myinf1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"

genes <- read.delim(myinf1,row.names=1)
colnames(genes) <- gsub("\\.","-",colnames(genes))

gene_sums <- apply(genes,1,sum)
sub_genes <- genes[which(gene_sums>0),]
gene_var <- apply(sub_genes,1,var)
gene_var <- gene_var[order(gene_var,decreasing=TRUE)]

#Load clinical data
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT al.aliquot_barcode, sa.sample_type, su.*, aliquot_batch
FROM biospecimen.aliquots al 
LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode
LEFT JOIN biospecimen.samples sa ON al.sample_barcode = sa.sample_barcode
WHERE aliquot_analyte_type = 'R'"

clin <- dbGetQuery(con,q)

log_genes <- log10(sub_genes+1)
log_sub_genes <- log_genes[names(gene_var[1:500]),]

#Annotate matrix and add heatmap annotations
big_annotation_table <- clin[match(colnames(log_genes),clin[,"aliquot_barcode"]),]

#Cluster differences
annotation_table <- big_annotation_table[,c("idh_codel_subtype","sample_type","aliquot_batch")]
annotation_table <- annotation_table %>% mutate(idh_codel_subtype = idh_codel_subtype,
										 		sample_type = recode(sample_type,"TP"="Primary","R1"="Recurrent","R2"="Recurrent","R3"="Recurrent"),
										 		aliquot_batch = aliquot_batch)

#Set colors for sidebars
sidebar_labels <- apply(annotation_table,2,function(x)unique(x[!is.na(x)])[order(unique(x[!is.na(x)]))])
idh_codel_subtype = c("#F8766D","#00BA38","#619CFF"); names(idh_codel_subtype) = sidebar_labels[[1]]
sample_type = c("black","white"); names(sample_type) = sidebar_labels[[2]]
aliquot_batch =  brewer.pal(12,"Set3")[1:11]; names(aliquot_batch) = sidebar_labels[[3]]
annotation_colors = list(idh_codel_subtype=idh_codel_subtype,
					sample_type=sample_type,
					aliquot_batch=aliquot_batch)
ha = HeatmapAnnotation(df = annotation_table,which="column",
	 col=annotation_colors,show_annotation_name=TRUE,show_legend=TRUE,na_col="#FFFFFF")


#Plot heatmap to check for batch effect
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/var_gene_heatmap.pdf",width=7,height=7)
hm = Heatmap(t(scale(t(log_sub_genes))),
		col = colorRamp2(c(-1,0,1),c("royalblue4","white","tomato3"),space="RGB"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=TRUE,
		show_column_names=FALSE,
		show_heatmap_legend=FALSE,
		bottom_annotation=ha)
hm
dev.off()

#----------------------

#Quantile normalize original gene expression matrix
myrk <- qn_genes <- matrix(0, nrow(genes), ncol(genes))
rownames(qn_genes) <- rownames(genes)
colnames(qn_genes) <- colnames(genes)
xx <- myrk
for(k in 1:ncol(genes))
{
	myrk[,k] <- rank(genes[,k])
	xx[,k] <- sort(genes[,k])
}
mymed = apply(xx, 1, median, na.rm=T)
for(k in 1:ncol(genes))
{
	qn_genes[,k] <- mymed[myrk[,k]]
}	

qn_gene_sums <- apply(qn_genes,1,sum)
sub_qn_genes <- qn_genes[which(qn_gene_sums>0),]
qn_gene_var <- apply(sub_qn_genes,1,var)
qn_gene_var <- qn_gene_var[order(qn_gene_var,decreasing=TRUE)]

log_qn_genes <- log10(sub_qn_genes+1)
log_sub_qn_genes <- log_qn_genes[names(qn_gene_var[1:500]),]

#Plot heatmap to see if quantile normalization diminished batch effect
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/var_gene_heatmap_qn.pdf",width=7,height=7)
hm = Heatmap(t(scale(t(log_sub_qn_genes))),
		col = colorRamp2(c(-1,0,1),c("royalblue4","white","tomato3"),space="RGB"),
		cluster_rows=TRUE,
		show_row_names=FALSE,
		cluster_columns=TRUE,
		show_column_names=FALSE,
		show_heatmap_legend=FALSE,
		bottom_annotation=ha)
hm
dev.off()
