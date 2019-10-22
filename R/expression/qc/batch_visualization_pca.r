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


#PCA plot
mypca <- prcomp(t(sub_genes),scale.=TRUE)
scores <- as.data.frame(mypca$x)
scores <- scores[clin[,1],]
scores[,"aliquot_batch"] <- as.factor(clin[,"aliquot_batch"])

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/pca_genes.pdf",width=7,height=7)
ggplot(data = scores, aes(x = PC1, y = PC2, colour=aliquot_batch)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  scale_colour_manual(values=brewer.pal(11,"Set3")) +
  ggtitle("PCA plot of GLASS RNAseq aliquots") +
  theme_classic()
dev.off()

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

mypca <- prcomp(t(sub_qn_genes),scale.=TRUE)
scores <- as.data.frame(mypca$x)
scores <- scores[clin[,1],]
scores[,"aliquot_batch"] <- as.factor(clin[,"aliquot_batch"])
scores[,"subtype"] <- as.factor(clin[,"idh_codel_subtype"])

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/pca_genes_qn.pdf",width=7,height=7)
ggplot(data = scores, aes(x = PC1, y = PC2, colour=aliquot_batch)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  scale_colour_manual(values=brewer.pal(11,"Set3")) +
  ggtitle("PCA plot of GLASS RNAseq aliquots") +
  theme_classic()
dev.off()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/pca_genes_qn_with_subtype.pdf",width=7,height=7)
ggplot(data = scores, aes(x = PC1, y = PC2, colour=aliquot_batch, shape = subtype)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  scale_colour_manual(values=brewer.pal(11,"Set3")) +
  ggtitle("PCA plot of GLASS RNAseq aliquots") +
  theme_classic()
dev.off()
