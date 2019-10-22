library(ggplot2)
library(odbc)
library(DBI)

rm(list=ls())

myinf1 <- "/projects/varnf/GLASS-III/GLASS-III/results/kallisto/kallisto/final/gene_tpm_matrix_all_samples.tsv"
myinf2 <- "/projects/varnf/SofWar/kallisto/Ensembl_v96_homo_sapiens/homo_sapiens/rRNA_annotation_table.txt"

data <- data.matrix(read.delim(myinf1,row.names=1))
colnames(data) <- gsub("\\.","-",colnames(data))
#data <- log10(data+1)

annotation <- read.delim(myinf2,stringsAsFactor=FALSE)
mygene <- unique(annotation[,"GeneSymbol"])

#Load clinical data
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT al.aliquot_barcode, sa.sample_type, su.*, aliquot_batch
FROM biospecimen.aliquots al 
LEFT JOIN clinical.surgeries su ON su.sample_barcode = al.sample_barcode
LEFT JOIN biospecimen.samples sa ON al.sample_barcode = sa.sample_barcode
WHERE aliquot_analyte_type = 'R'"

clin <- dbGetQuery(con,q)

mybatch <- unique(clin[,"aliquot_batch"])

res <- matrix(0,nrow=length(mygene),ncol=length(mybatch))
pval <- matrix(0,nrow=length(mygene),ncol=length(mybatch))
rownames(res) <- rownames(pval) <- mygene
colnames(res) <- colnames(pval) <- mybatch
for(i in 1:length(mybatch))
{
	myaliquots <- clin[which(clin[,"aliquot_batch"]==mybatch[i]),"aliquot_barcode"]
	sub_gene <- data[mygene,myaliquots]
	out_gene <- data[mygene,-which(colnames(data)%in%myaliquots)]
	for(j in 1:nrow(sub_gene))
	{
		pval[j,i] <- wilcox.test(sub_gene[j,],out_gene[j,])$p.value
	}
	res[,i] <- apply(sub_gene,1,mean)
	
}

apply(res,2,mean)

rrna_meta <- apply(data[mygene,],2,mean)
aliquot_batch <- gsub("-RNA","",clin[match(names(rrna_meta),clin[,"aliquot_barcode"]),"aliquot_batch"])
plot_res <- data.frame(rrna_meta,aliquot_batch)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/rrna_by_batch.pdf",width=7,height=4)
p1 <- ggplot(plot_res,aes(y = rrna_meta, x = aliquot_batch)) +
	geom_boxplot() +
	labs(x = "Aliquot batch", y = "rRNA TPM") +
	theme_bw() +
	theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
	axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.position="none") +
	coord_cartesian(ylim=c(0,200))
p1
dev.off()