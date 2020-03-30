library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT ds.*, ic.change
FROM analysis.driver_status_snv ds
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = ds.tumor_barcode_a AND gs.tumor_barcode_b = ds.tumor_barcode_b
JOIN analysis.immune_cluster_change ic ON ic.case_barcode = ds.case_barcode
ORDER BY change
"
dat <- dbGetQuery(con,q)

driver_genes <- dbReadTable(con, Id(schema="ref",table="driver_genes"))
driver_genes <- as.character(driver_genes[,"gene_symbol"])

res <- matrix(0,nrow=length(driver_genes),ncol=4)
rownames(res) <- driver_genes
colnames(res) <- c("inc_gain","inc_loss","dec_gain","dec_loss")
for(i in 1:length(driver_genes))
{
	#Increase association with gain
	sub1 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_a"]),]
	g1 <- nrow(sub1[which(sub1[,"change"]=="Increase"),])
	g2 <- nrow(sub1[which(sub1[,"change"]!="Increase"),])
	sub2 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_a"],invert=TRUE),]
	g3 <- nrow(sub2[which(sub2[,"change"]=="Increase"),])
	g4 <- nrow(sub2[which(sub2[,"change"]!="Increase"),])
	
	res[i,"inc_gain"] <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value
	
	#Increase association with loss
	sub1 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_b"]),]
	g1 <- nrow(sub1[which(sub1[,"change"]=="Increase"),])
	g2 <- nrow(sub1[which(sub1[,"change"]!="Increase"),])
	sub2 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_b"],invert=TRUE),]
	g3 <- nrow(sub2[which(sub2[,"change"]=="Increase"),])
	g4 <- nrow(sub2[which(sub2[,"change"]!="Increase"),])
	
	res[i,"inc_loss"] <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value

	#Decrease association with gain
	sub1 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_a"]),]
	g1 <- nrow(sub1[which(sub1[,"change"]=="Decrease"),])
	g2 <- nrow(sub1[which(sub1[,"change"]!="Decrease"),])
	sub2 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_a"],invert=TRUE),]
	g3 <- nrow(sub2[which(sub2[,"change"]=="Decrease"),])
	g4 <- nrow(sub2[which(sub2[,"change"]!="Decrease"),])
	
	res[i,"dec_gain"] <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value
	
	#Decrease association with loss
	sub1 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_b"]),]
	g1 <- nrow(sub1[which(sub1[,"change"]=="Decrease"),])
	g2 <- nrow(sub1[which(sub1[,"change"]!="Decrease"),])
	sub2 <- dat[grep(driver_genes[i],dat[,"snv_driver_change_b"],invert=TRUE),]
	g3 <- nrow(sub2[which(sub2[,"change"]=="Decrease"),])
	g4 <- nrow(sub2[which(sub2[,"change"]!="Decrease"),])
	
	res[i,"dec_loss"] <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value
}

#-----------------------------

rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT ds.*, im1.signature_name, im1.enrichment_score AS es_a, im2.enrichment_score AS es_b
FROM analysis.driver_status_snv ds
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = ds.tumor_barcode_a AND gs.tumor_barcode_b = ds.tumor_barcode_b
JOIN analysis.rna_dna_pairs rd ON rd.dna_barcode_a = ds.tumor_barcode_a AND rd.dna_barcode_b = ds.tumor_barcode_b
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
"
dat <- dbGetQuery(con,q)

driver_genes <- dbReadTable(con, Id(schema="ref",table="driver_genes"))
driver_genes <- as.character(driver_genes[,"gene_symbol"])

cells <- unique(dat[,"signature_name"])

gene_res <- list()
sig <- rep(0,length(driver_genes))
names(sig) <- driver_genes
for(i in 1:length(driver_genes))
{
	res <- matrix(0,nrow=length(cells),ncol=4)
	rownames(res) <- cells
	colnames(res) <- c("gain_pval","gain_eff","loss_pval","loss_eff")

	sub_gain <- dat[grep(driver_genes[i],dat[,"snv_driver_change_a"]),]
	sub_loss <- dat[grep(driver_genes[i],dat[,"snv_driver_change_b"]),]
	for(j in 1:length(cells))
	{
		gain_eff <- gain_pval <- loss_eff <- loss_pval <- NA
		sub1 <- sub_gain[which(sub_gain[,"signature_name"]==cells[j]),]
		sub2 <- sub_loss[which(sub_loss[,"signature_name"]==cells[j]),]
		
		if(nrow(sub1) > 0)
		{
			gain_pval <- wilcox.test(sub1[,"es_a"],sub1[,"es_b"],paired=TRUE)$p.value
			gain_eff <- median(sub1[,"es_b"]) - median(sub1[,"es_a"])
		}
		if(nrow(sub2) > 0)
		{			
			loss_pval <- wilcox.test(sub2[,"es_a"],sub2[,"es_b"],paired=TRUE)$p.value
			loss_eff <- median(sub2[,"es_b"]) - median(sub2[,"es_a"])
		}
		res[j,] <- c(gain_pval, gain_eff, loss_pval, loss_eff)
	}
	gene_res[[i]] <- res
	sig[i] <- as.numeric(sum(res[,"gain_pval"] < 0.05) > 0  | sum(res[,"loss_pval"] < 0.05) > 0)
}
names(gene_res) <- driver_genes