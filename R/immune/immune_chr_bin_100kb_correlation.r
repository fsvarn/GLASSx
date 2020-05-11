###################################################
# Test how the delta log2 copy ratio in each 1 Mb genome bin associates with changes in ESTIMATE immune score
# Updated: 2020.04.30
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(reshape)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT ps.*, an.chrom, an.bin, CONCAT(an.chrom, '_', an.bin) AS bin_id, an.delta_log2_copy_ratio, cnv_call,
im2.immune_score - im1.immune_score AS immune_dif,
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_bin_100kb_diff_call an ON an.tumor_barcode_a = ps.dna_barcode_a AND an.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.estimate im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate im2 ON im2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY chrom, bin, case_barcode
"
dat <- dbGetQuery(con,q)

chr_bins <- unique(dat[,"bin_id"])
subtypes <- unique(dat[,"idh_codel_subtype"])

##################################################
# Step 1: Correlate immune scores/purity score differences with % CNA difference
##################################################

#Get the number of cases for each subtype
max_subtype <- dat[,c("case_barcode","idh_codel_subtype")]
max_subtype <- max_subtype[-which(duplicated(max_subtype)),]
max_subtype <- table(max_subtype[,"idh_codel_subtype"])
max_subtype <- max_subtype[subtypes]

mycor <- cor_pval <- matrix(NA,nrow=length(chr_bins),ncol=length(subtypes))
rownames(mycor) <- rownames(cor_pval) <- chr_bins
colnames(mycor) <- colnames(cor_pval) <- subtypes

del_eff <- del_pval <- inc_eff <- inc_pval <- matrix(NA,nrow=length(chr_bins),ncol=length(subtypes))
rownames(del_eff) <- rownames(del_pval) <- rownames(inc_eff) <- rownames(inc_pval) <- chr_bins
colnames(del_eff) <- colnames(del_pval) <- colnames(inc_eff) <- colnames(inc_pval) <- subtypes

for(i in 1:length(chr_bins))
{
	cat("\r",i)
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"bin_id"] == chr_bins[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		if(nrow(sub_dat) < 0.7 * max_subtype[subtypes[j]] ){
			next}
		mycor[i,j] <- cor(sub_dat[,"immune_dif"],sub_dat[,"delta_log2_copy_ratio"],method="s")
		cor_pval[i,j] <- cor.test(sub_dat[,"immune_dif"],sub_dat[,"delta_log2_copy_ratio"],method="s")$p.value
		
		if(length(which(sub_dat[,"cnv_call"]==-1)) > 0 & length(which(sub_dat[,"cnv_call"] == 0)) > 0)
		{		
			del_eff[i,j] <- median(sub_dat[which(sub_dat[,"cnv_call"]==-1),"immune_dif"]) - median(sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])
			del_pval[i,j] <- wilcox.test(sub_dat[which(sub_dat[,"cnv_call"]==-1),"immune_dif"], sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])$p.value
		}		
		if(length(which(sub_dat[,"cnv_call"]==1)) > 0 & length(which(sub_dat[,"cnv_call"] == 0)) > 0)
		{
			inc_eff[i,j] <- median(sub_dat[which(sub_dat[,"cnv_call"]==1),"immune_dif"]) - median(sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])
			inc_pval[i,j] <- wilcox.test(sub_dat[which(sub_dat[,"cnv_call"]==1),"immune_dif"], sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])$p.value	
		}
	}
}

na_ind <-  which(apply(mycor,1,function(x)sum(is.na(x)))==3)
mycor <- mycor[-na_ind,]
cor_pval <- cor_pval[-na_ind,]
cor_pval_adj <- apply(cor_pval,2,function(x)p.adjust(x,"BH"))

na_ind <-  which(apply(del_pval,1,function(x)sum(is.na(x)))==3)
del_eff <- del_eff[-na_ind,]
del_pval <- del_pval[-na_ind,]
del_pval_adj <- apply(del_pval,2,function(x)p.adjust(x,"BH"))

na_ind <-  which(apply(inc_pval,1,function(x)sum(is.na(x)))==3)
inc_eff <- inc_eff[-na_ind,]
inc_pval <- inc_pval[-na_ind,]
inc_pval_adj <- apply(inc_pval,2,function(x)p.adjust(x,"BH"))


mycor[which(cor_pval_adj[,2] < 0.05),]
mycor[which(cor_pval[,2] < 0.01),]

inc_eff[which(inc_pval_adj[,2] < 0.1),]
del_eff[which(del_pval_adj[,2] < 0.1),]
del_eff[which(del_pval[,2] < 0.01),]

# Create plotting data

##################################################
# Step 2: Plot the differences for IDHwt only (only tumor with significance)
##################################################

# Create plotting data
chrom <- sapply(strsplit(rownames(mycor),"_"),function(x)x[1])
chrom <- factor(chrom, levels = unique(chrom)[order(as.numeric(unique(chrom)))])
bin <- sapply(strsplit(rownames(mycor),"_"),function(x)x[2])
bin <- as.numeric(gsub("\\[","",sapply(strsplit(bin,","),function(x)x[1])))
rho <- mycor[,2]
pval <- cor_pval_adj[,2]
sig <- as.numeric(pval < 0.05)
sig[which(sig==0)] <- NA
sig[which(sig==1)] <- rho[which(sig==1)]
eff <- as.factor(as.numeric(rho >= 0))
plot_cor <- data.frame(chrom, bin, rho, pval,sig)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cnv_100kb_bin_cor_bar.pdf",width=2,height=5.5)
ggplot(plot_cor, aes(x=bin,fill=eff)) +
geom_bar(stat="identity",aes(y=rho)) +
geom_point(aes(y=sig),colour="black",size=0.3,stroke=0, shape=16) +
facet_grid(chrom~., scales="free", space= "free") +
scale_fill_manual(values=c("royalblue4","tomato3")) +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.title=element_blank(),
axis.ticks = element_blank(),
axis.text.x=element_text(size=7),
axis.text.y=element_blank(),
strip.text = element_text(size=5, hjust=0.5),
strip.background = element_blank(),
panel.spacing = unit(0, "lines"),
legend.position="none") + 
coord_cartesian(ylim=c(-0.5,0.5)) +
coord_flip()
dev.off()


# Chromosome 4 region: PDGFRA at 55 Mb, TLRs at 38 Mb
# Chromosome 10 region: PTEN at 89 Mb


#---------------

q <- "
SELECT ps.*, cn1.gene_symbol, cn1.hlvl_call, cn2.hlvl_call, 
im1.immune_score AS immune_score_a, im2.immune_score AS immune_score_b,
cs.idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_cnv_by_gene cn1 ON cn1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.gatk_cnv_by_gene cn2 ON cn2.aliquot_barcode = ps.dna_barcode_b AND cn2.gene_symbol = cn1.gene_symbol
JOIN analysis.estimate im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate im2 ON im2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE cn1.gene_symbol = 'RB1' AND cn1.hlvl_call = 0 AND cn2.hlvl_call < 0  AND idh_codel_subtype = 'IDHwt' ORDER BY 1
"

dat <- dbGetQuery(con,q)

wilcox.test(dat[,"immune_score_a"], dat[,"immune_score_b"], paired=TRUE) # P = 0.125, n = 5

