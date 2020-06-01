###################################################
# Test how the delta log2 copy ratio in each gene associates with changes in ESTIMATE immune score
# Updated: 2020.05.06
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
SELECT ps.*, an.gene_symbol, an.chrom, pos, an.delta_log2_copy_ratio, cnv_call,
im2.immune_score - im1.immune_score AS immune_dif,
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_cnv_by_gene_diff_call an ON an.tumor_barcode_a = ps.dna_barcode_a AND an.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.estimate im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate im2 ON im2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY chrom, pos, case_barcode
"
dat <- dbGetQuery(con,q)

gene_symbol <- unique(dat[,"gene_symbol"])
subtypes <- unique(dat[,"idh_codel_subtype"])

##################################################
# Step 1: Correlate immune scores/purity score differences with % CNA difference
##################################################

#Get the number of cases for each subtype
max_subtype <- dat[,c("case_barcode","idh_codel_subtype")]
max_subtype <- max_subtype[-which(duplicated(max_subtype)),]
max_subtype <- table(max_subtype[,"idh_codel_subtype"])
max_subtype <- max_subtype[subtypes]

mycor <- cor_pval <- matrix(NA,nrow=length(gene_symbol),ncol=length(subtypes))
rownames(mycor) <- rownames(cor_pval) <- gene_symbol
colnames(mycor) <- colnames(cor_pval) <- subtypes

mean_change <- sd_change <- matrix(NA,nrow=length(gene_symbol),ncol=length(subtypes))
rownames(mean_change) <- rownames(sd_change) <- gene_symbol
colnames(mean_change) <- colnames(sd_change) <- subtypes

# del_eff <- del_pval <- inc_eff <- inc_pval <- matrix(NA,nrow=length(gene_symbol),ncol=length(subtypes))
# rownames(del_eff) <- rownames(del_pval) <- rownames(inc_eff) <- rownames(inc_pval) <- gene_symbol
# colnames(del_eff) <- colnames(del_pval) <- colnames(inc_eff) <- colnames(inc_pval) <- subtypes

for(i in 1:length(gene_symbol))
{
	cat("\r",i)
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"gene_symbol"] == gene_symbol[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		if(nrow(sub_dat) < 0.7 * max_subtype[subtypes[j]] ){
			next}
		mycor[i,j] <- cor(sub_dat[,"immune_dif"],sub_dat[,"delta_log2_copy_ratio"],method="s")
		cor_pval[i,j] <- cor.test(sub_dat[,"immune_dif"],sub_dat[,"delta_log2_copy_ratio"],method="s")$p.value
		
		mean_change[i,j] <- mean(sub_dat[,"delta_log2_copy_ratio"])
		sd_change[i,j] <- sd(sub_dat[,"delta_log2_copy_ratio"])

# 		if(length(which(sub_dat[,"cnv_call"]==-1)) > 0 & length(which(sub_dat[,"cnv_call"] == 0)) > 0)
# 		{		
# 			del_eff[i,j] <- median(sub_dat[which(sub_dat[,"cnv_call"]==-1),"immune_dif"]) - median(sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])
# 			del_pval[i,j] <- wilcox.test(sub_dat[which(sub_dat[,"cnv_call"]==-1),"immune_dif"], sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])$p.value
# 		}		
# 		if(length(which(sub_dat[,"cnv_call"]==1)) > 0 & length(which(sub_dat[,"cnv_call"] == 0)) > 0)
# 		{
# 			inc_eff[i,j] <- median(sub_dat[which(sub_dat[,"cnv_call"]==1),"immune_dif"]) - median(sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])
# 			inc_pval[i,j] <- wilcox.test(sub_dat[which(sub_dat[,"cnv_call"]==1),"immune_dif"], sub_dat[which(sub_dat[,"cnv_call"]==0),"immune_dif"])$p.value	
# 		}
	}
}

#na_ind <-  which(apply(mycor,1,function(x)sum(is.na(x)))==3)
#mycor <- mycor[-na_ind,]
#cor_pval <- cor_pval[-na_ind,]
cor_pval_adj <- apply(cor_pval,2,function(x)p.adjust(x,"BH"))

#na_ind <-  which(apply(del_pval,1,function(x)sum(is.na(x)))==3)
#del_eff <- del_eff[-na_ind,]
#del_pval <- del_pval[-na_ind,]
#del_pval_adj <- apply(del_pval,2,function(x)p.adjust(x,"BH"))

#na_ind <-  which(apply(inc_pval,1,function(x)sum(is.na(x)))==3)
#inc_eff <- inc_eff[-na_ind,]
#inc_pval <- inc_pval[-na_ind,]
#inc_pval_adj <- apply(inc_pval,2,function(x)p.adjust(x,"BH"))


mycor[which(cor_pval_adj[,2] < 0.05),]
mycor[which(cor_pval[,2] < 0.01),]
mycor[order(mycor[,2],decreasing=TRUE),]

#inc_eff[which(inc_pval_adj[,2] < 0.1),]
#del_eff[which(del_pval_adj[,2] < 0.1),]

#inc_eff[which(inc_pval[,2] < 0.01),]
#del_eff[which(del_pval[,2] < 0.01),]


##################################################
# Step 2: Plot the differences for IDHwt only (only tumor with significance)
##################################################

# Create position map 
pos_map <- dat[,c("gene_symbol","chrom","pos")]
pos_map <- aggregate(pos_map,by=list(pos_map[,"gene_symbol"]),function(x)x[1])
pos_map[,"pos_start"] <- as.numeric(gsub("\\[","",sapply(strsplit(pos_map[,"pos"],","),function(x)x[1])))
pos_map[,"pos_end"] <- as.numeric(gsub(")","",sapply(strsplit(pos_map[,"pos"],","),function(x)x[2])))
pos_map <- pos_map[order(pos_map[,"chrom"],pos_map[,"pos_start"],pos_map[,"pos_end"]),]

#For genes with same start site, add 1 to gene with later end site
#Doing it twice for HOXC genes which are clustered together
pos_map[which(duplicated(pos_map[,c("chrom","pos_start")])),"pos_start"] <- pos_map[which(duplicated(pos_map[,c("chrom","pos_start")])),"pos_start"] + 1
pos_map[which(duplicated(pos_map[,c("chrom","pos_start")])),"pos_start"] <- pos_map[which(duplicated(pos_map[,c("chrom","pos_start")])),"pos_start"] + 1

# Create plotting data for correlation barplot
rho <- mycor[pos_map[,"gene_symbol"],2]
pval <- cor_pval[,2]
adj_pval <- cor_pval_adj[,2]
sig <- as.numeric(adj_pval < 0.1)
sig[which(sig==0)] <- NA
sig[which(sig==1)] <- rho[which(sig==1)]
eff <- as.factor(as.numeric(rho >= 0))
plot_cor <- data.frame(pos_map, rho, pval, adj_pval, sig, eff)
plot_cor <- plot_cor[,-1]
plot_cor <- plot_cor[-which(is.na(plot_cor[,"rho"])),]

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cnv_gene_pos_cor_bar.pdf",width=2,height=5.5)
ggplot(plot_cor, aes(x=pos_start,fill=eff)) +
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
coord_flip(ylim=c(-0.5,0.5))
dev.off()

plot_cor[order(plot_cor[,"rho"]),c("gene_symbol","chrom","pos","rho")][1:10,]
plot_cor[order(plot_cor[,"rho"], decreasing=TRUE),c("gene_symbol","chrom","pos","rho")][1:10,]

sig_results <- plot_cor[which(plot_cor[,"adj_pval"] < 0.1),]
sig_results <- sig_results[order(sig_results[,"rho"]),c("gene_symbol","chrom","pos","rho","pval","adj_pval"),]
sig_results[,"mean_change"] <- mean_change[sig_results[,"gene_symbol"],"IDHwt"]

myoutf <- "/projects/verhaak-lab/GLASS-III/data/res/cnv_gene_immune_sig_results.txt"
write.table(sig_results, myoutf, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


##################################################
# Step 3: Plot the mean change in copy number for IDHwt only (only tumor with significance)
##################################################

# Not a very useful plot: Most copy number changes averages are right around 0

# Uses the position map in step 2
avg <- 2^mean_change[pos_map[,"gene_symbol"],2]
err <- sd_change[,2]
lower <- avg - err
upper <- avg + err
plot_avg <- data.frame(pos_map, avg, err, lower, upper)
plot_avg <- plot_avg[,-1]
plot_avg <- plot_avg[-which(is.na(plot_avg[,"avg"])),]

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cnv_gene_pos_avg_line.pdf",width=2,height=5.5)
ggplot(plot_avg, aes(x=pos_start)) +
geom_line(stat="identity",aes(y=avg)) +
geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.1) +
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
geom_hline(yintercept=0) +
#coord_flip()
coord_flip(ylim=c(0.5,1.5))
dev.off()

##################################################
# Step 4: Specific examples
##################################################

idhwt_dat <- dat[which(dat[,"idh_codel_subtype"] == "IDHwt"),]

#PDGFRA 
sub_gene1 <- idhwt_dat[which(idhwt_dat[,"gene_symbol"] == "PDGFRA"),]
g1 <- sub_gene1[which(sub_gene1[,"delta_log2_copy_ratio"] > 0), "immune_dif"]
g2 <- sub_gene1[which(sub_gene1[,"delta_log2_copy_ratio"] < 0), "immune_dif"]
wilcox.test(g1,g2)		#8e-5

#Top chromosome 13 hit: SLITRK5
sub_gene2 <- idhwt_dat[which(idhwt_dat[,"gene_symbol"] == "SLITRK5"),]
g1 <- sub_gene2[which(sub_gene2[,"delta_log2_copy_ratio"] > 0), "immune_dif"]
g2 <- sub_gene2[which(sub_gene2[,"delta_log2_copy_ratio"] < 0), "immune_dif"]
wilcox.test(g1,g2)		#8e-3



