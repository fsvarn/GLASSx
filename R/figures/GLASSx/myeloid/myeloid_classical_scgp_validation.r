###################################################
# Create boxplot examining how myeloid mesenchymal signature is distributed over time in GLASS
# Updated: 2021.01.14
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE cs.idh_codel_subtype = 'IDHwt' AND ts1.signature_name !='Mesenchymal' AND ts2.signature_name = 'Mesenchymal'"

dat <- dbGetQuery(con,q)


# Check how TCGA signature changes:
tcga_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v3.txt")
tcga_sig <- rownames(tcga_dat %>% filter(q.val < 0.05 & eff > log2(1.1)))

geps <- read.delim(myinf1[2], row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))

tcga_sig <- intersect(rownames(geps), tcga_sig)


g1 <- geps[tcga_sig,dat[,"tumor_barcode_a"]]
g2 <- geps[tcga_sig,dat[,"tumor_barcode_b"]]

g1 <- apply(g1,2,mean)
g2 <- apply(g2,2,mean)

wilcox.test(g1,g2,paired=TRUE)

# Boxplot of this result
aliquot_barcode <- c(names(g1), names(g2))
case_barcode <- sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[1:3],collapse="-"))
mes_score <- c(g1, g2)
timepoint <- c(rep("Initial",length(g1)), rep("Recurrent",length(g2)))
signature_name <- c(dat[,"subtype_a"], dat[,"subtype_b"])
mes_score <- data.frame(aliquot_barcode, case_barcode, mes_score, timepoint, signature_name)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_glass_mesenchymal_class_score.pdf",width=1.45,height=2.7)
ggplot(mes_score, aes(x=timepoint, y=mes_score)) + 
geom_boxplot(outlier.shape = NA)  +
geom_line(size=0.5,alpha=0.4, aes(group= case_barcode)) +
geom_point(size=1,aes(colour=signature_name)) +
theme_bw() +
labs(y = "Mesenchymal myeloid score") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_colour_manual(values=c("#008A22", "#8A0000","#00458A"))
dev.off()

diff_list <- cor_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	#geps <- log10(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]
	colnames(geps) <- gsub("\\.","-",colnames(geps))


	g1 <- geps[,dat[,"tumor_barcode_a"]]
	g2 <- geps[,dat[,"tumor_barcode_b"]]

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(g1[j,])
		group2 <- as.numeric(g2[j,])
	
		p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
		eff[j] <- log2(median(group2/group1))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	idhwt_res <- data.frame(p.val, q.val, eff)
	idhwt_res <- idhwt_res[order(p.val),]

	idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
	idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1
		
	mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])

	q <- 'SELECT es.*, gt.gene_symbol, gt.est_counts, gt.tpm
	FROM analysis.estimate es
	JOIN analysis.gene_tpm gt ON es.aliquot_barcode = gt.aliquot_barcode
	WHERE gene_symbol IN (mygenes)'

	mygenes_text <- paste("'", paste(mygenes, collapse="','"), "'", sep="")
	q <- gsub("mygenes", mygenes_text, q)

	pur_dat <- dbGetQuery(con,q)

	cor_res <- pur_dat %>%
	group_by(gene_symbol) %>%
	summarise(cor = cor(purity, tpm, method="s")) %>%
	inner_join(rownames_to_column(idhwt_res), by = c("gene_symbol" = "rowname")) %>%
	data.frame
	
	diff_list[[i]] <- idhwt_res
	cor_list[[i]] <- cor_res
}

names(diff_list) <- names(cor_list)  <- mytag

# Plot a volcano plot for myeloid cells

myeloid_diff <- diff_list[[2]]

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_myeloid_volcano_nonmes_mes.pdf",width=2.5,height=3)  
ggplot(myeloid_diff, aes(x=eff, y=logp)) + 
geom_point(size=0.5,aes(colour = sig)) +
labs(x = "log2(FC)", y = "-log10(adj p-val)") +
scale_colour_manual(values = c("black", "tomato3")) +
theme_bw() +
theme(
axis.title=element_text(size=7),
axis.text.x=element_text(size=7),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none") 
dev.off()
