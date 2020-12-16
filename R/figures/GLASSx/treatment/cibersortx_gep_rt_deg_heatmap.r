###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
# FROM analysis.rna_silver_set ps
# JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
# JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
# JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
# WHERE cs.idh_codel_subtype = 'IDHwt' AND ts1.signature_name ='Mesenchymal' AND ts2.signature_name = 'Classical'"

q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE cs.idh_codel_subtype = 'IDHwt' AND ts1.signature_name = ts2.signature_name AND received_rt"

dat <- dbGetQuery(con,q)

p <- se <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	geps <- log10(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]	
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	nrow(geps)

	g1 <- geps[,dat[,"tumor_barcode_a"]]
	g2 <- geps[,dat[,"tumor_barcode_b"]]

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(g1[j,])
		group2 <- as.numeric(g2[j,])
	
		p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
		eff[j] <- log2(mean(group2)/mean(group1))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	idhwt_res <- data.frame(p.val, q.val, eff)
	idhwt_res <- idhwt_res[order(eff),]

	idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
	idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1
	
	# Plot heatmaps
	mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
	sig_matrix <- geps[mygenes,c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"])]
	
	sig_matrix <- t(apply(sig_matrix, 1, function(x) x - median(x)))
	
	grp1 <- sig_matrix[,dat[,"tumor_barcode_a"]]
	grp1 <- apply(grp1, 1, mean)
	grp2 <- sig_matrix[,dat[,"tumor_barcode_b"]]
	grp2 <- apply(grp2, 1, mean)
	
	expr <- c(grp1, grp2)
	# Rescale for color
	expr[which(expr > quantile(expr,.99))] <- quantile(expr,.99)
	gene_symbol <- c(names(grp1), names(grp2))
	timepoint <- c(rep("Initial", length(grp1)), rep("Recurrent", length(grp2)))
	
	plot_hm <- data.frame(gene_symbol, expr, timepoint)
	lev <- names(grp1)[order(grp2 - grp1)]
	plot_hm$gene_symbol <- factor(plot_hm$gene_symbol, levels = lev)
	plot_hm[,"cell"] <- cell_state[i]
	
	p[[i]] <- plot_hm
	
	se[[i]] <- ggplot(data = plot_hm, aes(x = timepoint, y = gene_symbol)) +
	geom_tile(aes(fill=expr)) +
	scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
	space = "Lab", na.value = "grey50", guide = "colourbar",
	aesthetics = "fill") + 
	theme_void() +
	theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	legend.position = "none")	
}

plot_res <- do.call(rbind, p)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/radiation_pre_post_csx_heatmaps.pdf",width=7, height =3)
grid.arrange(se[[1]],se[[2]],se[[3]],se[[4]],nrow=1)
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/radiation_pre_post_csx_legends.pdf",width=7, height =3)
ggplot(data = plot_hm, aes(x = timepoint, y = gene_symbol)) +
	geom_tile(aes(fill=expr)) +
	scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
	space = "Lab", na.value = "grey50", guide = "colourbar",
	aesthetics = "fill") + 
	theme_void() +
	theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank())	
dev.off()
