###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Get signature genes
q <- "SELECT * 
FROM ref.immune_signatures
WHERE signature_set = 'Muller'"

sig_genes <- dbGetQuery(con,q)

q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b, tc.grade_change
FROM analysis.rna_silver_set ps
JOIN analysis.rna_dna_pairs rd ON rd.rna_barcode_a = ps.tumor_barcode_a AND rd.rna_barcode_b = ps.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
JOIN analysis.tumor_clinical_comparison tc ON tc.tumor_pair_barcode = rd.dna_pair_barcode
WHERE cs.idh_codel_subtype LIKE 'IDHmut%'" 
# AND grade_change = 'Grade up'"

dat <- dbGetQuery(con,q)

mean_log2 <- status <- cell_state <- neftel <- c()
sig_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	
	geps <- read.delim(myinf1[i], row.names=1)
	geps <- log2(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]

	colnames(geps) <- gsub("\\.","-",colnames(geps))
	
	g1 <- geps[,dat[,"tumor_barcode_a"]]
	g2 <- geps[,dat[,"tumor_barcode_b"]]

	sigs <- unique(sig_genes[,"signature_name"])
	p.value <- eff <- rep(0, length(sigs))
	for(j in 1:length(sigs))
	{
		mygenes <- sig_genes[which(sig_genes[,"signature_name"] == sigs[j]),"gene_symbol"]
	
		group1 <- apply(g1[mygenes,],2,function(x)mean(x,na.rm=TRUE)) 
		group2 <- apply(g2[mygenes,],2,function(x)mean(x,na.rm=TRUE)) 

		med_center <- median(c(group1, group2))
		
		mean_log2 <- c(mean_log2, c(group1-med_center, group2-med_center))
		status <- c(status, c(rep("Initial", length(group1)), rep("Recurrent", length(group2))))
		cell_state <- c(cell_state, rep(mytag[i], length(c(group1, group2))))
		neftel <- c(neftel, rep(sigs[j], length(c(group1, group2))))

		p.value[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
		eff[j] <- median(group2 - group1)
	}
	sig_res <- data.frame(sigs, p.value, eff)
	sig_res[,"cell"] <- mytag[i]
	sig_list[[i]] <- sig_res
}

aliquot_barcode <- names(mean_log2)
case_barcode <- substr(aliquot_barcode,1,12)

plot_res <- data.frame(case_barcode, aliquot_barcode, mean_log2, status, cell_state, neftel)
plot_res <- plot_res %>%
			mutate(cell_state = recode(cell_state, "differentiated_tumor" = "Differentiated tumor",
					"myeloid" = "Myeloid", "prolif_stemcell_tumor" = "Prolif. stem cell tumor",
					"stemcell_tumor" = "Stem cell tumor")) %>%
			mutate(cell_state = as_factor(cell_state)) %>%
			mutate(cell_state = fct_relevel(cell_state,"Differentiated tumor", "Stem cell tumor", "Prolif. stem cell tumor", "Myeloid")) 


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_neftel_states_idhmut_grade_increase_mg_ladder.pdf",width=2.5,height=5)  
ggplot(plot_res, aes(x=status, y=mean_log2)) + 
geom_boxplot(fill= "white", colour="black",outlier.shape=NA) +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1,aes(colour = status)) +
labs(y = "Normalized log2 expression") +
facet_grid(cell_state ~ neftel, scales="free") +
scale_colour_manual(values = c("#008A22", "#8A0000")) +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none") 
#coord_cartesian(ylim=c(-1.5,1.2))
dev.off()

# Zoom in on the significant differences in the differentiated tumors
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_neftel_differentiated_classical_mes_ladder.pdf",width=2,height=2)  
ggplot(plot_res %>% 
	filter(cell_state == "Differentiated tumor") %>%
	filter(neftel == "AC-like" | neftel == "MES-like"), 
	aes(x=status, y=mean_log2)) + 
geom_boxplot(fill= "white", colour="black",outlier.shape=NA) +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1,aes(colour = status)) +
labs(y = "Normalized log2 expression") +
facet_grid(~ neftel, scales="free") +
scale_colour_manual(values = c("#008A22", "#8A0000")) +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none") 
dev.off()
