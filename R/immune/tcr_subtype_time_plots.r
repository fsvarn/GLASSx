library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.shannon AS tcr_shannon_a, 
im2.shannon AS tcr_shannon_b, 
im1.evenness AS tcr_evenness_a,
im2.evenness AS tcr_evenness_b,
im1.richness AS tcr_richness_a,
im2.richness AS tcr_richness_b,
im1.total_tcr AS total_tcr_a,
im2.total_tcr AS total_tcr_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b,
cs1.idh_status AS idh_status
FROM analysis.rna_silver_set ps
JOIN analysis.tcr_stats im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.tcr_stats im2 ON im2.aliquot_barcode = ps.tumor_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
ORDER BY 1, 2, 6"

dat <- dbGetQuery(con, q)

wilcox.test(dat[,"total_tcr_a"], dat[,"total_tcr_b"],paired=TRUE)
wilcox.test(dat[,"tcr_evenness_a"], dat[,"tcr_evenness_b"],paired=TRUE)
wilcox.test(dat[,"tcr_richness_a"], dat[,"tcr_richness_b"],paired=TRUE)
wilcox.test(dat[,"tcr_shannon_a"], dat[,"tcr_shannon_b"],paired=TRUE)

time_subtype <- function(x, var_col, type_col)
{
	subtype <- unique(x[,type_col])
	p.val <- eff <- rep(0, length(unique(subtype)))
	
	for(i in 1:length(subtype))
	{
		sub_x <- x[which(x[,type_col] == subtype[i]),]
		
		var_col_a <- paste(var_col,"_a",sep="")
		var_col_b <- paste(var_col,"_b",sep="")
		
		p.val[i] <- wilcox.test(sub_x[,var_col_a], sub_x[,var_col_b], paired=TRUE)$p.value
		eff[i] <- median(sub_x[,var_col_b] - median(sub_x[,var_col_a]))
	}
	
	res <- data.frame(subtype,p.val,eff)
	return(res)
}

tcr_dif <- time_subtype(dat,"total_tcr","subtype_a")
richness_dif <- time_subtype(dat,"tcr_richness","subtype_a")
evenness_dif <- time_subtype(dat,"tcr_evenness","subtype_a")
shannon_dif <- time_subtype(dat,"tcr_shannon","subtype_a")

aliquot_barcode <- c(dat[,"tumor_barcode_a"],dat[,"tumor_barcode_b"])
case_barcode <- rep(dat[,"case_barcode"],2)
sample_type <- c(rep("Initial",nrow(dat)),rep("Recurrent",nrow(dat)))
total_tcr <- c(dat[,"total_tcr_a"],dat[,"total_tcr_b"])
richness <- c(dat[,"tcr_richness_a"],dat[,"tcr_richness_b"])
shannon <- c(dat[,"tcr_shannon_a"],dat[,"tcr_shannon_b"])
evenness <- c(dat[,"tcr_evenness_a"],dat[,"tcr_evenness_b"])
subtype <- c(dat[,"subtype_a"],dat[,"subtype_a"])

plot_dat <- data.frame(case_barcode,aliquot_barcode,sample_type,total_tcr,richness,shannon,evenness,subtype)


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/tcr_boxplots_initial_recurrent.pdf",width=6,height=2)
p1 <- ggplot(plot_dat,aes(x = sample_type, y = richness)) +
geom_boxplot(width=0.5,lwd=0.35,outlier.shape=NA) +
geom_line(aes(group=case_barcode,colour=subtype), size=0.45,alpha=0.4) +
geom_point(size=0.8,colour="black") +
labs(y = "TCR richness") +
facet_grid(.~subtype) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(0,250))

p2 <- ggplot(plot_dat,aes(x = sample_type, y = evenness)) +
geom_boxplot(width=0.5,lwd=0.35,outlier.shape=NA) +
geom_line(aes(group=case_barcode,colour=subtype), size=0.45,alpha=0.4) +
geom_point(size=0.8,colour="black") +
labs(y = "TCR evenness") +
facet_grid(.~subtype) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
scale_y_continuous(breaks = seq(0, 1.1, by = 25)) +
coord_cartesian(ylim=c(0,1.1))

grid.arrange(p1,p2,ncol=2)
dev.off()

