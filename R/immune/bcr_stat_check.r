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
im1.shannon AS bcr_shannon_a, 
im2.shannon AS bcr_shannon_b, 
im1.evenness AS bcr_evenness_a,
im2.evenness AS bcr_evenness_b,
im1.richness AS bcr_richness_a,
im2.richness AS bcr_richness_b,
im1.total_bcr AS total_bcr_a,
im2.total_bcr AS total_bcr_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b,
cs1.idh_status AS idh_status
FROM analysis.rna_silver_set ps
JOIN analysis.bcr_stats im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.bcr_stats im2 ON im2.aliquot_barcode = ps.tumor_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
WHERE cs1.idh_codel_subtype IS NOT NULL
ORDER BY 1, 2, 6"

dat <- dbGetQuery(con, q)

wilcox.test(dat[,"total_bcr_a"], dat[,"total_bcr_b"],paired=TRUE)
wilcox.test(dat[,"bcr_evenness_a"], dat[,"bcr_evenness_b"],paired=TRUE)
wilcox.test(dat[,"bcr_richness_a"], dat[,"bcr_richness_b"],paired=TRUE)
wilcox.test(dat[,"bcr_shannon_a"], dat[,"bcr_shannon_b"],paired=TRUE)

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

bcr_dif <- time_subtype(dat,"total_bcr","subtype_a")
richness_dif <- time_subtype(dat,"bcr_richness","subtype_a")
evenness_dif <- time_subtype(dat,"bcr_evenness","subtype_a")
shannon_dif <- time_subtype(dat,"bcr_shannon","subtype_a")