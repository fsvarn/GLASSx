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
WITH hla_loss AS
(
	SELECT
	lh.pair_barcode,
	lh.hla_type1,
	lh.hla_type2,
	lh.pval,
	CASE
		WHEN lh.pval > 0.05::double precision THEN 'FALSE'
		ELSE 'TRUE'
		END AS loss
	FROM variants.lohhla lh
	WHERE coverage_filter = 20 AND pval IS NOT NULL
),
any_loss AS
(
	SELECT pair_barcode, string_agg(loss,'; ') AS loss
	FROM hla_loss
	--WHERE pair_barcode LIKE '%-WXS'
	GROUP BY pair_barcode
)
SELECT al.pair_barcode, 
CASE WHEN loss LIKE '%TRUE%' THEN TRUE ELSE FALSE END AS any_loss,
im.signature_name,
im.enrichment_score,
idh_codel_subtype
FROM any_loss al
JOIN analysis.pairs pa ON pa.pair_barcode = al.pair_barcode
JOIN analysis.analyte_sets rd ON rd.dna_barcode = pa.tumor_barcode
JOIN analysis.davoli_immune_score im ON im.aliquot_barcode = rd.rna_barcode
JOIN analysis.platinum_set ps ON rd.dna_barcode = ps.dna_barcode_a --OR rd.dna_barcode = ps.dna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = rd.case_barcode
ORDER BY 1,3
"

dat <- dbGetQuery(con,q)
dat[,2] <- as.logical(as.numeric(dat[,2]))

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

res <- matrix(0,nrow=length(cells),ncol=length(subtypes))
rownames(res) <- cells
colnames(res) <- subtypes

for(j in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"idh_codel_subtype"]==subtypes[j]),]
	
	p.val <- eff <- rep(0,length(cells))
	names(p.val) <- names(eff) <- cells
	for(i in 1:length(cells))
	{
		c_dat <- sub_dat[which(sub_dat[,"signature_name"]==cells[i]),]
		g1 <- c_dat[which(c_dat[,"any_loss"]),"enrichment_score"]
		g2 <- c_dat[which(!c_dat[,"any_loss"]),"enrichment_score"]
	
		p.val[i] <- wilcox.test(g1,g2)$p.value
		eff[i] <- median(g1) - median(g2)
	}

	res[,j] <- p.val
}