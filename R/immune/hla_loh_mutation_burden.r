library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(ggbeeswarm)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "WITH collapse AS
(
SELECT *,
	CASE WHEN  (hla_type1_copy_number < 0.5 OR hla_type2_copy_number < 0.5) THEN 'loss' ELSE 'none' END AS loss_status
	FROM variants.lohhla
	WHERE coverage_filter=20
),
agg AS
(
	SELECT pair_barcode, string_agg(loss_status,'; ') AS any_loss
	FROM collapse
	GROUP BY pair_barcode
)
SELECT ag.*,  mf.coverage_adj_mut_freq, cs.idh_codel_subtype
FROM agg ag
JOIN analysis.pairs pa ON ag.pair_barcode = pa.pair_barcode
JOIN analysis.mut_freq mf ON mf.aliquot_barcode = pa.tumor_barcode
JOIN analysis.lohhla_set gs ON mf.aliquot_barcode = gs.tumor_barcode_a OR mf.aliquot_barcode = gs.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = gs.case_barcode
ORDER BY 3 DESC
"

dat <- dbGetQuery(con, q)

subtypes <- unique(dat[,"idh_codel_subtype"])

p.val <- eff <- rep(0,length(subtypes))
names(p.val) <- subtypes
for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"idh_codel_subtype"]==subtypes[i]),]
	
	g1 <- sub_dat[grep("loss",sub_dat[,"any_loss"]),"coverage_adj_mut_freq"]
	g2 <- sub_dat[grep("loss",sub_dat[,"any_loss"],invert=TRUE),"coverage_adj_mut_freq"]

	eff[i] <- median(g1) - median(g2)
	p.val[i] <- wilcox.test(g1,g2)$p.value
}

data.frame(eff,p.val)


#Neoantigen and mutation count
#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "WITH collapse AS
(
	SELECT *,
	CASE WHEN pval < 0.1 AND (hla_type1_copy_number < 0.5 OR hla_type2_copy_number < 0.5) THEN 'loss' ELSE 'none' END AS loss_status
	FROM variants.lohhla
	WHERE coverage_filter=20
),
agg AS
(
	SELECT pair_barcode, string_agg(loss_status,'; ') AS any_loss
	FROM collapse
	GROUP BY pair_barcode
)
SELECT ag.*,  mf.mt_count, mf.neoag_count, cs.idh_codel_subtype
FROM agg ag
JOIN analysis.pairs pa ON ag.pair_barcode = pa.pair_barcode
JOIN analysis.neoag_freq mf ON mf.aliquot_barcode = pa.tumor_barcode
JOIN analysis.lohhla_set gs ON mf.aliquot_barcode = gs.tumor_barcode_a OR mf.aliquot_barcode = gs.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = gs.case_barcode
ORDER BY 3 DESC
"

dat <- dbGetQuery(con, q)

subtypes <- unique(dat[,"idh_codel_subtype"])

#Exonic mutation count
p.val <- eff <- rep(0,length(subtypes))
names(p.val) <- subtypes
for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"idh_codel_subtype"]==subtypes[i]),]
	
	g1 <- sub_dat[grep("loss",sub_dat[,"any_loss"]),"mt_count"]
	g2 <- sub_dat[grep("loss",sub_dat[,"any_loss"],invert=TRUE),"mt_count"]

	eff[i] <- median(g1) - median(g2)
	p.val[i] <- wilcox.test(g1,g2)$p.value
}

data.frame(eff,p.val)

#                  eff     p.val
# IDHwt           -4.5 0.1433124
# IDHmut-noncodel  9.0 0.0252692
# IDHmut-codel    -5.5 0.5827564

#Neoantigen count
p.val <- eff <- rep(0,length(subtypes))
names(p.val) <- subtypes
for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"idh_codel_subtype"]==subtypes[i]),]
	
	g1 <- sub_dat[grep("loss",sub_dat[,"any_loss"]),"neoag_count"]
	g2 <- sub_dat[grep("loss",sub_dat[,"any_loss"],invert=TRUE),"neoag_count"]

	eff[i] <- median(g1) - median(g2)
	p.val[i] <- wilcox.test(g1,g2)$p.value
}

data.frame(eff,p.val)

#                  eff      p.val
# IDHwt           -3.0 0.31761287
# IDHmut-noncodel  7.0 0.01523149
# IDHmut-codel    -2.5 0.93080872