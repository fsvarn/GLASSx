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
WITH collapse AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, string_agg(hla_loh_change,'; ') AS hla_change
	FROM variants.hla_loh_change 
	GROUP BY case_barcode, tumor_barcode_a, tumor_barcode_b
	ORDER BY 1
)
SELECT co.*, sp1.cellularity AS cellularity_a, sp2.cellularity AS cellularity_b, sp2.cellularity - sp1.cellularity AS cellularity_change
FROM collapse co
JOIN analysis.pairs tp1 ON tp1.tumor_barcode = co.tumor_barcode_a
JOIN analysis.pairs tp2 ON tp2.tumor_barcode = co.tumor_barcode_b
JOIN variants.seqz_params sp1 ON sp1.pair_barcode = tp1.pair_barcode
JOIN variants.seqz_params sp2 ON sp2.pair_barcode = tp2.pair_barcode
"

dat <- dbGetQuery(con,q)

loh_dat <- dat[grep("gain",dat[,"hla_change"]),]
wilcox.test(loh_dat[,"cellularity_a"], loh_dat[,"cellularity_b"],paired=TRUE)

noloh_dat <- dat[grep("loss",dat[,"hla_change"]),]
wilcox.test(noloh_dat[,"cellularity_a"], noloh_dat[,"cellularity_b"],paired=TRUE)


g1 <- dat[grep("gain",dat[,"hla_change"]),"cellularity_change"]
g2 <- dat[grep("gain",dat[,"hla_change"],invert=TRUE),"cellularity_change"]

wilcox.test(g1,g2)		#P = 0.07


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
WITH collapse AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, string_agg(hla_loh_change,'; ') AS hla_change
	FROM variants.hla_loh_change 
	GROUP BY case_barcode, tumor_barcode_a, tumor_barcode_b
	ORDER BY 1
)
SELECT co.*, tp1.purity AS purity_a, tp2.purity AS purity_b, tp2.purity - tp1.purity AS purity_change
FROM collapse co
JOIN analysis.pairs pa1 ON pa1.tumor_barcode = co.tumor_barcode_a
JOIN analysis.pairs pa2 ON pa2.tumor_barcode = co.tumor_barcode_b
JOIN variants.titan_params tp1 ON tp1.pair_barcode = pa1.pair_barcode
JOIN variants.titan_params tp2 ON tp2.pair_barcode = pa2.pair_barcode
"

dat <- dbGetQuery(con,q)

loh_dat <- dat[grep("gain",dat[,"hla_change"]),]
wilcox.test(loh_dat[,"purity_a"], loh_dat[,"purity_b"],paired=TRUE)		#P = 3e-3

noloh_dat <- dat[grep("loss",dat[,"hla_change"]),]
wilcox.test(noloh_dat[,"purity_a"], noloh_dat[,"purity_b"],paired=TRUE)


g1 <- dat[grep("gain",dat[,"hla_change"]),"purity_change"]
g2 <- dat[grep("gain",dat[,"hla_change"],invert=TRUE),"purity_change"]

wilcox.test(g1,g2)		#P = 5e-4