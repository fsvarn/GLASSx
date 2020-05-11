###################################################
# Test % of genome with amp or deletion associates with immune levels
# Updated: 2020.04.28
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
WITH arrange_segments_by_arm AS 
(
	SELECT gs.aliquot_barcode,
	gs.chrom,
	ca.arm,
	ca.pos * gs.pos AS pos,
	upper(ca.pos * gs.pos)::numeric - lower(ca.pos * gs.pos)::numeric - 1::numeric AS seg_size,
	sum(upper(ca.pos * gs.pos)::numeric - lower(ca.pos * gs.pos)::numeric - 1::numeric) OVER w AS arm_size,
	2::numeric ^ gs.log2_copy_ratio AS cr,
	gs.cnv_call AS cnv,
	row_number() OVER w2 - row_number() OVER w3 AS grp
	FROM analysis.gatk_seg_call gs
	JOIN ref.chr_arms ca ON ca.chrom = gs.chrom AND ca.pos && gs.pos
	WHERE (gs.chrom <> ALL (ARRAY[23, 24])) AND ca.arm <> '21p'::text
	WINDOW w AS (PARTITION BY gs.aliquot_barcode, ca.arm), w2 AS (PARTITION BY gs.aliquot_barcode, ca.arm ORDER BY (ca.pos * gs.pos)), w3 AS (PARTITION BY gs.aliquot_barcode, ca.arm, gs.cnv_call ORDER BY (ca.pos * gs.pos))
), 
join_adjacent_segments AS 
(
	SELECT t1.aliquot_barcode,
	t1.chrom,
	t1.arm,
	t1.grp,
	int4range(min(lower(t1.pos)), max(upper(t1.pos))) AS pos,
	t1.cnv,
	count(*) AS num_seg,
	sum(t1.seg_size * t1.cr) / sum(t1.seg_size) AS wcr,
	CASE
	WHEN count(*) > 1 AND sum(t1.seg_size * (t1.cr ^ 2::numeric)) > ((sum(t1.seg_size * t1.cr) ^ 2::numeric) / sum(t1.seg_size)) THEN sqrt((sum(t1.seg_size * (t1.cr ^ 2::numeric)) - (sum(t1.seg_size * t1.cr) ^ 2::numeric) / sum(t1.seg_size)) / (sum(t1.seg_size) - 1::numeric))
	ELSE 0::numeric
	END AS wsd,
	max(upper(t1.pos))::numeric - min(lower(t1.pos))::numeric - 1::numeric AS join_seg_size
	FROM arrange_segments_by_arm t1
	WHERE t1.cr > 0::numeric AND t1.seg_size > 0::numeric
	GROUP BY t1.aliquot_barcode, t1.chrom, t1.arm, t1.grp, t1.cnv
	ORDER BY t1.aliquot_barcode, t1.chrom, t1.arm, (int4range(min(lower(t1.pos)), max(upper(t1.pos))))
),
cnv_size AS 
(
	SELECT ad.aliquot_barcode,
	ad.chrom,
	ad.arm,
	ad.pos,
	join_seg_size AS seg_size,
	sum(
	CASE
	WHEN ad.cnv = 0 THEN upper(ad.pos)::numeric - lower(ad.pos)::numeric - 1::numeric
	ELSE 0
	END) AS het_size,
	sum(
	CASE
	WHEN ad.cnv = 1 THEN upper(ad.pos)::numeric - lower(ad.pos)::numeric - 1::numeric
	ELSE 0
	END) AS amp_size,
	sum(
	CASE
	WHEN ad.cnv = -1 THEN upper(ad.pos)::numeric - lower(ad.pos)::numeric - 1::numeric
	ELSE 0
	END) AS del_size
	FROM join_adjacent_segments ad
	WHERE ad.chrom <> ALL (ARRAY[23, 24])
	GROUP BY ad.aliquot_barcode, ad.chrom, ad.arm, ad.pos, ad.join_seg_size
	ORDER BY 1, 2, 3, 4
),
prop_amp_del AS
(
	SELECT 
	aliquot_barcode, 
	chrom,
	arm,
	ROUND(1.0 - (SUM(het_size)::numeric/SUM(seg_size)::numeric), 4) AS prop_aneuploidy,
	ROUND(SUM(amp_size)::numeric/SUM(seg_size)::numeric, 4) AS prop_amp,
	ROUND(SUM(del_size)::numeric/SUM(seg_size)::numeric, 4) AS prop_del
	FROM cnv_size
	GROUP BY aliquot_barcode, chrom, arm
	ORDER BY 1, 2, 3
)
SELECT ps.*, 
an1.chrom,
an1.arm,
an1.prop_aneuploidy AS prop_aneuploidy_a, 
an2.prop_aneuploidy AS prop_aneuploidy_b,
an2.prop_aneuploidy - an1.prop_aneuploidy AS aneuploidy_dif,
an1.prop_amp AS prop_amp_a, 
an2.prop_amp AS prop_amp_b,
an2.prop_amp - an1.prop_amp AS amp_dif,
an1.prop_del AS prop_del_a, 
an2.prop_del AS prop_del_b,
an2.prop_del - an1.prop_del AS del_dif,
im1.immune_score AS immune_score_a,
im1.immune_score AS immune_score_b,
im2.immune_score - im1.immune_score AS immune_dif,
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN prop_amp_del an1 ON an1.aliquot_barcode = ps.dna_barcode_a 
JOIN prop_amp_del an2 ON an2.aliquot_barcode = ps.dna_barcode_b AND an1.chrom = an2.chrom AND an1.arm = an2.arm
JOIN analysis.estimate im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate im2 ON im2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY 1,2,3
"
dat <- dbGetQuery(con,q)

arms <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

##################################################
# Step 1: Correlate immune scores/purity score differences with % CNA difference
##################################################

#Immune score


mycor1 <- pval1 <- mycor2 <- pval2 <- difcor <- difpval <- matrix(0,nrow=length(cells),ncol=length(subtypes))
rownames(mycor1) <- rownames(pval1) <- rownames(mycor2) <- rownames(pval2) <- rownames(difcor) <- rownames(difpval) <- cells
colnames(mycor1) <- colnames(pval1) <- colnames(mycor2) <- colnames(pval2) <- colnames(difcor) <- colnames(difpval) <- subtypes
for(i in 1:length(arms))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"arm"] == arms[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		mycor1[i,j] <- cor(sub_dat[,"immune_score_a"],sub_dat[,"prop_aneuploidy_a"],method="s")
		pval1[i,j] <- cor.test(sub_dat[,"immune_score_a"],sub_dat[,"prop_aneuploidy_a"],method="s")$p.value
		mycor2[i,j] <- cor(sub_dat[,"immune_score_b"],sub_dat[,"prop_aneuploidy_b"],method="s")
		pval2[i,j] <- cor.test(sub_dat[,"immune_score_b"],sub_dat[,"prop_aneuploidy_b"],method="s")$p.value

		difcor[i,j] <- cor(sub_dat[,"immune_dif"] , sub_dat[,"aneuploidy_dif"],method="s")
		difpval[i,j] <- cor.test(sub_dat[,"es_dif"],sub_dat[,"aneuploidy_dif"],method="s")$p.value
	}
}

