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
WITH ranked_neoag AS
(
	SELECT aliquot_barcode, case_barcode, variant_id, transcript, pvacseq_variant_type, gene_name, hla_allele, netmhcpan_mt_score,
	rank() OVER (PARTITION BY variant_id ORDER BY netmhcpan_mt_score) AS affinity
	FROM analysis.neoantigens_by_aliquot
	WHERE ssm2_pass_call IS true --AND pyclone_ccf >= 0.1 AND  pyclone_ccf < 0.5
),
top_neo AS
(
	SELECT *
	FROM ranked_neoag
	WHERE affinity = 1
	ORDER BY 1
),
hla_loss AS
(
	SELECT aliquot_barcode, case_barcode, variant_id, transcript, pvacseq_variant_type, gene_name, netmhcpan_mt_score, hla_allele, loss_allele, kept_allele, pval,
	CASE WHEN hla_allele = SUBSTRING(lh.loss_allele, 1, 11) THEN 1 ELSE 0 END AS hla_loss,
	CASE WHEN hla_allele = SUBSTRING(lh.kept_allele, 1, 11) THEN 1 ELSE 0 END AS hla_kept
	FROM top_neo tn
	JOIN analysis.pairs pa ON pa.tumor_barcode = tn.aliquot_barcode
	JOIN variants.lohhla lh ON lh.pair_barcode = pa.pair_barcode AND (tn.hla_allele = SUBSTRING(lh.hla_type1,1,11) OR tn.hla_allele = SUBSTRING(lh.hla_type2,1,11))
	WHERE lh.coverage_filter = 20 AND (pval < 0.1 AND (hla_type1_copy_number < 0.5 OR hla_type2_copy_number < 0.5))
)
SELECT aliquot_barcode, SUM(hla_loss) AS loss_count, SUM(hla_kept) AS kept_count, cs.idh_codel_subtype
FROM hla_loss hl
JOIN analysis.diamond_set ds ON ds.tumor_barcode_a = hl.aliquot_barcode OR ds.tumor_barcode_b = hl.aliquot_barcode
JOIN clinical.subtypes cs ON ds.case_barcode = cs.case_barcode
GROUP BY aliquot_barcode, idh_codel_subtype
ORDER BY 1
"

dat <- dbGetQuery(con,q)

wilcox.test(as.numeric(dat[,2]),as.numeric(dat[,3]),paired=TRUE)	#P = 0.57

status <- c(rep("Kept",nrow(dat)),rep("Lost",nrow(dat)))
count <- as.numeric(c(dat[,"kept_count"],dat[,"loss_count"]))
aliquot_barcode <- rep(dat[,"aliquot_barcode"],2)
plot_dat <- data.frame(aliquot_barcode,status,count)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/hla_neoag_loss_v_kept.pdf",width=2,height=2)
ggplot(data = plot_dat, aes(x = status, y = count, colour= status)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=aliquot_barcode),colour= "black") +
geom_point(size=1) +
scale_colour_manual(values=c("royalblue4","tomato3")) +
labs(y = "Neoantigen count") +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim=c(1,40))
dev.off()

#Subquery to account for expression:
# expr_neo AS
# (
# 	SELECT tn.*, tr.tpm
# 	FROM top_neo tn
# 	JOIN analysis.analyte_sets an ON an.dna_barcode = tn.aliquot_barcode
# 	JOIN analysis.transcript_tpm tr ON tr.target_id = tn.transcript AND tr.aliquot_barcode = an.rna_barcode
# 	WHERE tpm > 0
# 	ORDER BY 1
# ),