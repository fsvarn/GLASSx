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
WITH neoag AS
(
	SELECT aliquot_barcode, case_barcode, variant_id, transcript, pvacseq_variant_type, gene_name, hla_allele, netmhcpan_mt_score
	FROM analysis.neoantigens_by_aliquot
	WHERE ssm2_pass_call IS true
),
hla_loss AS
(
	SELECT aliquot_barcode, case_barcode, variant_id, transcript, pvacseq_variant_type, gene_name, netmhcpan_mt_score, hla_allele, loss_allele, kept_allele, pval,
	CASE WHEN hla_allele = SUBSTRING(lh.loss_allele, 1, 11) THEN netmhcpan_mt_score ELSE NULL END AS hla_loss,
	CASE WHEN hla_allele = SUBSTRING(lh.kept_allele, 1, 11) THEN netmhcpan_mt_score ELSE NULL END AS hla_kept
	FROM neoag tn
	JOIN analysis.pairs pa ON pa.tumor_barcode = tn.aliquot_barcode
	JOIN variants.lohhla lh ON lh.pair_barcode = pa.pair_barcode AND (tn.hla_allele = SUBSTRING(lh.hla_type1,1,11) OR tn.hla_allele = SUBSTRING(lh.hla_type2,1,11))
	WHERE lh.coverage_filter = 10 AND (pval < 0.05 AND (hla_type1_copy_number < 0.5 OR hla_type2_copy_number < 0.5))
),
together AS
(
SELECT aliquot_barcode, COUNT(hla_loss) AS loss_count, COUNT(hla_kept) AS kept_count,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM hla_loss hl
JOIN analysis.gold_set gs ON gs.tumor_barcode_a = hl.aliquot_barcode OR gs.tumor_barcode_b = hl.aliquot_barcode
JOIN clinical.subtypes cs ON gs.case_barcode = cs.case_barcode
GROUP BY aliquot_barcode, idh_codel_subtype
ORDER BY 1
)
SELECT * FROM together 
"

dat <- dbGetQuery(con,q)

dat %>%
group_by(idh_status) %>%
summarise(pval = wilcox.test(as.numeric(loss_count),as.numeric(kept_count),paired=TRUE)$p.value,
		  eff = median(loss_count) - median(kept_count))

ss <- dbReadTable(con, Id(schema="analysis",table="gold_set"))
timepoint <- rep("",nrow(dat))
timepoint[which(dat$aliquot_barcode %in% ss$tumor_barcode_a)] <- "Initial"
timepoint[which(dat$aliquot_barcode %in% ss$tumor_barcode_b)] <- "Recurrent"

dat$timepoint <- timepoint

status <- c(rep("Kept",nrow(dat)),rep("Lost",nrow(dat)))
count <- as.numeric(c(dat[,"kept_count"],dat[,"loss_count"]))
aliquot_barcode <- rep(dat[,"aliquot_barcode"],2)
timepoint <- rep(dat[,"timepoint"],2)
idh_status <- rep(dat[,"idh_status"],2)
plot_dat <- data.frame(aliquot_barcode,status,timepoint,count,idh_status)

plot_dat <- plot_dat %>%
mutate(idh_status = fct_relevel(idh_status, "IDHwt", "IDHmut"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/hla_neoag_count_loss_v_kept.pdf",width=2,height=1.8)
ggplot(data = plot_dat, aes(x = status, y = count, colour= timepoint)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.6,alpha=0.4,aes(group=aliquot_barcode),colour= "black") +
geom_point(size=1) +
scale_colour_manual(values=c("#a6611a","#018571")) +
labs(y = "Neoantigen count") +
facet_grid(.~idh_status) +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim = c(0,400))
dev.off()

