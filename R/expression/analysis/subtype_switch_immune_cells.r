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
q <- "
WITH subtype_label AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	ss.tumor_barcode_a,
	ss.tumor_barcode_b,
	ts1.signature_name, 
	ts1.p_value AS pval_a,
	ts2.p_value AS pval_b,
	CASE WHEN ts1.p_value < 0.05 THEN ts1.signature_name ELSE NULL END AS subtype_a,
	CASE WHEN ts2.p_value < 0.05 THEN ts2.signature_name ELSE NULL END AS subtype_b,
	su.idh_codel_subtype
	FROM analysis.rna_silver_set ss
	JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
	JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
	JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
	JOIN clinical.surgeries su ON su.sample_barcode = al1.sample_barcode
	WHERE idh_codel_subtype IS NOT NULL AND
	ts1.p_value < 0.05 OR ts2.p_value < 0.05
),
subtype_collapse AS
(
	SELECT sl.tumor_pair_barcode,
	sl.case_barcode,
	sl.tumor_barcode_a,
	sl.tumor_barcode_b,
	string_agg(subtype_a,',') AS subtype_a,
	string_agg(subtype_b,',') AS subtype_b,
	idh_codel_subtype
	FROM subtype_label sl
	GROUP BY tumor_pair_barcode, sl.case_barcode, sl.tumor_barcode_a, sl.tumor_barcode_b, idh_codel_subtype
	ORDER BY 1
)
SELECT sc.*,
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
im2.enrichment_score - im1.enrichment_score AS es_dif
FROM subtype_collapse sc
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = sc.tumor_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = sc.tumor_barcode_b AND im2.signature_name = im1.signature_name
"
dat <- dbGetQuery(con, q)

dat <- dat[which(dat[,"signature_name"] %in% c("CD4.mature","CD8.effector","NK.cells","B.cells","T.reg","Dendritic","Macrophages")),]

dat[which(is.na(dat[,"subtype_b"])),"subtype_b"] <- "Mixed"
dat[which(is.na(dat[,"subtype_a"])),"subtype_a"] <- "Mixed"

mesenchymal_status <- rep("none",nrow(dat))
mesenchymal_status[which(!(str_detect("Mesenchymal",dat[,"subtype_a"])) & str_detect("Mesenchymal",dat[,"subtype_b"]) |
	str_detect(",",dat[,"subtype_a"]) & str_detect("Mesenchymal",dat[,"subtype_b"]))] <- "gain"
mesenchymal_status[which(str_detect("Mesenchymal",dat[,"subtype_a"]) & !(str_detect("Mesenchymal",dat[,"subtype_b"])))] <- "loss"
dat[,"mesenchymal_status"] <- mesenchymal_status

cells <- unique(dat[,"signature_name"])

gain_pval <- gain_eff <- loss_pval <- loss_eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"] == cells[i]),]
	g1 <- sub_dat[which(sub_dat[,"mesenchymal_status"]=="gain"),"es_dif"]
	g2 <- sub_dat[which(sub_dat[,"mesenchymal_status"]=="none"),"es_dif"]
	g3 <- sub_dat[which(sub_dat[,"mesenchymal_status"]=="loss"),"es_dif"]
	gain_pval[i] <- wilcox.test(g1,g2)$p.value
	gain_eff[i] <- median(g1,na.rm=TRUE) - median(g2,na.rm=TRUE)
	loss_pval[i] <- wilcox.test(g3,g2)$p.value
	loss_eff[i] <- median(g3,na.rm=TRUE) - median(g2,na.rm=TRUE)
}
data.frame(cells,gain_pval,gain_eff,loss_pval,loss_eff)

dat[,"mesenchymal_status"] <- factor(dat[,"mesenchymal_status"],levels=c("none","gain","loss"))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/mesenchymal_change_immune_cell.pdf",width=2.6,height=1.2)
ggplot(dat,aes(x = signature_name, y = es_dif, fill=mesenchymal_status)) +
geom_hline(yintercept=0, color="gray50") +
geom_boxplot(outlier.size=0.25) +
labs(y = "Change in enrichment score") +
scale_fill_manual(values=c("gray50","tomato3","royalblue4")) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(-0.5,0.5))
dev.off()

