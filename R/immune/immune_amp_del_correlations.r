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
WITH sizes AS
(
	SELECT gs.aliquot_barcode,
	sum(upper(gs.pos) - lower(gs.pos) - 1) AS seg_size,
	sum(
	CASE
	WHEN gs.cnv_call = 0 THEN upper(gs.pos) - lower(gs.pos) - 1
	ELSE 0
	END) AS het_size,
	sum(
	CASE
	WHEN gs.cnv_call = 1 THEN upper(gs.pos) - lower(gs.pos) - 1
	ELSE 0
	END) AS amp_size,
	sum(
	CASE
	WHEN gs.cnv_call = -1 THEN upper(gs.pos) - lower(gs.pos) - 1
	ELSE 0
	END) AS del_size
	FROM analysis.gatk_seg_call gs
	WHERE gs.chrom <> ALL (ARRAY[23, 24])
	GROUP BY gs.aliquot_barcode
),
prop_amp_del AS
(
	SELECT 
	aliquot_barcode, 
	ROUND(1.0 - (het_size::numeric/seg_size::numeric), 4) AS prop_aneuploidy,
	ROUND(amp_size::numeric/seg_size::numeric, 4) AS prop_amp,
	ROUND(del_size::numeric/seg_size::numeric, 4) AS prop_del
	FROM sizes
),
signatures AS
(
	SELECT * FROM analysis.muller_tam_score
	UNION
	SELECT * FROM analysis.davoli_immune_score
	WHERE signature_name NOT LIKE 'Macrophages%' AND signature_name != 'CD8.effector.NK.cells'
	UNION 
	SELECT aliquot_barcode, 'Purity' AS signature_name, purity AS enrichment_score
	FROM analysis.estimate
)
SELECT ps.*, 
an1.prop_aneuploidy AS prop_aneuploidy_a, 
an2.prop_aneuploidy AS prop_aneuploidy_b,
an2.prop_aneuploidy - an1.prop_aneuploidy AS aneuploidy_dif,
an1.prop_amp AS prop_amp_a, 
an2.prop_amp AS prop_amp_b,
an2.prop_amp - an1.prop_amp AS amp_dif,
an1.prop_del AS prop_del_a, 
an2.prop_del AS prop_del_b,
an2.prop_del - an1.prop_del AS del_dif,
im1.signature_name,
im1.enrichment_score AS enrichment_score_a,
im2.enrichment_score AS enrichment_score_b,
im2.enrichment_score - im1.enrichment_score AS es_dif,
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN prop_amp_del an1 ON an1.aliquot_barcode = ps.dna_barcode_a 
JOIN prop_amp_del an2 ON an2.aliquot_barcode = ps.dna_barcode_b 
JOIN signatures im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN signatures im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im1.signature_name = im2.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY 1,9
"
dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

##################################################
# Step 1: Correlate immune scores/purity score differences with % CNA difference
##################################################

#Immune score

subtype <- rep(subtypes, each=81)
signature <- rep(rep(cells,each=3),3)
cnv_class <- rep(c("Amp","Del","Any"),81)
rho <- rep(0, length(subtype))
pval <- rep(0, length(subtype))
res <- data.frame(subtype, signature, cnv_class, rho, pval)

for(i in 1:length(cells))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"signature_name"] == cells[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		
		res[which(res[,"signature"] == cells[i] & res[,"subtype"]==subtypes[j] & res[,"cnv_class"] == "Any"),"rho"] <- 
			cor(sub_dat[,"es_dif"],sub_dat[,"aneuploidy_dif"],method="s")
		res[which(res[,"signature"] == cells[i] & res[,"subtype"]==subtypes[j] & res[,"cnv_class"] == "Any"),"pval"] <- 
			cor.test(sub_dat[,"es_dif"],sub_dat[,"aneuploidy_dif"],method="s")$p.value
		
		res[which(res[,"signature"] == cells[i] & res[,"subtype"]==subtypes[j] & res[,"cnv_class"] == "Amp"),"rho"] <- 
			cor(sub_dat[,"es_dif"],sub_dat[,"amp_dif"],method="s")
		res[which(res[,"signature"] == cells[i] & res[,"subtype"]==subtypes[j] & res[,"cnv_class"] == "Amp"),"pval"] <- 
			cor.test(sub_dat[,"es_dif"],sub_dat[,"amp_dif"],method="s")$p.value

		res[which(res[,"signature"] == cells[i] & res[,"subtype"]==subtypes[j] & res[,"cnv_class"] == "Del"),"rho"] <- 
			cor(sub_dat[,"es_dif"],sub_dat[,"del_dif"],method="s")
		res[which(res[,"signature"] == cells[i] & res[,"subtype"]==subtypes[j] & res[,"cnv_class"] == "Del"),"pval"] <- 
			cor.test(sub_dat[,"es_dif"],sub_dat[,"del_dif"],method="s")$p.value

	}
}

marker <- ifelse(res[,"pval"] < 0.1, "·","")
marker[which(res[,"pval"] < 0.05)] <- "*"
res <- cbind(res, marker)

res[,"signature"] <- factor(res[,"signature"],levels=c("B.cells","CD4.mature","CD8.effector","Dendritic","Macrophages","Microglia","NK.cells","T.reg","Purity"))
res[,"cnv_class"] <- factor(res[,"cnv_class"],levels=rev(c("Any","Amp","Del")))
res[,"subtype"] <- factor(res[,"subtype"],levels=c("IDHwt","IDHmut-noncodel","IDHmut-codel"))

plot_res <- res[which(res[,"signature"] != "Purity"),]

##################################################
# Step 2: Plot the results
##################################################

# All subtypes (does not look good)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cna_amp_del_barplot.pdf",width=4.4,height=4)
ggplot(plot_res, aes(x=signature,y=rho,fill=cnv_class)) +
geom_bar(stat="identity",position = position_dodge()) +
facet_wrap(subtype~.,ncol=3) +
scale_fill_manual(values=rev(c("#33a02c","#e31a1c", "#1f78b4"))) +
theme_bw() +
theme(axis.title=element_blank(),
axis.text.x=element_text(size=7,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "bottom") +
guides(color = guide_legend(override.aes = list(size = 0.1))) +
guides(shape = guide_legend(override.aes = list(size = 0.1))) +
coord_flip()
dev.off()


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cna_amp_del_heatmap.pdf",width=3.2,height=3)
ggplot(plot_res, aes(x=signature,y=cnv_class,fill=rho)) +
geom_tile() +
facet_wrap(subtype~.,nrow=3) +
geom_text(aes(label=marker),colour="black") +
scale_fill_gradient2(low="royalblue4", mid="#ffffff", high="tomato3") +
guides(fill = guide_colourbar(barwidth = 0.75, barheight = 3)) +
theme_bw() +
theme(axis.title=element_blank(),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "right")
dev.off()
