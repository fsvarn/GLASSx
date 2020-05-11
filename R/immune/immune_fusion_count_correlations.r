###################################################
# Test how fusion counts associate with immune levels
# Updated: 2020.04.22
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
WITH fusion_count AS
(
	SELECT aliquot_barcode, COUNT(*) AS fusion_count
	FROM variants.prada_fusions
	WHERE discordant_n >= 2 AND jsr_n >= 1 AND NOT called_in_normal
	GROUP BY aliquot_barcode
),
full_count AS
(
	SELECT al.aliquot_barcode, COALESCE(fusion_count, 0) AS fusion_count
	FROM biospecimen.aliquots al
	LEFT JOIN fusion_count fu ON al.aliquot_barcode = fu.aliquot_barcode
	WHERE aliquot_analysis_type = 'RNA'
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
im1.signature_name,
im1.enrichment_score AS enrichment_score_a,
im2.enrichment_score AS enrichment_score_b,
im2.enrichment_score - im1.enrichment_score AS es_dif,
fc1.fusion_count AS fusion_count_a,
fc2.fusion_count AS fusion_count_b,
fc2.fusion_count - fc1.fusion_count AS fusion_dif,
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a 
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b 
JOIN signatures im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN signatures im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im1.signature_name = im2.signature_name
JOIN full_count fc1 ON fc1.aliquot_barcode = ps.rna_barcode_a
JOIN full_count fc2 ON fc2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY 1,9
"
dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

dat[,"fusion_count_a"] <- as.numeric(dat[,"fusion_count_a"])
dat[,"fusion_count_b"] <- as.numeric(dat[,"fusion_count_b"])
dat[,"fusion_dif"] <- as.numeric(dat[,"fusion_dif"])

##################################################
# Step 1: Correlate immune scores/purity score with % CNAs
##################################################

#Immune score

mycor1 <- pval1 <- mycor2 <- pval2 <- mycor3 <- mycor4 <- difcor <- difpval <- matrix(0,nrow=length(cells),ncol=length(subtypes))
rownames(mycor1) <- rownames(pval1) <- rownames(mycor2) <- rownames(pval2) <- rownames(mycor3) <- rownames(mycor4) <- rownames(difcor) <- rownames(difpval) <- cells
colnames(mycor1) <- colnames(pval1) <- colnames(mycor2) <- colnames(pval2) <- colnames(mycor3) <- colnames(mycor4) <- colnames(difcor) <- colnames(difpval) <- subtypes
for(i in 1:length(cells))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"signature_name"] == cells[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		mycor1[i,j] <- cor(sub_dat[,"enrichment_score_a"],sub_dat[,"fusion_count_a"],method="s")
		pval1[i,j] <- cor.test(sub_dat[,"enrichment_score_a"],sub_dat[,"fusion_count_a"],method="s")$p.value
		mycor2[i,j] <- cor(sub_dat[,"enrichment_score_b"],sub_dat[,"fusion_count_b"],method="s")
		pval2[i,j] <- cor.test(sub_dat[,"enrichment_score_b"],sub_dat[,"fusion_count_b"],method="s")$p.value

		difcor[i,j] <- cor(sub_dat[,"es_dif"] , sub_dat[,"fusion_dif"], method="s")
		difpval[i,j] <- cor.test(sub_dat[,"es_dif"] , sub_dat[,"fusion_dif"], method="s")$p.value
	}
}

ini_cor <- melt(mycor1)
rec_cor <- melt(mycor2)
dif_cor <- melt(difcor)

ini_pval <- melt(pval1)
rec_pval <- melt(pval2)
dif_pval <- melt(difpval)
all_pval <- rbind(ini_pval, rec_pval, dif_pval)

timepoint <- c(rep("Initial",nrow(ini_cor)), rep("Recurrent",nrow(rec_cor)), rep("Difference",nrow(dif_cor)))
p_val <- all_pval[,"value"]
marker <- ifelse(p_val < 0.1, "·","")
marker[which(p_val < 0.05)] <- "*"

plot_res <- rbind(ini_cor, rec_cor, dif_cor)
plot_res <- data.frame(plot_res, timepoint, p_val, marker)
colnames(plot_res) <- c("signature","subtype","correlation","timepoint","p_val","marker")
plot_res[,"signature"] <- factor(plot_res[,"signature"],levels=c("B.cells","CD4.mature","CD8.effector","Dendritic","Macrophages","Microglia","NK.cells","T.reg","Purity"))
plot_res[,"timepoint"] <- factor(plot_res[,"timepoint"],levels=c("Initial","Recurrent","Difference"))

##################################################
# Step 2: Plot the results
##################################################

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_fusion_heatmap.pdf",width=3.2,height=3)
ggplot(plot_res, aes(x=signature,y=subtype,fill=correlation)) +
geom_tile() +
facet_wrap(timepoint~.,nrow=3) +
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

##################################################
# Step 3: Correlate fusion count with chromosomal instability
##################################################


init_cin_cor <- init_cin_pval <- rec_cin_cor <- rec_cin_pval <- dif_cin_cor <- dif_cin_pval <- rep(0, length(subtypes))
for(j in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"signature_name"] == cells[1] & dat[,"idh_codel_subtype"]==subtypes[j]),]
	
	init_cin_cor[j] <- cor(sub_dat[,"prop_aneuploidy_a"], sub_dat[,"fusion_count_a"], method="s")
	init_cin_pval[j] <- cor.test(sub_dat[,"prop_aneuploidy_a"], sub_dat[,"fusion_count_a"], method="s")$p.value

	rec_cin_cor[j] <- cor(sub_dat[,"prop_aneuploidy_b"], sub_dat[,"fusion_count_b"], method="s")
	rec_cin_pval[j] <- cor.test(sub_dat[,"prop_aneuploidy_b"], sub_dat[,"fusion_count_b"], method="s")$p.value

	dif_cin_cor[j] <- cor(sub_dat[,"aneuploidy_dif"], sub_dat[,"fusion_dif"], method="s")
	dif_cin_pval[j] <- cor.test(sub_dat[,"aneuploidy_dif"], sub_dat[,"fusion_dif"], method="s")$p.value
}

