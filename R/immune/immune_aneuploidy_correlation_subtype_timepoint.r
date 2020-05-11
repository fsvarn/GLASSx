###################################################
# Test % of genome with CNV associates with immune levels
# Updated: 2020.04.16
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
WITH signatures AS
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
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a 
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b 
JOIN signatures im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN signatures im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im1.signature_name = im2.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY 1,9
"
dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

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
		mycor1[i,j] <- cor(sub_dat[,"enrichment_score_a"],sub_dat[,"prop_aneuploidy_a"],method="s")
		pval1[i,j] <- cor.test(sub_dat[,"enrichment_score_a"],sub_dat[,"prop_aneuploidy_a"],method="s")$p.value
		mycor2[i,j] <- cor(sub_dat[,"enrichment_score_b"],sub_dat[,"prop_aneuploidy_b"],method="s")
		pval2[i,j] <- cor.test(sub_dat[,"enrichment_score_b"],sub_dat[,"prop_aneuploidy_b"],method="s")$p.value
		mycor3[i,j] <- cor(sub_dat[,"enrichment_score_a"],sub_dat[,"prop_aneuploidy_b"],method="s")
		mycor4[i,j] <- cor(sub_dat[,"enrichment_score_b"],sub_dat[,"prop_aneuploidy_a"],method="s")
	
		difcor[i,j] <- cor(sub_dat[,"es_dif"] , sub_dat[,"aneuploidy_dif"],method="s")
		difpval[i,j] <- cor.test(sub_dat[,"es_dif"],sub_dat[,"aneuploidy_dif"],method="s")$p.value
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

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cna_heatmap.pdf",width=3.2,height=3)
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
# Step 3: GLM to examine whether this association exists after adjusting for purity
##################################################

q <- "
WITH signatures AS
(
	SELECT * FROM analysis.muller_tam_score
	UNION
	SELECT * FROM analysis.davoli_immune_score
	WHERE signature_name NOT LIKE 'Macrophages%' AND signature_name != 'CD8.effector.NK.cells'
)

SELECT ps.*, 
an1.prop_aneuploidy AS prop_aneuploidy_a, 
an2.prop_aneuploidy AS prop_aneuploidy_b,
an2.prop_aneuploidy - an1.prop_aneuploidy AS aneuploidy_dif,
im1.signature_name,
im1.enrichment_score AS enrichment_score_a,
im2.enrichment_score AS enrichment_score_b,
im2.enrichment_score - im1.enrichment_score AS es_dif,
pu1.purity AS purity_a,
pu2.purity AS purity_b,
pu2.purity - pu1.purity AS purity_dif,
idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a 
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b 
JOIN signatures im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN signatures im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im1.signature_name = im2.signature_name
JOIN analysis.estimate pu1 ON pu1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.estimate pu2 ON pu2.aliquot_barcode = ps.rna_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY 1,9
"
dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"idh_codel_subtype"])

eff1 <- pval1 <- eff2 <- pval2 <- difeff <- difpval <- matrix(0,nrow=length(cells),ncol=length(subtypes))
rownames(eff1) <- rownames(pval1) <- rownames(eff2) <- rownames(pval2) <- rownames(difeff) <- rownames(difpval) <- cells
colnames(eff1) <- colnames(pval1) <- colnames(eff2) <- colnames(pval2) <- colnames(difeff) <- colnames(difpval) <- subtypes
for(i in 1:length(cells))
{
	for(j in 1:length(subtypes))
	{
		sub_dat <- dat[which(dat[,"signature_name"] == cells[i] & dat[,"idh_codel_subtype"]==subtypes[j]),]
		eff1[i,j] <- summary(glm(enrichment_score_a ~ prop_aneuploidy_a + purity_a, family ="gaussian", data = sub_dat))[["coefficients"]][2]
		pval1[i,j] <- summary(glm(enrichment_score_a ~ prop_aneuploidy_a + purity_a, family ="gaussian", data = sub_dat))[["coefficients"]][11]
		eff2[i,j] <- summary(glm(enrichment_score_b ~ prop_aneuploidy_b + purity_b, family ="gaussian", data = sub_dat))[["coefficients"]][2]
		pval2[i,j] <- summary(glm(enrichment_score_b ~ prop_aneuploidy_b + purity_b, family ="gaussian", data = sub_dat))[["coefficients"]][11]

		difeff[i,j] <- summary(glm(es_dif ~ aneuploidy_dif + purity_dif, family ="gaussian", data = sub_dat))[["coefficients"]][2]
		difpval[i,j] <- summary(glm(es_dif ~ aneuploidy_dif + purity_dif, family ="gaussian", data = sub_dat))[["coefficients"]][11]
	}
}


##################################################
# Step 4: Unified GLM that adjusts for subtype
##################################################

uni_eff1 <- uni_pval1 <- uni_eff2 <- uni_pval2 <- uni_difeff <- uni_difpval <- rep(0,nrow=length(cells))
for(i in 1:length(cells))
{

	sub_dat <- dat[which(dat[,"signature_name"] == cells[i]),]
	uni_eff1[i] <- summary(glm(enrichment_score_a ~ prop_aneuploidy_a + purity_a + idh_codel_subtype, family ="gaussian", data = sub_dat))[["coefficients"]][2]
	uni_pval1[i] <- summary(glm(enrichment_score_a ~ prop_aneuploidy_a + purity_a + idh_codel_subtype, family ="gaussian", data = sub_dat))[["coefficients"]][17]
	uni_eff2[i] <- summary(glm(enrichment_score_b ~ prop_aneuploidy_b + purity_b + idh_codel_subtype, family ="gaussian", data = sub_dat))[["coefficients"]][2]
	uni_pval2[i] <- summary(glm(enrichment_score_b ~ prop_aneuploidy_b + purity_b + idh_codel_subtype, family ="gaussian", data = sub_dat))[["coefficients"]][17]

	uni_difeff[i] <- summary(glm(es_dif ~ aneuploidy_dif + purity_dif + idh_codel_subtype, family ="gaussian", data = sub_dat))[["coefficients"]][2]
	uni_difpval[i] <- summary(glm(es_dif ~ aneuploidy_dif + purity_dif + idh_codel_subtype, family ="gaussian", data = sub_dat))[["coefficients"]][17]
}
names(uni_eff1) <- names(uni_pval1) <- names(uni_eff2) <- names(uni_pval2) <- names(uni_difeff) <- names(uni_difpval) <- cells

