library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggbeeswarm)


#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
WITH neoag_by_ali AS
(
	SELECT aliquot_barcode, variant_id, gene_name, mutation, pvacseq_protein_position, peptide_length, sub_peptide_position, mt_epitope_seq 
	FROM analysis.neoantigens_by_aliquot neo
	WHERE ssm2_pass_call = TRUE 
	GROUP BY aliquot_barcode, variant_id, gene_name, mutation, pvacseq_protein_position, peptide_length, sub_peptide_position, mt_epitope_seq
),
idh_counts AS
(
	SELECT aliquot_barcode, 
	CASE 
		WHEN gene_name = 'IDH1' AND mutation = 'R/H' AND pvacseq_protein_position = '132' THEN 1
		ELSE 0
	END AS idh1_neoag
	FROM neoag_by_ali
),
idh_neoag AS
(
	SELECT aliquot_barcode, SUM(idh1_neoag) 
	FROM idh_counts
	GROUP BY aliquot_barcode
)
SELECT ps.case_barcode, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b, 
na1.sum AS idh_neoag_a, 
na2.sum AS idh_neoag_b,
bs1.sample_type AS sample_type_a
FROM analysis.platinum_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN idh_neoag na1 ON na1.aliquot_barcode = ps.dna_barcode_a
JOIN idh_neoag na2 ON na2.aliquot_barcode = ps.dna_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
WHERE cs1.idh_status = 'IDHmut' AND bs1.sample_type = 'TP'
ORDER BY 2
"

dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])

pval_a <- pval_b <- pval_idh <- pval_no <- eff_a <- eff_b <- eff_idh <- eff_no <- rep(0,length(cells))
for(i in 1:length(cells))
{

	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	g1a <- sub_dat[which(sub_dat[,"idh_neoag_a"] > 0),"es_a"]
	g2a <- sub_dat[which(sub_dat[,"idh_neoag_a"] == 0),"es_a"]

	pval_a[i] <- wilcox.test(g1a,g2a)$p.value
	eff_a[i] <- median(g1a) - median(g2a)
	
	g1b <- sub_dat[which(sub_dat[,"idh_neoag_b"] > 0),"es_b"]
	g2b <- sub_dat[which(sub_dat[,"idh_neoag_b"] == 0),"es_b"]

	pval_b[i] <- wilcox.test(g1b,g2b)$p.value
	eff_b[i] <- median(g1b) - median(g2b)
	
	pval_idh[i] <- wilcox.test(g1a,g1b,paired=TRUE)$p.value
	eff_idh[i] <- median(g1b)-median(g1a)
	
	pval_no[i] <- wilcox.test(g2a,g2b,paired=TRUE)$p.value
	eff_no[i] <- median(g2b)-median(g2a)
}

res <- data.frame(eff_a, pval_a, eff_b, pval_b, pval_idh, eff_idh, pval_no, eff_no)
rownames(res) <- cells


es <- c(dat[,"es_a"],dat[,"es_b"])
idh_neoag <- as.factor(as.numeric(c(dat[,"idh_neoag_a"],dat[,"idh_neoag_b"]) > 0))
signature_name <- rep(dat[,"signature_name"],2)
status <- c(rep("Initial",nrow(dat)), rep("Recurrent",nrow(dat)))
case_barcode <- rep(dat[,"case_barcode"],2)
plot_res <- data.frame(case_barcode, signature_name, es, idh_neoag, status)

plot_res <- plot_res[which(plot_res[,"signature_name"] %in% c("B.cells","CD4.mature","CD8.effector","CD8.effector.NK.cells","Dendritic","Macrophages","NK.cells","T.reg")),]

#Plot showing ladder plots between initial and recurrence for all subtypes
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/idh1_neoantigen_idhmut.pdf",width=7,height=3)
ggplot(plot_res,aes(x = idh_neoag, y = es, fill= idh_neoag)) +
geom_boxplot(lwd=0.4,outlier.colour=NA,fatten=2) +
geom_jitter(width = 0.3, size = 0.5) +
#stat_summary(fun.y=median,geom="line",colour="red") +
labs(y = "Enrichment score") +
scale_fill_manual(values= c("royalblue4","tomato3")) +
facet_grid(status~signature_name) +
theme_bw() +
theme(axis.text.x= element_text(size=7,hjust=0.5),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none")
dev.off()

#Plot showing ladder plots between initial and recurrence for all subtypes
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/idh1_neoantigen_idhmut_time.pdf",width=7,height=3)
ggplot(plot_res,aes(x = status, y = es, group= case_barcode, colour=idh_neoag)) +
geom_line(size=0.45,alpha=0.4) +
geom_point(size=1,colour="black") +
#stat_summary(fun.y=median,geom="line",colour="red") +
labs(y = "Enrichment score") +
scale_colour_manual(values= c("royalblue4","tomato3")) +
facet_grid(.~signature_name) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none")
dev.off()
