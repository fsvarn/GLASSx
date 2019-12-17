library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggbeeswarm)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "SELECT ps.case_barcode, 
ps.rna_barcode_a, 
ps.rna_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.shannon AS tcr_shannon_a, 
im2.shannon AS tcr_shannon_b, 
im1.evenness AS tcr_evenness_a,
im2.evenness AS tcr_evenness_b,
im1.richness AS tcr_richness_a,
im2.richness AS tcr_richness_b,
im1.total_tcr AS total_tcr_a,
im2.total_tcr AS total_tcr_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b,
CASE WHEN mf1.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf1.coverage_adj_mut_freq < 10 THEN 0 END AS hm_a,
CASE WHEN mf2.coverage_adj_mut_freq >= 10 THEN 1 WHEN mf2.coverage_adj_mut_freq < 10 THEN 0 END AS hm_b
FROM analysis.platinum_set ps
JOIN analysis.tcr_stats im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.tcr_stats im2 ON im2.aliquot_barcode = ps.rna_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
JOIN analysis.mut_freq mf1 ON mf1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.mut_freq mf2 ON mf2.aliquot_barcode = ps.dna_barcode_b
ORDER BY 1, 2, 6"

hm <- dbGetQuery(con, q)

#See if there are significant changes across samples over time

#Where are those missing TCGA samples???? (Microarray??)

hm_only <- hm[which(hm[,"hm_b"] ==1 & hm[,"hm_a"]==0),]

eff <- median(hm_only[,"total_tcr_a"]) - median(hm_only[,"total_tcr_b"])
p.value <- wilcox.test(hm_only[,"total_tcr_a"], hm_only[,"total_tcr_b"],paired=TRUE)$p.value

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/scatterplot_tcr_hypermutator.pdf",width=2,height=2)
ggplot(hm_only,aes(x = total_tcr_a, y = total_tcr_b,colour=subtype_b)) +
geom_point(size=0.5) +
geom_abline(intercept=0,slope=1,colour="gray") +
scale_colour_manual(values=c("#00BA38","#619CFF")) +
labs(x = "TCR abundance (initial)", y = "TCR abundance (recurrent)", title="TCR abundance in hypermutators") +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")  +
coord_cartesian(xlim=c(0,50),ylim=c(0,50))
dev.off()
