library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Query to select all IDHwt samples that received radiation
q <- "
SELECT ps.*, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs.idh_codel_subtype
FROM analysis.tumor_clinical_comparison ps
JOIN analysis.rna_dna_pairs rd ON rd.dna_pair_barcode = ps.tumor_pair_barcode
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = rd.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = rd.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
WHERE received_rt IS true AND idh_codel_subtype ='IDHwt'
ORDER BY 2, 1, 16
"


dat <- dbGetQuery(con, q)

cells <- unique(dat[,"signature_name"])

#for loop testing which cells differ pre- and post-radiation therapy
eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}

sig <- rep("",length(p.value))
sig[which(p.value < 0.05)] <- "*"
sig[which(p.value < 0.01)] <- "**"
time_res <- data.frame(cells,eff,p.value,sig)

#Plot the result
aliquot_barcode <- c(dat[,"tumor_barcode_a"],dat[,"tumor_barcode_b"])
case_barcode <- rep(dat[,"case_barcode"],2)
signature_name <- rep(dat[,"signature_name"],2)
enrichment_score <- c(dat[,"es_a"],dat[,"es_b"])
subtype <- rep(dat[,"idh_codel_subtype"],2)
timepoint <- c(rep("Initial",nrow(dat)),rep("Recurrent",nrow(dat)))
plot_res <- data.frame(aliquot_barcode, case_barcode, signature_name, enrichment_score, timepoint, subtype)
cd4_res <- plot_res[which(plot_res[,"signature_name"]=="CD4.mature"),]
treg_res <- plot_res[which(plot_res[,"signature_name"]=="T.reg"),]

pval1 <- signif(time_res[which(time_res[,"cells"]=="CD4.mature"),"p.value"],digits=1)
pval1 <- formatC(pval1, format = "e", digits = 0)
pval2 <- signif(time_res[which(time_res[,"cells"]=="T.reg"),"p.value"],digits=1)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/idhwt_cd4_ladderplot.pdf",width=3,height=2)
cd4 <- ggplot(cd4_res,aes(x = timepoint, y = enrichment_score)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode,colour=subtype)) +
geom_point(size=1,colour="black") +
scale_colour_manual(values="#619CFF") +
labs(title="CD4+ T cells", y = "Enrichment score") +
annotate("text", x=1.1, y = min(cd4_res[,"enrichment_score"])-(0.15 * (max(cd4_res[,"enrichment_score"]) - min(cd4_res[,"enrichment_score"]))) + 
	0.08 * abs(min(cd4_res[,"enrichment_score"])-(0.15 * (max(cd4_res[,"enrichment_score"]) - min(cd4_res[,"enrichment_score"])))), 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval1)))), parse=TRUE, size=2.5) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(min(cd4_res[,"enrichment_score"])-(0.15 * (max(cd4_res[,"enrichment_score"]) - min(cd4_res[,"enrichment_score"]))),
max(cd4_res[,"enrichment_score"])))

tregs <- ggplot(treg_res,aes(x = timepoint, y = enrichment_score)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode,colour=subtype)) +
geom_point(size=1,colour="black") +
scale_colour_manual(values="#619CFF") +
labs(title="T-regs", y = "Enrichment score") +
annotate("text", x=1.1, y = min(treg_res[,"enrichment_score"])-(0.15 * (max(treg_res[,"enrichment_score"]) - min(treg_res[,"enrichment_score"]))) + 
	0.04 * abs(min(treg_res[,"enrichment_score"])-(0.15 * (max(treg_res[,"enrichment_score"]) - min(treg_res[,"enrichment_score"])))), 
	label=deparse((bquote(italic("P") ~" = " ~ .(pval2)))), parse=TRUE, size=2.5) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(min(treg_res[,"enrichment_score"])-(0.15 * (max(treg_res[,"enrichment_score"]) - min(treg_res[,"enrichment_score"]))),
max(treg_res[,"enrichment_score"])))

grid.arrange(cd4,tregs,nrow=1)
dev.off()
