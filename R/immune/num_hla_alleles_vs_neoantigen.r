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
grouped_neo AS
(
	SELECT aliquot_barcode, hla_allele, COUNT(*) AS neoag_ct
	FROM top_neo
	GROUP BY aliquot_barcode, hla_allele
)
SELECT aliquot_barcode, COUNT(*) AS hla_num, SUM(neoag_ct) AS neoag_ct
FROM grouped_neo
GROUP BY aliquot_barcode

"

dat <- dbGetQuery(con,q)
dat[which(dat[,"hla_num"] < 3),] <- 3
dat[,"neoag_ct"] <- log10(dat[,"neoag_ct"])

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/num_hla_alleles_vs_neoag.pdf",width=3,height=3)
se1 <- ggplot(data = dat, aes(x = factor(hla_num), y = neoag_ct)) +
#geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
#geom_point(size=1) +
geom_boxplot() +
labs(x = "# unique HLA genes", y = "log10(neoantigen count)") +
# annotate(geom="text", size=2.5, 
# 	x=max(plot_res[,"lohhla_minor_cn"]) - (0.80*(max(plot_res[,"lohhla_minor_cn"]) - min(plot_res[,"lohhla_minor_cn"]))), 
# 	y=max(plot_res[,"minor_cn"]) - (0.97*(max(plot_res[,"minor_cn"]) - min(plot_res[,"minor_cn"]))), 
# 	label=deparse(bquote(italic("R")~" = "~.(cor1))),parse=TRUE) +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")
se1
dev.off()

