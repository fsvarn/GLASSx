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
WITH hla_overlap AS
(
	SELECT ts.pair_barcode, ts.chrom, ts.pos, 
	CASE WHEN hla_type1 LIKE 'HLA-A%' THEN pos * int4range(29909037, 29913661)
		WHEN hla_type1 LIKE 'HLA-B%' THEN pos * int4range(31321649,31324964)
		WHEN hla_type1 LIKE 'HLA-C%' THEN pos * int4range(31236526,31239869)
	END AS overlap,
	ts.corrected_call, ts.copy_number, ts.major_cn, ts.minor_cn,
	hla_type1, hla_type2, hla_type1_copy_number, hla_type2_copy_number, pval, loss_allele, kept_allele
	FROM variants.titan_seg ts
	JOIN variants.lohhla lh ON lh.pair_barcode = ts.pair_barcode
	WHERE chrom = 6 AND pos && (
	CASE WHEN hla_type1 LIKE 'HLA-A%' THEN int4range(29909037, 29913661)
	WHEN hla_type1 LIKE 'HLA-B%' THEN int4range(31321649,31324964)
	WHEN hla_type1 LIKE 'HLA-C%' THEN int4range(31236526,31239869)
	END) AND
	coverage_filter = covfilt AND
	num_snp >= 100
	ORDER BY 1
),
weight AS
(
	SELECT pair_barcode, chrom, pos, 
	CAST(upper(ts.overlap) - lower(ts.overlap) AS numeric) AS overlap_range,
	CAST(upper(ts.overlap) - lower(ts.overlap) AS numeric) * CAST(copy_number AS numeric) AS wt_copy_number,
	CAST(upper(ts.overlap) - lower(ts.overlap) AS numeric) * CAST(major_cn AS numeric) AS wt_major_cn,
	CAST(upper(ts.overlap) - lower(ts.overlap) AS numeric) * CAST(minor_cn AS numeric) AS wt_minor_cn,
	corrected_call, copy_number, major_cn, minor_cn,
	hla_type1, hla_type2, hla_type1_copy_number, hla_type2_copy_number, pval, loss_allele, kept_allele
	FROM hla_overlap ts
	ORDER BY 1
),
weighted_hla AS
(
	SELECT pair_barcode, 
	SUM(wt_copy_number)/SUM(overlap_range) AS weighted_copy_number,
	SUM(wt_major_cn)/SUM(overlap_range) AS weighted_major_cn,
	SUM(wt_minor_cn)/SUM(overlap_range) AS weighted_minor_cn,
	hla_type1, hla_type2, hla_type1_copy_number, hla_type2_copy_number, pval, loss_allele, kept_allele
	FROM weight 
	GROUP BY pair_barcode, hla_type1, hla_type2, hla_type1_copy_number, hla_type2_copy_number, pval, loss_allele, kept_allele
	ORDER BY 1
)
SELECT * 
FROM weighted_hla wh
JOIN analysis.pairs pa ON pa.pair_barcode = wh.pair_barcode
JOIN analysis.diamond_set gs ON pa.tumor_barcode = gs.tumor_barcode_a OR pa.tumor_barcode = gs.tumor_barcode_b
"


cov_filt <- c(5,10,20,30)
res <- matrix(0,nrow=length(cov_filt),ncol=4)
rownames(res) <- cov_filt
colnames(res) <- c("major_cn_cor","minor_cn_cor","all_cor","num")
for(i in 1:length(cov_filt))
{
	newq <- gsub("covfilt",cov_filt[i],q)
	
	dat <- dbGetQuery(con,newq)

	major_cn <- minor_cn <- rep(0,nrow(dat))

	hla1_maj <- dat[,"hla_type1"] == dat[,"kept_allele"]
	major_cn[which(hla1_maj)] <- dat[which(hla1_maj),"hla_type1_copy_number"]
	major_cn[which(!hla1_maj)] <- dat[which(!hla1_maj),"hla_type2_copy_number"]
	minor_cn[which(!hla1_maj)] <- dat[which(!hla1_maj),"hla_type1_copy_number"]
	minor_cn[which(hla1_maj)] <- dat[which(hla1_maj),"hla_type2_copy_number"]


	res[i,"major_cn_cor"] <- cor(major_cn,dat[,"weighted_major_cn"],method="s")
	res[i,"minor_cn_cor"] <- cor(minor_cn,dat[,"weighted_minor_cn"],method="s")
	res[i,"all_cor"] <- cor(c(major_cn,minor_cn),c(dat[,"weighted_major_cn"],dat[,"weighted_minor_cn"]),method="s")
	
	res[i,"num"] <- nrow(dat)
	
}

#Diamond set
#    major_cn_cor minor_cn_cor   all_cor lohhla_in_titan titan_in_lohhla num
# 5  0.0723043568    0.5468275 0.3134111       0.4761905       0.8000000 394
# 10 0.0571462346    0.5460195 0.3055258       0.4651163       0.8166667 395
# 20 0.0004122647    0.5755360 0.3101671       0.4700000       0.8297872 301
# 30 0.1455366278    0.5259812 0.3666969       0.6065574       0.8378378 215

#Gold set
#    major_cn_cor minor_cn_cor   all_cor lohhla_in_titan titan_in_lohhla num
# 5    0.06237273    0.4983777 0.2945715       0.5125000       0.7439024 517
# 10   0.03728792    0.4937882 0.2802155       0.4846626       0.7594937 516
# 20  -0.02029421    0.4934318 0.2634977       0.4285714       0.7894737 414
# 30   0.16583864    0.4559259 0.3421834       0.5921053       0.8444444 284

#LOHHLA set
#    major_cn_cor minor_cn_cor   all_cor lohhla_in_titan titan_in_lohhla num
# 5    0.04908109    0.5682941 0.3081891       0.4234234       0.7966102 381
# 10   0.03092913    0.5672325 0.2956315       0.4210526       0.8275862 383
# 20   0.02558996    0.5935105 0.3252208       0.4000000       0.8163265 327
# 30   0.18656408    0.5471850 0.3897928       0.5322581       0.8461538 250

#Plot figures:

newq <- gsub("covfilt",20,q)

dat <- dbGetQuery(con,newq)

major_cn <- minor_cn <- rep(0,nrow(dat))

hla1_maj <- dat[,"hla_type1"] == dat[,"kept_allele"]
lohhla_major_cn <- lohhla_minor_cn <- rep(0,nrow(dat))
lohhla_major_cn[which(hla1_maj)] <- dat[which(hla1_maj),"hla_type1_copy_number"]
lohhla_major_cn[which(!hla1_maj)] <- dat[which(!hla1_maj),"hla_type2_copy_number"]
lohhla_minor_cn[which(!hla1_maj)] <- dat[which(!hla1_maj),"hla_type1_copy_number"]
lohhla_minor_cn[which(hla1_maj)] <- dat[which(hla1_maj),"hla_type2_copy_number"]

cor1 <- cor(lohhla_minor_cn,dat[,"major_cn"],method="s")
cor2 <- cor(lohhla_major_cn,dat[,"minor_cn"],method="s")

plot_res <- cbind(dat,lohhla_minor_cn,lohhla_major_cn)

se1 <- ggplot(data = plot_res, aes(x = factor(minor_cn), y = lohhla_minor_cn)) +
#geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
#geom_point(size=1) +
geom_boxplot() +
labs(x = "TITAN minor copy number", y = "LOHHLA minor copy number") +
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

se2 <- ggplot(data = plot_res, aes(x = factor(major_cn), y = lohhla_major_cn)) +
# geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
# geom_point(size=1) +
geom_boxplot() +
labs(x = "TITAN major copy number", y = "LOHHLA major copy number") +
# annotate(geom="text", size=2.5, 
# 	x=max(plot_res[,"lohhla_major_cn"]) - (0.80*(max(plot_res[,"lohhla_major_cn"]) - min(plot_res[,"lohhla_major_cn"]))), 
# 	y=max(plot_res[,"major_cn"]) - (0.97*(max(plot_res[,"major_cn"]) - min(plot_res[,"major_cn"]))), 
# 	label=deparse(bquote(italic("R")~" = "~.(cor2))),parse=TRUE) +
theme_bw() +
theme(axis.text = element_text(size=7),
	axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/lohhla_titan_scatterplot.pdf",width=4,height=2)
grid.arrange(se1,se2,nrow=1)
dev.off()

#Get numbers for Venn diagram

titan_loss <- dat[grep("LOH",dat[,"corrected_call"]),]
sum(titan_loss[,"pval"]<0.1 & (titan_loss[,"hla_type1_copy_number"] < 0.5 | titan_loss[,"hla_type2_copy_number"] < 0.5)) #40
sum(titan_loss[,"hla_type1_copy_number"] < 0.5 | titan_loss[,"hla_type2_copy_number"] < 0.5) #40
nrow(titan_loss) #100

hla_loss <- dat[which(dat[,"pval"]<0.1 & (dat[,"hla_type1_copy_number"] < 0.5 | dat[,"hla_type2_copy_number"] < 0.5)),]
nrow(hla_loss[grep("LOH",hla_loss[,"corrected_call"]),]) #40
nrow(hla_loss)
	
hla_loss_no_pval <- dat[which(dat[,"hla_type1_copy_number"] < 0.5 | dat[,"hla_type2_copy_number"] < 0.5),]
nrow(hla_loss_no_pval[grep("LOH",hla_loss_no_pval[,"corrected_call"]),]) #40
nrow(hla_loss_no_pval)

#Alternative method where any HLA loss in a region is compared to an LOH call for that region (as outlined in methods):

q <- "
WITH hla AS
(
	SELECT ts.pair_barcode, ts.chrom, ts.pos, ts.median_ratio, ts.median_logr, ts.corrected_call, ts.copy_number, ts.major_cn, ts.minor_cn,
	hla_type1, hla_type2, hla_type1_copy_number, hla_type2_copy_number, pval, loss_allele, kept_allele
	FROM variants.titan_seg ts
	JOIN variants.lohhla lh ON lh.pair_barcode = ts.pair_barcode
	WHERE chrom = 6 AND pos @> int4range(29909037, 31239869)
	ORDER BY 1
),
sig AS
(
	SELECT pair_barcode, chrom, pos, median_ratio, median_logr, corrected_call, copy_number, major_cn, minor_cn,
	SUM(CASE WHEN pval < 0.05 THEN 1 ELSE 0 END) AS hla_loss
	FROM hla
	GROUP BY pair_barcode, chrom, pos, median_ratio, median_logr, corrected_call, copy_number, major_cn, minor_cn
	ORDER BY 1
)
SELECT * 
FROM sig s
JOIN analysis.pairs pa ON pa.pair_barcode = s.pair_barcode
JOIN analysis.diamond_set gs ON pa.tumor_barcode = gs.tumor_barcode_a OR pa.tumor_barcode = gs.tumor_barcode_b
ORDER BY 1
"
dat <- dbGetQuery(con,q)
