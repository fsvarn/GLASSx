###################################################
# Examine HLA LOH in initial vs recurrent samples across subtypes 
# Author: Frederick Varn
# Date: 2021.12.16
# Figure 6A
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(gridExtra)
library(ggalluvial)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "WITH lohhla_pairs AS
(
	SELECT ps.*,
	lh1.hla_type1, lh1.hla_type2, 
	lh1.hla_type1_copy_number AS hla_type1_copy_number_a, lh2.hla_type1_copy_number AS hla_type1_copy_number_b,
	lh1.hla_type2_copy_number AS hla_type2_copy_number_a, lh2.hla_type2_copy_number AS hla_type2_copy_number_b,
	lh1.pval AS pval_a, lh1.loss_allele AS loss_allele_a, 
	lh2.pval AS pval_b, lh2.loss_allele AS loss_allele_b
	FROM analysis.gold_set ps
	JOIN analysis.pairs pa1 ON pa1.tumor_barcode = ps.tumor_barcode_a
	JOIN analysis.pairs pa2 ON pa2.tumor_barcode = ps.tumor_barcode_b
	JOIN variants.lohhla lh1 ON lh1.pair_barcode = pa1.pair_barcode
	JOIN variants.lohhla lh2 ON lh2.pair_barcode = pa2.pair_barcode AND lh2.hla_type1 = lh1.hla_type1
	WHERE lh1.coverage_filter = 10 AND lh2.coverage_filter = 10
),
any_loss AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, 
	SUM(CASE WHEN pval_a < 0.05 AND (hla_type1_copy_number_a < 0.5 OR hla_type2_copy_number_a < 0.5) THEN 1 ELSE 0 END) AS loss_a, 
	SUM(CASE WHEN pval_b < 0.05 AND (hla_type1_copy_number_b < 0.5 OR hla_type2_copy_number_b < 0.5) THEN 1 ELSE 0 END) AS loss_b
	FROM lohhla_pairs lp
	GROUP BY case_barcode, tumor_barcode_a, tumor_barcode_b
),
full_data AS
(
	SELECT al.case_barcode, tumor_barcode_a, tumor_barcode_b,
	CASE WHEN loss_a > 0 THEN 'HLA LOH' ELSE 'No HLA LOH' END AS loss_a,
	CASE WHEN loss_b > 0 THEN 'HLA LOH' ELSE 'No HLA LOH' END AS loss_b,
	ic.idh_codel_subtype
	FROM any_loss al
	JOIN clinical.subtypes ic ON al.case_barcode = ic.case_barcode
	ORDER BY 1
)
SELECT idh_codel_subtype, loss_a, loss_b, COUNT(*)
FROM full_data
GROUP BY loss_a, loss_b, idh_codel_subtype
ORDER BY idh_codel_subtype
"

dat <- dbGetQuery(con, q)

# Number of patients experiencing LOH HLA at some time:
as.numeric(sum(dat %>% filter(loss_a == "HLA LOH" | loss_b == "HLA LOH") %>% .$count))/sum(dat$count)


dat$count <- as.numeric(dat$count)
dat$colour <- dat$idh_codel_subtype
dat[which(dat$loss_a == "No HLA LOH"),"colour"] <- "None"
dat[which(dat$loss_a == "No HLA LOH"),"colour"] <- "None"
dat[which(dat$loss_a == "No HLA LOH" & dat$loss_b == "HLA LOH"),"colour"] <- "Gain"


# Test to see if IDHmut-noncodels have a larger transition fraction
g1 <- as.numeric(sum(dat[which(dat$loss_a == "No HLA LOH" & dat$loss_b == "HLA LOH" & dat$idh_codel_subtype == "IDHmut-noncodel"),"count"]))
g2 <- as.numeric(sum(dat[which(dat$idh_codel_subtype == "IDHmut-noncodel"),"count"])) - g1
g3 <- as.numeric(sum(dat[which(dat$loss_a == "No HLA LOH" & dat$loss_b == "HLA LOH" & dat$idh_codel_subtype == "IDHwt"),"count"]))
g4 <- as.numeric(sum(dat[which(dat$idh_codel_subtype != "IDHmut-noncodel"),"count"])) - g3

ct <- matrix(c(g1,g2,g3,g4),nrow=2)
fisher.test(ct)

# Plot the fraction

bar_dat <- dat %>%
#filter(idh_codel_subtype !=  "IDHmut-codel") %>%
mutate(colour = recode(colour, "IDHmut-noncodel" = "None", "IDHwt" = "None", "IDHmut-codel" ="None")) %>%
group_by(idh_codel_subtype, colour) %>%
summarise(count = sum(count)) %>%
ungroup() %>%
group_by(idh_codel_subtype) %>%
summarise(fraction = (count/sum(count)) * 100, colour = colour) %>%
ungroup() %>%
#mutate(colour = paste(idh_codel_subtype, colour, sep ="_")) %>%
#mutate(idh_codel_subtype = fct_relevel(idh_codel_subtype, "IDHwt","IDHmut-noncodel")) %>%
#mutate(colour = fct_relevel(colour, "IDHwt_None", "IDHmut-noncodel_None", "IDHwt_Gain","IDHmut-noncodel_Gain")) %>%
mutate(idh_codel_subtype = fct_relevel(idh_codel_subtype, "IDHwt","IDHmut-noncodel", "IDHmut-codel")) %>%
mutate(colour = fct_relevel(colour, "None","Gain")) %>%
mutate(test = factor("test"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/acq_loh_idh_barplot.pdf",width=1.25,height=1.7)
ggplot(bar_dat, aes(fill=colour, y=fraction, x=idh_codel_subtype)) + 
geom_bar(position="stack", stat="identity") +
#scale_fill_manual(values=c("IDHmut-noncodel_Gain" = "#8B8B8B", "IDHmut-noncodel_None" = "#ECECEC", "IDHwt_Gain"="#8B8B8B", "IDHwt_None"="#ECECEC")) +
scale_fill_manual(values=c("Gain" = "#8B8B8B", "None" = "#ECECEC")) +
labs(y = "Proportion acquired HLA LOH (%)") +
facet_grid(.~test) +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,hjust=0.5),
axis.text.y=element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "none")
dev.off()


dat <- dat %>%
mutate(idh_codel_subtype = fct_relevel(idh_codel_subtype,"IDHwt","IDHmut-noncodel","IDHmut-codel"))

p1 <- ggplot(dat %>% filter(idh_codel_subtype == 'IDHwt'), aes(axis1 = loss_a, axis2 = loss_b, y=count)) +
geom_alluvium(aes(fill = colour), width = 1/12, knot.pos = 0.4) +
geom_stratum(width = 1/3, fill = "white", color = "grey") +
scale_x_discrete(limits = c("Initial", "Recurrent"),expand = c(0,0)) +
scale_y_continuous(expand = c(0,0)) +
scale_fill_manual(values=c("IDHwt" = "#619CFF","None" = "#D3D3D3", "Gain" = "Black")) +
ggtitle("IDHwt") +
geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
theme_bw() +
theme(axis.text.x= element_text(size=7), axis.text.y=element_blank(),
axis.ticks = element_blank(),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust = 0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 

p2 <- ggplot(dat %>% filter(idh_codel_subtype != 'IDHwt'), aes(axis1 = loss_a, axis2 = loss_b, y=count)) +
geom_alluvium(aes(fill = colour), width = 1/12, knot.pos = 0.4) +
geom_stratum(width = 1/3, fill = "white", color = "grey") +
scale_x_discrete(limits = c("Initial", "Recurrent"),expand = c(0,0)) +
scale_y_continuous(expand = c(0,0)) +
scale_fill_manual(values=c("IDHmut-noncodel" = "#00BA38","IDHmut-codel" = "#F8766D", "None" = "#D3D3D3", "Gain" = "Black")) +
facet_grid(idh_codel_subtype~.,scale="free_y",space = "free_y") +
ggtitle("IDHmut") +
geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
theme_bw() +
theme(axis.text.x= element_text(size=7), axis.text.y=element_blank(),
axis.ticks = element_blank(),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust = 0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_blank(),
strip.background = element_blank(),
legend.position="none") 


# p3 <- ggplot(dat %>% filter(idh_codel_subtype == 'IDHmut-codel'), aes(axis1 = loss_a, axis2 = loss_b, y=count)) +
# geom_alluvium(aes(fill = loss_a), width = 1/12, knot.pos = 0.4) +
# geom_stratum(width = 1/4, fill = "white", color = "grey") +
# scale_x_discrete(limits = c("Initial", "Recurrent"),expand = c(0,0)) +
# scale_y_continuous(expand = c(0,0)) +
# scale_fill_manual(values=c("#F8766D","#D3D3D3")) +
# ggtitle("IDHmut-codel") +
# geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
# theme_bw() +
# theme(axis.text.x= element_text(size=7), axis.text.y=element_blank(),
# axis.ticks = element_blank(),
# axis.title = element_blank(),
# plot.title = element_text(size=7, hjust = 0.5),
# panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
# legend.position="none") 


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/lohhla_sankey_plot_facet.pdf",width=3,height=1.7)
grid.arrange(p1, p2, ncol = 2)
dev.off()

# 
# pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/lohhla_sankey_plot_facet.pdf",width=3.5,height=2.5)
# grid.arrange(p1, p2, p3, ncol = 2,
#   layout_matrix = rbind(c(1, 2),
#                         c(1, 3))
# )
# dev.off()


# SQL script to count number of samples where LOHHLA occurs:

WITH lohhla_pairs AS
(
	SELECT ps.*,
	lh1.hla_type1, lh1.hla_type2, 
	lh1.hla_type1_copy_number AS hla_type1_copy_number_a, lh2.hla_type1_copy_number AS hla_type1_copy_number_b,
	lh1.hla_type2_copy_number AS hla_type2_copy_number_a, lh2.hla_type2_copy_number AS hla_type2_copy_number_b,
	lh1.pval AS pval_a, lh1.loss_allele AS loss_allele_a, 
	lh2.pval AS pval_b, lh2.loss_allele AS loss_allele_b
	FROM analysis.gold_set ps
	JOIN analysis.pairs pa1 ON pa1.tumor_barcode = ps.tumor_barcode_a
	JOIN analysis.pairs pa2 ON pa2.tumor_barcode = ps.tumor_barcode_b
	JOIN variants.lohhla lh1 ON lh1.pair_barcode = pa1.pair_barcode
	JOIN variants.lohhla lh2 ON lh2.pair_barcode = pa2.pair_barcode AND lh2.hla_type1 = lh1.hla_type1
	WHERE lh1.coverage_filter = 10 AND lh2.coverage_filter = 10
),
any_loss AS
(
	SELECT case_barcode, tumor_barcode_a, tumor_barcode_b, 
	SUM(CASE WHEN pval_a < 0.05 AND (hla_type1_copy_number_a < 0.5 OR hla_type2_copy_number_a < 0.5) THEN 1 ELSE 0 END) AS loss_a, 
	SUM(CASE WHEN pval_b < 0.05 AND (hla_type1_copy_number_b < 0.5 OR hla_type2_copy_number_b < 0.5) THEN 1 ELSE 0 END) AS loss_b
	FROM lohhla_pairs lp
	GROUP BY case_barcode, tumor_barcode_a, tumor_barcode_b
),
full_data AS
(
	SELECT al.case_barcode, tumor_barcode_a, tumor_barcode_b,
	CASE WHEN loss_a > 0 THEN 1 ELSE 0 END AS loss_a,
	CASE WHEN loss_b > 0 THEN 1 ELSE 0 END AS loss_b,
	ic.idh_codel_subtype
	FROM any_loss al
	JOIN clinical.subtypes ic ON al.case_barcode = ic.case_barcode
	ORDER BY 1
)
SELECT SUM(loss_a), SUM(loss_b), COUNT(*) FROM full_data

# Add loss_a sum + loss_b sum and then divide by count to get the number