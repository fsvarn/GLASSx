###################################################
# Test frequency at which HLA LOH is copy neutral
# Author: Frederick Varn
# Date: 2021.10.29
# Revision comment
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

q <- "
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
loss_status AS
(
	SELECT lp.case_barcode, tumor_barcode_a, tumor_barcode_b, 
	CASE WHEN pval_a < 0.05 AND (((hla_type1_copy_number_a < 0.5 AND hla_type2_copy_number_a > 1.5) OR(hla_type2_copy_number_a < 0.5 AND hla_type1_copy_number_a > 1.5))) THEN 'copy neutral'
			 WHEN pval_a < 0.05 AND (((hla_type1_copy_number_a < 0.5 AND hla_type2_copy_number_a <= 1.5) OR(hla_type2_copy_number_a < 0.5 AND hla_type1_copy_number_a <= 1.5))) THEN 'deletion'
			 ELSE 'None' END AS loss_a, 
	CASE WHEN pval_a < 0.05 AND (((hla_type1_copy_number_b < 0.5 AND hla_type2_copy_number_b > 1.5) OR(hla_type2_copy_number_b < 0.5 AND hla_type1_copy_number_b > 1.5))) THEN 'copy neutral'
			 WHEN pval_a < 0.05 AND (((hla_type1_copy_number_b < 0.5 AND hla_type2_copy_number_b <= 1.5) OR(hla_type2_copy_number_b < 0.5 AND hla_type1_copy_number_b <= 1.5))) THEN 'deletion'
			 ELSE 'None' END AS loss_b,
	ic.idh_codel_subtype
	FROM lohhla_pairs lp
    JOIN clinical.subtypes ic ON lp.case_barcode = ic.case_barcode
	
)
SELECT 	CASE WHEN idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
loss_a, loss_b, COUNT(*) FROM loss_status GROUP BY loss_a, loss_b, idh_status
"
dat <- dbGetQuery(con, q)

dat_a <- dat %>%
		 filter(!(loss_a == "None")) %>%
		 group_by(idh_status, loss_a) %>%
		 summarise(total = sum(count)) %>%
		 ungroup() %>%
		 group_by(idh_status) %>%
		 summarise(loss = loss_a, total = total, fraction = total/sum(total))%>%
		 mutate(timepoint = "Initial")
		 
dat_b <- dat %>%
		 filter(!(loss_b == "None")) %>%
		 group_by(idh_status, loss_b) %>%
		 summarise(total = sum(count)) %>%
		 ungroup() %>%
		 group_by(idh_status) %>%
		 summarise(loss = loss_b, total = total, fraction = total/sum(total)) %>%
		 mutate(timepoint = "Recurrent")
		 
		
bar_dat <- rbind(dat_a, dat_b)
bar_dat <- bar_dat %>%
mutate(idh_status = fct_relevel(idh_status, "IDHwt", "IDHmut"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/loh_copy_neutral_del.pdf",width=1.5,height=1.7)
ggplot(bar_dat, aes(fill=loss, y=fraction, x=timepoint)) + 
geom_bar(position="stack", stat="identity") +
#scale_fill_manual(values=c("IDHmut-noncodel_Gain" = "#8B8B8B", "IDHmut-noncodel_None" = "#ECECEC", "IDHwt_Gain"="#8B8B8B", "IDHwt_None"="#ECECEC")) +
scale_fill_manual(values=c("deletion" = "#8B8B8B", "copy neutral" = "#ECECEC")) +
labs(y = "Proportion (%)") +
facet_grid(.~idh_status) +
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
