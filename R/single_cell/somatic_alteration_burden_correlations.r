###################################################
# Create stacked barplots (transcriptional classifier/simplicity score/CIBERSORTx) of each GLASS sample
# Updated: 2020.07.06
# Author: Frederick Varn
##################################################
library(tidyverse)
library(odbc)
library(DBI)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")


#######################################################
# Step 1: Mutation burden changes
#######################################################

q <- "
SELECT ps.*, mf1.coverage_adj_mut_freq AS mut_freq_a, mf2.coverage_adj_mut_freq AS mut_freq_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b, 
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' WHEN cs.idh_codel_subtype = 'IDHmut-noncodel' THEN 'IDHmut' WHEN cs.idh_codel_subtype = 'IDHmut-codel' THEN 'IDHmut' END AS idh_status
FROM analysis.platinum_set ps
JOIN analysis.mut_freq mf1 ON ps.dna_barcode_a = mf1.aliquot_barcode
JOIN analysis.mut_freq mf2 ON ps.dna_barcode_b = mf2.aliquot_barcode
JOIN analysis.cibersortx_scgp cs1 ON ps.rna_barcode_a = cs1.aliquot_barcode
JOIN analysis.cibersortx_scgp cs2 ON ps.rna_barcode_b = cs2.aliquot_barcode AND cs1.cell_state = cs2.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con, q)

cor_res <- dat %>% 
group_by(cell_state, idh_status) %>%
summarise(cor_init = cor(fraction_a, mut_freq_a, method = "s"), 
		  cor_rec = cor(fraction_b, mut_freq_b, method = "s"), 
		  cor_diff = cor(fraction_b - fraction_a, mut_freq_b - mut_freq_a, method ="s"),
		  cor_diff_pval = cor.test(fraction_b - fraction_a, mut_freq_b - mut_freq_a, method ="s")$p.value) %>%
data.frame()

hm_pair <- dat %>%
filter(mut_freq_b > 100) %>%
group_by(cell_state) %>%
summarise(hm_pval = wilcox.test(fraction_b, fraction_a, paired=TRUE)$p.value,
		  hm_eff = median(fraction_b - fraction_a))

# Recurrent hypermutants vs recurrent non-hypermutants
cells <- unique(dat[,"cell_state"])
hm_p_val <- hm_eff <- rep(0, length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat %>% filter(cell_state == cells[i])
	hm_p_val[i] <- wilcox.test(sub_dat[which(sub_dat[,"mut_freq_b"] > 100), "fraction_b"], sub_dat[which(sub_dat[,"mut_freq_b"] > 100), "fraction_a"])$p.value
	hm_eff[i] <- median(sub_dat[which(sub_dat[,"mut_freq_b"] > 100), "fraction_b"]) - median(sub_dat[which(sub_dat[,"mut_freq_b"] > 100), "fraction_a"])
}
hm_rec <- data.frame(cells, hm_p_val, hm_eff)

cells <- unique(dat[,"cell_state"])
conc_p_val <- inc_med <- dec_med <- rep(0, length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat %>% 
	filter(cell_state == cells[i], idh_status == 'IDHwt') %>%
	mutate(mut_diff = mut_freq_b - mut_freq_a) %>%
	mutate(fraction_diff = fraction_b - fraction_a)
	
	g1 <- sub_dat[which(sub_dat[,"mut_diff"] >= 0), "fraction_diff"]
	g2 <- sub_dat[which(sub_dat[,"mut_diff"] < 0), "fraction_diff"]
	
	conc_p_val[i] <- wilcox.test(g1,g2)$p.value
	inc_med[i] <- median(g1)
	dec_med[i] <- median(g2)
}
conc_res <- data.frame(cells, conc_p_val, inc_med, dec_med)


# Plot the one significant association in the difference (ranked barplot)
plot_res <- dat %>% 
			filter(idh_status == 'IDHwt', cell_state == 'prolif_stemcell_tumor') %>%
			arrange(mut_freq_b - mut_freq_a) %>%
			mutate(fraction_dif = fraction_b - fraction_a)
plot_res <- plot_res %>%
			mutate(mut_rank = 1:nrow(plot_res))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/mutation_diff_prolif_stemcell_barplot.pdf",width=5,height=3.5)  
ggplot(plot_res, aes(x=mut_rank,y=fraction_dif)) +
	geom_bar(stat = "identity") +
	labs(title = "Proliferating stem cell", x = "Mutation difference rank", y = "Change in fraction") +
	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.title=element_text(size=7,hjust=0.5),
	axis.title=element_text(size=7),
	axis.text.x=element_text(size=7),
	axis.text.y=element_text(size=7),
	legend.position="none") 			
dev.off()		

# Plot the one significant association in the difference (boxplot)
plot_res <- dat %>% 
			filter(idh_status == 'IDHwt', cell_state == 'prolif_stemcell_tumor') %>%
			mutate(mut_freq_change = ifelse((mut_freq_b - mut_freq_a) > 0, "Increase", "Decrease")) %>%
			mutate(fraction_diff = (fraction_b - fraction_a) * 100) 

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/mutation_diff_prolif_stemcell_boxplot.pdf",width=1.5,height=2)  
ggplot(plot_res, aes(x=mut_freq_change,y=fraction_diff)) +
	geom_boxplot() +
	labs(title = "Proliferating stem cell", y = "Change in fraction (%)") +
	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.title=element_text(size=7,hjust=0.5),
	axis.title.y=element_text(size=7),
	axis.title.x=element_blank(),
	axis.text.x=element_text(size=7),
	axis.text.y=element_text(size=7),
	legend.position="none") 			
dev.off()		



#######################################################
# Step 2: Copy number load changes
#######################################################

#Read in data for analysis
q <- "
SELECT ps.*, an1.prop_aneuploidy AS prop_aneuploidy_a, an2.prop_aneuploidy AS prop_aneuploidy_b, cs1.cell_state, cs1.fraction AS fraction_a, cs2.fraction AS fraction_b, 
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' WHEN cs.idh_codel_subtype = 'IDHmut-noncodel' THEN 'IDHmut' WHEN cs.idh_codel_subtype = 'IDHmut-codel' THEN 'IDHmut' END AS idh_status
FROM analysis.platinum_set ps
JOIN analysis.gatk_aneuploidy an1 ON ps.dna_barcode_a = an1.aliquot_barcode
JOIN analysis.gatk_aneuploidy an2 ON ps.dna_barcode_b = an2.aliquot_barcode
JOIN analysis.cibersortx_scgp cs1 ON ps.rna_barcode_a = cs1.aliquot_barcode
JOIN analysis.cibersortx_scgp cs2 ON ps.rna_barcode_b = cs2.aliquot_barcode AND cs1.cell_state = cs2.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con, q)


cor_res <- dat %>% 
group_by(cell_state, idh_status) %>%
summarise(cor_init = cor(fraction_a, prop_aneuploidy_a, method = "s"), 
		  cor_rec = cor(fraction_b, prop_aneuploidy_b, method = "s"), 
		  cor_diff = cor(fraction_b - fraction_a, prop_aneuploidy_b - prop_aneuploidy_a, method ="s"),
		  cor_diff_pval = cor.test(fraction_b - fraction_a, prop_aneuploidy_b - prop_aneuploidy_a, method ="s")$p.value) %>%
data.frame()

# Plot the one significant association in the difference (boxplot)
plot_res <- dat %>% 
			filter(idh_status == 'IDHwt', cell_state == 'myeloid') %>%
			mutate(prop_aneuploidy_change = ifelse((prop_aneuploidy_b - prop_aneuploidy_a) > 0, "Increase", "Decrease")) %>%
			mutate(fraction_diff = (fraction_b - fraction_a) * 100) 

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/aneuploidy_diff_myeloid_boxplot.pdf",width=1.5,height=2)  
ggplot(plot_res, aes(x=prop_aneuploidy_change,y=fraction_diff)) +
	geom_boxplot() +
	labs(title = "Myeloid", y = "Change in fraction (%)") +
	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.title=element_text(size=7,hjust=0.5),
	axis.title.y=element_text(size=7),
	axis.title.x=element_blank(),
	axis.text.x=element_text(size=7),
	axis.text.y=element_text(size=7),
	legend.position="none") 			
dev.off()		

