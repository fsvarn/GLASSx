library(tidyverse)
library(odbc)
library(DBI)
library(gridExtra)
library(ggalluvial)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# LOHHLA analysis with cell state changes
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
	CASE WHEN ic.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
	ic.idh_codel_subtype
	FROM any_loss al
	JOIN clinical.subtypes ic ON al.case_barcode = ic.case_barcode
	ORDER BY 1
)
SELECT ps.*, lh.loss_a, lh.loss_b, lh.idh_status, lh.idh_codel_subtype, sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS fraction_diff
FROM analysis.platinum_set ps
JOIN full_data lh ON lh.tumor_barcode_a = ps.dna_barcode_a AND lh.tumor_barcode_b = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ps.rna_barcode_b AND sc1.cell_state = sc2.cell_state
"

dat <- dbGetQuery(con, q)


cell_states <- unique(dat$cell_state)
wt_p.val_a <- wt_p.val_b <- mut_p.val_a <- mut_p.val_b <- rep(0, length(cell_states))
for(i in 1:length(cell_states))
{
	g1 <- dat %>% filter(loss_a == "No HLA LOH", cell_state == cell_states[i], idh_status == "IDHwt")  %>% .$fraction_a
	g2 <- dat %>% filter(loss_a == "HLA LOH", cell_state == cell_states[i], idh_status == "IDHwt")  %>% .$fraction_a
	wt_p.val_a[i] <- wilcox.test(g1,g2)$p.value

	g1 <- dat %>% filter(loss_b == "No HLA LOH", cell_state == cell_states[i], idh_status == "IDHwt")  %>% .$fraction_b
	g2 <- dat %>% filter(loss_b == "HLA LOH", cell_state == cell_states[i], idh_status == "IDHwt")  %>% .$fraction_b
	wt_p.val_b[i] <- wilcox.test(g1,g2)$p.value

	g1 <- dat %>% filter(loss_a == "No HLA LOH", cell_state == cell_states[i], idh_status == "IDHmut")  %>% .$fraction_a
	g2 <- dat %>% filter(loss_a == "HLA LOH", cell_state == cell_states[i], idh_status == "IDHmut")  %>% .$fraction_a
	mut_p.val_a[i] <- wilcox.test(g1,g2)$p.value

	g1 <- dat %>% filter(loss_b == "No HLA LOH", cell_state == cell_states[i], idh_status == "IDHmut")  %>% .$fraction_b
	g2 <- dat %>% filter(loss_b == "HLA LOH", cell_state == cell_states[i], idh_status == "IDHmut")  %>% .$fraction_b
	mut_p.val_b[i] <- wilcox.test(g1,g2)$p.value
}

pvals <- data.frame(cell_states, wt_p.val_a, wt_p.val_b, mut_p.val_a, mut_p.val_b)


loh_status <- c(dat$loss_a, dat$loss_b)
cell_state <- rep(dat$cell_state, 2)
fraction <- c(dat$fraction_a, dat$fraction_b)
timepoint <- c(rep("Initial", nrow(dat)), rep("Recurrent", nrow(dat)))
idh_status <- rep(dat$idh_status, 2)
idh_codel_subtype <- rep(dat$idh_codel_subtype, 2)

plot_dat <- data.frame(loh_status, cell_state, fraction, timepoint, idh_status, idh_codel_subtype)
plot_dat <-  plot_dat %>% 
			 mutate(idh_status = fct_relevel(idh_status, "IDHwt", "IDHmut")) %>%
			 mutate(loh_status = fct_relevel(loh_status, "No HLA LOH", "HLA LOH"))

test_dat1 <- plot_dat %>% filter(cell_state == "t_cell", idh_status == "IDHwt")
g1 <- test_dat1[which(test_dat1$loh_status == "HLA LOH"),"fraction"]
g2 <- test_dat1[which(test_dat1$loh_status == "No HLA LOH"),"fraction"]
wilcox.test(g1,g2)

test_dat2 <- plot_dat %>% filter(cell_state == "t_cell", idh_status == "IDHmut")
g1 <- test_dat2[which(test_dat2$loh_status == "HLA LOH"),"fraction"]
g2 <- test_dat2[which(test_dat2$loh_status == "No HLA LOH"),"fraction"]
wilcox.test(g1,g2)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/hla_loss_t_cells.pdf",width=2,height=1.8)
ggplot(data = plot_dat %>% filter(cell_state=="t_cell"), aes(x = loh_status, y = fraction * 100)) +
geom_violin(aes(fill = loh_status)) +
geom_boxplot(width=0.1,outlier.size = 0.1, fatten = 1.5) +
scale_fill_manual(values = c("#27408B", "#CD4F39")) +
labs(y = "T cell fraction (%)") +
facet_grid(.~idh_status, scale = "free_y") +
theme_classic() +
theme(
	axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim = c(0,7))
dev.off()




# Supplemental examination of neoantigen depletion with cell state analysis

q <- "
SELECT ps.*, nd1.rneo AS rneo_a, nd2.rneo AS rneo_b, nd2.rneo - nd1.rneo AS rneo_diff, sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS fraction_diff,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.platinum_set ps
JOIN analysis.neoantigen_depletion nd1 ON nd1.aliquot_barcode = ps.dna_barcode_a 
JOIN analysis.neoantigen_depletion nd2 ON nd2.aliquot_barcode = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ps.rna_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"


dat <- dbGetQuery(con, q)

res <- dat %>%
group_by(idh_status, cell_state) %>%
summarise(cor = cor(fraction_diff, rneo_diff,method="s")) %>%
data.frame()

res <- dat %>%
group_by(idh_status, cell_state) %>%
summarise(cor = cor(fraction_b, rneo_b, method="s")) %>%
data.frame()
