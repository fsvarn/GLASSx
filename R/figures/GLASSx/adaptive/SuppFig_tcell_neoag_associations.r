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
SELECT ps.*, nd1.rneo AS rneo_a, nd2.rneo AS rneo_b, nd2.rneo - nd1.rneo AS rneo_diff, 
nd1.bobs AS neoag_a, nd2.bobs AS neoag_b,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS fraction_diff,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.platinum_set ps
JOIN analysis.neoantigen_depletion nd1 ON nd1.aliquot_barcode = ps.dna_barcode_a 
JOIN analysis.neoantigen_depletion nd2 ON nd2.aliquot_barcode = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ps.rna_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"


dat <- dbGetQuery(con, q)

dat %>%
filter(cell_state == "t_cell") %>%
group_by(idh_status) %>%
summarise(cor_a = cor(fraction_a, rneo_a, method="p"), pval_a = cor.test(fraction_a, rneo_a, method="p")$p.value,
cor_b = cor(fraction_b, rneo_b, method="p"),  pval_b = cor.test(fraction_b, rneo_b, method="p")$p.value) %>%
data.frame()

dat <- dat %>% filter(cell_state == "t_cell")

rneo <- c(dat$rneo_a, dat$rneo_b)
fraction <- c(dat$fraction_a, dat$fraction_b)
timepoint <- c(rep("Initial",nrow(dat)), rep("Recurrent", nrow(dat)))
idh_status <- rep(dat$idh_status, 2)
case_barcode <- rep(dat$case_barcode, 2)
plot_dat <- data.frame(case_barcode, timepoint, rneo, fraction, idh_status)

plot_dat <- plot_dat %>%
			mutate(idh_status = fct_relevel(idh_status, "IDHwt", "IDHmut"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/neoag_depletion_t_cells.pdf",width=2.3,height=2.2)
ggplot(plot_dat, aes(x = fraction*100, y = rneo)) + 
geom_point()  +
geom_smooth(method="lm", se = FALSE, fullrange=TRUE) +
facet_grid(timepoint ~ idh_status, scale = "free_x") +
theme_bw() +
labs(x="T cell proportion (%)", y = "Neoantigen depletion rate") +
#scale_x_continuous(breaks = breaks_fun, limits = c(0, NA)) +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") 
dev.off()

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT ps.*, nd1.coverage_adj_mut_freq AS freq_a, nd2.coverage_adj_mut_freq AS freq_b,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS fraction_diff,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.platinum_set ps
JOIN analysis.mut_freq nd1 ON nd1.aliquot_barcode = ps.dna_barcode_a 
JOIN analysis.mut_freq nd2 ON nd2.aliquot_barcode = ps.dna_barcode_b
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ps.rna_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"


dat <- dbGetQuery(con, q)

dat %>%
filter(cell_state == "t_cell") %>%
group_by(idh_status) %>%
summarise(cor_a = cor(fraction_a, freq_a, method="s"), pval_a = cor.test(fraction_a, freq_a, method="s")$p.value,
cor_b = cor(fraction_b, freq_b, method="s"),  pval_b = cor.test(fraction_b, freq_b, method="s")$p.value) %>%
data.frame()