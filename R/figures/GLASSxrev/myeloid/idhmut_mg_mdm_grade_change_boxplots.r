##################################################
# Examine how microglia and macrophage gene expression changes in IDH-mutant samples that increase grade
# Author: Frederick Varn
# Date: 2021.12.15
# Figure 5B
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)
library(survival)
library(topGO)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype LIKE 'IDHmut%'"


dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))

geps <- read.delim(myinf1, row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]	
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))
nrow(geps)

# Read in microglia signature
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
mg_sig <- sigs %>% filter(signature_set == "Muller", signature_name == "Microglia") %>% .$gene_symbol
mdm_sig <- sigs %>% filter(signature_set == "Muller", signature_name == "Macrophages") %>% .$gene_symbol

mg_sig <- intersect(mg_sig, rownames(geps))
mdm_sig <- intersect(mdm_sig, rownames(geps))

mg_score <- apply(geps[mg_sig,], 2, mean)
mdm_score <- apply(geps[mdm_sig,], 2, mean)

up_dat <- dat %>% filter(grade_change == "Grade up")
g1 <- up_dat$tumor_barcode_a
g2 <- up_dat$tumor_barcode_b

wilcox.test(mg_score[g1],mg_score[g2], paired=TRUE)
wilcox.test(mdm_score[g1],mdm_score[g2], paired=TRUE)
t.test(mg_score[g1],mg_score[g2], paired=TRUE)
t.test(mdm_score[g1],mdm_score[g2], paired=TRUE)


st_dat <- dat %>% filter(grade_change == "Grade stable")
g1 <- st_dat$tumor_barcode_a
g2 <- st_dat$tumor_barcode_b

wilcox.test(mg_score[g1],mg_score[g2], paired=TRUE)
wilcox.test(mdm_score[g1],mdm_score[g2], paired=TRUE)
t.test(mg_score[g1],mg_score[g2], paired=TRUE)
t.test(mdm_score[g1],mdm_score[g2], paired=TRUE)

case_barcode <- rep(dat$case_barcode, 4)
aliquot_barcode <- rep(c(dat$tumor_barcode_a, dat$tumor_barcode_b), 2)
grade_change <- rep(dat$grade_change, 4)
idh_codel_subtype <- rep(dat$idh_codel_subtype, 4) 
status <- rep(rep(c("Initial", "Recurrent"), each = nrow(dat)), 2)
score <- c(mg_score[dat$tumor_barcode_a], mg_score[dat$tumor_barcode_b], mdm_score[dat$tumor_barcode_a], mdm_score[dat$tumor_barcode_b])
cell_state <- c(rep("Microglia", nrow(dat)*2), rep("Macrophages", nrow(dat)*2))
plot_dat <- data.frame(case_barcode, aliquot_barcode, status, grade_change, idh_codel_subtype, cell_state, score)

plot_dat <- plot_dat %>%
			mutate(grade_change = fct_relevel(grade_change, "Grade up", "Grade stable"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_glass_idhmut_mdm_score.pdf",width=2.5,height=1.4)
ggplot(plot_dat %>% filter(cell_state == "Macrophages"), aes(x=status, y=score)) + 
geom_boxplot(outlier.shape = NA)  +
geom_line(size=0.5,alpha=0.4, aes(group= case_barcode)) +
geom_point(size=1,aes(colour=idh_codel_subtype)) +
facet_wrap(. ~ grade_change) + 
scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks = c(1.6, 1.7, 1.8, 1.9)) +
theme_bw() +
labs(y = "Mean expression") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7,),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") +
scale_colour_manual(values=c("#F8766D", "#00BA38"))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_glass_idhmut_mg_score.pdf",width=2.5,height=1.4)
ggplot(plot_dat %>% filter(cell_state == "Microglia"), aes(x=status, y=score)) + 
geom_boxplot(outlier.shape = NA)  +
geom_line(size=0.5,alpha=0.4, aes(group= case_barcode)) +
geom_point(size=1,aes(colour=idh_codel_subtype)) +
facet_wrap(. ~ grade_change) + 
scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
theme_bw() +
labs(y = "Mean expression") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_text(size=7,),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") +
scale_colour_manual(values=c("#F8766D", "#00BA38")) +
coord_cartesian(ylim=c(1.67, 2.1))
dev.off()
