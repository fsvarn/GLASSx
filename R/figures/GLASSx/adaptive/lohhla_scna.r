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
	an1.prop_aneuploidy AS scna_a,
	an2.prop_aneuploidy AS scna_b,
	CASE WHEN ic.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
	ic.idh_codel_subtype
	FROM any_loss al
	JOIN clinical.subtypes ic ON al.case_barcode = ic.case_barcode
	JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = al.tumor_barcode_a
	JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = al.tumor_barcode_b
	ORDER BY 1
)
SELECT * FROM full_data
"

dat <- dbGetQuery(con, q)


# Check IDHwt at both timepoints
dat1 <- dat %>%
filter(idh_status == "IDHwt")

wt_init_pval <- wilcox.test(dat1 %>% filter(loss_a == "HLA LOH") %>% .$scna_a, dat1 %>% filter(loss_a == "No HLA LOH") %>% .$scna_a)$p.value
wt_init_eff <- median(dat1 %>% filter(loss_a == "HLA LOH") %>% .$scna_a) - median(dat1 %>% filter(loss_a == "No HLA LOH") %>% .$scna_a)

wt_rec_pval <- wilcox.test(dat1 %>% filter(loss_b == "HLA LOH") %>% .$scna_b, dat1 %>% filter(loss_b == "No HLA LOH") %>% .$scna_b)$p.value
wt_rec_eff <-  median(dat1 %>% filter(loss_b == "HLA LOH") %>% .$scna_b) - median(dat1 %>% filter(loss_b == "No HLA LOH") %>% .$scna_b)

# Plot these results

timepoint <- c(rep("Initial", nrow(dat1)), rep("Recurrent", nrow(dat1)))
loss <- c(dat1$loss_a, dat1$loss_b)
scna <- c(dat1$scna_a, dat1$scna_b)

plot_dat <- data.frame(timepoint, loss, scna)
plot_dat <- plot_dat %>% mutate(loss = fct_relevel(loss, "No HLA LOH", "HLA LOH"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/hla_loss_scna_idhwt.pdf",width=2,height=1.8)
ggplot(data = plot_dat, aes(x = loss, y = scna)) +
geom_violin(aes(fill = loss)) +
geom_boxplot(width=0.1,outlier.size = 0.1, fatten = 1.5) +
scale_fill_manual(values = c("#27408B", "#CD4F39")) +
labs(y = "SCNA burden") +
facet_grid(.~timepoint, scale = "free_y") +
theme_classic() +
theme(
	axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") +
coord_cartesian(ylim = c(0,1))
dev.off()



# Check IDHmut at both timepoints
dat2 <- dat %>%
filter(idh_status == "IDHmut")

mut_init_pval <- wilcox.test(dat2 %>% filter(loss_a == "HLA LOH") %>% .$scna_a, dat2 %>% filter(loss_a == "No HLA LOH") %>% .$scna_a)$p.value
mut_init_eff <- median(dat2 %>% filter(loss_a == "HLA LOH") %>% .$scna_a) - median(dat2 %>% filter(loss_a == "No HLA LOH") %>% .$scna_a)

mut_rec_pval <- wilcox.test(dat2 %>% filter(loss_b == "HLA LOH") %>% .$scna_b, dat2 %>% filter(loss_b == "No HLA LOH") %>% .$scna_b)$p.value
mut_rec_eff <-  median(dat2 %>% filter(loss_b == "HLA LOH") %>% .$scna_b) - median(dat2 %>% filter(loss_b == "No HLA LOH") %>% .$scna_b)

# Check change over time in IDHmut
dat3 <- dat %>%
filter(loss_a == "No HLA LOH", loss_b == "HLA LOH", idh_codel_subtype == "IDHmut-noncodel") #idh_status == "IDHmut")

wilcox.test(dat3$scna_a, dat3$scna_b, paired=TRUE)
median(dat3$scna_b - dat3$scna_a)

l1 <- dat3$scna_b - dat3$scna_a

dat4 <- dat %>%
filter((loss_a == "HLA LOH" & loss_b == "No HLA LOH") | loss_a == loss_b,  idh_codel_subtype == "IDHmut-noncodel") #idh_status == "IDHmut")

wilcox.test(dat4$scna_a, dat4$scna_b,paired=TRUE)
median(dat4$scna_b - dat4$scna_a)

l2 <- dat4$scna_b - dat4$scna_a

wilcox.test(l1,l2)

# Plot the results

status <- rep("Not acquired", nrow(dat))
status[which(dat[,"loss_a"] == "No HLA LOH" & dat["loss_b"] == "HLA LOH")] <- "Acquired"
dat[,"status"] <- status
dat <- dat %>% mutate(status = fct_relevel(status, "Not acquired", "Acquired"))

plot_dat <- dat %>%
filter(idh_codel_subtype == "IDHmut-noncodel") 

plot_dat <- plot_dat %>%
pivot_longer(cols = starts_with("scna_"), names_to = "timepoint", values_to = "scna_burden") %>%
mutate(timepoint = recode(timepoint, "scna_a" = "Initial", "scna_b" = "Recurrent")) %>%
data.frame()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_hlaloh_acq.pdf",width=2,height=1.8)
ggplot(data = plot_dat, aes(x = timepoint, y = scna_burden)) +
#geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.6,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1, colour = "#00BA38") +
#scale_colour_manual(values=c("royalblue4","tomato3")) +
labs(y = "SCNA burden") +
facet_grid(.~status) +
theme_classic() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")  +
coord_cartesian(ylim=c(0,1.00))
dev.off()

# Plot differences

plot_dat2 <- dat %>%
filter(idh_codel_subtype == "IDHmut-noncodel")
plot_dat2$diff <- plot_dat2[,"scna_b"] - plot_dat2[,"scna_a"]

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_hlaloh_acq_diff.pdf",width=1.2,height=1.8)
ggplot(data = plot_dat2, aes(x = status, y = diff)) +
geom_hline(yintercept=0, colour="gray") +
geom_boxplot(outlier.size=0,colour="black") +
#geom_line(size=0.6,alpha=0.4,aes(group=case_barcode),colour= "black") +
geom_point(size=1, colour = "#00BA38") +
#scale_colour_manual(values=c("royalblue4","tomato3")) +
labs(y = "SCNA burden difference") +
facet_grid(.~idh_status) +
theme_classic() +
theme(axis.text = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")  +
coord_cartesian(ylim=c(-0.2,0.8))
dev.off()

