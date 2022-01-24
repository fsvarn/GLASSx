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
myinf1 <- myinf1[grep("stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)

#  Read in post-treatment stem cell signature
myDir2 <- "data/res/CIBERSORTx/analysis/"
myinf2 <- dir(myDir2)
myinf2 <- myinf2[grep("idhmut", myinf2)]
myinf2 <- myinf2[grep("_postreatment_result", myinf2)]
mytag <- myinf2[grep("stemcell|differentiated",myinf2)]		# These are the signatures we have some confidence in
myinf2 <- paste(myDir2, mytag, sep = "/")
mytag <- gsub("_postreatment_result.txt","",mytag)
mytag <- gsub("GLASS_idhmut_","",mytag)


#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype LIKE 'IDHmut%' AND subtype_a = subtype_b AND received_treatment"


dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))


p_gu_up <- p_gs_up <- eff_gu_up <- eff_gs_up <- p_gu_dn <- p_gs_dn <- eff_gu_dn <- eff_gs_dn <- rep(0, length(myinf1))
plot_res <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	geps <- log10(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]	
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	nrow(geps)

	sigtable <- read.delim(myinf2[i])	
	sigup <- rownames(sigtable %>% filter(sig, eff > 0))
	sigdn <- rownames(sigtable %>% filter(sig, eff < 0))

	# Calculate the up and down signature scores
	sigup_score <- apply(geps[sigup,], 2, mean)
	sigdn_score <- apply(geps[sigdn,], 2, mean)
	
	# Test changes in grade up pairs
	g1 <- dat %>% filter(grade_change == "Grade up") %>% .$tumor_barcode_a
	g2 <- dat %>% filter(grade_change == "Grade up") %>% .$tumor_barcode_b
			
	p_gu_up[i] <- wilcox.test(sigup_score[g1], sigup_score[g2], paired=TRUE)$p.value
	eff_gu_up[i] <- median(sigup_score[g2] - sigup_score[g1])	

	p_gu_dn[i] <- wilcox.test(sigdn_score[g1], sigdn_score[g2], paired=TRUE)$p.value
	eff_gu_dn[i] <- median(sigdn_score[g2] - sigdn_score[g1])

	# Collect data for plotting grade up
	aliquot_barcode <- rep(c(g1, g2), 2)
	sig <- c(sigup_score[c(g1,g2)], sigdn_score[c(g1,g2)])
	direction <- c(rep("Up",length(c(g1,g2))), rep("Down", length(c(g1,g2))))
	cell_state = mytag[i]
	grade_change <- "Grade up"
	timepoint <- rep(c(rep("Initial",length(g1)), rep("Recurrent",length(g2))), 2)
	
	plot_gu <- data.frame(aliquot_barcode, sig, direction, cell_state, grade_change, timepoint)

	# Test changes in grade stable pairs
	g1 <- dat %>% filter(grade_change == "Grade stable") %>% .$tumor_barcode_a
	g2 <- dat %>% filter(grade_change == "Grade stable") %>% .$tumor_barcode_b
	
	p_gs_up[i] <- wilcox.test(sigup_score[g1], sigup_score[g2], paired=TRUE)$p.value
	eff_gs_up[i] <- median(sigup_score[g2] - sigup_score[g1])
	
	p_gs_dn[i] <- wilcox.test(sigdn_score[g1], sigdn_score[g2], paired=TRUE)$p.value
	eff_gs_dn[i] <- median(sigdn_score[g2] - sigdn_score[g1])
	
	# Collect data for plotting grade stable
	aliquot_barcode <- rep(c(g1, g2), 2)
	sig <- c(sigup_score[c(g1,g2)], sigdn_score[c(g1,g2)])
	direction <- c(rep("Up",length(c(g1,g2))), rep("Down", length(c(g1,g2))))
	cell_state = mytag[i]
	grade_change <- "Grade stable"
	timepoint <- rep(c(rep("Initial",length(g1)), rep("Recurrent",length(g2))), 2)

	plot_gs <- data.frame(aliquot_barcode, sig, direction, cell_state, grade_change, timepoint)

	# Merge data into a single plotting table
	plot_res[[i]] <- rbind(plot_gu, plot_gs)
}

res_up <- data.frame(mytag, p_gu_up, eff_gu_up, p_gs_up, eff_gs_up)
res_dn <- data.frame(mytag, p_gu_dn, eff_gu_dn, p_gs_dn, eff_gs_dn)

plot_res <- do.call(rbind, plot_res)

plot_res <- plot_res %>%
			filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor")) %>%
			mutate(cell_state = recode(cell_state, "differentiated_tumor" = "Diff.-like", "stemcell_tumor" = "Stem-like")) %>%
			mutate(status = fct_relevel(timepoint, "Initial","Recurrent")) %>%
			mutate(direction = fct_relevel(direction, "Up","Down")) %>%
			filter(direction == "Up")

# Plot it
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cibersortx_gep_idhmut_grade_change_boxplot.pdf", width=1.8,height=1.8)
ggplot(data = plot_res, aes(x = cell_state, y = sig, fill = timepoint)) +
geom_boxplot(outlier.size=1,colour="black") +
scale_fill_manual(values=c("#a6611a", "#018571")) +
facet_grid(grade_change ~ .) +
labs(y = "Signature score") +
theme_classic() +
theme(plot.title = element_text(size= 7, hjust = 0.5),
	axis.text.x = element_text(size=7, hjust = 0.5),
	axis.text.y = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none")
dev.off()


