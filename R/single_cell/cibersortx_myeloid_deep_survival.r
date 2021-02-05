library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(survival)
library(topGO)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT tc.*, cc.case_age_diagnosis_years,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
WHERE idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

#Mark surgery (all happened)
dat$event <- 1
n = 4
fraction_res <- dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ diff + case_age_diagnosis_years + received_rt + received_tmz))$coefficient[(1*n) + 1])

myeloid_dat <- dat %>%
filter(cell_state == "myeloid")
summary(coxph(Surv(surgical_interval, event) ~ subtype_a, data = subtype_dat))
summary(coxph(Surv(surgical_interval, event) ~ subtype_b, data = subtype_dat))

# Make Kaplan Meier plot
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/createSurvivalFrame.r")
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/qplot_survival.r")

myeloid_dat[,"tag"] <- as.numeric(myeloid_dat$fraction_b > mean(myeloid_dat$fraction_b))

diff = survdiff(Surv(surgical_interval, event)~tag, data=myeloid_dat)
p.value = 1-pchisq(diff$chisq, length(diff$n)-1)

fit <- survfit(Surv(surgical_interval, event)~tag, data=myeloid_dat)
frame <- createSurvivalFrame(fit)
p.value <- signif(p.value, digits = 1)

se1 <- qplot_survival(frame, FALSE) +
scale_colour_manual(values = c("tomato3","royalblue4"), 
	labels = c(paste("High (n=",sum(myeloid_dat$tag==1),")",sep=""),paste("Low (n=",sum(myeloid_dat$tag==0),")",sep=""))) +
annotate("text", x=0.83*max(myeloid_dat[,"surgical_interval"],na.rm=T), y =0.93, 
	label=deparse((bquote(italic("P") ~" = " ~ .(p.value)))), parse=TRUE, size=2.5) +
labs(x = "Surgical interval (months)", y = "Probability of recurrence") +
theme_bw() +
theme(axis.text= element_text(size=7),axis.title = element_text(size=7),
	plot.title = element_text(size=7,hjust=0.5),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	legend.justification=c(1,1), legend.position=c(x=1.03,y=0.96), 
	legend.text=element_text(size=7), legend.title=element_blank(),
	legend.key.size=unit(0.5, "line"),legend.background = element_rect(fill = "transparent",colour=NA) ) +
coord_cartesian(ylim=c(0,1))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/myeloid_fraction_rec_km_plot.pdf",width=2,height=1.5)
se1
dev.off()

# See how the myeloid group changes in the bad survival group and how it doesn't in the not bad group
#--------------------------------
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

gep_list <- good_list <- bad_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	geps <- log10(geps+1)
	gep_list[[i]] <- geps

	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]

	# Comparison 1: Examine how bad survival myeloid group changes
	g1 <- myeloid_dat %>% filter(tag == 1) %>% .$tumor_barcode_a
	g2 <- myeloid_dat %>% filter(tag == 1) %>% .$tumor_barcode_b

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(geps[j,g1])
		group2 <- as.numeric(geps[j,g2])
	
		p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
		eff[j] <- log2(median(group1)/median(group2))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	bad_res <- data.frame(p.val, q.val, eff)
	bad_res <- bad_res[order(bad_res$p.val),]

	# Comparison 2: Examine how bad survival myeloid group changes
	g1 <- myeloid_dat %>% filter(tag == 0) %>% .$tumor_barcode_a
	g2 <- myeloid_dat %>% filter(tag == 0) %>% .$tumor_barcode_b

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(geps[j,g1])
		group2 <- as.numeric(geps[j,g2])
	
		p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
		eff[j] <- log2(median(group1)/median(group2))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	good_res <- data.frame(p.val, q.val, eff)
	good_res <- good_res[order(good_res$p.val),]

	sig1 <- rownames(bad_res %>% filter(q.val < 0.05 & eff > log2(1.1)))
	sig2 <- rownames(good_res %>% filter(q.val < 0.05 & eff > log2(1.1)))

	bad_list[[i]] <- bad_res
	good_list[[i]] <- good_res
}

bad_diff <- bad_list[[4]]
good_diff <- good_list[[4]]

sig1 <- rownames(bad_diff %>% filter(q.val < 0.05 & eff > log2(1.1)))
sig2 <- rownames(good_diff %>% filter(q.val < 0.05 & eff > log2(1.1)))

bad_diff <- bad_list[[2]]
good_diff <- good_list[[2]]

sig3 <- rownames(bad_diff %>% filter(q.val < 0.05 & eff < log2(0.9)))
sig4 <- rownames(good_diff %>% filter(q.val < 0.05 & eff < log2(0.9)))


##################################################
# Step 4: Perform functional enrichment analysis on the signature
##################################################

#all_genes <- rownames(gep_list[[2]])
#bg_genes <- as.numeric(all_genes %in% sig3)

all_genes <- rownames(gep_list[[4]])
bg_genes <- as.numeric(all_genes %in% sig1)
names(bg_genes) <- all_genes

diffGenes <- function(bg_genes) {
	return(bg_genes == 1)}

# Functional enrichment of signature
sampleGOdata <- new("topGOdata",
				description = "Simple session", 
				ontology = "BP",
			    allGenes = bg_genes,
			    geneSelectionFun = diffGenes, 
			    nodeSize = 10,
			    annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

# Fishers test
# resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
# fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
# fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
# 
# #Parentchild
# resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
# pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
# pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
# pcRes[which(pcRes$q.value < 0.05 & pcRes$Significant > pcRes$Expected),]

# Use this one: elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
elimRes[which(elimRes$raw.p.value < 0.05 & elimRes$Significant > elimRes$Expected),]
elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

term_level <- rev(as.character(elimRes$Term))
plotRes <- elimRes %>% 
		   mutate(Term = as_factor(Term)) %>%
		   mutate(Term = fct_relevel(Term, term_level)) %>%
		   filter(q.value < 0.05 & elimRes$Significant > elimRes$Expected)








gep_list <- diff_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	geps <- log10(geps+1)
	gep_list[[i]] <- geps

	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]

	g1 <- myeloid_dat %>% filter(tag == 1) %>% .$tumor_barcode_b
	g2 <- myeloid_dat %>% filter(tag == 0) %>% .$tumor_barcode_b

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(geps[j,g1])
		group2 <- as.numeric(geps[j,g2])
		batch <- substr()
	
		p.val[j] <- wilcox.test(group1,group2)$p.value
		eff[j] <- log2(median(group1)/median(group2))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	res <- data.frame(p.val, q.val, eff)
	res <- res[order(res$p.val),]
	
	diff_list[[i]] <- res
}

diff_list[[4]] %>% filter(q.val < 0.05)
test <- rownames(diff_list[[4]] %>% filter(q.val < 0.05 & eff > log2(1.1)))

all_genes <- rownames(gep_list[[4]])
bg_genes <- as.numeric(all_genes %in% test)
names(bg_genes) <- all_genes

diffGenes <- function(bg_genes) {
	return(bg_genes == 1)}

# Functional enrichment of signature
sampleGOdata <- new("topGOdata",
				description = "Simple session", 
				ontology = "BP",
			    allGenes = bg_genes,
			    geneSelectionFun = diffGenes, 
			    nodeSize = 10,
			    annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

# Use this one: elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
elimRes[which(elimRes$raw.p.value < 0.05 & elimRes$Significant > elimRes$Expected),]
elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

term_level <- rev(as.character(elimRes$Term))
plotRes <- elimRes %>% 
		   mutate(Term = as_factor(Term)) %>%
		   mutate(Term = fct_relevel(Term, term_level)) %>%
		   filter(q.value < 0.05 & elimRes$Significant > elimRes$Expected)




q <- "SELECT tc.*, 
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
"

dat <- dbGetQuery(con,q)


dat <- dat %>% mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligo.",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		mutate(cell_state = fct_relevel(cell_state, "Prolif. stem-like", "Stem-like","Diff.-like",
										"Fibroblast", "Pericyte","Endothelial", "Oligo.",
										"Myeloid", "Dendritic cell", "T cell", "Granulocyte","B cell"))

# T/myeloid ratio
t_my <- dat %>%
	   group_by(case_barcode) %>%
	   summarise(
	   surgical_interval = surgical_interval[cell_state=="T cell"], idh_status = idh_status[cell_state == "T cell"],
	   case_age_diagnosis_years = case_age_diagnosis_years[cell_state=="T cell"],
	   case_overall_survival_mo = case_overall_survival_mo[cell_state=="T cell"],
	   case_vital_status = case_vital_status[cell_state=="T cell"],
	   received_treatment = received_treatment[cell_state == "T cell"],
	   ratio_a = (fraction_a[cell_state == "T cell"]) / (fraction_a[cell_state == "Myeloid"] + 0.01),
	   ratio_b = (fraction_b[cell_state == "T cell"]) / (fraction_b[cell_state == "Myeloid"] +0.01))	%>%
	   mutate(ratio_change = ratio_b - ratio_a, event = 1, case_vital_status = recode(case_vital_status, 'dead' = 1, 'alive' = 0))

summary(coxph(Surv(surgical_interval, event) ~ ratio_change + idh_status + case_age_diagnosis_years + received_treatment, data = t_my))		# P = 0.1
summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ ratio_change + idh_status + case_age_diagnosis_years, data = t_my))	# P = 0.9


cor.test(t_my$surgical_interval, t_my$ratio_change, method="p")
cor.test(t_my$surgical_interval, t_my$ratio_change, method="s")

cor.test(t_my$surgical_interval, t_my$ratio_a, method="s")
cor.test(t_my$surgical_interval, t_my$ratio_b, method="s")


normal_dat <- dat %>% 
		filter(!grepl("-like", cell_state)) %>%
		filter(!is.na(surgical_interval)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, 
			idh_status, surgical_interval, received_treatment, case_age_diagnosis_years) %>%
	   	mutate(fraction_diff = fraction_b - fraction_a, event = 1) %>%
	   	ungroup()
   	
normal_res <- normal_dat %>%
		filter(idh_status == "IDHwt") %>%
		group_by(cell_state) %>%
		summarise(cor_a = cor(fraction_a, surgical_interval, method="s"),
				  cor_b = cor(fraction_b, surgical_interval, method="s"),
				  diff = cor(fraction_diff, surgical_interval, method="s"))

t_my <- normal_dat %>%
	   group_by(case_barcode) %>%
	   filter(idh_status == "IDHwt") %>%	   
	   summarise(
	   surgical_interval = surgical_interval[cell_state=="T cell"], idh_status = idh_status[cell_state == "T cell"],
	   case_age_diagnosis_years = case_age_diagnosis_years[cell_state=="T cell"],
	   received_treatment = received_treatment[cell_state == "T cell"],
	   ratio_a = (fraction_a[cell_state == "T cell"]) / (fraction_a[cell_state == "Myeloid"] + 0.01),
	   ratio_b = (fraction_b[cell_state == "T cell"]) / (fraction_b[cell_state == "Myeloid"] +0.01))	%>%
	   mutate(ratio_change = ratio_b - ratio_a, event = 1)

n = 4
normal_cox <- normal_dat %>%
group_by(cell_state) %>%
summarise(init_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  init_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_a + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  rec_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  rec_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_b + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1],
		  diff_pval = summary(coxph(Surv(surgical_interval, event) ~ fraction_diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(4*n) + 1],
		  diff_hr = summary(coxph(Surv(surgical_interval, event) ~ fraction_diff + idh_status + case_age_diagnosis_years + received_treatment))$coefficient[(1*n) + 1])

summary(coxph(Surv(surgical_interval, event) ~ subtype_a, data = subtype_dat))
summary(coxph(Surv(surgical_interval, event) ~ subtype_b, data = subtype_dat))
