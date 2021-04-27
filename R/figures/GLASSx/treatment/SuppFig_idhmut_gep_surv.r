###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHmut tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)
library(survival)
library(topGO)

#######################################################rm(list=ls())
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

sigup_score <- sigdn_score <- list()
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
	sigup_score[[i]] <- apply(geps[sigup,], 2, mean)
	sigdn_score[[i]] <- apply(geps[sigdn,], 2, mean)
}

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype LIKE 'IDHmut%' AND subtype_a = subtype_b AND received_treatment"


dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))
		
dat$recur_status <- 1
	   
# Diff-like plot
#Kaplan-Meier plots
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/createSurvivalFrame.r")
source("/projects/verhaak-lab/USERS/varnf/SofWar/R/createSurvivalFrame.r")

diff_up <- sigup_score[[1]]
diff_dn <- sigdn_score[[1]]

uptag <- diff_up[dat$tumor_barcode_b] - diff_up[dat$tumor_barcode_a]
dntag <- diff_dn[dat$tumor_barcode_b] - diff_dn[dat$tumor_barcode_a]
fulltag <- uptag - dntag

thr <- median(fulltag)
fulltag <- fulltag > thr

diff_dat <- dat
diff_dat$tag <- fulltag

fit <- survfit(Surv(surgical_interval, recur_status)~tag, data=diff_dat)
frame <- createSurvivalFrame(fit)

diff = survdiff(Surv(surgical_interval, recur_status)~tag, data=diff_dat)

summary(coxph(Surv(surgical_interval, recur_status)~tag, data=diff_dat))





