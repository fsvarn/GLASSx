###################################################
# Test whether the stem-like recurrence signature associates with survival
# Author: Frederick Varn
# Date: 2022.01.21
# Reported in text, no figure
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(survival)

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

#  Post-treatment signature
myinf2 <- paste("data/res/CIBERSORTx/analysis/GLASS_idhwt_", cell_state,"_postreatment_result.txt",sep="")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in IDHwt data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ag1.signature_name AS subtype_a, 
ag2.signature_name AS subtype_b,
cs.idh_codel_subtype,
ca.case_overall_survival_mo,
ca.case_vital_status,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ss
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
JOIN clinical.cases ca ON ca.case_barcode = ss.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con, q)


sig_score <- list()
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


	mysig <- read.delim(myinf2[i])
	mysig <- rownames(mysig %>% filter(sig, eff > 0))
		
	sig_score[[i]] <- apply(geps[mysig,], 2, mean)
}

aliquot_barcode <- names(sig_score[[1]])
sig_mat <- data.frame(sig_score[[1]], sig_score[[2]], sig_score[[3]])
colnames(sig_mat) <- cell_state

# Make a table with these scores
#write.table(sig_mat, "data/res/CIBERSORTx/analysis/idhwt_signature_scores.txt", sep = "\t", row.names=TRUE, quote=FALSE)

sub_dat <- dat %>% filter()
sig_mat_a <- sig_mat[dat[,"tumor_barcode_a"],]
colnames(sig_mat_a) <- paste(colnames(sig_mat_a), "_a", sep="")
dat <- cbind(dat,sig_mat_a)

sig_mat_b <- sig_mat[dat[,"tumor_barcode_b"],]
colnames(sig_mat_b) <- paste(colnames(sig_mat_b), "_b", sep="")
dat <- cbind(dat,sig_mat_b)

dat <- dat %>%
mutate(case_vital_status = recode(case_vital_status, 'alive' = 0, 'dead' = 1))

dat$stemcell_diff <- dat$stemcell_tumor_b - dat$stemcell_tumor_a
mytag <- rep(0, nrow(dat))
mytag[which(dat$stemcell_diff <= 0)] <- 1
mytag[which(dat$stemcell_diff > 0)] <- 2
dat$mytag <- mytag

summary(coxph(Surv(case_overall_survival_mo, case_vital_status)~mytag + , data=dat)) # P = 0.8
