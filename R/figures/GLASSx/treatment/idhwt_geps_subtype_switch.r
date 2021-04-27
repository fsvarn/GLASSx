###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
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

myinf2 <- paste("data/res/CIBERSORTx/analysis/GLASS_idhwt_", cell_state,"_postreatment_result.txt",sep="")


#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype = 'IDHwt' AND subtype_a != subtype_b AND received_treatment"

dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))

ss <- dat[,c("subtype_a","subtype_b")]
ss <- ss[-which(duplicated(ss)),]
ss <- ss[order(ss[,1],ss[,2]),]

# Start with stem cell tumor

geps <- read.delim(myinf1[3], row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]	
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))

mysig <- read.delim(myinf2[3])
mysig <- rownames(mysig[which(mysig$sig),])

subtype_switch <- paste(ss[,1],ss[,2],sep="_")
p.val <- eff <- rep(NA, length(subtype_switch))
for(i in 1:nrow(ss))
{
	sub_dat <- dat %>% filter(subtype_a == ss[i,"subtype_a"], subtype_b == ss[i,"subtype_b"])
	if(nrow(sub_dat) == 1){
		next}
	group1 <- sub_dat$tumor_barcode_a
	group2 <- sub_dat$tumor_barcode_b
	
	sub_gep1 <- geps[mysig,group1]
	sub_gep2 <- geps[mysig,group2]
	
	score1 <- apply(sub_gep1, 2, mean)
	score2 <- apply(sub_gep2, 2, mean)

	p.val[i] <- wilcox.test(score1, score2, paired=TRUE)$p.value
	eff[i] <- median(score2 - score1)
}

res <- data.frame(subtype_switch, p.val, eff)