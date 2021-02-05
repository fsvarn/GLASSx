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
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
WHERE sc1.cell_state = 'myeloid' AND idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con, q)

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-combat_20200106_scgp/"
#myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")


gep_list <- good_list <- bad_list <- list()

geps <- read.delim(myinf1[2], row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
geps <- log10(geps+1)

rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
all_genes <- rownames(geps)

myinf2 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-combat_20200106_scgp/CIBERSORTxGEP_GLASS_Fractions-Adjusted.txt"
#myinf2 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/CIBERSORTxGEP_GLASS_Fractions-Adjusted.txt"
fract <- read.delim(myinf2,row.names=1)
rownames(fract) <- gsub("\\.","-",rownames(fract))

geps <- t(geps)

geps <- geps[dat$tumor_barcode_b,]
fract <- fract[dat$tumor_barcode_b,]

cors <- apply(geps, 2, function(x)cor(x, fract$myeloid, method="s"))
p.val <- apply(geps, 2, function(x)cor.test(x, fract$myeloid, method="s")$p.value)
q.val <- p.adjust(p.val,"BH")

res <- data.frame(cors,p.val,q.val)
res <- res[order(cors,decreasing=TRUE),]

# Read in macrophage modules
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"
modules <- read.delim(myinf3,stringsAsFactors=FALSE,row.names=1)
modules <- as.list(modules)
modules <- sapply(modules, function(x) x[-which(x=="")])

mod_res <- sapply(modules, function(x) mean(res[x,"cors"],na.rm=TRUE))
mod_res <- mod_res[order(mod_res, decreasing=TRUE)]

il4 <- c("X10","X11","X12","X13","X14","X15","X37","X45","X46","X48")

#how to visualize?