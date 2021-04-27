
library(tidyverse)
library(odbc)
library(DBI)

rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in clinical information
q <- "SELECT ps.*,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
"

dat <- dbGetQuery(con,q)
idh_glass <- rep(dat$idh_status, 2)
names(idh_glass) <- c(dat$tumor_barcode_a, dat$tumor_barcode_b)


myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("stemcell|differentiated|myeloid",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"


counts <- read.csv(myinf2,row.names=1)
cpm <- apply(counts,2, function(x) (x/sum(x))*1000000) 
cpm <- log10(cpm + 1)

info <- read.csv(myinf3,row.names=1,stringsAsFactor=FALSE)

# Glioma samples (Klemm):
glioma <- cpm[,grep("_glioma_",colnames(cpm))]
glioma <- glioma[,grep("_mdm|_mg|_cd45n",colnames(glioma))] 			####
glioma <- t(scale(t(glioma)))

# IDHwt
#idhwt_klemm <- rownames(info[which(info$IDH.Status == "wild type"),])
#idhwt_pull <- sapply(strsplit(colnames(glioma),"_"),function(x)paste(x[1:3],collapse="_"))
#idhwt_glioma <- glioma[,which(idhwt_pull %in% idhwt_klemm)]
idhwt_glioma <- glioma

# Create average myeloid and microglia signatures from Klemm
mdm <- idhwt_glioma[,grep("_mdm", colnames(idhwt_glioma))]
mdm <- apply(mdm,1,mean)
mg <- glioma[,grep("_mg", colnames(glioma))]
mg <- apply(mg,1,mean)

# Create average tumor signatures from Klemm
tumor_sub <- idhwt_glioma[,grep("_cd45n", colnames(idhwt_glioma))]
tumor <- apply(tumor_sub,1,mean)

# Extract uniquely imputed genes from CIBERSORTxGEP
gep1 <- read.delim(myinf1[1], row.names=1)
gep2 <- read.delim(myinf1[2], row.names=1)
gep3 <- read.delim(myinf1[3], row.names=1)
gep4 <- read.delim(myinf1[4], row.names=1)
colnames(gep1) <- paste("differentiated_tumor",colnames(gep1),sep="__")
colnames(gep2) <- paste("myeloid",colnames(gep2),sep="__")
#colnames(gep3) <- paste("stemcell_tumor",colnames(gep3),sep="__")
#colnames(gep4) <- paste("prolif_stemcell_tumor",colnames(gep4),sep="__")
full_gep <- data.frame(gep1,gep2)#,gep3,gep4)
full_gep <- log10(full_gep + 1)

exgenes <- apply(full_gep, 1, function(x)sum(is.na(x)))
outgenes <- names(exgenes[which(exgenes==0)])
outgep <- full_gep[-which(rownames(full_gep) %in% outgenes),]


# Convert NA to 0(?)
#full_gep <- full_gep %>%  mutate_if(is.numeric, ~replace(., is.na(.), 0))
#vars <- apply(full_gep, 1, var)
#scale_gep <- full_gep[-which(vars==0),]

scale_gep <- t(scale(t(full_gep)))

mean_cor <- low_cor <- high_cor <- flow_cell <- cell_state <- c()
for(i in 1:length(myinf1))
{
	mycell <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_","",myinf1[i])
	mycell <- gsub("_Window48.txt","",mycell)
	geps <- scale_gep[,grep(mycell,colnames(scale_gep))]
	colnames(geps) <- sapply(strsplit(colnames(geps),"__"),function(x)x[2])
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	
	geps <- geps[,names(idh_glass[which(idh_glass=="IDHwt")])]

	comxx <- intersect(rownames(idhwt_glioma), rownames(geps))
	geps <- geps[comxx,]

	sub_mdm <- mdm[rownames(geps)]
	sub_mg <- mg[rownames(geps)]
	sub_tumor <- tumor[rownames(geps)]

	mdm_cor <- apply(geps, 2, function(x)cor(x,sub_mdm,method="s",use="pairwise.complete.obs"))
	mg_cor <- apply(geps, 2, function(x)cor(x,sub_mg,method="s",use="pairwise.complete.obs"))
	tumor_cor <- apply(geps, 2, function(x)cor(x,sub_tumor,method="s",use="pairwise.complete.obs"))

	cell_state <- c(cell_state, rep(mycell, 3))
	flow_cell <- c(flow_cell, "mdm","mg","tumor")
	mean_cor <- c(mean_cor, mean(mdm_cor), mean(mg_cor), mean(tumor_cor))
	low_cor <- c(low_cor, min(mdm_cor), min(mg_cor), min(tumor_cor))
	high_cor <- c(high_cor, max(mdm_cor), max(mg_cor), max(tumor_cor))

}
idhwt_res <- data.frame(cell_state, flow_cell, mean_cor, low_cor, high_cor)
idhwt_res