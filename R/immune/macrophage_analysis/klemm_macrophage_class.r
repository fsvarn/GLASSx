###################################################
# Calculate TPM for purified macrophage profiles
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.05.29
# Author: Frederick Varn
##################################################

library(tidyverse)
library(DBI)
library(ssgsea.GBM.classification)
library(DESeq2)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

#Establish connection to db
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
ref <- dbReadTable(con, Id(schema="ref",table="genes"))

expr <- read.csv(myinf1,row.names=1)
info <- read.csv(myinf2,stringsAsFactor=FALSE)

##################################################
# Helper functions
##################################################

# TPM function
# Michael's version
# https://support.bioconductor.org/p/91218/
#-------------------------------------------

calc_tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}


# GCT creator function (for transcriptional classifier)
#-------------------------------------------

gct <- function(data, gct_out)
{
	Description <- rep(NA,nrow(data))
	data2 <- cbind(rownames(data),Description,data[,1:ncol(data)])
	colnames(data2)[1] <- "NAME"
	write.table(data2, gct_out, sep="\t", quote=FALSE, row.names=FALSE)

	conIn <- file(gct_out, "r")
	rawfile = readLines(conIn)
	close(conIn)

	mytext <- c("#1.2", paste(nrow(data2),"\t",(ncol(data)-1),sep=""),rawfile)
	conOut = file(gct_out, "w")
	writeLines(mytext, conOut)
	close(conOut)
}	
	
##################################################
# Step 1: Calculate (rough) TPM
##################################################
 
lengths <- ref[,"cds_size"]
names(lengths) <- ref[,"gene_symbol"]

overlap <- intersect(rownames(expr), names(lengths))
lengths <- lengths[overlap]
counts <- expr[overlap,]

tpm <- calc_tpm(counts,lengths)

##################################################
# Step 2: Get subtype of each CD45- cells from glioma
##################################################

glioma <- tpm[,grep("glioma",colnames(tpm))]
glioma <- glioma[,grep("cd45n",colnames(glioma))]
	
# Save glioma-specific TPM file in GCT format for classifier

cd45n_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_tpm_cd45n.gct"
gct(glioma, cd45n_path)
	
# Run Qianghu's transcriptional classifier
runSsGSEAwithPermutation(cd45n_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/p_result_BrainTIME_tpm_cd45n.gct.txt")

aliquot_barcode <- rep(rownames(subtype_ssgsea),3)
signature_name <- c(rep("Proneural",nrow(subtype_ssgsea)),rep("Classical",nrow(subtype_ssgsea)),rep("Mesenchymal",nrow(subtype_ssgsea)))
enrichment_score <- c(subtype_ssgsea[,"Proneural"],subtype_ssgsea[,"Classical"],subtype_ssgsea[,"Mesenchymal"])
p_value <- c(subtype_ssgsea[,"Proneural_pval"],subtype_ssgsea[,"Classical_pval"],subtype_ssgsea[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(aliquot_barcode,signature_name,enrichment_score,p_value)
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"aliquot_barcode"]),]

#Calculate simplicity scores and assign subtypes
aliquots <- unique(as.character(transcriptional_subtype[,"aliquot_barcode"]))

simplicity_score <- subtype_class <- rep(NA,length(aliquots))
for(i in 1:length(aliquots))
{
	sub_dat <- transcriptional_subtype[which(transcriptional_subtype[,"aliquot_barcode"] == aliquots[i]),]
	sub_dat[,"p_rank"] <- rank(sub_dat[,"p_value",],ties.method="min")
	
	subtype_class[i] <- paste(sub_dat[which(sub_dat[,"p_value"] == min(sub_dat[,"p_value"])),"signature_name"],collapse=",")
	
	r0 <- sub_dat[which(sub_dat[,"p_rank"] ==1),"p_value"][1]
	ri <- sub_dat[which(sub_dat[,"p_rank"] > 1),"p_value"]
	ri <- ri[order(ri)]
	
	adds <- sum(ri - r0)
	
	d <- abs(outer(ri,ri,"-"))
	diag(d) <- NA
	d[lower.tri(d)] <- NA
	adns <- sum(d,na.rm=TRUE)
	
	rn1 <- sub_dat[which(sub_dat[,"p_rank"] == max(sub_dat[,"p_rank"])),"p_value"][1]
	n1 <- 2	#Number of unique subtypes - 1
	simplicity_score[i] <- (adds - adns) * (rn1 - r0)/n1
}
transcriptional_subtype <- data.frame(aliquots, subtype_class, simplicity_score,stringsAsFactors=FALSE)
colnames(transcriptional_subtype) <- c("aliquot_barcode","transcriptional_subtype","simplicity_score")


##################################################
# Step 3: Identify the differentially expressed genes in macrophages from each subtype (using DESeq2)
##################################################

# Focusing on IDHwt only:
idhwt <- info[which(info[,"Tumor.Type"]=="Glioma" & info[,"IDH.Status"]=="wild type"),"Patient"]
idhwt_subtype <- transcriptional_subtype[which(gsub("_cd45n","",transcriptional_subtype[,"aliquot_barcode"]) %in% idhwt),]

subtypes <- c("Proneural","Classical","Mesenchymal")
for(i in 1:length(subtypes))
{
	sub_pt <- idhwt_subtype %>%
		filter(transcriptional_subtype == "Mesenchymal") %>%
		pull(aliquot_barcode) %>%
		gsub("_cd45n","",.)
		
	# Use DESeq2 to compare subtype versus no subtype
	dds <- DESeqDataSetFromMatrix(countData = cts,
								  colData = coldata,
								  design= ~ batch + condition)
	dds <- DESeq(dds)
	resultsNames(dds) # lists the coefficients
	res <- results(dds, name="condition_trt_vs_untrt")
	# or to shrink log fold changes association with condition:
	res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")

}