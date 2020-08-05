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
library(GSVA)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

#Establish connection to db
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
ref <- dbReadTable(con, Id(schema="ref",table="genes"))

expr <- read.csv(myinf1,row.names=1)
info <- read.csv(myinf2,row.names=1,stringsAsFactor=FALSE)

# Add a batch column to info
info[,"batch"] <- sapply(strsplit(rownames(info),"_"),function(x)x[1])

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

# overrepresentation test of macrophage modules
#-------------------------------------------

or_test <- function(dif, sig, bg)
{
	# Differentially expressed and in signature
	g1 <- length(bg[which(bg %in% dif & bg %in% sig)])
	# Differentially expressed and not in signature
	g2 <- length(bg[which(bg %in% dif & !(bg %in% sig))])
	# Not differentially expressed, in signature
	g3 <- length(bg[which(!(bg %in% dif) & bg %in% sig)])
	# Not differentially expressed, not in signature
	g4 <- length(bg[which(!(bg %in% dif) & !(bg %in% sig))])
	
	p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2))$p.value
	
	return(p.val)
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

#only use genes expressed in at least half of the cohort
sums <- apply(tpm,1,function(x)sum(x>0))
#means <- apply(expr,1,mean)
tpm <- tpm[which(sums>(0.5*max(sums))),]


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

# Store full results in a data frame
transcriptional_subtype <- data.frame(aliquots, subtype_class, simplicity_score,stringsAsFactors=FALSE)
colnames(transcriptional_subtype) <- c("aliquot_barcode","transcriptional_subtype","simplicity_score")

# Add results to clinical info table
names(subtype_class) <- gsub("_cd45n","",aliquots)
info[,"transcriptional_subtype"] <- subtype_class[rownames(info)]

##################################################
# Step 3: Examine whether macrophages or microglia activation modules are enriched in different subtypes (ssGSEA)
##################################################

# Load macrophage modules and convert to signature lists:
module_inf <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"

mac_module_df <- read.delim(module_inf,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- gsub("X","module",colnames(mac_module_df))

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)

myeloid_tpm <- cbind(tpm[,grep("glioma_\\d\\d\\d\\d_mdm",colnames(tpm))], tpm[,grep("glioma_\\d\\d\\d\\d_mg",colnames(tpm))])

# Run GSVA with macrophage modules on myeloid cell tpm:
myeloid_res <- t(gsva(myeloid_tpm, mac_modules, method="ssgsea",parallel.sz=1))

# Focusing on IDHwt only:
idhwt <- info[which(info[,"Tumor.Type"]=="Glioma" & info[,"IDH.Status"]=="wild type"),]

subtypes <- c("Proneural","Classical","Mesenchymal")
modules <- colnames(myeloid_res)
mac_pval <- mg_pval <- mac_eff <- mg_eff <- matrix(0,nrow= length(modules), ncol = length(subtypes))
rownames(mac_es) <- rownames(mg_es) <- rownames(mac_eff) <- rownames(mg_eff) <- modules
colnames(mac_es) <- colnames(mg_es) <- colnames(mac_eff) <- colnames(mg_eff) subtypes
for(i in 1:length(subtypes))
{
	cat("\r",i)
	
	# Macrophages
	#--------------------
	mac_res <- myeloid_res[grep("glioma_\\d\\d\\d\\d_mdm",rownames(myeloid_res)),]
	rownames(mac_res) <- gsub("_mdm","",rownames(mac_res))
	
	# Match the info and expression names
	comxx <- intersect(rownames(mac_res),rownames(idhwt))
	mac_res <- mac_res[comxx,]
	mac_info <- idhwt[comxx,]
	
	subtype_res <- c()
	for(j in 1:length(modules))
	{
		g1 <- mac_res[rownames(mac_info[which(mac_info[,"transcriptional_subtype"]==subtypes[i]),]),j]
		g2 <- mac_res[rownames(mac_info[which(mac_info[,"transcriptional_subtype"]!=subtypes[i]),]),j]
		
		mac_es[j,i] <- wilcox.test(g1,g2)$p.value
		mac_eff[j,i] <- median(g1) - median(g2)
	}
	
	# Microglia
	#--------------------
	mg_res <- myeloid_res[grep("glioma_\\d\\d\\d\\d_mg",rownames(myeloid_res)),]
	rownames(mg_res) <- gsub("_mg","",rownames(mg_res))
	
	# Match the info and expression names
	comxx <- intersect(rownames(mg_res),rownames(idhwt))
	mg_res <- mg_res[comxx,]
	mg_info <- idhwt[comxx,]
	
	subtype_res <- c()
	for(j in 1:length(modules))
	{
		g1 <- mg_res[rownames(mg_info[which(mg_info[,"transcriptional_subtype"]==subtypes[i]),]),j]
		g2 <- mg_res[rownames(mg_info[which(mg_info[,"transcriptional_subtype"]!=subtypes[i]),]),j]
		
		mg_es[j,i] <- wilcox.test(g1,g2)$p.value
		mg_eff[j,i] <- median(g1) - median(g2)
	}
}


##################################################
# Step 4: Identify the differentially expressed genes in macrophages from each subtype (using DESeq2)
##################################################

subtypes <- c("Proneural","Classical","Mesenchymal")
mac_up <- mg_up <- list()
for(i in 1:length(subtypes))
{
	cat("\r",i)
	
	# Macrophages
	#--------------------
	mac_cts <- expr[,grep("glioma_\\d\\d\\d\\d_mdm",colnames(expr))]
	colnames(mac_cts) <- gsub("_mdm","",colnames(mac_cts))
	
	# Match the info and expression names
	comxx <- intersect(colnames(mac_cts),rownames(idhwt))
	mac_cts <- mac_cts[,comxx]
	mac_info <- idhwt[comxx,]
	
	# Add subtype of interest indicator
	mac_info[,"mysubtype"] <- mac_info[,"transcriptional_subtype"] == subtypes[i]

	# Use DESeq2 to compare subtype versus no subtype
	mac_dds <- DESeqDataSetFromMatrix(countData = mac_cts,
								  colData = mac_info,
								  design= ~ batch + mysubtype)
	dds <- DESeq(mac_dds)
	res <- results(dds, name="mysubtypeTRUE")
	mac_up[[i]] <- rownames(res[which(res[,"padj"]< 0.05 & res[,"log2FoldChange"]> 3),])
	
	
	# Microglia
	#--------------------
	mg_cts <- expr[,grep("glioma_\\d\\d\\d\\d_mg",colnames(expr))]
	colnames(mg_cts) <- gsub("_mg","",colnames(mg_cts))
	
	# Match the info and expression names
	comxx <- intersect(colnames(mg_cts),rownames(idhwt))
	mg_cts <- mg_cts[,comxx]
	mg_info <- idhwt[comxx,]
	
	# Add subtype of interest indicator
	mg_info[,"mysubtype"] <- mg_info[,"transcriptional_subtype"] == subtypes[i]

	# Use DESeq2 to compare subtype versus no subtype
	mg_dds <- DESeqDataSetFromMatrix(countData = mg_cts,
								  colData = mg_info,
								  design= ~ batch + mysubtype)
	dds <- DESeq(mg_dds)
	res <- results(dds, name="mysubtypeTRUE")
	mg_up[[i]] <- rownames(res[which(res[,"padj"]< 0.05 & res[,"log2FoldChange"]> 3),])
}
names(mac_up) <- names(mg_up) <- subtypes

##################################################
# Step 5: Overrepresentation analysis of different macrophage activation states
##################################################

module_inf <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"

mac_module_df <- read.delim(module_inf,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- gsub("X","module",colnames(mac_module_df))

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)


measured_genes <- rownames(expr)
mac_pval <- mg_pval <- matrix(0, nrow = length(mac_modules),ncol = length(mac_up))
rownames(mac_pval) <- rownames(mg_pval) <- names(mac_modules)
colnames(mac_pval) <- colnames(mg_pval) <- names(mac_up)
for(i in 1:length(mac_modules))
{
	mymodule <- mac_modules[[i]]
	mymodule <- mymodule[which(mymodule %in% measured_genes)]
	
	mac_pval[i,] <- sapply(mac_up, function(x) or_test(x, mymodule, measured_genes))
	mg_pval[i,] <- sapply(mg_up, function(x) or_test(x, mymodule, measured_genes))
}

mac_qval <- apply(mac_pval, 2, function(x) p.adjust(x,"BH"))
mg_qval <- apply(mg_pval, 2, function(x) p.adjust(x,"BH"))

which(mac_qval[,3] < 0.1)
# module41 
#       41 
which(mg_qval[,3] < 0.1)
#  module3  module6 module15 
#        3        6       15 