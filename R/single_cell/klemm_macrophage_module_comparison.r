###################################################
# Compare macrophage activation profiles in Klemm between transcriptional subtypes
# RNAseq count data was variance stabilized using DESeq2 and then batch corrected to get final expression values
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.07.22
# Author: Frederick Varn
##################################################

library(tidyverse)
library(DESeq2)
library(ssgsea.GBM.classification)
library(GSVA)
library(limma)
library(RColorBrewer)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

cts <- read.csv(myinf1,row.names=1)
info <- read.csv(myinf2,row.names=1,stringsAsFactor=FALSE)

# Add a batch column to info
info[,"batch"] <- sapply(strsplit(rownames(info),"_"),function(x)x[1])

sort_name <- colnames(cts)
sample_name <- sapply(strsplit(sort_name, "_"),function(x)paste(x[1:3],collapse="_"))
sub_info <- info[sample_name,]
rownames(sub_info) <- NULL
coldata <- data.frame(sort_name, sample_name, sub_info)

# Create DESeq2 object and perform variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ batch)                     
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)

# Batch correction
expr <- limma::removeBatchEffect(mat, vsd$batch)

# Test to ensure batch effect was removed
prevars <- apply(t(cts),2,var)
postvars <- apply(t(expr),2,var)

pcacts <- t(cts)[,-which(prevars==0)]
pcaexpr <- t(expr)[,-which(postvars==0)]

prepca <- prcomp(pcacts, scale.=TRUE)
postpca <- prcomp(pcaexpr, scale.=TRUE)


prescores <- as.data.frame(prepca$x)
prescores[,"batch"] <- coldata[,"batch"]
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/klemm_pre_correction_pca.pdf",width=7,height=7)
ggplot(data = prescores, aes(x = PC1, y = PC2, colour=batch)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  theme_classic()
dev.off()

postscores <- as.data.frame(postpca$x)
postscores[,"batch"] <- coldata[,"batch"]
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/qc/klemm_post_correction_pca.pdf",width=7,height=7)
ggplot(data = postscores, aes(x = PC1, y = PC2, colour=batch)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  theme_classic()
dev.off()

##################################################
# Helper functions
##################################################

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
# Step 1: Get subtype of each CD45- cells from glioma
##################################################

glioma <- expr[,grep("glioma",colnames(expr))]
glioma <- glioma[,grep("cd45n",colnames(glioma))]
	
# Save glioma-specific expr file in GCT format for classifier

cd45n_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_vst_cd45n_deseq.gct"
gct(glioma, cd45n_path)
	
# Run Qianghu's transcriptional classifier (only needed once)
#runSsGSEAwithPermutation(cd45n_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/p_result_BrainTIME_vst_cd45n_deseq.gct.txt")

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

# If ties, give mesenchymal the benefit of the doubt due to low macrophages during dissociation
transcriptional_subtype[grep(",Mesenchymal",transcriptional_subtype[,"transcriptional_subtype"]),"transcriptional_subtype"] <- "Mesenchymal"
subtype_class[grep(",Mesenchymal",subtype_class)] <- "Mesenchymal"

# Add results to clinical info table
names(subtype_class) <- gsub("_cd45n","",aliquots)
info[,"transcriptional_subtype"] <- subtype_class[rownames(info)]
info[,"sampleID"] <- rownames(info)

##################################################
# Step 2: Examine whether macrophages or microglia activation modules are enriched in different subtypes (ssGSEA)
##################################################

# Load macrophage modules and convert to signature lists:
module_inf <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"

mac_module_df <- read.delim(module_inf,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- gsub("X","module_",colnames(mac_module_df))

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)

myeloid_vst <- cbind(expr[,grep("glioma_\\d\\d\\d\\d_mdm",colnames(expr))], expr[,grep("glioma_\\d\\d\\d\\d_mg",colnames(expr))])

# Run GSVA with macrophage modules on myeloid cell vst:
res <- t(gsva(myeloid_vst, mac_modules, method="ssgsea",parallel.sz=1))

subtype <- unique(info[,"transcriptional_subtype"])
subtype <- subtype[-which(is.na(subtype))]

eff <- p.val <- q.val <- matrix(0, nrow = ncol(res), ncol = length(subtype))
rownames(eff) <- rownames(p.val) <- rownames(q.val) <- colnames(res)
colnames(eff) <- colnames(p.val) <- colnames(q.val) <- subtype
for(i in 1:length(subtype))
{
	mes_info <- info %>%
			filter(Tissue.Type == "Glioma", IDH.Status == "wild type", transcriptional_subtype == subtype[i])
	cn1 <- c(paste(mes_info[,"sampleID"],"_mdm",sep=""), paste(mes_info[,"sampleID"],"_mg",sep=""))
	cn1 <- intersect(cn1, rownames(res))
	#cn1 <- cn1[grep("_mg",cn1)]
	g1 <- res[cn1,]		

	nonmes_info <- info %>%
				filter(Tissue.Type == "Glioma", IDH.Status == "wild type", transcriptional_subtype != subtype[i])
	cn2 <- c(paste(nonmes_info[,"sampleID"],"_mdm",sep=""), paste(nonmes_info[,"sampleID"],"_mg",sep=""))
	cn2 <- intersect(cn2, rownames(res))
	#cn2 <- cn2[grep("_mg",cn2)]
	g2 <- res[cn2,]		

	test_data <- rbind(g1,g2)
	
	eff[,i] <- apply(test_data, 2, function(x)median(x[cn1]) - median(x[cn2]))
	mypval <- apply(test_data, 2, function(x)wilcox.test(x[cn1],x[cn2])$p.value)
	p.val[,i] <- mypval
	q.val[,i] <- p.adjust(mypval, "BH")	
}

# Write results to table
module <- rep(rownames(p.val), 3)
subtype <- rep(subtype, each = nrow(p.val))
p.value <- c(p.val[,1], p.val[,2], p.val[,3])
q.value <- c(q.val[,1], q.val[,2], q.val[,3])
median_effect <- c(eff[,1], eff[,2], eff[,3])
full_result <- data.frame(module, subtype, p.value, q.value, median_effect)

write.table(full_result, "/projects/verhaak-lab/GLASS-III/results/cibersortx/downstream/myeloid/klemm_transcriptional_subtype_comparisons.txt",sep="\t",quote=FALSE,row.names=FALSE)



# Look at all samples

mes_info <- info %>%
			filter(Tissue.Type == "Glioma", transcriptional_subtype == "Mesenchymal")
cn1 <- c(paste(mes_info[,"sampleID"],"_mdm",sep=""), paste(mes_info[,"sampleID"],"_mg",sep=""))
cn1 <- intersect(cn1, rownames(res))
cn1 <- cn1[grep("_mg",cn1)]
g1 <- res[cn1,]		

nonmes_info <- info %>%
			filter(Tissue.Type == "Glioma", transcriptional_subtype != "Mesenchymal")
cn2 <- c(paste(nonmes_info[,"sampleID"],"_mdm",sep=""), paste(nonmes_info[,"sampleID"],"_mg",sep=""))
cn2 <- intersect(cn2, rownames(res))
cn2 <- cn2[grep("_mg",cn2)]
g2 <- res[cn2,]		

test_data <- rbind(g1,g2)

p.val <- apply(test_data, 2, function(x)wilcox.test(x[cn1],x[cn2])$p.value)
eff <- apply(test_data, 2, function(x)median(x[cn1]) - median(x[cn2]))

q.val <- p.adjust(p.val,"BH")
q.val[which(q.val<0.1)]

