###################################################
# Calculate macrophage module enrichment in purified myeloid cell profiles
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.05.29
# Author: Frederick Varn
##################################################

library(tidyverse)
library(DBI)
library(ssgsea.GBM.classification)
library(DESeq2)
library(GSVA)
library(ggdendro)

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

# Single batch analysis (nyc because its bigger)
info <- info %>%
	filter(batch == "nyc")

##################################################
# Helper functions
##################################################

# TPM calculation function
# Source: https://www.biostars.org/p/307603/
# Modified from: https://gist.github.com/slowkow/c6ab0348747f86e2748b
#-------------------------------------------

Counts_to_tpm <- function(counts, featureLength) {

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))

  # Compute effective lengths of features in each library.
  effLen <- featureLength

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
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
# Step 1: Calculate TPM from count data 
##################################################
	
# Get gene lengths
geneLengths <- read.delim("/projects/verhaak-lab/verhaak_ref/star/GRCh38/GRCh38.d1.vd1/gencode.gene.info.v22.tsv",stringsAsFactor=FALSE)

featureLengths <- geneLengths %>%
	filter(gene_name %in% rownames(expr)) %>%
	filter(gene_type == "protein_coding") %>%
	group_by(gene_name) %>%
	slice_min(order_by = exon_num) %>%
	slice_min(order_by = exon_length) %>%
	slice_min(order_by = start) %>%
	dplyr::select(gene_name, exon_length)
	
featureVector <- featureLengths$exon_length
names(featureVector) <- featureLengths$gene_name
featureVector <- featureVector[rownames(expr)]

tpm <- Counts_to_tpm(expr, featureVector)

# Test with Klemm preprocessed data on the webserver (Test 1: ITGA4, Test 2: HLA-A)
# tpm_comp <- tpm["ITGA4",test[,"Sample.ID"]]
# test <- read.csv("/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/normCounts_ITGA4.csv",stringsAsFactor=FALSE)
# cor(tpm_comp, test[,"normCounts"])
# 0.9399454
# cor(tpm_comp, test[,"normCounts"],method="s")
# 
# test <- read.csv("/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/normCounts_HLA-A.csv",stringsAsFactor=FALSE)
# tpm_comp <- tpm["HLA-A",test[,"Sample.ID"]]
# cor(tpm_comp, test[,"normCounts"])
# 0.6070107
# cor(tpm_comp, test[,"normCounts"],method="s")
# 0.8384406


##################################################
# Step 2: Get subtype of each CD45- cells from glioma
##################################################

glioma <- tpm[,grep("glioma",colnames(tpm))]
glioma <- glioma[,grep("cd45n",colnames(glioma))]
	
# Save glioma-specific TPM file in GCT format for classifier

cd45n_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_tpm_cd45n.gct"
gct(glioma, cd45n_path)
	
# Run Qianghu's transcriptional classifier (only needed once)
# runSsGSEAwithPermutation(cd45n_path,100)

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

# If ties, give mesenchymal the benefit of the doubt due to low macrophages during dissociation
transcriptional_subtype[grep(",Mesenchymal",transcriptional_subtype[,"transcriptional_subtype"]),"transcriptional_subtype"] <- "Mesenchymal"
subtype_class[grep(",Mesenchymal",subtype_class)] <- "Mesenchymal"

# Add results to clinical info table
names(subtype_class) <- gsub("_cd45n","",aliquots)
info[,"transcriptional_subtype"] <- subtype_class[rownames(info)]
info[,"sampleID"] <- rownames(info)

##################################################
# Step 3: Examine whether macrophages or microglia activation modules are enriched in different subtypes (ssGSEA)
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

myeloid_tpm <- cbind(tpm[,grep("glioma_\\d\\d\\d\\d_mdm",colnames(tpm))], tpm[,grep("glioma_\\d\\d\\d\\d_mg",colnames(tpm))])

# Run GSVA with macrophage modules on myeloid cell tpm:
res <- gsva(myeloid_tpm, mac_modules, method="ssgsea",parallel.sz=1)

module <- rownames(res)
res <- data.frame(module, res,stringsAsFactors=FALSE)
colnames(res) <- gsub("\\.","-",colnames(res))


plot_res <- res %>%
 			pivot_longer(-module, names_to = "aliquot_barcode", values_to = "es") %>%
 			mutate(sampleID = sapply(strsplit(aliquot_barcode,"_"),function(x)paste(x[1:3],collapse="_"))) %>%
 			inner_join(info, "sampleID") %>%
 			as.data.frame()

write.table(plot_res, "/projects/verhaak-lab/GLASS-III/data/res/macrophage_modules/klemm_facs_myeloid_xue_modules.txt",sep="\t",quote=FALSE,row.names=FALSE)

# Rescale each measurement (independently) to have a mean of 0 and variance of 1
plot_scaled <- plot_res
for(i in 1:length(module))
{
	plot_scaled[which(plot_res[,"module"]==module[i]),"es"] <- scale(plot_scaled[which(plot_res[,"module"]==module[i]),"es"])
}

plot_scaled <- plot_scaled %>%
			   mutate(idhmut = recode(IDH.Status, 'wild type' = 'IDHwt', 'mutated' = 'IDHmut')) %>%
			   mutate(codel = recode(X1p.19q.co.deletion, "Yes" = "codel")) %>%
			   mutate(codel = replace_na(codel, "noncodel")) %>%
			   unite("idh_codel_status", idhmut:codel, sep="-") %>%
			   mutate(idh_codel_status = recode(idh_codel_status, "IDHwt-codel" = "IDHwt", "IDHwt-noncodel" = "IDHwt")) %>%
			   mutate(module = gsub("module_","",module))
	
# Macrophage module order for purified macrophages (see macrophage module script)
module_levels <- as.character(c(23, 25, 1, 26, 34, 19, 36, 11, 24, 2, 27, 49, 10, 35, 17, 37, 13, 48, 31, 40, 6, 5, 22, 38, 16, 21, 39, 3, 4, 7, 8, 9, 18, 28, 43, 20, 14, 46, 41, 44, 42, 47, 45, 12, 15, 30, 32, 29, 33))
plot_scaled[,"module"] <- factor(plot_scaled[,"module"], levels = module_levels, ordered=TRUE)
	
			   
original_scaled <- plot_scaled

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Make heatmap
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# IDHwt

plot_scaled <- original_scaled %>%
	filter(idh_codel_status == "IDHwt")

# plot_res <- plot_res %>% filter(IDH.Status == "wild type")
# p.val <- rep(0, length(module))
# names(p.val) <- module
# for(i in 1:length(module))
# {
# 	sub_res <- plot_res[which(plot_res[,"module"] == module[i]),]
# 	g1 <- sub_res[which(sub_res[,"Status"] == "initial"),"es"]
# 	g2 <- sub_res[which(sub_res[,"Status"] == "recurrent"),"es"]
# 
# 	p.val[i] <- wilcox.test(g1,g2 )$p.value
# }
# p.val[which(p.val < 0.05)]

#Make a matrix for clustering
plot_mat <- matrix(0,nrow=length(unique(plot_scaled[,"module"])),ncol=length(unique(plot_scaled[,"sampleID"])))
rownames(plot_mat) <- unique(plot_scaled[,"module"])
colnames(plot_mat) <- unique(plot_scaled[,"sampleID"])
for(i in 1:nrow(plot_scaled))
{
	mysigname <- as.character(plot_scaled[i,"module"])
	myalicode <- as.character(plot_scaled[i,"sampleID"])
	plot_mat[mysigname,myalicode] <- plot_scaled[i,"es"]
}
plot_mat <- t(plot_mat)

# Run clustering
dendro <- as.dendrogram(hclust(d = dist(x = plot_mat)))
#Reverse to get high on the right and low on the left
dendro <- rev(dendro)

# Create dendrogram plot
dendro.plot <- ggdendrogram(data = dendro) + 
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
dendro.plot

# Heatmap
# Data wrangling
# Extract the order of the samples in the dendrogram
plot.order <- order.dendrogram(dendro)
# Order the levels according to their position in the cluster
plot_scaled[,"aliquot_barcode"] <- factor(plot_scaled[,"aliquot_barcode"],
                               levels = rownames(plot_mat)[plot.order], 
                               ordered = TRUE)
# Height for sidebar                          
plot_scaled[,"status"] <- 1

#Create sidebars with sample information

# Initial/recurrence status
gg_timepoint <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = Status)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("white", "black")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
  
# IDHmut status
gg_idhmut <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = idh_codel_status)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("#619CFF")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Transcriptional subtype
gg_ts <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = transcriptional_subtype)) +
	geom_bar(stat="identity") +
	theme_void() +
  	scale_fill_manual(values=c("#008A22", "#8A0000", "#00458A")) +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Batch
gg_batch <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = batch)) +
	geom_bar(stat="identity") +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
        
# Create heatmap plot

heatmap.plot <- ggplot(data = plot_scaled, aes(x = sampleID, y = module)) +
  geom_tile(aes(fill=es)) +
  scale_fill_gradient2(low ="royalblue4", mid="white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  theme_void() +
  theme(axis.text.y = element_text(size=7),
  		axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(gg_batch)
gb3 <- ggplot_build(gg_timepoint)
gb4 <- ggplot_build(gg_idhmut)
gb5 <- ggplot_build(gg_ts)
gb6 <- ggplot_build(heatmap.plot)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)
n4 <- length(gb4$layout$panel_scales_y[[1]]$range$range)
n5 <- length(gb5$layout$panel_scales_y[[1]]$range$range)
n6 <- length(gb6$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)
gF <- ggplot_gtable(gb6)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")
g <- gtable:::rbind_gtable(g, gF, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n6/4, "null")
g$heights[panels[2]] <- unit(n6/16, "null")
g$heights[panels[3]] <- unit(n6/16,"null")
g$heights[panels[4]] <- unit(n6/16,"null")
g$heights[panels[5]] <- unit(n6/16,"null")
g$heights[panels[6]] <- unit(n6,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/macrophage_heatmap_klemm_idhwt.pdf",width=2,height=7, useDingbats = FALSE)
grid.draw(g)
dev.off()

#write.table(plot_res,"/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/macrophage_subtype_res.txt",sep="\t",quote=FALSE,row.names=FALSE)

# IDHmut

plot_scaled <- original_scaled %>%
	filter(idh_codel_status != "IDHwt")

#Make a matrix for clustering
plot_mat <- matrix(0,nrow=length(unique(plot_scaled[,"module"])),ncol=length(unique(plot_scaled[,"sampleID"])))
rownames(plot_mat) <- unique(plot_scaled[,"module"])
colnames(plot_mat) <- unique(plot_scaled[,"sampleID"])
for(i in 1:nrow(plot_scaled))
{
	mysigname <- as.character(plot_scaled[i,"module"])
	myalicode <- as.character(plot_scaled[i,"sampleID"])
	plot_mat[mysigname,myalicode] <- plot_scaled[i,"es"]
}
plot_mat <- t(plot_mat)

# Run clustering
dendro <- as.dendrogram(hclust(d = dist(x = plot_mat)))
#Reverse to get high on the right and low on the left
dendro <- rev(dendro)

# Create dendrogram plot
dendro.plot <- ggdendrogram(data = dendro) + 
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
dendro.plot

# Heatmap
# Data wrangling
# Extract the order of the samples in the dendrogram
plot.order <- order.dendrogram(dendro)
# Order the levels according to their position in the cluster
plot_scaled[,"aliquot_barcode"] <- factor(plot_scaled[,"aliquot_barcode"],
                               levels = rownames(plot_mat)[plot.order], 
                               ordered = TRUE)
# Height for sidebar                          
plot_scaled[,"status"] <- 1

#Create sidebars with sample information

# Initial/recurrence status
gg_timepoint <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = Status)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("white", "black")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
  
# IDHmut status
gg_idhmut <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = idh_codel_status)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("#F8766D", "#00BA38")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Transcriptional subtype
gg_ts <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = transcriptional_subtype)) +
	geom_bar(stat="identity") +
	theme_void() +
  	scale_fill_manual(values=c("#008A22", "#8A0000", "#00458A")) +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Batch
gg_batch <- ggplot(data = plot_scaled, aes(x = sampleID, y = status, fill = batch)) +
	geom_bar(stat="identity") +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
        
# Create heatmap plot

heatmap.plot <- ggplot(data = plot_scaled, aes(x = sampleID, y = module)) +
  geom_tile(aes(fill=es)) +
  scale_fill_gradient2(low ="royalblue4", mid="white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  theme_void() +
  theme(axis.text.y = element_text(size=7),
  		axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(gg_batch)
gb3 <- ggplot_build(gg_timepoint)
gb4 <- ggplot_build(gg_idhmut)
gb5 <- ggplot_build(gg_ts)
gb6 <- ggplot_build(heatmap.plot)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)
n4 <- length(gb4$layout$panel_scales_y[[1]]$range$range)
n5 <- length(gb5$layout$panel_scales_y[[1]]$range$range)
n6 <- length(gb6$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)
gF <- ggplot_gtable(gb6)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")
g <- gtable:::rbind_gtable(g, gF, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n6/4, "null")
g$heights[panels[2]] <- unit(n6/16, "null")
g$heights[panels[3]] <- unit(n6/16,"null")
g$heights[panels[4]] <- unit(n6/16,"null")
g$heights[panels[5]] <- unit(n6/16,"null")
g$heights[panels[6]] <- unit(n6,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/macrophage_heatmap_klemm_idhmut.pdf",width=1,height=7, useDingbats = FALSE)
grid.draw(g)
dev.off()

# Create average profiles? Correlate clusters?? How do macrophages vs microglia affect this may need to subset.