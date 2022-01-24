##################################################
# Compare the TCGA GEPs (transcriptional profiles) between IDHwt tumors across transcriptional subtypes
# Author: Kevin Johnson
# Date: 2022.01.08
# Revision comment
#######################################################


# Necessary packages:
library(tidyverse)
library(ssgsea.GBM.classification)
library(GSVA)
library(RColorBrewer)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)
library(topGO)

## Prompted by the following reviewer comment:
# It would be interesting to show the subtype specific malignant cell differences in a non-paired manner 
# (i.e. compare mesenchymal vs. pro-neural tumors) and associated subtype specific pathways that may underly 
# the related histological or immunological features.

## Standard stripped down plotting theme.
plot_theme    <- theme_bw(base_size = 12) + theme(axis.title = element_text(size = 12),
                                                  axis.text = element_text(size = 12),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank(),
                                                  strip.background = element_blank(),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(), 
                                                  panel.border = element_blank(),
                                                  axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                                                  axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))

## Provide output folder for temporary files and figures.
main_dir <- "/fastscratch"
sub_dir <- "varnf"
output_dir <- file.path(main_dir, sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

# Apply the RNAseq classifier to the TCGA data.
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/p_result_TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct.txt")

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

# Store full results in a data frame.
transcriptional_subtype <- data.frame(aliquots, subtype_class, simplicity_score,stringsAsFactors=FALSE)
colnames(transcriptional_subtype) <- c("aliquot_barcode","transcriptional_subtype","simplicity_score")
transcriptional_subtype[,"case_barcode"] <- sapply(strsplit(transcriptional_subtype[,"aliquot_barcode"],"\\."),function(x)paste(x[1:3],collapse="-"))

## TCGA clinical data and combine with the transcriptional subtype.
## Restrict to "pure" transcriptional subtypes. That is, removed mixed subtypes.
clin_dat <- read.delim("/projects/verhaak-lab/USERS/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt")
clin_dat_filt <- clin_dat %>% 
  filter(RNAseq=="Yes") %>% 
  dplyr::select(case_barcode = Case, idh_status = IDH.status) %>% 
  inner_join(transcriptional_subtype, by="case_barcode") %>% 
  filter(transcriptional_subtype%in%c("Classical", "Mesenchymal", "Proneural"))

## Finalize the two datasets for IDHwt and IDHmut.
clin_dat_filt_wt <- clin_dat_filt %>%
  filter(idh_status=="WT") %>% 
  mutate(gep_barcode = substr(aliquot_barcode, 1, 15))
clin_dat_filt_mut <- clin_dat_filt %>% filter(idh_status=="Mutant")

## Load in the GEPs and perform pairwise differential gene expression testing (stem-like & diff-like: pro vs. others, clas vs. others, mes vs. others)
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in, but will also incorporate prolif. stem if interested.
mytag <- gsub("CIBERSORTxHiRes_toil_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/CIBERSORTxHiRes_toil_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)

## Set the specific data subset.
## Load in example geps:
geps <- read.delim(myinf1[1], row.names=1)
which(!clin_dat_filt_wt$gep_barcode[clin_dat_filt_wt$transcriptional_subtype=="Mesenchymal"]%in%colnames(geps))
clin_dat_filt_wt$gep_barcode[clin_dat_filt_wt$transcriptional_subtype=="Mesenchymal"][22]

## Remove the object:
rm(geps)

## Remove the sample that does not have a GEP.
dat = clin_dat_filt_wt %>%
  filter(gep_barcode!="TCGA.02.2486.01")


###### Proneural vs. Others 
## Set up the transcriptional comparisons, we are interested in Proneural vs. not Proneural.
p <- se <- sigs <- list()
for(i in 1:length(myinf1))
{
  cat("\r", i)
  geps <- read.delim(myinf1[i], row.names=1)
  geps <- log10(geps+1)
  rem <- apply(geps,1,function(x)sum(is.na(x)))
  geps <- geps[-which(rem==ncol(geps)),]	
  vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
  geps <- geps[which(vars > 0),]
  nrow(geps)
  
  g1 <- geps[,dat$gep_barcode[dat$transcriptional_subtype=="Proneural"]]
  g2 <- geps[,dat$gep_barcode[dat$transcriptional_subtype!="Proneural"]]
  
  p.val <- eff <- rep(0, nrow(geps))
  names(p.val) <- names(eff) <- rownames(geps)
  for(j in 1:nrow(geps))
  {
    group1 <- as.numeric(g1[j,])
    group2 <- as.numeric(g2[j,])
    
    p.val[j] <- wilcox.test(group1,group2,paired=FALSE)$p.value
    eff[j] <- log2(mean(group2)/mean(group1))
    #eff[i] <- median(group2) - median(group1)
  }
  q.val <- p.adjust(p.val,"BH")
  idhwt_res <- data.frame(p.val, q.val, eff)
  idhwt_res <- idhwt_res[order(eff),]
  
  idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
  idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.05
  
  idhwt_res <- idhwt_res[order(idhwt_res$p.val),]
  
  myoutf <- paste(output_dir, "/tcga_idhwt_", cell_state[i],"_proneural_vs_other_result.txt",sep="")
  
  write.table(idhwt_res, myoutf, sep="\t",quote=FALSE,row.names=TRUE) 
  
  sigs[[i]] <- idhwt_res
  
  # Plot heatmaps (if desired), irrelevant to final image.
  mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
  sig_matrix <- geps[mygenes,c(dat$gep_barcode[dat$transcriptional_subtype=="Proneural"], dat$gep_barcode[dat$transcriptional_subtype!="Proneural"])]
  
  sig_matrix <- t(apply(sig_matrix, 1, function(x) x - median(x)))
  
  grp1 <- sig_matrix[,dat$gep_barcode[dat$transcriptional_subtype=="Proneural"]]
  grp1 <- apply(grp1, 1, mean)
  grp2 <- sig_matrix[,dat$gep_barcode[dat$transcriptional_subtype!="Proneural"]]
  grp2 <- apply(grp2, 1, mean)
  
  expr <- c(grp1, grp2)
  # Rescale for color
  expr[which(expr > 0.33)] <- 0.33 # Plug the number in from the whole dataset (obtained below) back in for scaling
  gene_symbol <- c(names(grp1), names(grp2))
  subtype <- c(rep("Proneural", length(grp1)), rep("Other", length(grp2)))
  
  plot_hm <- data.frame(gene_symbol, expr, subtype)
  lev <- names(grp1)[order(grp2 - grp1)]
  plot_hm$gene_symbol <- factor(plot_hm$gene_symbol, levels = lev)
  plot_hm[,"cell"] <- cell_state[i]
  
  p[[i]] <- plot_hm
  
  se[[i]] <- ggplot(data = plot_hm, aes(x = subtype, y = gene_symbol)) +
    geom_tile(aes(fill=expr)) +
    scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
                         space = "Lab", na.value = "grey50", guide = "colourbar",
                         aesthetics = "fill", limits = c(-0.33, 0.33)) + 
    theme_void() +
    theme(axis.text = element_blank(),
          axis.title= element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.title = element_blank(),
          legend.position = "none")	
}


###### GO enrichment  ########
tumor_sigs <- sigs
tumor_tag <- mytag


#### Up-regulated in Proneural vs. Others ###
# Diff gene function for GO enrichment
diffGenes <- function(bg_genes) {
  return(bg_genes == 1)}

goList <- list()
for(i in 1:length(tumor_sigs))
{
  cat("\r", i)
  mysig <- tumor_sigs[[i]]
  up <- rownames(mysig %>% filter(sig, eff < 0))
  
  all_genes <- rownames(mysig)
  bg_genes <- as.numeric(all_genes %in% up)
  names(bg_genes) <- all_genes
  
  # Functional enrichment of signature
  sampleGOdata <- new("topGOdata",
                      description = "Simple session", 
                      ontology = "BP",
                      allGenes = bg_genes,
                      geneSelectionFun = diffGenes, 		# Call function here
                      nodeSize = 10,
                      annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
  
  # Fishers test
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
  fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
  fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
  
  fishRes[,"cell_state"] <- tumor_tag[i]
  goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)

## Output GO results for Proneural vs. others
write.table(goRes, paste0(output_dir, "/tcga_idhwt_proneural_vs_others_sig_GO_result.txt"), sep="\t",quote=FALSE,row.names=TRUE) 

############ Classical vs. Others   #############
## Repeat analysis for Classical vs. non-classical tumors.
p <- se <- sigs <- list()
for(i in 1:length(myinf1))
{
  cat("\r", i)
  geps <- read.delim(myinf1[i], row.names=1)
  geps <- log10(geps+1)
  rem <- apply(geps,1,function(x)sum(is.na(x)))
  geps <- geps[-which(rem==ncol(geps)),]	
  vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
  geps <- geps[which(vars > 0),]
  nrow(geps)
  
  g1 <- geps[,dat$gep_barcode[dat$transcriptional_subtype=="Classical"]]
  g2 <- geps[,dat$gep_barcode[dat$transcriptional_subtype!="Classical"]]
  
  p.val <- eff <- rep(0, nrow(geps))
  names(p.val) <- names(eff) <- rownames(geps)
  for(j in 1:nrow(geps))
  {
    group1 <- as.numeric(g1[j,])
    group2 <- as.numeric(g2[j,])
    
    p.val[j] <- wilcox.test(group1,group2,paired=FALSE)$p.value
    eff[j] <- log2(mean(group2)/mean(group1))
    #eff[i] <- median(group2) - median(group1)
  }
  q.val <- p.adjust(p.val,"BH")
  idhwt_res <- data.frame(p.val, q.val, eff)
  idhwt_res <- idhwt_res[order(eff),]
  
  idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
  idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.05
  
  idhwt_res <- idhwt_res[order(idhwt_res$p.val),]
  
  myoutf <- paste(output_dir, "/tcga_idhwt_", cell_state[i],"_classical_vs_other_result.txt",sep="")
  
  write.table(idhwt_res, myoutf, sep="\t",quote=FALSE,row.names=TRUE) 
  
  sigs[[i]] <- idhwt_res
  
  # Plot heatmaps (unnecessary).
  mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
  sig_matrix <- geps[mygenes,c(dat$gep_barcode[dat$transcriptional_subtype=="Classical"], dat$gep_barcode[dat$transcriptional_subtype!="Classical"])]
  
  sig_matrix <- t(apply(sig_matrix, 1, function(x) x - median(x)))
  
  grp1 <- sig_matrix[,dat$gep_barcode[dat$transcriptional_subtype=="Classical"]]
  grp1 <- apply(grp1, 1, mean)
  grp2 <- sig_matrix[,dat$gep_barcode[dat$transcriptional_subtype!="Classical"]]
  grp2 <- apply(grp2, 1, mean)
  
  expr <- c(grp1, grp2)
  # Rescale for color
  expr[which(expr > 0.33)] <- 0.33 # Plug the number in from the whole dataset (obtained below) back in for scaling
  gene_symbol <- c(names(grp1), names(grp2))
  subtype <- c(rep("Classical", length(grp1)), rep("Other", length(grp2)))
  
  plot_hm <- data.frame(gene_symbol, expr, subtype)
  lev <- names(grp1)[order(grp2 - grp1)]
  plot_hm$gene_symbol <- factor(plot_hm$gene_symbol, levels = lev)
  plot_hm[,"cell"] <- cell_state[i]
  
  p[[i]] <- plot_hm
  
  se[[i]] <- ggplot(data = plot_hm, aes(x = subtype, y = gene_symbol)) +
    geom_tile(aes(fill=expr)) +
    scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
                         space = "Lab", na.value = "grey50", guide = "colourbar",
                         aesthetics = "fill", limits = c(-0.33, 0.33)) + 
    theme_void() +
    theme(axis.text = element_blank(),
          axis.title= element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.title = element_blank(),
          legend.position = "none")	
}


###### GO enrichment  ########
tumor_sigs <- sigs
tumor_tag <- mytag

#### Up-regulated in Classical vs. Others ###
# Diff gene function for GO enrichment
diffGenes <- function(bg_genes) {
  return(bg_genes == 1)}

goList <- list()
for(i in 1:length(tumor_sigs))
{
  cat("\r", i)
  mysig <- tumor_sigs[[i]]
  up <- rownames(mysig %>% filter(sig, eff < 0))
  
  all_genes <- rownames(mysig)
  bg_genes <- as.numeric(all_genes %in% up)
  names(bg_genes) <- all_genes
  
  # Functional enrichment of signature
  sampleGOdata <- new("topGOdata",
                      description = "Simple session", 
                      ontology = "BP",
                      allGenes = bg_genes,
                      geneSelectionFun = diffGenes, 		# Call function here
                      nodeSize = 10,
                      annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
  
  # Fishers test
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
  fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
  fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
  
  fishRes[,"cell_state"] <- tumor_tag[i]
  goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)

## Output Classical vs. others GO results
write.table(goRes, paste0(output_dir, "/tcga_idhwt_classical_vs_others_sig_GO_result.txt"), sep="\t",quote=FALSE,row.names=TRUE) 


############ Mesenchymal vs. Others   #############
## Apply differential expression to Mesenchymal vs. non-Mesenchymal
p <- se <- sigs <- list()
for(i in 1:length(myinf1))
{
  cat("\r", i)
  geps <- read.delim(myinf1[i], row.names=1)
  geps <- log10(geps+1)
  rem <- apply(geps,1,function(x)sum(is.na(x)))
  geps <- geps[-which(rem==ncol(geps)),]	
  vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
  geps <- geps[which(vars > 0),]
  nrow(geps)
  
  g1 <- geps[,dat$gep_barcode[dat$transcriptional_subtype=="Mesenchymal"]]
  g2 <- geps[,dat$gep_barcode[dat$transcriptional_subtype!="Mesenchymal"]]
  
  p.val <- eff <- rep(0, nrow(geps))
  names(p.val) <- names(eff) <- rownames(geps)
  for(j in 1:nrow(geps))
  {
    group1 <- as.numeric(g1[j,])
    group2 <- as.numeric(g2[j,])
    
    p.val[j] <- wilcox.test(group1,group2,paired=FALSE)$p.value
    eff[j] <- log2(mean(group2)/mean(group1))
  }
  q.val <- p.adjust(p.val,"BH")
  idhwt_res <- data.frame(p.val, q.val, eff)
  idhwt_res <- idhwt_res[order(eff),]
  
  idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
  idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.05
  
  idhwt_res <- idhwt_res[order(idhwt_res$p.val),]
  
  myoutf <- paste(output_dir, "/tcga_idhwt_", cell_state[i],"_mesenchymal_vs_other_result.txt",sep="")
  
  write.table(idhwt_res, myoutf, sep="\t",quote=FALSE,row.names=TRUE) 
  
  sigs[[i]] <- idhwt_res
  
  # Plot heatmaps (unnecessary)
  mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
  sig_matrix <- geps[mygenes,c(dat$gep_barcode[dat$transcriptional_subtype=="Mesenchymal"], dat$gep_barcode[dat$transcriptional_subtype!="Mesenchymal"])]
  
  sig_matrix <- t(apply(sig_matrix, 1, function(x) x - median(x)))
  
  grp1 <- sig_matrix[,dat$gep_barcode[dat$transcriptional_subtype=="Mesenchymal"]]
  grp1 <- apply(grp1, 1, mean)
  grp2 <- sig_matrix[,dat$gep_barcode[dat$transcriptional_subtype!="Mesenchymal"]]
  grp2 <- apply(grp2, 1, mean)
  
  expr <- c(grp1, grp2)
  # Rescale for color
  expr[which(expr > 0.33)] <- 0.33 # Plug the number in from the whole dataset (obtained below) back in for scaling
  gene_symbol <- c(names(grp1), names(grp2))
  subtype <- c(rep("Mesenchymal", length(grp1)), rep("Other", length(grp2)))
  
  plot_hm <- data.frame(gene_symbol, expr, subtype)
  lev <- names(grp1)[order(grp2 - grp1)]
  plot_hm$gene_symbol <- factor(plot_hm$gene_symbol, levels = lev)
  plot_hm[,"cell"] <- cell_state[i]
  
  p[[i]] <- plot_hm
  
  se[[i]] <- ggplot(data = plot_hm, aes(x = subtype, y = gene_symbol)) +
    geom_tile(aes(fill=expr)) +
    scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
                         space = "Lab", na.value = "grey50", guide = "colourbar",
                         aesthetics = "fill", limits = c(-0.33, 0.33)) + 
    theme_void() +
    theme(axis.text = element_blank(),
          axis.title= element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.title = element_blank(),
          legend.position = "none")	
}


###### GO enrichment  ########
tumor_sigs <- sigs
tumor_tag <- mytag


#### Up-regulated in Mesenchymal vs. Others ###
# Diff gene function for GO enrichment
diffGenes <- function(bg_genes) {
  return(bg_genes == 1)}

goList <- list()
for(i in 1:length(tumor_sigs))
{
  cat("\r", i)
  mysig <- tumor_sigs[[i]]
  up <- rownames(mysig %>% filter(sig, eff < 0))
  
  all_genes <- rownames(mysig)
  bg_genes <- as.numeric(all_genes %in% up)
  names(bg_genes) <- all_genes
  
  
  
  # Functional enrichment of signature
  sampleGOdata <- new("topGOdata",
                      description = "Simple session", 
                      ontology = "BP",
                      allGenes = bg_genes,
                      geneSelectionFun = diffGenes, 		# Call function here
                      nodeSize = 10,
                      annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")
  
  # Fishers test
  resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
  fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
  fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
  
  fishRes[,"cell_state"] <- tumor_tag[i]
  goList[[i]] <- fishRes
}


goRes <- do.call(rbind, goList)

## Output Mesenchymal vs. others GO results.
write.table(goRes, paste0(output_dir, "/tcga_idhwt_mesenchymal_vs_others_sig_GO_result.txt"), sep="\t",quote=FALSE,row.names=TRUE) 

#####################
### Read in all the GO results and generate a tiled heatmap with -log10 p-values for comparisons
#####################
proneural <- read.table(paste0(output_dir, "/tcga_idhwt_proneural_vs_others_sig_GO_result.txt"), sep = "\t", header = T)
proneural$group <- paste(proneural$cell_state, "up_proneural", sep="_")
proneural$raw.p.value <- as.numeric(as.character(proneural$raw.p.value))
proneural$go_analysis <- "Proneural"

classical <- read.table(paste0(output_dir, "/tcga_idhwt_classical_vs_others_sig_GO_result.txt"), sep = "\t", header = T)
classical$group <- paste(classical$cell_state, "up_classical", sep="_")
classical$raw.p.value <- as.numeric(as.character(classical$raw.p.value))
classical$go_analysis <- "Classical"

mesenchymal <- read.table(paste0(output_dir, "/tcga_idhwt_mesenchymal_vs_others_sig_GO_result.txt"), sep = "\t", header = T)
mesenchymal$group <- paste(mesenchymal$cell_state, "up_mesenchymal", sep="_")
mesenchymal$raw.p.value <- as.numeric(as.character(mesenchymal$raw.p.value))
mesenchymal$go_analysis <- "Mesenchymal"

all_comparisons <- bind_rows(proneural, 
                             classical,
                             mesenchymal) 

# To make figure interpretable, use the top 2 ontologies for each comparison.
go_terms = all_comparisons %>% 
  group_by(go_analysis) %>%
  top_n(2, -log10(q.value)) 

all_comparisons_filt <- all_comparisons %>% 
  filter(Term%in%go_terms$Term) %>% 
  filter(cell_state!="prolif_stemcell_tumor") %>% 
  mutate(cell_state = recode(cell_state, `stemcell_tumor` = "Stem-like",
                             `differentiated_tumor` = "Diff.-like"))
all_comparisons_filt$Term <- droplevels(all_comparisons_filt$Term)

## Set order of terms so that the terms align with the transcriptional subtype.
term_levels = c("modulation of chemical synaptic transmission", "regulation of trans-synaptic signaling", "cell-cell signaling",
                "regulation of transcription by RNA polymerase II", "regulation of RNA metabolic process",
                "immune system process", "establishment of protein localization to membrane")

all_comparisons_filt$Term <- factor(x = all_comparisons_filt$Term, levels = rev(term_levels))

all_comparisons_filt$GO.ID <- droplevels(all_comparisons_filt$GO.ID)
id_levels = c("GO:0050804", "GO:0099177", "GO:0007267",
              "GO:0006357", "GO:0051252",
              "GO:0002376", "GO:0090150")
all_comparisons_filt$GO.ID <- factor(x = all_comparisons_filt$GO.ID, levels = rev(id_levels))

## Set the levels so that they match with the ordering of the subtypes elsewhere in the manuscript.
compare_levels = c("Proneural", 
                   "Classical", 
                   "Mesenchymal")
all_comparisons_filt$go_analysis <- factor(x = all_comparisons_filt$go_analysis, levels = compare_levels)

## Create final plot.
go_plot <- ggplot(all_comparisons_filt, aes(x=go_analysis, y = GO.ID, fill=-log10(q.value))) +
  geom_tile() + 
  scale_fill_continuous(low = "white",
                         high = "#CD4F39",
                         na.value = "grey50") +
  theme_classic() +
  facet_grid(.~cell_state, scales = "free_y", space = "free_y") +
  theme(axis.text.x = element_text(size=7, angle = 45, hjust=1),
        axis.text.y = element_text(size=7),
        axis.title = element_blank(),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_blank(), 
        strip.text.y = element_blank()) +
  guides(fill=FALSE) +
  labs(x = "") 

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/suppfig4D_enriched_go_id_heatmap.pdf", width=2, height = 2.5)
go_plot
dev.off()

### END ###