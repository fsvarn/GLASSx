library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

#######################################################
rm(list=ls())

# Core vs periphery analysis (includes neurons/astrocytes)
myinf1a <- "/projects/verhaak-lab/GLASS-III/data/dataset/darmanis_2017/GBM_normalized_gene_counts.csv"
myinf1b <- "/projects/verhaak-lab/GLASS-III/data/dataset/darmanis_2017/GBM_metadata.csv"

expr <- read.csv(myinf1a,sep=" ")
colnames(expr) <- substr(colnames(expr),2,nchar(colnames(expr)))
expr <- expr

info <- read.csv(myinf1b,sep=" ")

g1 <- rownames(info %>% filter(Selection == "Unpanned"))
g2 <- rownames(info %>% filter(Selection == "Neurons(Thy1)"))
g3 <- rownames(info %>% filter(Selection == "Astrocytes(HEPACAM)"))

#######################################################


cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)


#  Read in post-treatment stem cell signature
myDir2 <- "data/res/CIBERSORTx/analysis/"
myinf2 <- dir(myDir2)
myinf2 <- myinf2[grep("_postreatment_result", myinf2)]
myinf2 <- myinf2[grep("myeloid|stemcell|differentiated",myinf2)]		# These are the signatures we have some confidence in
myinf2 <- paste(myDir2, myinf2, sep = "/")

sigtable <- read.delim(myinf2[4])	# 4 = stem cells
mysig <- rownames(sigtable %>% filter(sig, eff > 0))


p.val <- apply(expr, 1, function(x)wilcox.test(x[g1], x[g2])$p.value)
eff <- apply(expr, 1, function(x)log2(mean(x[g1])/(mean(x[g2]+0.01))))
q.val <- p.adjust(p.val, "BH")

res <- data.frame("gene_symbol" = rownames(expr), p.val, q.val, eff)
res <- res[order(res$p.val),]

neursig <- res %>% filter(q.val < 0.1)

overlap1 <- intersect(mysig, rownames(neursig))

neures <- res[overlap1,]

p.val <- apply(expr, 1, function(x)wilcox.test(x[g1], x[g3])$p.value)
eff <- apply(expr, 1, function(x)log2(mean(x[g1])/(mean(x[g3]+0.01))))
q.val <- p.adjust(p.val, "BH")

res <- data.frame("gene_symbol" = rownames(expr), p.val, q.val, eff)
res <- res[order(res$p.val),]

astrosig <- res %>% filter(q.val < 0.1)

overlap2<- intersect(mysig, rownames(astrosig))

astrores <- res[overlap2,]

# Normal genes
sum(neures$eff < 0)			#130
sum(astrores$eff < 0)		#90
length(intersect(rownames(astrores[which(astrores$eff<0),]), rownames(neures[which(neures$eff<0),]))) #73

# Tumor genes
sum(neures$eff > 0)			#7
sum(astrores$eff > 0)		#50
length(intersect(rownames(astrores[which(astrores$eff>0),]), rownames(neures[which(neures$eff>0),]))) #5

# Within 404 gene stem cell signature:
# 147 genes neuron/astro-specific
# 52 genes tumor-specific
# 205 genes unknown

newsig <-  mysig[-which(mysig %in% rownames(astrores) | mysig %in% rownames(neures))] 
all_genes <- rownames(sigtable)
bg_genes <- as.numeric(all_genes %in% newsig)
names(bg_genes) <- all_genes

diffGenes <- function(bg_genes) {
	return(bg_genes == 1)}

# Functional enrichment of signature
sampleGOdata <- new("topGOdata",
				description = "Simple session", 
				ontology = "BP",
				allGenes = bg_genes,
				geneSelectionFun = diffGenes, 
				nodeSize = 10,
				annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

# Fishers test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")
fishRes[which(fishRes$q.value < 0.05 & fishRes$Significant > fishRes$Expected),]

#######################################################

#  Read in post-treatment differentiated cell signature
myinf2 <- paste("data/res/CIBERSORTx/analysis/GLASS_", cell_state,"_postreatment_result.txt",sep="")
sigtable <- read.delim(myinf2[1])	# 4 = differentiated cells
mysig <- rownames(sigtable %>% filter(sig, eff > 0))

p.val <- apply(expr, 1, function(x)wilcox.test(x[g1], x[g2])$p.value)
eff <- apply(expr, 1, function(x)log2(mean(x[g1])/(mean(x[g2]+0.01))))
q.val <- p.adjust(p.val, "BH")

res <- data.frame("gene_symbol" = rownames(expr), p.val, q.val, eff)
res <- res[order(res$p.val),]

neursig <- res %>% filter(q.val < 0.05)

overlap <- intersect(mysig, rownames(neursig))

neures <- res[overlap,]

p.val <- apply(expr, 1, function(x)wilcox.test(x[g1], x[g3])$p.value)
eff <- apply(expr, 1, function(x)log2(mean(x[g1])/(mean(x[g3]+0.01))))
q.val <- p.adjust(p.val, "BH")

res <- data.frame("gene_symbol" = rownames(expr), p.val, q.val, eff)
res <- res[order(res$p.val),]

astrosig <- res %>% filter(q.val < 0.05)

overlap <- intersect(mysig, rownames(astrosig))

astrores <- res[overlap,]

# Normal genes
sum(neures$eff < 0)			#55
sum(astrores$eff < 0)		#48
length(intersect(rownames(astrores[which(astrores$eff<0),]), rownames(neures[which(neures$eff<0),]))) #31

# Tumor genes
sum(neures$eff > 0)			#8
sum(astrores$eff > 0)		#42
length(intersect(rownames(astrores[which(astrores$eff>0),]), rownames(neures[which(neures$eff>0),]))) #6

# Within 232 gene differentiated cell signature:
# 72 genes neuron/astro-specific
# 44 genes tumor-specific
# 116 genes unknown

# Overlap in unknown between each signature?





# Infiltrating glioma signature
presig <- c("ATP1A2","FGFR3","LMO3","NCAN","FXYD1","PSD2","PRODH","HIF3A","HRSP12","KCNN3","PPM1K","KCNJ10","ADCYAP1R1","BMP7","KAT2B","CNTN1","SAMD9L","SLC7A11","ECHDC2","FAM181B","SALL2","SASH1")


mysig <- rownames(sigtable %>% filter(sig, eff > 0))

intersect(presig,mysig)

# Minimal overlap