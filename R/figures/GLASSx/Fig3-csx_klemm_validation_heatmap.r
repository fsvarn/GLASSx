###################################################
# Create a macrophage expression profile correlation heatmap
# Updated: 2020.07.31
# Author: Frederick Varn
##################################################

library(tidyverse)
library(DBI)
library(odbc)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-v3/"

mytag<- dir(myDir1)
mytag <- mytag[grep("Window40.txt",mytag)]
myinf1 <- paste(myDir1, mytag, sep="/")
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", mytag)
mytag <- gsub("_Window40.txt", "", mytag)

##################################################
# Prepare Klemm dataset
##################################################

myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_RawCounts.csv"
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/BrainTIME_ClinicalAnnotation.csv"

counts <- read.csv(myinf2,row.names=1)
info <- read.csv(myinf3,row.names=1,stringsAsFactor=FALSE)

# Glioma or normal samples:
counts <- counts[,c(grep("_glioma_|_nonTumor_",colnames(counts)))]
cpm <- apply(counts,2, function(x) (x/sum(x))*1000000) 
cpm <- log2(cpm + 1)

##################################################
# Correlation matrices
##################################################

cor_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	myfile <- read.delim(myinf1[i],row.names=1)
	myfile <- log2(myfile)
	
	comxx <- intersect(rownames(myfile), rownames(cpm))
	myfile <- myfile[comxx,]
	mycpm <- cpm[comxx,]
	
	mycor <- cor(myfile, mycpm, method="p", use = "complete.obs")
	rownames(mycor) <- paste(rownames(mycor), mytag[i], sep="_")
	
	cell_state <- mytag[i]
	mycor <- data.frame(mycor, cell_state)
	cor_list[[i]] <- mycor
}

full_cor <- do.call(rbind, cor_list)
rownames(full_cor) <- gsub("\\.","-",rownames(full_cor))

##################################################
# Build the IDHwt matrix
##################################################

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps1.rna_barcode_a
FROM analysis.platinum_set ps1
JOIN clinical.subtypes cs ON cs.case_barcode = ps1.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
UNION
SELECT ps2.rna_barcode_b
FROM analysis.platinum_set ps2
JOIN clinical.subtypes cs ON cs.case_barcode = ps2.case_barcode
WHERE cs.idh_codel_subtype = 'IDHwt'
"
aliquots <- dbGetQuery(con, q)[,1]

idhwt_klemm <- info %>%
			   filter(IDH.Status == "wild type")

indices <- sapply(strsplit(colnames(full_cor),"_"),function(x) paste(x[1:3],collapse="_"))

idhwt_pull <- which(sapply(strsplit(rownames(full_cor),"_"),function(x)x[1]) %in% aliquots)
idhwt_cor <- full_cor[idhwt_pull, which(indices %in% rownames(idhwt_klemm) | colnames(full_cor) == "cell_state")]

cells <- unique(sapply(strsplit(colnames(idhwt_cor[,1:(ncol(idhwt_cor)-1)]),"_"),function(x) paste(x[4])))

idhwt_final <- matrix(0, nrow = length(unique(idhwt_cor[,"cell_state"])), ncol = length(cells))
colnames(idhwt_final) <- cells 
for(i in 1:length(cells))
{
	sub_idhwt_cor <- idhwt_cor[,which(grepl(cells[i], colnames(idhwt_cor)) | colnames(idhwt_cor) == "cell_state")]
	
	mean_res <- sub_idhwt_cor %>%
				group_by(cell_state) %>%
				summarise_all(mean) %>%
				data.frame() %>%
				column_to_rownames(var = "cell_state")
				
	meanomeans <- apply(mean_res, 1, mean)
	idhwt_final[,i] <- meanomeans
}
rownames(idhwt_final) <- rownames(mean_res)


##################################################
# Build the IDHmut matrix
##################################################

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps1.rna_barcode_a
FROM analysis.platinum_set ps1
JOIN clinical.subtypes cs ON cs.case_barcode = ps1.case_barcode
WHERE cs.idh_codel_subtype LIKE 'IDHmut%'
UNION
SELECT ps2.rna_barcode_b
FROM analysis.platinum_set ps2
JOIN clinical.subtypes cs ON cs.case_barcode = ps2.case_barcode
WHERE cs.idh_codel_subtype LIKE 'IDHmut%'
"
aliquots <- dbGetQuery(con, q)[,1]

idhmut_klemm <- info %>%
			   filter(IDH.Status == "mutated")

indices <- sapply(strsplit(colnames(full_cor),"_"),function(x) paste(x[1:3],collapse="_"))

idhmut_pull <- which(sapply(strsplit(rownames(full_cor),"_"),function(x)x[1]) %in% aliquots)
idhmut_cor <- full_cor[idhmut_pull, which(indices %in% rownames(idhmut_klemm) | colnames(full_cor) == "cell_state")]

cells <- unique(sapply(strsplit(colnames(idhmut_cor[,1:(ncol(idhmut_cor)-1)]),"_"),function(x) paste(x[4])))

idhmut_final <- matrix(0, nrow = length(unique(idhmut_cor[,"cell_state"])), ncol = length(cells))
colnames(idhmut_final) <- cells 
for(i in 1:length(cells))
{
	sub_idhmut_cor <- idhmut_cor[,which(grepl(cells[i], colnames(idhmut_cor)) | colnames(idhmut_cor) == "cell_state")]
	
	mean_res <- sub_idhmut_cor %>%
				group_by(cell_state) %>%
				summarise_all(mean) %>%
				data.frame() %>%
				column_to_rownames(var = "cell_state")
				
	meanomeans <- apply(mean_res, 1, mean)
	idhmut_final[,i] <- meanomeans
}
rownames(idhmut_final) <- rownames(mean_res)


##################################################
# Plot the results
##################################################

plot_idhwt <- idhwt_final %>%
			  data.frame() %>%
			  rownames_to_column(var = "cell_state") %>%
			  pivot_longer(-cell_state) %>%
			  add_column(idh_codel_subtype = "IDHwt")
			  
plot_idhmut <- idhmut_final %>%
			  data.frame() %>%
			  rownames_to_column(var = "cell_state") %>%
			  pivot_longer(-cell_state) %>%
			  add_column(idh_codel_subtype = "IDHmut")
			  
plot_res <- rbind(plot_idhwt, plot_idhmut) %>%
			mutate(cell_state = recode(cell_state, "tumor" = "Tumor", "myeloid" = "Myeloid", "t_cell" = "T cell")) %>%
			mutate(name = recode(name, "cd45n" = "CD45-", "mg" = "Microglia", "mdm" = "Macrophage", "neutrophils" = "Neutrophil", "cd4" = "CD4+ T cell", "cd8" = "CD8+ T cell")) %>%
			data.frame()
			
plot_res[,"cell_state"] <- factor(plot_res[,"cell_state"], levels = rev(c("Tumor", "Myeloid", "T cell", "b_cell", "dendritic_cell", "endothelial", "fibroblast", "granulocyte", "oligodendrocyte", "pericyte")))
plot_res[,"name"] <- factor(plot_res[,"name"], levels = rev(c("Neutrophil","CD4+ T cell", "CD8+ T cell", "Microglia", "Macrophage", "CD45-")))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/csx_klemm_validation_heatmap.pdf",width=2.4,height=0.95)
ggplot(data = plot_res %>% filter(cell_state %in% c("Tumor", "Myeloid")), aes(x = name, y = cell_state)) +
  geom_tile(aes(fill=value)) +
  scale_fill_gradient2(low ="white", high = "#260000", limits = c(0,0.85),
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  facet_grid(.~idh_codel_subtype, scales = "free", space = "free") + 
  theme_void() +
  theme(axis.text.x = element_text(size=7, angle = 90, hjust = 1),
  	axis.text.y = element_text(size=7),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.title = element_blank(),
 	legend.key.size = unit(0.25, "cm"))
dev.off()

