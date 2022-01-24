###################################################
# Compare CIBERSORTx profiles across transcriptional subtypes across TCGA
# Author: Frederick Varn
# Date: 2021.12.29
# Figure S1E
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(gridExtra)


#######################################################
rm(list=ls())
set.seed(11)
##################################################
# Step 1: Subtype each TCGA sample
##################################################

gct_path <- "/projects/verhaak-lab/varnf/data/xenahub/toil/TcgaTargetGtex_GBMLGGBRAIN_rsem_gene_symbol_tpm_antilog.gct"
clin_path <- "/projects/verhaak-lab/USERS/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt"

info <- read.delim(clin_path)

# Run Qianghu's transcriptional classifier (only needed once)
#runSsGSEAwithPermutation(gct_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/varnf/data/xenahub/toil/p_result_TcgaTargetGtex_GBMLGGBRAIN_rsem_gene_symbol_tpm_antilog.gct.txt")
rownames(subtype_ssgsea) <- gsub("\\.","-",rownames(subtype_ssgsea))

aliquot_barcode <- rep(rownames(subtype_ssgsea),3)
signature_name <- c(rep("Proneural",nrow(subtype_ssgsea)),rep("Classical",nrow(subtype_ssgsea)),rep("Mesenchymal",nrow(subtype_ssgsea)))
enrichment_score <- c(subtype_ssgsea[,"Proneural"],subtype_ssgsea[,"Classical"],subtype_ssgsea[,"Mesenchymal"])
p_value <- c(subtype_ssgsea[,"Proneural_pval"],subtype_ssgsea[,"Classical_pval"],subtype_ssgsea[,"Mesenchymal_pval"])

transcriptional_subtype <- data.frame(aliquot_barcode,signature_name,enrichment_score,p_value)
transcriptional_subtype <- transcriptional_subtype[order(transcriptional_subtype[,"aliquot_barcode"]),]

transcriptional_subtype[,"signif"] <- transcriptional_subtype[,"p_value"] < 0.05

sig_sub <- transcriptional_subtype %>%
		   group_by(aliquot_barcode, signif) %>%
		   summarise(signature_name = paste(signature_name, collapse=",")) %>%
		   data.frame()
# Pull out significant subtypes first
sig_sub1 <- sig_sub %>%
		    filter(signif)
signif_ali <- sig_sub1[,"aliquot_barcode"]

sig_sub2 <- sig_sub %>%
			filter(!(aliquot_barcode %in% signif_ali) & !signif)
sig_sub <- rbind(sig_sub1, sig_sub2)
  
sig_sub[which(sig_sub[,"signature_name"] == "Proneural,Classical,Mesenchymal"),"signature_name"] <- "Mixed"


# TCGA-specific modifications
tcga_sub <- sig_sub %>% filter(grepl("TCGA",aliquot_barcode))
tcga_sub[,"case_barcode"] <- sapply(strsplit(as.character(tcga_sub[,"aliquot_barcode"]),"-"),function(x)paste(x[1:3],collapse="-"))
tcga_sub <- tcga_sub %>%
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case")) %>%
		   inner_join((transcriptional_subtype %>% filter(signature_name=="Mesenchymal"))[,c("aliquot_barcode","enrichment_score")], by = "aliquot_barcode")

# Combine everything together
sig_sub <- tcga_sub

##################################################
# Step 2: Link to fractions
##################################################

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/CIBERSORTxGEP_toil_Fractions-Adjusted.txt"

frac <- read.delim(myinf1)

dat <- frac %>%
	   select(-c("P.value", "Correlation", "RMSE")) %>%
	   pivot_longer(-Mixture, names_to = "cell_state", values_to = "fraction") %>%
	   filter(grepl("TCGA", Mixture)) %>%
	   mutate(Mixture = gsub("\\.","-", Mixture)) %>%
	   inner_join(sig_sub, by = c("Mixture" = "aliquot_barcode")) %>%
	   filter(!grepl(",",signature_name)) %>%
	   filter(signature_name != "Mixed") %>%
	   group_by(signature_name, cell_state) %>%
	   summarise(fraction = mean(fraction))  %>% 
	   mutate(fraction = fraction*100) %>%
	   mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
			  		   				"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
			  		   				"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
			  		   				"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
			  		   				"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
			  		   				"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
	   mutate(cell_state = as_factor(cell_state)) %>%
	   mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
												 "Oligodendrocyte", 
												 "Endothelial", "Pericyte",
												 "Fibroblast", 
												 "Diff.-like", "Stem-like", "Prolif. stem-like"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cibersortx_stacked_barplot_ts_subtype_tcga.pdf",width=1,height=1.82)
ggplot(dat, aes(fill=cell_state, y=fraction, x=signature_name)) + 
geom_bar(position="stack", stat="identity") +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
						 "Fibroblast" = "#feb24c",
						 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
labs(y = "Proportion") +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "none")
dev.off()


