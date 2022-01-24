#######################################################
# Compare the GEPs between TCGA IDHwt tumors across transcriptional subtypes
# Author: Kevin Johnson
# Date: 2021.10.21
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
library(broom)

## Address the following question using TCGAxGEP:
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
sub_dir <- "glassx"
output_dir <- file.path(main_dir, sub_dir)

if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Directory already exists!")
}

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

## Load the toil batch corrected fractions.
tcga_toil_prop <- read.delim("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/CIBERSORTxGEP_toil_Fractions-Adjusted.txt", header = TRUE)

tcga_toil_prop_long <- tcga_toil_prop %>% 
  filter(grepl("TCGA", Mixture)) %>% 
  dplyr::select(-c(P.value, Correlation, RMSE)) %>% 
  pivot_longer(cols = stemcell_tumor:b_cell, 
             names_to = "cell_state",
             values_to = "fraction") %>% 
  mutate(aliquot_barcode = gsub("\\.", "-", Mixture)) %>% 
  dplyr::select(aliquot_barcode, cell_state, fraction) %>% 
  inner_join(sig_sub, by="aliquot_barcode") %>% 
  dplyr::rename(transcriptional_subtype = signature_name) %>%
  dplyr::rename(sample_barcode = aliquot_barcode) %>%
  filter(transcriptional_subtype %in% c("Mesenchymal", "Classical", "Proneural")) %>%
  mutate(idh_status = recode(IDH.status, `WT` = "IDHwt",
                             `Mutant` = "IDHmut")) %>% 
  mutate(cell_state = recode(cell_state, `stemcell_tumor` = "Stem-like",
                             `differentiated_tumor` = "Diff.-like",
                             `prolif_stemcell_tumor` = "Prolif. stem-like",
                             `endothelial` = "Endothelial",
                             `pericyte` = "Pericyte",
                             `granulocyte` = "Granulocyte",
                             `myeloid` = "Myeloid",
                             `t_cell` = "T cell",
                             `fibroblast` = "Fibroblast",
                             `oligodendrocyte` = "Oligodendrocyte",
                             `b_cell` = "B cell",
                             `dendritic_cell` = "Dendritic cell"))

my_levels = c("Prolif. stem-like", "Stem-like", "Diff.-like", 
              "Oligodendrocyte",
              "Fibroblast",
              "Endothelial", "Pericyte",
              "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid")
tcga_toil_prop_long$cell_state <- factor(x = tcga_toil_prop_long$cell_state, levels = rev(my_levels))

tcga_toil_prop_long %>% 
  dplyr::select(idh_status, transcriptional_subtype, sample_barcode) %>% 
  distinct() %>% 
  group_by(idh_status, transcriptional_subtype) %>% 
  summarise(sample_size = n())

tcga_toil_prop_long_prop <- tcga_toil_prop_long %>% 
  group_by(transcriptional_subtype, idh_status, cell_state) %>% 
  summarise(avg_fractions = mean(fraction)*100,
            sample_count = n()) %>% 
  filter(!is.na(idh_status))
tcga_toil_prop_long_prop$transcriptional_subtype <- factor(x = tcga_toil_prop_long_prop$transcriptional_subtype, levels = c("Classical", "Mesenchymal", "Proneural"))
tcga_toil_prop_long_prop$idh_status <- factor(x = tcga_toil_prop_long_prop$idh_status, levels = c("IDHwt", "IDHmut"))

kw_results <- tcga_toil_prop_long %>% 
  filter(!is.na(idh_status)) %>% 
  group_by(idh_status, cell_state) %>% 
  summarise(fit = list(kruskal.test(fraction ~ transcriptional_subtype) %>% tidy)) %>% 
  unnest_wider(fit) %>% 
  mutate(q.value = p.adjust(p.value,"BH")) %>% 
  mutate(sig = ifelse(q.value < 0.05, TRUE, FALSE))

tcga_toil_prop_long_prop_wt <- tcga_toil_prop_long_prop %>% filter(idh_status=="IDHwt")
tcga_toil_prop_long_prop_mut <- tcga_toil_prop_long_prop %>% filter(idh_status=="IDHmut")


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/figure-1.5.1A.pdf", width = 4, height = 5, useDingbats = FALSE)
ggplot(tcga_toil_prop_long_prop_wt, aes(x = transcriptional_subtype, y = avg_fractions, fill=cell_state)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                             "Oligodendrocyte" = "#2ca25f",
                             "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                             "Fibroblast" = "#feb24c",
                             "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  labs(y = "Proportion (%)", x = "TCGA bulk RNA subtype (IDHwt)", fill = "Pan-glioma\ncell state") +
  plot_theme + theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/figure-1.5.1B.pdf", width = 4, height = 5, useDingbats = FALSE)
ggplot(tcga_toil_prop_long_prop_mut, aes(x = transcriptional_subtype, y = avg_fractions, fill=cell_state)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
                             "Oligodendrocyte" = "#2ca25f",
                             "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
                             "Fibroblast" = "#feb24c",
                             "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
  labs(y = "Proportion (%)", x = "TCGA bulk RNA subtype (IDHmut)", fill = "Pan-glioma\ncell state") +
  plot_theme + theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1))
dev.off()
### END ###