library(tidyverse)
library(data.table)


#######################################################
rm(list=ls())

gct_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct"
clin_path <- "/projects/verhaak-lab/USERS/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt"

info <- read.delim(clin_path)
mycase <- info %>% 
		filter(IDH.codel.subtype == "IDHwt") %>%
		dplyr::select(Case) %>%
		.$Case %>%
		as.character()

# Run Qianghu's transcriptional classifier (only needed once)
#runSsGSEAwithPermutation(gct_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/p_result_TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct.txt")
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

sig_sub[,"case_barcode"] <- sapply(strsplit(as.character(sig_sub[,"aliquot_barcode"]),"-"),function(x)paste(x[1:3],collapse="-"))
sig_sub <- sig_sub %>%
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case")) %>%
		   inner_join((transcriptional_subtype %>% filter(signature_name=="Mesenchymal"))[,c("aliquot_barcode","enrichment_score")], by = "aliquot_barcode")
#----------------------------

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/CIBERSORTxGEP_toil_Fractions-Adjusted.txt"
myinf2 <- "/projects/verhaak-lab/varnf/data/xenahub/toil/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"

frac <- read.delim(myinf1,stringsAsFactor=FALSE)
frac <- frac[,1:13]

frac <- frac %>% 
pivot_longer(-Mixture, names_to = "cell_state") %>%
data.frame()
frac[,"Mixture"] <- gsub("\\.","-",frac[,"Mixture"])

gtex_anno <- read.delim(myinf2,stringsAsFactors=FALSE)

# Get brain samples (cortex)
brain <- gtex_anno %>% 
		filter(SMTSD == "Brain - Frontal Cortex (BA9)" | SMTSD == "Brain - Cortex")
brain_map <- brain[,"SMTSD"]
names(brain_map) <- brain$SAMPID

idhmap <- as.character(sig_sub$IDH.codel.subtype)
names(idhmap) <- sig_sub$aliquot_barcode
names(idhmap) <- sapply(strsplit(names(idhmap),"-"), function(x) paste(x[1:4],collapse="-"))
names(idhmap) <- substr(names(idhmap),1,15)

map <- c(brain_map, idhmap)
tissue <- map[frac$Mixture]

frac$tissue <- tissue

plot_res <- frac %>%
filter(!is.na(tissue)) %>%
group_by(cell_state, tissue) %>%
summarise(fraction = mean(value)) %>%
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

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tumor_tissue_fractions.pdf",width=4,height=4)
ggplot(plot_res, aes(x=tissue, y = fraction, fill = cell_state)) +
geom_bar(position="stack", stat="identity") +
theme_bw() +
theme(axis.text.x = element_text(size=7, angle = 45, hjust = 1),
	axis.text.y = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y = element_text(size=7)) +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) 
dev.off()
