###################################################
# Correlate Ivy GAP feature abundance with myeloid cell expression signatures
# Author: Frederick Varn
# Date: 2022.01.22
# Figures S5F, S5K
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(gridExtra)
library(ssgsea.GBM.classification)
library(umap)

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

tcga_sub[,"case_barcode"] <- sapply(strsplit(as.character(tcga_sub[,"aliquot_barcode"]),"-"),function(x)paste(x[1:3],collapse="-"))
tcga_sub <- tcga_sub %>%
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case")) %>%
		   inner_join((transcriptional_subtype %>% filter(signature_name=="Mesenchymal"))[,c("aliquot_barcode","enrichment_score")], by = "aliquot_barcode")


# GTEx-specific modifications
gtex_sub <- sig_sub %>% filter(grepl("GTEX",aliquot_barcode))
gtex_case <- gtex_sub$aliquot_barcode

# Add GTEx info to this table
gtex_anno_path <- "/projects/verhaak-lab/varnf/data/xenahub/toil/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
gtex_anno <- read.delim(gtex_anno_path,stringsAsFactors=FALSE)

# Get brain samples (cortex)
brain <- gtex_anno %>% 
		filter(SMTSD == "Brain - Frontal Cortex (BA9)" | SMTSD == "Brain - Cortex")
brain_map <- brain[,"SMTSD"]
names(brain_map) <- brain$SAMPID

# Change tumor-specific variables to tissue source site because we don't care about glioma subtypes in normal brain
gtex_sub[,"signature_name"] <- brain_map[gtex_sub$aliquot_barcode]
gtex_sub[,"case_barcode"] <- gtex_case
gtex_sub[,"IDH.codel.subtype"] <- brain_map[gtex_sub$aliquot_barcode]
gtex_sub[,"IDH.status"] <- brain_map[gtex_sub$aliquot_barcode]
gtex_sub[,"enrichment_score"] <-NA

# Combine everything together
sig_sub <- rbind(tcga_sub, gtex_sub)

##################################################
# Step 2: Calculate macrophage and microglia score and mesenchymal score (for revision)
##################################################

# Read in CIBERSORTx TCGA gene expression profiles
myeloid_gep <- read.delim("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/CIBERSORTxHiRes_toil_myeloid_Window48.txt",row.names=1)

colnames(myeloid_gep) <- gsub("\\.","-",colnames(myeloid_gep))
myeloid_gep <- log10(myeloid_gep+1)
myeloid_gep <- t(myeloid_gep) %>%
			   as.data.frame() %>%
			   rownames_to_column(var = "aliquot_barcode")


# Calculate average score for macrophages and microglia 

# Establish connection to db and get macrophage and microglia score
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
sigs <- sigs %>%
		filter(signature_set == "Muller")

#MDM score
mdm_sig <- sigs %>%
		   filter(signature_name == "Macrophages") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mdm_gep <- myeloid_gep[,c("aliquot_barcode",mdm_sig)]
mdm_rem <- apply(mdm_gep, 2, function(x)sum(is.na(x)))
#mdm_gep <- mdm_gep[,-which(mdm_rem == nrow(mdm_gep))]
mdm_score <- apply(mdm_gep[,2:ncol(mdm_gep)], 1, mean)
names(mdm_score) <- myeloid_gep$aliquot_barcode

#MG score
mg_sig <- sigs %>%
		   filter(signature_name == "Microglia") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mg_gep <- myeloid_gep[,c("aliquot_barcode",mg_sig)]
mg_rem <- apply(mg_gep, 2, function(x)sum(is.na(x)))
#mg_gep <- mg_gep[,-which(mg_rem == nrow(mg_gep))]
mg_score <- apply(mg_gep[,2:ncol(mg_gep)], 1, mean)
names(mg_score) <- myeloid_gep$aliquot_barcode

#Mesenchymal score (revision)
tcga_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v4.txt")
mes_sig <- rownames(tcga_dat %>% filter(q.val < 0.05 & eff > log2(1.5)))

mes_gep <- myeloid_gep[,c("aliquot_barcode",mes_sig)]
mes_rem <- apply(mes_gep, 2, function(x)sum(is.na(x)))	#No NAs
#mes_gep <- mes_gep[,-which(mes_rem == nrow(mes_gep))]
mes_score <- apply(mes_gep[,2:ncol(mes_gep)], 1, mean)
names(mes_score) <- myeloid_gep$aliquot_barcode

plot_dat <- data.frame(sig_sub, "mdm_score" = mdm_score[sig_sub$aliquot_barcode], "mg_score" = mg_score[sig_sub$aliquot_barcode], "mes_score" = mes_score[sig_sub$aliquot_barcode])

##################################################
# Step 3: Test how MDM/MG/MES correlates with IvyGAP structures
##################################################

fraction <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/toil/CIBERSORT_toil_ivygap_fractions.txt")
fraction$Mixture <- gsub("\\.","-",fraction$Mixture)
fraction <- fraction[,1:5]

fraction_res <- fraction %>%
			pivot_longer(-Mixture, names_to = "cell_state", values_to = "fraction") %>%
			inner_join(plot_dat,  by = c("Mixture" = "aliquot_barcode")) %>%
			#rename(aliquot_barcode = Mixture) %>%
			group_by(IDH.status, cell_state) %>%
			summarise(mac_cor = cor(fraction, mdm_score, method="s"),
					  mg_cor = cor(fraction, mg_score, method="s"),
					  mes_cor = cor(fraction, mes_score, method="s")) %>%
			pivot_longer(ends_with("_cor"), names_to = "cell", values_to = "cor") %>%
			mutate(cell = recode(cell, "mac_cor" = "Macrophages", "mg_cor" = "Microglia","mes_cor" = "Mesenchymal")) %>%
			mutate(IDH.status = recode(IDH.status, "WT" = "IDHwt", "Mutant" = "IDHmut")) %>%
			mutate(IDH.status = fct_relevel(IDH.status, "IDHwt","IDHmut")) %>%
			mutate(cell_state = fct_relevel(cell_state, "LE","CT","CTmvp","CTpan")) %>%
			mutate(cell = fct_relevel(cell, "Macrophages", "Microglia", "Mesenchymal"))
			
# Make a barplot of this finding
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/toil_ivygap_mdm_mg_mes_cor_bar.pdf",width=2.4,height=2.7)  
ggplot(fraction_res %>% filter(!is.na(IDH.status), IDH.status %in% c("IDHmut","IDHwt")), aes(x=cell_state,y=cor)) +
geom_bar(stat = "identity") +
geom_hline(yintercept = 0, colour = "black") +
labs(y = "Spearman correlation") +
facet_grid(cell ~IDH.status) +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.title=element_text(size=7,hjust=0.5),
axis.title.y=element_text(size=7),
axis.title.x = element_blank(),
axis.text.x=element_text(size=7, angle = 45, hjust = 1),
axis.text.y=element_text(size=7),
strip.background=element_blank(),
strip.text = element_text(size=7),
legend.position="none")			
dev.off()		

