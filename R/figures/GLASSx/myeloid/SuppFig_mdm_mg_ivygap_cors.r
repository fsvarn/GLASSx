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

gct_path <- "/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct"
clin_path <- "/projects/verhaak-lab/USERS/varnf/Data/Firehose/TCGA_GBMLGG_clinical_info_2016.txt"

info <- read.delim(clin_path)

# Run Qianghu's transcriptional classifier (only needed once)
#runSsGSEAwithPermutation(gct_path,100)

# Read output of classifier
subtype_ssgsea <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/p_result_TCGA_Firehose_RNASeqV2_LGGGBM_for_CIBERSORTx.gct.txt",stringsAsFactor=FALSE)
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
		   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case"))

##################################################
# Step 2: Calcualte macrophage and microglia score
##################################################

# Read in CIBERSORTx TCGA gene expression profiles
myeloid_gep <- read.delim("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/CIBERSORTxHiRes_TCGA_myeloid_Window48.txt",row.names=1)

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
mdm_gep <- mdm_gep[,-which(mdm_rem == nrow(mdm_gep))]
mdm_score <- apply(mdm_gep[,2:ncol(mdm_gep)], 1, mean)
names(mdm_score) <- myeloid_gep$aliquot_barcode

#MG score
mg_sig <- sigs %>%
		   filter(signature_name == "Microglia") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mg_gep <- myeloid_gep[,c("aliquot_barcode",mg_sig)]
mg_rem <- apply(mg_gep, 2, function(x)sum(is.na(x)))
mg_gep <- mg_gep[,-which(mg_rem == nrow(mg_gep))]
mg_score <- apply(mg_gep[,2:ncol(mg_gep)], 1, mean)
names(mg_score) <- myeloid_gep$aliquot_barcode

plot_dat <- data.frame(sig_sub, "mdm_score" = mdm_score[sig_sub$aliquot_barcode], "mg_score" = mg_score[sig_sub$aliquot_barcode])

##################################################
# Step 3: Test how MDM/MG correlates with IvyGAP structures
##################################################

fraction <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/TCGA_Firehose_RNASeqV2_LGGGBM_ivygap_CIBERSORT_fractions.txt")
fraction$Mixture <- gsub("\\.","-",fraction$Mixture)
fraction <- fraction[,1:5]

fraction_res <- fraction %>%
			pivot_longer(-Mixture, names_to = "cell_state", values_to = "fraction") %>%
			inner_join(plot_dat,  by = c("Mixture" = "aliquot_barcode")) %>%
			rename(aliquot_barcode = Mixture) %>%
			group_by(IDH.status, cell_state) %>%
			summarise(mac_cor = cor(fraction, mdm_score, method="s"),
					  mg_cor = cor(fraction, mg_score, method="s")) %>%
			pivot_longer(ends_with("_cor"), names_to = "cell", values_to = "cor") %>%
			mutate(cell = recode(cell, "mac_cor" = "Macrophages", "mg_cor" = "Microglia")) %>%
			mutate(IDH.status = recode(IDH.status, "WT" = "IDHwt", "Mutant" = "IDHmut")) %>%
			mutate(IDH.status = fct_relevel(IDH.status, "IDHwt","IDHmut")) %>%
			mutate(cell_state = fct_relevel(cell_state, "LE","CT","CTmvp","CTpan"))
			
# Make a barplot of this finding
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tcga_ivygap_mdm_mg_cor_bar.pdf",width=3,height=2)  
ggplot(fraction_res %>% filter(!is.na(IDH.status)), aes(x=cell_state,y=cor)) +
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

