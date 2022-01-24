###################################################
# Create UMAP/boxplots displaying distribution of blood-derived macrophage and microglia gene expression in TCGA myeloid profiles
# Author: Frederick Varn
# Updated: 2022.01.06
# Figures 5A, S5D
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)
library(UpSetR)
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
# Step 2: UMAP of each CIBERSORTx GEP
##################################################

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_TCGA_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

# Plot the UMAP
geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
geps <- log10(geps+1)
gep_list <- t(geps)

rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
geps <- geps[,which(colnames(geps) %in% sig_sub[,"aliquot_barcode"])]

embedding <- umap(t(geps))$layout

plot_res <- embedding
colnames(plot_res) <- c("UMAP1", "UMAP2")
plot_res <- plot_res %>%
	data.frame() %>%
	rownames_to_column(var = "aliquot_barcode") %>%
	inner_join(sig_sub, by = "aliquot_barcode")


# UMAP plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_umap.pdf",width=4,height=2.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
geom_point()  +
theme_bw() +
labs(title = "Myeloid") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.text = element_text(size=7),
legend.title = element_text(size=7)) +
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B")) +
scale_y_continuous(breaks = seq(-4, 4, by = 2))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_umap_legend.pdf",width=4,height=8)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
geom_point()  +
theme_bw() +
labs(title = "Myeloid") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.text = element_text(size=7),
legend.title = element_text(size=7)) +
scale_colour_manual(values=c("#000000", "#A9A9A9", "#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B")) +
scale_y_continuous(breaks = seq(-4, 4, by = 2))
dev.off()

##################################################
# Step 3: Macrophage and microglia plot
##################################################

# Calculate average score for macrophages and microglia 

# Establish connection to db and get macrophage and microglia score
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
sigs <- sigs %>%
		filter(signature_set == "Muller")

myeloid_gep <- t(geps)
case_barcode <- sapply(strsplit(rownames(myeloid_gep),"-"),function(x)paste(x[1:3],collapse="-"))
myeloid_gep <- data.frame(case_barcode, myeloid_gep)

#MDM score
mdm_sig <- sigs %>%
		   filter(signature_name == "Macrophages") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mdm_gep <- myeloid_gep[,c("case_barcode",mdm_sig)]
#mdm_rem <- apply(mdm_gep, 2, function(x)sum(is.na(x)))
#mdm_gep <- mdm_gep[,-which(mdm_rem == nrow(mdm_gep))]
mdm_score <- apply(mdm_gep[,2:ncol(mdm_gep)], 1, mean)
mdm_score <- data.frame(case_barcode, mdm_score)

#MG score
mg_sig <- sigs %>%
		   filter(signature_name == "Microglia") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mg_gep <- myeloid_gep[,c("case_barcode",mg_sig)]
#mg_rem <- apply(mg_gep, 2, function(x)sum(is.na(x)))
#mg_gep <- mg_gep[,-which(mg_rem == nrow(mg_gep))]
mg_score <- apply(mg_gep[,2:ncol(mg_gep)], 1, mean)
mg_score <- data.frame(case_barcode, mg_score)

plot_mdm <- plot_res %>%
			filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])) %>%
		    inner_join(mdm_score, by = "case_barcode") %>%
		    dplyr::rename(myeloid_score = mdm_score) %>%
		    add_column(cell_type = "Macrophages")
plot_mg <- plot_res %>%
			filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])) %>%
		    inner_join(mg_score, by = "case_barcode") %>%
		    dplyr::rename(myeloid_score = mg_score) %>%
		    add_column(cell_type = "Microglia")
plot_score <- rbind(plot_mdm, plot_mg)

colors <- c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B")
names(colors) <- levels(factor(plot_score$signature_name))

# Microglia score and macrophage score
p1 <- ggplot(plot_mdm, aes(UMAP1, UMAP2, color = myeloid_score, shape = IDH.codel.subtype)) + 
geom_point(size=0.5)  +
labs(title = "Macrophages") +
theme_bw() +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_color_gradient(low = "#D3D3D3", high = "#a50000", na.value = NA) +
scale_y_continuous(breaks = seq(-4, 4, by = 2))


p2 <- ggplot(plot_mg, aes(UMAP1, UMAP2, color = myeloid_score, shape = IDH.codel.subtype)) + 
geom_point(size=0.5)  +
labs(title = "Microglia") +
theme_bw() +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_color_gradient(low = "#D3D3D3", high = "#a50000", na.value = NA) +
scale_y_continuous(breaks = seq(-4, 4, by = 2))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_microglia_mdm_umap.pdf",width=1.25,height=2.75)
grid.arrange(p1,p2,nrow=2)
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_microglia_mdm_umap_legend.pdf",width=1.25,height=2.75)
ggplot(plot_mg, aes(UMAP1, UMAP2, color = myeloid_score, shape = IDH.codel.subtype)) + 
geom_point(size=0.5)  +
labs(title = "Microglia") +
theme_bw() +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5)) +
scale_color_gradient(low = "#D3D3D3", high = "#a50000", na.value = NA)
dev.off()


##################################################
# Step 4: Create boxplot showing this more explicitly
##################################################

plot_score <- plot_score %>%
mutate(signature_name = fct_relevel(signature_name, "Proneural", "Proneural,Classical","Mixed","Classical","Proneural,Mesenchymal","Classical,Mesenchymal","Mesenchymal"))
reord_colors <- colors[levels(plot_score$signature_name)]

# Transcriptional subtype
p1 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Macrophages"), aes(signature_name, myeloid_score, fill = signature_name)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
labs(title = "Macrophages") +
#scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=reord_colors) +
coord_cartesian(ylim = c(1.48, 2))

p2 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Microglia"), aes(signature_name, myeloid_score, fill = signature_name)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
labs(title = "Microglia") +
#scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=reord_colors) +
coord_cartesian(ylim = c(1.59, 2.1))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_ts_boxplot.pdf",width=1.75, height=2.56)
grid.arrange(p1,p2,nrow=2)
dev.off()

p1 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Macrophages"), aes(IDH.codel.subtype, myeloid_score, fill = IDH.codel.subtype)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
#scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
labs(title = "Macrophages") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
coord_cartesian(ylim = c(1.48, 2))

p2 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Microglia"), aes(IDH.codel.subtype, myeloid_score, fill = IDH.codel.subtype)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
#scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
labs(title = "Microglia") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none")  +
coord_cartesian(ylim = c(1.59, 2.1))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_idh_boxplot.pdf",width=0.8, height=2.56)
grid.arrange(p1,p2,nrow=2)
dev.off()

# Get p-values
g1 <- plot_score %>% filter(grepl("Mesenchymal",signature_name) & cell_type == "Macrophages") %>% .$myeloid_score
g2 <- plot_score %>% filter(!grepl("Mesenchymal",signature_name) & cell_type == "Macrophages") %>% .$myeloid_score
wilcox.test(g1,g2)

g1 <- plot_score %>% filter(grepl("Mesenchymal",signature_name) & cell_type == "Microglia") %>% .$myeloid_score
g2 <- plot_score %>% filter(grepl("Mixed",signature_name) & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHmut-non-codel" & cell_type == "Microglia") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype == "IDHmut-codel" & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHwt" & cell_type == "Microglia") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype == "IDHmut-non-codel" & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHwt" & cell_type == "Macrophages") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype != "IDHwt" & cell_type == "Macrophages") %>% .$myeloid_score
wilcox.test(g1,g2)

