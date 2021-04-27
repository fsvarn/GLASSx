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
# Step 2: UMAP of each CIBERSORTx GEP
##################################################

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated|t_cell|oligodendrocyte",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_TCGA_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

diff_list <- embedding <- gep_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	colnames(geps) <- paste(mytag[i], colnames(geps), sep="__")
	geps <- log10(geps+1)
	gep_list[[i]] <- t(geps)

	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]


	embedding[[i]] <- umap(t(geps))$layout

}

full_embedding <- do.call(rbind, embedding)

plot_res <- full_embedding
aliquot_barcode <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[2])
cell_state <- sapply(strsplit(rownames(plot_res),"__"),function(x)x[1])
plot_res <- data.frame(aliquot_barcode, plot_res, cell_state)
colnames(plot_res) <- c("aliquot_barcode","UMAP1","UMAP2", "cell_state")
plot_res <- plot_res %>%
	inner_join(sig_sub, by = "aliquot_barcode")
rownames(plot_res) <- paste(plot_res[,"cell_state"], plot_res[,"aliquot_barcode"], sep = "__")

# UMAP plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_facet_umap.pdf",width=7,height=3.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
geom_point()  +
theme_bw() +
facet_wrap(~ cell_state, scales = 'free', nrow = 2) +
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
#scale_colour_manual(values=c("#008A22", "#8A0000", "#00458A"))
dev.off()


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_umap.pdf",width=4,height=2.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid"), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
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
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_idhwt_myeloid_umap.pdf",width=4,height=2.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid", IDH.codel.subtype=="IDHwt"), aes(UMAP1, UMAP2, color = signature_name)) + 
geom_point(shape = "square")  +
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
coord_cartesian()
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_umap_legend.pdf",width=4,height=8)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid"), aes(UMAP1, UMAP2, color = signature_name,shape = IDH.codel.subtype)) + 
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
scale_colour_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
dev.off()

##################################################
# Step 3: Calcualte macrophage and microglia score
##################################################

# Calculate average score for macrophages and microglia 

# Establish connection to db and get macrophage and microglia score
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
sigs <- sigs %>%
		filter(signature_set == "Muller")

myeloid_gep <- gep_list[[2]]
case_barcode <- sapply(strsplit(sapply(strsplit(rownames(myeloid_gep),"__"),function(x)x[2]),"-"),function(x)paste(x[1:3],collapse="-"))
myeloid_gep <- data.frame(case_barcode, myeloid_gep)

#MDM score
mdm_sig <- sigs %>%
		   filter(signature_name == "Macrophages") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mdm_gep <- myeloid_gep[,c("case_barcode",mdm_sig)]
mdm_rem <- apply(mdm_gep, 2, function(x)sum(is.na(x)))
mdm_gep <- mdm_gep[,-which(mdm_rem == nrow(mdm_gep))]
mdm_score <- apply(mdm_gep[,2:ncol(mdm_gep)], 1, mean)
mdm_score <- data.frame(case_barcode, mdm_score)

#MG score
mg_sig <- sigs %>%
		   filter(signature_name == "Microglia") %>%
		   .$gene_symbol %>%
		   intersect(colnames(myeloid_gep))
mg_gep <- myeloid_gep[,c("case_barcode",mg_sig)]
mg_rem <- apply(mg_gep, 2, function(x)sum(is.na(x)))
mg_gep <- mg_gep[,-which(mg_rem == nrow(mg_gep))]
mg_score <- apply(mg_gep[,2:ncol(mg_gep)], 1, mean)
mg_score <- data.frame(case_barcode, mg_score)

plot_mdm <- plot_res %>%
			filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid") %>%
		    inner_join(mdm_score, by = "case_barcode") %>%
		    rename(myeloid_score = mdm_score) %>%
		    add_column(cell_type = "Macrophages")
plot_mg <- plot_res %>%
			filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"]),cell_state == "myeloid") %>%
		    inner_join(mg_score, by = "case_barcode") %>%
		    rename(myeloid_score = mg_score) %>%
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
scale_color_gradient(low = "#D3D3D3", high = "#a50000", na.value = NA)

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
scale_color_gradient(low = "#D3D3D3", high = "#a50000", na.value = NA)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_microglia_mdm_umap.pdf",width=1.25,height=2.75)
grid.arrange(p1,p2,nrow=2)
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_microglia_mdm_umap_legend.pdf",width=1.25,height=2.75)
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
scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=reord_colors) +
coord_cartesian(ylim = c(2.97, 3.3))

p2 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Microglia"), aes(signature_name, myeloid_score, fill = signature_name)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
labs(title = "Microglia") +
scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_fill_manual(values=reord_colors) +
coord_cartesian(ylim = c(3.07, 3.45))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_ts_boxplot.pdf",width=1.75, height=2.56)
grid.arrange(p1,p2,nrow=2)
dev.off()

p1 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Macrophages"), aes(IDH.codel.subtype, myeloid_score, fill = IDH.codel.subtype)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
labs(title = "Macrophages") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") + 
coord_cartesian(ylim = c(2.97, 3.3))

p2 <- ggplot(plot_score %>% filter(!is.na(IDH.codel.subtype) & cell_type == "Microglia"), aes(IDH.codel.subtype, myeloid_score, fill = IDH.codel.subtype)) + 
geom_boxplot(outlier.size = 0.1)  +
theme_bw() +
scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", x))) +
labs(title = "Microglia") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text.x = element_blank(),
axis.text.y = element_text(size=7),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
coord_cartesian(ylim = c(3.07, 3.45))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_tcga_myeloid_idh_boxplot.pdf",width=0.8, height=2.56)
grid.arrange(p1,p2,nrow=2)
dev.off()

# Get p-values
g1 <- plot_score %>% filter(grepl("Mesenchymal",signature_name) & cell_type == "Macrophages") %>% .$myeloid_score
g2 <- plot_score %>% filter(!grepl("Mesenchymal",signature_name) & cell_type == "Macrophages") %>% .$myeloid_score
wilcox.test(g1,g2)

g1 <- plot_score %>% filter(grepl("Proneural",signature_name) & cell_type == "Microglia") %>% .$myeloid_score
g2 <- plot_score %>% filter(grepl("Mixed",signature_name) & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHmut-non-codel" & cell_type == "Microglia") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype == "IDHmut-codel" & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHwt" & cell_type == "Microglia") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype == "IDHmut-non-codel" & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)	# P = 0.08

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHwt" & cell_type == "Macrophages") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype != "IDHwt" & cell_type == "Macrophages") %>% .$myeloid_score
wilcox.test(g1,g2)

##################################################
# Step 5: Test how MDM/MG correlates with endothelial cell fraction
##################################################

myinf2 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/CIBERSORTxGEP_TCGA_Fractions-Adjusted.txt"

fraction <- read.delim(myinf2)
fraction$Mixture <- gsub("\\.","-",fraction$Mixture)
fraction <- fraction[,1:13]

fraction_res <- fraction %>%
			pivot_longer(-Mixture, names_to = "cell_state", values_to = "fraction") %>%
			inner_join((plot_score %>% 
						filter(cell_type == "Macrophages") %>% 
						select(aliquot_barcode, myeloid_score, IDH.status)), 
			by = c("Mixture" = "aliquot_barcode")) %>%
			inner_join((plot_score %>% 
						filter(cell_type == "Microglia") %>% 
						select(aliquot_barcode, myeloid_score)), 
			by = c("Mixture" = "aliquot_barcode")) %>%
			rename(aliquot_barcode = Mixture, macrophage_score = myeloid_score.x, microglia_score = myeloid_score.y) %>%
			filter(!grepl("tumor",cell_state)) %>%
			group_by(aliquot_barcode) %>%
			summarise(cell_state, fraction = fraction/sum(fraction), IDH.status, macrophage_score, microglia_score) %>%
			ungroup() %>%
			group_by(IDH.status, cell_state) %>%
			summarise(mac_cor = cor(fraction, macrophage_score, method="s"),
					  mg_cor = cor(fraction, microglia_score, method="s")) %>%
			data.frame()
# IvyGAP fraction?			
