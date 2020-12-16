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
g2 <- plot_score %>% filter(IDH.codel.subtype == "IDHmut-codel" & cell_type == "Microglia") %>% .$myeloid_score
wilcox.test(g1,g2)	# P = 0.08

g1 <- plot_score %>% filter(IDH.codel.subtype == "IDHwt" & cell_type == "Macrophages") %>% .$myeloid_score
g2 <- plot_score %>% filter(IDH.codel.subtype != "IDHwt" & cell_type == "Macrophages") %>% .$myeloid_score
wilcox.test(g1,g2)


##################################################
# Step 5: Validate with single cell profiles (not the ones in the original signature)
##################################################

myeloid_gep <- gep_list[[2]]
case_barcode <- sapply(strsplit(sapply(strsplit(rownames(myeloid_gep),"__"),function(x)x[2]),"-"),function(x)paste(x[1:3],collapse="-"))
myeloid_gep <- data.frame(case_barcode, myeloid_gep)
myeloid_gep <- myeloid_gep %>% 
			   inner_join(info[,c("Case", "IDH.codel.subtype", "IDH.status")], by = c("case_barcode"="Case"))

rem <- apply(myeloid_gep, 2, function(x)sum(is.na(x)))
myeloid_gep <- myeloid_gep[,-which(rem == nrow(myeloid_gep))]

# Get differentially expressed genes in IDHmut vs IDHwt
diff_idh <- myeloid_gep %>%
			group_by(IDH.status) %>%
			select(2:(ncol(.)-2)) %>%
			data.frame()

idh_status <- diff_idh[,"IDH.status"]
idh_pval <- apply(diff_idh[,2:ncol(diff_idh)], 2, function(x) wilcox.test(x[which(diff_idh[,"IDH.status"]=="WT")], x[which(diff_idh[,"IDH.status"]=="Mutant")])$p.value)
idh_qval <- p.adjust(idh_pval, "BH")

# Take top 500 since there are > 3000 significant genes
idh_sig <- idh_qval[which(idh_qval < 0.05)]
mygenes <- names(idh_sig)
idh_ord <- idh_qval[order(idh_qval)]
idh_ord <- idh_ord[1:500]


# Location of the SCGP single cell data
myinf1 <- "/projects/verhaak-lab/scgp/data/10x/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds"

load(myinf1)

# Remove QC genes:
qc_genes <- c("ENSGGENES","ENSGUMI","ENSGMITO", "ENSGSEQSAT","ENSGSAMP") 
log2cpm <- log2cpm[-which(rownames(log2cpm) %in% qc_genes),]
featuredata <- featuredata[-which(rownames(featuredata) %in% qc_genes),]

# Annotate clusters using previous definitions
clust_annot = tsne.data %>%
	rownames_to_column('cell') %>%
	mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
                              `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
                              `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) %>%
	column_to_rownames('cell')

# clust_annot file and log2cpm files are in the same order so no need to match the order between them
cell_type <- clust_annot[,"cell_type"]

# Assign each cell a subtype
log2cpm_annot <- log2cpm
colnames(log2cpm_annot) <- cell_type

# Get sample names
sample_id <- sapply(strsplit(rownames(clust_annot), "-"), "[[", 3)
sample_id <- recode(sample_id, "0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006",
					"5" = "SM011", "6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018")
clust_annot[,"sample_id"] <- sample_id
names(sample_id) <- rownames(clust_annot)

# Convert ensembl ID to gene symbol
featuredata <- featuredata[rownames(log2cpm_annot),]

gene <- featuredata[,"Associated.Gene.Name"]
names(gene) <- rownames(featuredata)

rownames(log2cpm_annot) <- gene

idhwt_scgp <- c("SM006","SM012, SM017", "SM018", "SM011")
idhmut_scgp <- c("SM004", "SM001", "SM015", "UC917", "SM002", "SM008")

# Use scran to find markers between these subtypes

# Check percent of significant genes that are significant in single cell data
#g1 <- log2cpm_annot[mygenes, which(sample_id %in% idhwt_scgp & grepl("myeloid", colnames(log2cpm_annot)))]	# IDHwt
#g2 <- log2cpm_annot[mygenes, which(sample_id  %in% idhmut_scgp & grepl("myeloid", colnames(log2cpm_annot)))]	# IDHmut
g1 <- sample_id[which(sample_id %in% idhwt_scgp & grepl("myeloid", colnames(log2cpm_annot)))]
g2 <- sample_id[which(sample_id %in% idhmut_scgp & grepl("myeloid", colnames(log2cpm_annot)))]

tmp <- names(log2cpm_annot[,c(g1,g2)])

scgp_p.val <- apply(log2cpm_annot[mygenes,], 1, function(x) wilcox.test(x[g1], x[g2])$p.value)
p.val <- eff <- rep(0, length(mygenes))
for(i in 1:length(mygenes))
{
	cat("\r",i)
	p.val[i] <- wilcox.test(as.numeric(g1[i,]), as.numeric(g2[i,]))$p.value
	eff[i] <- mean(as.numeric(g1[i,])) - mean(as.numeric(g2[i,]))
}
scres <- data.frame(mygenes, p.val, eff)

# Plot a heatmap

genedata <- featuredata[rownames(log2cpm),]
gene_symbol <- genedata[,"Associated.Gene.Name"]
gene_cpm <- data.frame(gene_symbol, log2cpm)

clust_annot_join <- clust_annot
clust_annot_join[,"cell_id"] <- gsub("-","\\.",rownames(clust_annot_join))
mycells <- clust_annot_join %>%
		   filter(cell_type == "myeloid" & (sample_id == "SM006" | sample_id == "SM018")) %>%
		   dplyr::select(cell_id) %>%
		   .$cell_id
			
plot_cpm <- gene_cpm[,c("gene_symbol",mycells)] %>%
			filter(gene_symbol %in% mygenes) %>%
			pivot_longer(-gene_symbol, names_to = "cell_id", values_to = "cpm") %>%
			inner_join(clust_annot_join, by = "cell_id") %>%
			group_by(gene_symbol) %>% 
			mutate(z = scale(cpm))


# Order cells based on their expression of each gene
gene_means <- plot_cpm %>%
			  group_by(cell_id) %>%
			  summarise(avg = mean(cpm)) %>%
			  arrange(desc(avg))

plot_cpm$cell_id = factor(plot_cpm$cell_id, levels = gene_means$cell_id)

# 10x data so sparse expression. Include only genes expressed in 20% of cells
gene_include <- gene_cpm[,c("gene_symbol",mycells)] %>%
		  		filter(gene_symbol %in% mygenes) %>%
		  	 	remove_rownames() %>%
		   		column_to_rownames("gene_symbol")	   
include_sums <- apply(gene_include, 1, function(x)sum(x > 0))
include_crit <- names(include_sums[which(include_sums > 0.2 * ncol(gene_include))])
	
# Rescale the values for visualization
quant <- 0.95
plot_cpm[which(plot_cpm[,"z"] > quantile(plot_cpm$z,probs=quant,na.rm=TRUE)),"z"] <- quantile(plot_cpm$z,probs=quant,na.rm=TRUE)
plot_cpm[which(plot_cpm[,"z"] < (-1) * quantile(plot_cpm$z,probs=quant,na.rm=TRUE)),"z"] <- quantile(plot_cpm$z,probs=quant,na.rm=TRUE) * -1


# Order the genes (rows) based on their fold-change in the CIBERSORTx data
gene_order <- myeloid_diff %>%
			  filter(q.val < 0.1) %>%
			  arrange(desc(eff))

plot_cpm$gene_symbol = factor(plot_cpm$gene_symbol, levels = rev(rownames(gene_order)))

p1 <- ggplot(data = plot_cpm %>% filter(gene_symbol %in% include_crit), aes(x = cell_id, y = gene_symbol)) +
  geom_tile(aes(fill=z)) +
  scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") + 
  theme_void() +
  theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	legend.position = "none")
# 	legend.key.size = unit(0.25, "cm"))

p2 <- ggplot(data = plot_cpm %>% filter(gene_symbol == "OSM"), aes(x = cell_id, y = gene_symbol)) +
  geom_tile(aes(fill=sample_id)) +
  scale_fill_manual(values=c("#008A22","#8A0000")) + 
  theme_void() +
  theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	legend.position = "none")
# 	legend.key.size = unit(0.25, "cm"))

#Align figures for printing
gb1 <- ggplot_build(p1)
gb2 <- ggplot_build(p2)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)

g <- gtable:::rbind_gtable(gA, gB, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n1, "null")
g$heights[panels[2]] <- unit(n1/30, "null")

grid.newpage()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_tcga_mes_vs_class_heatmap.pdf",width=3.5,height=3)
grid.draw(g)
dev.off()


# Create a barplot with p-values for this cohort
plot_merge <- plot_cpm %>%
			 inner_join(scres, by = c("gene_symbol"="mygenes")) %>%
			 filter(cell_id == "AGAGAATAGTCAGCCC.1.10")

plot_pval <- myeloid_diff %>%
			 rownames_to_column("gene_symbol") %>%
			 inner_join(plot_merge, by = "gene_symbol")

plot_pval[,"gene_symbol"] <- factor(plot_pval[,"gene_symbol"], levels=levels(plot_cpm$gene_symbol))

plot_pval[,"agree"] <- factor(plot_pval[,"p.val.y"] < 0.05 & ((plot_pval[,"eff.x"] > 0 & plot_pval[,"eff.y"] > 0) | (plot_pval[,"eff.x"] < 0 & plot_pval[,"eff.y"] < 0) ))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_tcga_mes_vs_class_signif_bar.pdf",width=0.1,height=2.9)
p3 <- ggplot(data = plot_pval %>% filter(cell_id == "AGAGAATAGTCAGCCC.1.10" & gene_symbol %in% include_crit), aes(x = cell_id, y = gene_symbol)) +
  geom_tile(aes(fill=agree)) +
  scale_fill_manual(values=c("white","black")) + 
  theme_void() +
  theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	legend.position = "none")
p3
dev.off()

# Get legends
p1 <- ggplot(data = plot_cpm %>% filter(gene_symbol %in% include_crit), aes(x = cell_id, y = gene_symbol)) +
  geom_tile(aes(fill=z)) +
  scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") + 
  theme_void() +
  theme(axis.text.x = element_blank(),
  	axis.text.y = element_text(size=5, hjust = 1),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	legend.key.size = unit(0.25, "cm"))

p2 <- ggplot(data = plot_cpm %>% filter(gene_symbol == "RASGEF1B"), aes(x = cell_id, y = gene_symbol)) +
  geom_tile(aes(fill=sample_id)) +
  scale_fill_manual(values=c("#008A22","#8A0000")) + 
  theme_void() +
  theme(axis.text.x = element_blank(),
  	axis.text.y = element_text(size=5),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	legend.key.size = unit(0.25, "cm"))

#Align figures for printing
gb1 <- ggplot_build(p1)
gb2 <- ggplot_build(p2)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)

g <- gtable:::rbind_gtable(gA, gB, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n1, "null")
g$heights[panels[2]] <- unit(n1/30, "null")

grid.newpage()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_mes_vs_class_legends.pdf",width=3.5,height=3)
grid.draw(g)
dev.off()
