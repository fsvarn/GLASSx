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
scale_colour_manual(values=c("#000000", "#A9A9A9", "#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
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
scale_colour_manual(values=c("#000000", "#A9A9A9", "#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B"))
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

myeloid_gep <- gep_list
case_barcode1 <- rownames(myeloid_gep)[grep("GTEX-",rownames(myeloid_gep))]
case_barcode1 <- sapply(strsplit(case_barcode1,"-"),function(x)paste(x[1:2],collapse="-"))
case_barcode2 <- rownames(myeloid_gep)[grep("TCGA-",rownames(myeloid_gep))]
case_barcode2 <- sapply(strsplit(case_barcode2,"-"),function(x)paste(x[1:3],collapse="-"))
case_barcode <- c(case_barcode1, case_barcode2)
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

colors <- c("#000000", "#A9A9A9", "#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B")
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

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_microglia_mdm_umap.pdf",width=1.25,height=2.75)
grid.arrange(p1,p2,nrow=2)
dev.off()


g1 <- sig_sub[which(sig_sub$IDH.status == "WT"),"aliquot_barcode"]
g2 <- sig_sub[which(sig_sub$IDH.status == "Brain - Cortex"),"aliquot_barcode"]

wilcox.test(mg_score[g1,2], mg_score[g2,2])

##################################################
# Step 4: Identify differentially-expressed genes
##################################################
jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

mynorm <- "Brain - Cortex"
mysubtype <- c("Classical","Mesenchymal","Proneural")
diff_list <- list()
for(i in 1:length(mysubtype))
{
	cat("\r", i)
	g1 <- geps[,sig_sub[which(sig_sub[,"signature_name"]==mysubtype[i] & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"]]
	g2 <- geps[,sig_sub[which(sig_sub[,"signature_name"]=="Brain - Cortex"),"aliquot_barcode"]]

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(g1[j,])
		group2 <- as.numeric(g2[j,])
	
		p.val[j] <- wilcox.test(group1,group2)$p.value
		eff[j] <- log2(median(group1)/median(group2))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	idhwt_res <- data.frame(p.val, q.val, eff)
	idhwt_res <- idhwt_res[order(p.val),]

	idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
	idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1
		
	mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
	
	diff_list[[i]] <- idhwt_res
}
names(diff_list) <- mysubtype

# Identify gene sets
cla_genes <- rownames(diff_list[["Classical"]] %>% filter(eff > log2(1.1), sig))
mes_genes <- rownames(diff_list[["Mesenchymal"]] %>% filter(eff > log2(1.1), sig))
pro_genes <- rownames(diff_list[["Proneural"]] %>% filter(eff > log2(1.1), sig))

upset_list <- list(cla_genes, mes_genes, pro_genes)
names(upset_list) <- mysubtype
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_upset.pdf")  
upset(fromList(upset_list), keep.order = TRUE)
dev.off()

basal_sig <- intersect(cla_genes, mes_genes)
basal_sig <- intersect(basal_sig, pro_genes)

clames <- intersect(mes_genes, cla_genes)


##################

mynorm <- "Brain - Cortex"
mysubtype <- c("IDHwt","IDHmut-non-codel","IDHmut-codel")
diff_list <- list()
for(i in 1:length(mysubtype))
{
	cat("\r", i)
	g1 <- geps[,sig_sub[which(sig_sub[,"IDH.codel.subtype"]==mysubtype[i] ),"aliquot_barcode"]]
	g2 <- geps[,sig_sub[which(sig_sub[,"IDH.codel.subtype"]=="Brain - Cortex"),"aliquot_barcode"]]

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(g1[j,])
		group2 <- as.numeric(g2[j,])
	
		p.val[j] <- wilcox.test(group1,group2)$p.value
		eff[j] <- log2(median(group1)/median(group2))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	idhwt_res <- data.frame(p.val, q.val, eff)
	idhwt_res <- idhwt_res[order(p.val),]

	idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
	idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1
		
	mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
	
	diff_list[[i]] <- idhwt_res
}
names(diff_list) <- mysubtype

# Identify gene sets
idhwt_genes <- rownames(diff_list[["IDHwt"]] %>% filter(eff > log2(1.1), sig))
noncodel_genes <- rownames(diff_list[["IDHmut-non-codel"]] %>% filter(eff > log2(1.1), sig))
codel_genes <- rownames(diff_list[["IDHmut-codel"]] %>% filter(eff > log2(1.1), sig))

upset_list <- list(idhwt_genes, noncodel_genes, codel_genes)
names(upset_list) <- mysubtype

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idh_codel_upset.pdf")  
upset(fromList(upset_list), keep.order = TRUE)
dev.off()


pan_sig <- intersect(idhwt_genes, noncodel_genes)
pan_sig <- intersect(pan_sig, codel_genes)
##################################################
# Step 5: Compare the gene sets with the genes identified in Klemm et al
##################################################

klemm_mdm_file <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/klemm_mdm_glioma_vs_in_vitro.txt"
klemm_mg_file <- "/projects/verhaak-lab/GLASS-III/data/dataset/klemm_tme_2020/klemm_mg_glioma_vs_nontumor.txt"

klemm_mdm <- read.delim(klemm_mdm_file)
klemm_mg <- read.delim(klemm_mg_file)

ensembl_map <- dbReadTable(con, Id(schema = "ref", table = "ensembl_gene_mapping"))
ensembl_key <- ensembl_map[,"gene_symbol"]
names(ensembl_key) <- ensembl_map[,"ensembl_gene_id"]

klemm_mdm[,"gene_symbol"] <- ensembl_key[klemm_mdm[,"ensembl_ID"]]
