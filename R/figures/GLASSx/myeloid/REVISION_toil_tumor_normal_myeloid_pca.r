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
# Step 2: PCA of the myeloid CIBERSORTx GEPs
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

pca <- prcomp(t(geps))$x

plot_res <- pca[,c("PC1", "PC2")]
plot_res <- plot_res %>%
	data.frame() %>%
	rownames_to_column(var = "aliquot_barcode") %>%
	inner_join(sig_sub, by = "aliquot_barcode") %>%
	mutate(IDH.codel.subtype = recode(IDH.codel.subtype, "Brain - Frontal Cortex (BA9)" = "Normal frontal cortex", "Brain - Cortex" = "Normal cortex")) %>%
	mutate(IDH.codel.subtype = fct_relevel(IDH.codel.subtype, "IDHmut-codel", "IDHmut-non-codel", "IDHwt", "Normal frontal cortex", "Normal cortex"))	%>%
	mutate(signature_name = recode(signature_name, "Brain - Frontal Cortex (BA9)" = "Normal frontal cortex", "Brain - Cortex" = "Normal cortex")) %>%
	mutate(signature_name = fct_relevel(signature_name, "Normal cortex", "Normal frontal cortex", "Proneural", "Proneural,Classical", "Mixed", "Classical", "Proneural,Mesenchymal","Classical,Mesenchymal","Mesenchymal"))



# PCA plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_pca_new.pdf",width=4,height=2.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(PC1, PC2, color = signature_name,shape = IDH.codel.subtype)) + 
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
scale_colour_manual(values=c("Normal cortex" = "#000000", "Normal frontal cortex" = "#A9A9A9", 
						   "Classical" = "#008A22", "Classical,Mesenchymal" = "#4B3F0F", "Mesenchymal" = "#8A0000", "Mixed" = "#D3D3D3", 
						   "Proneural" = "#00458A", "Proneural,Classical" = "#00645B", "Proneural,Mesenchymal" = "#3F264B"))
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_pca_new_legend.pdf",width=4,height=4)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(PC1, PC2, color = signature_name,shape = IDH.codel.subtype)) + 
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
scale_colour_manual(values=c("Normal cortex" = "#000000", "Normal frontal cortex" = "#A9A9A9", 
						   "Classical" = "#008A22", "Classical,Mesenchymal" = "#4B3F0F", "Mesenchymal" = "#8A0000", "Mixed" = "#D3D3D3", 
						   "Proneural" = "#00458A", "Proneural,Classical" = "#00645B", "Proneural,Mesenchymal" = "#3F264B"))
dev.off()

# PCA density plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_toil_myeloid_pca_dens.pdf",width=2.5,height=2.5)
ggplot(plot_res %>% filter(!grepl("-02A-|-02B-",plot_res[,"aliquot_barcode"])), aes(PC1, fill = signature_name)) + 
geom_density()  +
theme_bw() +
labs(title = "Myeloid") +
facet_grid(signature_name~., scales="free_y", switch = "y") +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title.x = element_text(size=7), 
axis.title.y = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_blank(),
strip.text.y.left = element_text(size=7,angle = 0, hjust = 1),
strip.background = element_blank(),
plot.title = element_text(size=7),
strip.placement = "outside",
legend.position="none")+
scale_fill_manual(values=c("Normal cortex" = "#000000", "Normal frontal cortex" = "#A9A9A9", 
						   "Classical" = "#008A22", "Classical,Mesenchymal" = "#4B3F0F", "Mesenchymal" = "#8A0000", "Mixed" = "#D3D3D3", 
						   "Proneural" = "#00458A", "Proneural,Classical" = "#00645B", "Proneural,Mesenchymal" = "#3F264B"))
dev.off()

##################################################
# Step 3: upSet plot to identify transcriptional overlap
##################################################

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
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_upset.pdf", width = 3.25, height = 2.5)  
upset(fromList(upset_list), keep.order = TRUE, text.scale = 0.875)
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

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idh_codel_upset.pdf", width = 3.25, height = 2.5)  
upset(fromList(upset_list), keep.order = TRUE, text.scale = 0.875)
dev.off()

##################################################
# Step 4: Ligand-receptor and prior signature analysis
##################################################

# Test the intersection of the microglia and MDM signatures:

# Establish connection to db and get macrophage and microglia score
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
sigs <- sigs %>%
		filter(signature_set == "Muller")

intersect(sigs[,"gene_symbol"], basal_sig)

# Test ligand-receptor pairs

lig_rec <- read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/ramilowski_2015/receptor_ligand_pairs_fantom5_ramilowski_ncomms_2015.txt",stringsAsFactor=FALSE)
lig_rec <- lig_rec %>%
filter(lig_rec[,"Pair.Evidence"] == "literature supported")
lig_rec <- lig_rec[,c("Ligand.ApprovedSymbol", "Receptor.ApprovedSymbol")]

pair_list <- c(lig_rec$Ligand.ApprovedSymbol, lig_rec$Receptor.ApprovedSymbol)
intersect(pair_list, cla_genes)
intersect(pair_list, pro_genes)
intersect(pair_list, mes_genes)

length(intersect(pair_list, cla_genes))/length(cla_genes)
length(intersect(pair_list, pro_genes))/length(pro_genes)
length(intersect(pair_list, mes_genes))/length(mes_genes)

##################################################
# Step 5: Create final tables
##################################################

# Overall final table
gene_list <- unique(c(cla_genes, pro_genes, mes_genes))
in_classical <- gene_list %in% cla_genes
in_proneural <- gene_list %in% pro_genes
in_mesenchymal <- gene_list %in% mes_genes

final_table <- data.frame(gene_list, in_proneural, in_classical, in_mesenchymal)
ordvar <- rep(NA, nrow(final_table))
ordvar[which(in_proneural & in_classical & in_mesenchymal)] <- "a"
ordvar[which(in_proneural & in_classical & !(in_mesenchymal))] <- "b"
ordvar[which(in_proneural & !(in_classical) & in_mesenchymal)] <- "c"
ordvar[which(!(in_proneural) & in_classical & in_mesenchymal)] <- "d"
ordvar[which(in_proneural & !(in_classical) & !(in_mesenchymal))] <- "e"
ordvar[which(!(in_proneural) & in_classical & !(in_mesenchymal))] <- "f"
ordvar[which(!(in_proneural) & !(in_classical) & in_mesenchymal)] <- "g"


final_table <- final_table[order(ordvar),]
final_table[,2] <- as.numeric(final_table[,2])
final_table[,3] <- as.numeric(final_table[,3])
final_table[,4] <- as.numeric(final_table[,4])
colnames(final_table) <- c("Gene symbol" ,"Proneural", "Classical", "Mesenchymal")

write.table(final_table, "/projects/verhaak-lab/GLASS-III/data/res/myeloid_overlap/idhwt_subtype_overlap_tab.txt", sep = "\t", quote=FALSE, row.names=FALSE)



# Ligand-receptor table
lig_coll <- lig_rec %>%
 group_by(Ligand.ApprovedSymbol) %>%
 summarize(Receptor = str_c(Receptor.ApprovedSymbol, collapse = ", "))
 
rec_coll <- lig_rec %>%
 group_by(Receptor.ApprovedSymbol) %>%
 summarize(Ligand = str_c(Ligand.ApprovedSymbol, collapse = ", "))
 
lig_vect <- lig_coll$Receptor
names(lig_vect) <- lig_coll$Ligand.ApprovedSymbol
rec_vect <- rec_coll$Ligand
names(rec_vect) <- rec_coll$Receptor.ApprovedSymbol

ligand_table <- final_table[which(final_table[,"Gene symbol"] %in% lig_rec[,"Ligand.ApprovedSymbol"]),]
ligand_table[,"Receptors"] <- lig_vect[as.character(ligand_table[,"Gene symbol"])]
receptor_table <- final_table[which(final_table[,"Gene symbol"] %in% lig_rec[,"Receptor.ApprovedSymbol"]),]
receptor_table[,"Ligands"] <- rec_vect[as.character(receptor_table[,"Gene symbol"])]

write.table(ligand_table, "/projects/verhaak-lab/GLASS-III/data/res/myeloid_overlap/idhwt_subtype_overlap_ligand_tab.txt", sep = "\t", quote=FALSE, row.names=FALSE)
write.table(receptor_table, "/projects/verhaak-lab/GLASS-III/data/res/myeloid_overlap/idhwt_subtype_overlap_receptor_tab.txt", sep = "\t", quote=FALSE, row.names=FALSE)


write.table(final_table, "/projects/verhaak-lab/GLASS-III/data/res/myeloid_overlap/idhwt_subtype_overlap_tab.txt", sep = "\t", quote=FALSE, row.names=FALSE)

