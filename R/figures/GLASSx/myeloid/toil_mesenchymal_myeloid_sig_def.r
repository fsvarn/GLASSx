library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(ssgsea.GBM.classification)
library(umap)
library(topGO)

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
# Step 2: Differential expression analysis of myeloid GEP
##################################################

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/toil/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_toil_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
colnames(geps) <- paste(mytag, colnames(geps), sep="__")
#geps <- log10(geps+1) No need because we are taking log2 fold-change later
gep_list <- t(geps)

rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]

colnames(geps) <- gsub(paste(mytag,"__",sep=""), "", colnames(geps))

g1 <- geps[,as.character(sig_sub[which(sig_sub[,"signature_name"]=="Mesenchymal" & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"])]
g2 <- geps[,as.character(sig_sub[which(!grepl("Mesenchymal",sig_sub[,"signature_name"]) & sig_sub[,"IDH.codel.subtype"]=="IDHwt"),"aliquot_barcode"])]

p.val <- eff <- rep(0, nrow(geps))
names(p.val) <- names(eff) <- rownames(geps)
for(j in 1:nrow(geps))
{
	group1 <- as.numeric(g1[j,])
	group2 <- as.numeric(g2[j,])

	p.val[j] <- wilcox.test(group1,group2)$p.value
	eff[j] <- log2(mean(group1)/mean(group2))
	#eff[i] <- median(group2) - median(group1)
}
q.val <- p.adjust(p.val,"BH")
idhwt_res <- data.frame(p.val, q.val, eff)
idhwt_res <- idhwt_res[order(p.val),]

idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1
	
mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])

diff_list <- idhwt_res

# Save the myeloid results to apply to GLASS
mes_myel_gene_sig <- diff_list
write.table(mes_myel_gene_sig, "data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v4.txt",sep="\t",quote=FALSE,row.names=TRUE) 
#v2 is with the directionality so that mes is the "change" and non-mes is the reference
#v3 is the signature developed from IDHwt only
#v4 is the signature developed from IDHwt in TOIL

##################################################
# Step 3: Define the myeloid signature and plot a heatmap
##################################################

# Make a TCGA heatmap of mesenchymal
mes_geps <- gep_list
mes_genes <- rownames(mes_myel_gene_sig[which(mes_myel_gene_sig[,"sig"] & mes_myel_gene_sig[,"eff"] > log2(1.5)),])

# Normalize the profile
mes_sig_prof <- log10(mes_geps[,mes_genes]+1)
#mes_sig_norm <- apply(mes_sig_prof, 2, function(x) ((x - mean(x)))/sd(x))

heatmap_res <- mes_sig_prof %>% 
			   as.data.frame() %>%
			   rownames_to_column(var="aliquot_barcode") %>%
			   mutate(aliquot_barcode=gsub("myeloid__","",aliquot_barcode)) %>%
			   pivot_longer(!aliquot_barcode,names_to= "gene" ,values_to= "z") %>%
			   inner_join(sig_sub, by = "aliquot_barcode") %>%
			   filter(IDH.codel.subtype == "IDHwt") %>%
			   mutate(enrichment_score = (enrichment_score - abs(min(enrichment_score)))/1000) %>%
			   group_by(gene) %>%
			   summarise(aliquot_barcode, z = (z-mean(z))/sd(z),signature_name, IDH.codel.subtype,enrichment_score)
heatmap_res[,"logp"] <- mes_myel_gene_sig[heatmap_res$gene, "logp"]
			   
level <- heatmap_res %>% filter(gene == "OSM") %>% arrange(desc(enrichment_score)) %>% .$aliquot_barcode
heatmap_res <- heatmap_res %>%
			   mutate(aliquot_barcode = as_factor(aliquot_barcode)) %>%
			   mutate(aliquot_barcode = fct_relevel(aliquot_barcode, level)) %>%
			   mutate(gene = as_factor(gene)) %>%
			   mutate(gene = fct_relevel(gene, rev(mes_genes)))

# Set limits for plotting
colmax <- quantile(heatmap_res$z,.99)
colmin <- -1 * colmax

heatmap_res <- heatmap_res %>%
mutate(z = ifelse(z > colmax, colmax, z)) %>%
mutate(z = ifelse(z < colmin, colmin, z))

# Barplot of mesenchymal score
p1 <- ggplot(data = heatmap_res %>% filter(gene == "TMEM150B"), aes(x = aliquot_barcode, y = enrichment_score)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=c(0, 12)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.title= element_blank(),
	axis.ticks.x = element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	plot.margin = unit(c(-1.2, -1.2, -1.5, 5.5), "pt"),
	legend.position = "none") 

# Heatmap
p2 <- ggplot(data = heatmap_res, aes(x = aliquot_barcode, y = gene)) +
  geom_tile(aes(fill=z)) +
  scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0,
  	  space = "Lab", na.value = "grey50", guide = "colourbar",
	  aesthetics = "fill") + 
  theme_void() +
  theme(axis.text = element_blank(),
	axis.title= element_blank(),
	axis.ticks = element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	plot.margin = unit(c(-1.2, -1.2, -1.5, 5.5), "pt"),
	legend.position = "none")
# 	legend.key.size = unit(0.25, "cm"))

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tcga_heatmap_legend_v2.pdf",width=0.2907,height=2)
ggplot(data = heatmap_res, aes(x = aliquot_barcode, y = gene)) +
geom_tile(aes(fill=z)) +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, 
space = "Lab", na.value = "grey50", guide = "colourbar",
aesthetics = "fill") + 
theme_void() +
theme(axis.text = element_blank(),
axis.title= element_blank(),
axis.ticks = element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank())
dev.off()

# Subtype classification tile
p3 <- ggplot(data = heatmap_res %>% filter(gene == "TMEM150B"), aes(x = aliquot_barcode, y = gene)) +
  geom_tile(aes(fill=signature_name)) +
  scale_fill_manual(values=c("#008A22", "#4B3F0F", "#8A0000", "#D3D3D3", "#00458A", "#00645B", "#3F264B")) +
  theme_void() +
  theme(axis.text = element_blank(),
	axis.title= element_blank(),
	axis.ticks = element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	plot.margin = unit(c(-1.2, -1.2, -1.5, 5.5), "pt"),
	legend.position = "none")
# 	legend.key.size = unit(0.25, "cm"))

# Barplot of log pval
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tcga_sig_logp_test_v2.pdf",width=0.2907,height=2)
ggplot(data = heatmap_res %>% filter(aliquot_barcode=="TCGA-02-0047-01A-01R-1849-01"), aes(x = logp, y = gene)) +
geom_bar(stat="identity") +
geom_vline(xintercept=1) + 
theme_bw() +
theme(axis.text.x = element_text(size=7),
axis.text.y = element_blank(),
axis.ticks.y = element_blank(),
axis.title = element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
plot.margin = unit(c(-1.2, -1.2, -1.5, -1.5), "pt"),
legend.position = "none") 
dev.off()

#Align figures for printing
gb1 <- ggplot_build(p1)
gb2 <- ggplot_build(p2)
gb3 <- ggplot_build(p3)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n2/6, "null")
g$heights[panels[2]] <- unit(n2, "null")
g$heights[panels[3]] <- unit(n2/15, "null")

grid.newpage()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tcga_sig_heatmap_v2.pdf",width=2,height=2.2)
grid.draw(g)
dev.off()


##################################################
# Step 4: Perform functional enrichment analysis on the signature
##################################################

all_genes <- rownames(diff_list)
bg_genes <- as.numeric(all_genes %in% mes_genes)
names(bg_genes) <- all_genes

diffGenes <- function(bg_genes) {
	return(bg_genes == 1)}

# Functional enrichment of signature
sampleGOdata <- new("topGOdata",
				description = "Simple session", 
				ontology = "BP",
			    allGenes = bg_genes,
			    geneSelectionFun = diffGenes, 
			    nodeSize = 10,
			    annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

# Fishers test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")

#Parentchild
resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
pcRes[which(pcRes$q.value < 0.05 & pcRes$Significant > pcRes$Expected),]

# Use this one: elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
elimRes[which(elimRes$q.value < 0.05 & elimRes$Significant > elimRes$Expected),]
elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

term_level <- rev(as.character(elimRes$Term))
plotRes <- elimRes %>% 
		   mutate(Term = as_factor(Term)) %>%
		   mutate(Term = fct_relevel(Term, term_level)) %>%
		   filter(q.value < 0.05 & elimRes$Significant > elimRes$Expected)

		   
# Plot the results
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tcga_mes_sig_functional_enrichment.pdf",width=4,height=2.5)
ggplot(data = plotRes, aes(x = logp, y = Term)) +
geom_bar(stat="identity") +
geom_vline(xintercept=-log10(0.05),linetype=2) +
labs(x="-log10(P)") +
theme_classic() +
theme(axis.text = element_text(size=7),
axis.title.x= element_text(size=7),
axis.title.y= element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
legend.position = "none") 
dev.off()