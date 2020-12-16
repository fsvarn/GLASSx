###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT ps.*, cs.idh_codel_subtype, ts1.signature_name AS subtype_a, ts2.signature_name AS subtype_b
FROM analysis.rna_silver_set ps
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ps.tumor_barcode_b
WHERE cs.idh_codel_subtype = 'IDHwt' AND ts1.signature_name !='Mesenchymal' AND ts2.signature_name = 'Mesenchymal'"

dat <- dbGetQuery(con,q)


# Check how TCGA signature changes:
tcga_sig <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result.txt")
tcga_sig <- rownames(tcga_sig %>% filter(q.val < 0.05 & eff < 0))

geps <- read.delim(myinf1[2], row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))

tcga_sig <- intersect(rownames(geps), tcga_sig)


g1 <- geps[tcga_sig,dat[,"tumor_barcode_a"]]
g2 <- geps[tcga_sig,dat[,"tumor_barcode_b"]]

g1 <- apply(g1,2,mean)
g2 <- apply(g2,2,mean)

wilcox.test(g1,g2,paired=TRUE)

# Boxplot of this result
aliquot_barcode <- c(names(g1), names(g2))
case_barcode <- sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[1:3],collapse="-"))
mes_score <- c(g1, g2)
timepoint <- c(rep("Initial",length(g1)), rep("Recurrent",length(g2)))
signature_name <- c(dat[,"subtype_a"], dat[,"subtype_b"])
mes_score <- data.frame(aliquot_barcode, case_barcode, mes_score, timepoint, signature_name)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/CIBERSORTx_glass_mesenchymal_class_score.pdf",width=1.5,height=2.5)
ggplot(mes_score, aes(x=timepoint, y=mes_score)) + 
geom_boxplot(outlier.size = NA)  +
geom_line(size=0.5,alpha=0.4, aes(group= case_barcode)) +
geom_point(size=1,aes(colour=signature_name)) +
theme_bw() +
labs(title = "Mesenchymal signature") +
theme(
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7),
plot.title = element_text(size=7, hjust=0.5),
legend.position="none") +
scale_colour_manual(values=c("#008A22", "#8A0000","#00458A"))
dev.off()

diff_list <- cor_list <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	#geps <- log10(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]
	colnames(geps) <- gsub("\\.","-",colnames(geps))


	g1 <- geps[,dat[,"tumor_barcode_a"]]
	g2 <- geps[,dat[,"tumor_barcode_b"]]

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(g1[j,])
		group2 <- as.numeric(g2[j,])
	
		p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
		eff[j] <- log2(median(group2/group1))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	idhwt_res <- data.frame(p.val, q.val, eff)
	idhwt_res <- idhwt_res[order(p.val),]

	idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
	idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.1
		
	mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])

	q <- 'SELECT es.*, gt.gene_symbol, gt.est_counts, gt.tpm
	FROM analysis.estimate es
	JOIN analysis.gene_tpm gt ON es.aliquot_barcode = gt.aliquot_barcode
	WHERE gene_symbol IN (mygenes)'

	mygenes_text <- paste("'", paste(mygenes, collapse="','"), "'", sep="")
	q <- gsub("mygenes", mygenes_text, q)

	pur_dat <- dbGetQuery(con,q)

	cor_res <- pur_dat %>%
	group_by(gene_symbol) %>%
	summarise(cor = cor(purity, tpm, method="s")) %>%
	inner_join(rownames_to_column(idhwt_res), by = c("gene_symbol" = "rowname")) %>%
	data.frame
	
	diff_list[[i]] <- idhwt_res
	cor_list[[i]] <- cor_res
}

names(diff_list) <- names(cor_list)  <- mytag

# Plot a volcano plot for myeloid cells

myeloid_diff <- diff_list[[2]]

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_myeloid_volcano_nonmes_mes.pdf",width=2.5,height=3)  
ggplot(myeloid_diff, aes(x=eff, y=logp)) + 
geom_point(size=0.5,aes(colour = sig)) +
labs(x = "log2(FC)", y = "-log10(adj p-val)") +
scale_colour_manual(values = c("black", "tomato3")) +
theme_bw() +
theme(
axis.title=element_text(size=7),
axis.text.x=element_text(size=7),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.position = "none") 
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Location of the SCGP single cell data
myinf1 <- "/projects/verhaak-lab/scgp/data/10x/10X-aggregation-20190905/RV18001-2-3-RV19001-2-4-5-6-7-8-9_20190827.Rds"

load(myinf1)

# Remove all cells with no expression
#sums <- apply(log2cpm, 2, sum)	
#range(sums)
# [1]  2453.769 45315.686
# No cells with 0 expression

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

g1 <- log2cpm_annot[,which(sample_id == "SM006" & grepl("myeloid", colnames(log2cpm_annot)))]	# AC-like
g2 <- log2cpm_annot[,which(sample_id == "SM018" & grepl("myeloid", colnames(log2cpm_annot)))]	# MES-like

mygenes <- rownames(myeloid_diff[which(myeloid_diff[,"q.val"] < 0.1),])
mygenes <- intersect(mygenes, rownames(log2cpm_annot))

p.val <- eff <- rep(0, length(mygenes))
for(i in 1:length(mygenes))
{
	p.val[i] <- wilcox.test(as.numeric(g1[mygenes[i],]), as.numeric(g2[mygenes[i],]))$p.value
	eff[i] <- mean(as.numeric(g2[mygenes[i],])) - mean(as.numeric(g1[mygenes[i],]))
}
scres <- data.frame(mygenes, p.val, eff)
scres[,"q.val"] <- p.adjust(scres[,"p.val"],"BH")

# Plot a heatmap

genedata <- featuredata[rownames(log2cpm),]
gene_symbol <- genedata[,"Associated.Gene.Name"]
gene_cpm <- data.frame(gene_symbol, log2cpm)

clust_annot_join <- clust_annot
clust_annot_join[,"cell_id"] <- gsub("-","\\.",rownames(clust_annot_join))
mycells <- clust_annot_join %>%
		   filter(cell_type == "myeloid" & (sample_id == "SM006" | sample_id == "SM018")) %>%
		   select(cell_id) %>%
		   .$cell_id
			
plot_cpm <- gene_cpm[,c("gene_symbol",mycells)] %>%
			filter(gene_symbol %in% mygenes) %>%
			pivot_longer(-gene_symbol, names_to = "cell_id", values_to = "cpm") %>%
			inner_join(clust_annot_join, by = "cell_id") %>%
			group_by(gene_symbol) %>% 
			mutate(z = scale(cpm))

# Rescale the values for visualization
quant <- 0.9
plot_cpm[which(plot_cpm[,"z"] > quantile(plot_cpm$z,probs=quant,na.rm=TRUE)),"z"] <- quantile(plot_cpm$z,probs=quant,na.rm=TRUE)
plot_cpm[which(plot_cpm[,"z"] < (-1) * quantile(plot_cpm$z,probs=quant,na.rm=TRUE)),"z"] <- quantile(plot_cpm$z,probs=quant,na.rm=TRUE) * -1


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
	
# Order the genes (rows) based on their fold-change in the CIBERSORTx data
gene_order <- diff_list[[2]] %>%
			  filter(q.val < 0.1) %>%
			  arrange(desc(eff))

plot_cpm$gene_symbol = factor(plot_cpm$gene_symbol, levels = rev(rownames(gene_order)))

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
	legend.position = "none")
# 	legend.key.size = unit(0.25, "cm"))

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

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/scgp_mes_vs_class_heatmap.pdf",width=3.5,height=3)
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

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/scgp_mes_vs_class_signif_bar.pdf",width=0.1,height=2.9)
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

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/scgp_mes_vs_class_legends.pdf",width=3.5,height=3)
grid.draw(g)
dev.off()