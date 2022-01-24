library(tidyverse)
library(Seurat)
library(grid)
library(limma)
library(topGO)

#######################################################
rm(list=ls())

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Step 1: Identify the cells to compare
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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


		
# Get sample names
sample_id <- sapply(strsplit(rownames(clust_annot), "-"), "[[", 3)
sample_id <- recode(sample_id, "0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006",
					"5" = "SM011", "6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018")
clust_annot[,"sample_id"] <- sample_id
names(sample_id) <- rownames(clust_annot)


# Convert ensembl ID to gene symbol
featuredata <- featuredata[rownames(log2cpm),]

gene <- featuredata[,"Associated.Gene.Name"]
names(gene) <- rownames(featuredata)

rownames(log2cpm) <- gene

#pat_keep <- c("SM006", "SM012", "SM017", "SM018")
cell_keep <- clust_annot %>% filter(cell_type == "myeloid")# & sample_id %in% pat_keep)
sub_cpm <- log2cpm[,which(colnames(log2cpm) %in% rownames(cell_keep))]

sub_cpm <- t(sub_cpm[c("HMMR","ITGAX","CD68"),])

cd68 <- sub_cpm %>% data.frame() %>% filter(CD68 > 0, HMMR == 0, ITGAX == 0) %>% rownames()					#2146			7276
hmmr <- sub_cpm %>% data.frame() %>% filter(CD68 == 0, HMMR > 0, ITGAX == 0) %>% rownames()					#2				9
itgax <- sub_cpm %>% data.frame() %>% filter(CD68 == 0, HMMR == 0, ITGAX > 0) %>% rownames()				#286			847
cd68_hmmr <- sub_cpm %>% data.frame() %>% filter(CD68 > 0, HMMR > 0, ITGAX == 0) %>% rownames()				#82				114
cd68_itgax <- sub_cpm %>% data.frame() %>% filter(CD68 > 0, HMMR == 0, ITGAX > 0) %>% rownames()			#2542			5443
hmmr_itgax <- sub_cpm %>% data.frame() %>% filter(CD68 == 0, HMMR > 0, ITGAX > 0) %>% rownames()			#3				4
cd68_hmmr_itgax <- sub_cpm %>% data.frame() %>% filter(CD68 > 0, HMMR > 0, ITGAX > 0) %>% rownames()		#76				93
#5137/5534 (93%)
trip_neg <- sub_cpm %>% data.frame() %>% filter(CD68 == 0, HMMR == 0, ITGAX == 0) %>% rownames()			#397			1881

identity_vector <- rep(0, ncol(log2cpm))
names(identity_vector) <- colnames(log2cpm)
identity_vector[cd68] <- 1
identity_vector[hmmr] <- 2
identity_vector[itgax] <- 3
identity_vector[cd68_hmmr] <- 4
identity_vector[cd68_itgax] <- 5
identity_vector[hmmr_itgax] <- 6
identity_vector[cd68_hmmr_itgax] <- 7
identity_vector[trip_neg] <- 8



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Step 2: Read in and annotate the raw counts data
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Location of raw counts data
count_file <- "/projects/verhaak-lab/GLASS-III/data/scgp/analysis_scRNAseq_tumor_counts.h5ad"

counts <- ReadH5AD(count_file)

# Set identities for differential expression analysis
ident_count <- counts
Idents(ident_count) <- identity_vector

markers <- list()
for(i in 1:8)
{
	cat("\r", i)
	out_group <- 1:8
	out_group <- out_group[-which(out_group == i)] 
	
	if(sum(identity_vector==i) < 3){
		markers[[i]] <- NA
		next}
	
	tmp_markers <- FindMarkers(ident_count, ident.1=i, ident.2 = out_group)
	tmp_markers$cell_state <- i
	markers[[i]] <-  tmp_markers
	
}

final_markers <- do.call(rbind, markers)
final_markers <- final_markers %>%
				 mutate(cell_state = recode(cell_state, "1" = "CD68+/CD163-/CD11c-","2" = "CD68-/CD163+/CD11c-",
				 										"3" = "CD68-/CD163-/CD11c+","4" = "CD68+/CD163+/CD11c-",
				 										"5" = "CD68+/CD163-/CD11c+","6" = "CD68-/CD163+/CD11c+",
				 										"7" = "CD68+/CD163+/CD11c+","8" = "CD68-/CD163-/CD11c-"))
#write.table(final_markers, "/projects/verhaak-lab/GLASS-III/data/scgp/heimberger/myeloid_degs.txt",sep="\t", quote=FALSE, row.names=TRUE)
write.table(final_markers, "/projects/verhaak-lab/GLASS-III/data/scgp/heimberger/myeloid_degs_all_samps.txt",sep="\t", quote=FALSE, row.names=TRUE)
write.table(data.frame(rownames(counts)), "/projects/verhaak-lab/GLASS-III/data/scgp/heimberger/bg_genes.txt", quote=FALSE, header=FALSE, row.names=FALSE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Step 3: Perform the differential expression analysis
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Loading results in (picking up where Step 2 left off)
mydegs <- read.delim("/projects/verhaak-lab/GLASS-III/data/scgp/heimberger/myeloid_degs.txt")
mybg <- read.delim("/projects/verhaak-lab/GLASS-III/data/scgp/heimberger/bg_genes.txt",stringsAsFactors=FALSE)
mybg <- mybg[,1]

# GO enrichment
all_genes <- mybg

tumor_sigs <- mydegs
tumor_tag <- unique(mydegs$cell_state)

# Diff gene function for GO enrichment
diffGenes <- function(bg_genes) {
	return(bg_genes == 1)}

goList <- list()
for(i in 1:length(tumor_tag))
{
	cat("\r", i)
	mysig <- mydegs %>% 
			 filter(cell_state == tumor_tag[i]) %>%
			 filter(p_val_adj < 0.05 & avg_logFC > 0)
	
	up <- rownames(mysig) 
	
	bg_genes <- as.numeric(all_genes %in% up)
	names(bg_genes) <- all_genes



	# Functional enrichment of signature
	sampleGOdata <- new("topGOdata",
					description = "Simple session", 
					ontology = "BP",
					allGenes = bg_genes,
					geneSelectionFun = diffGenes, 		# Call function here
					nodeSize = 10,
					annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

	# Fishers test
# 	resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
# 	fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
# 	fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")

	# Other options
	#Parentchild
# 	resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
# 	pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
# 	pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
# 	pcRes[which(pcRes$q.value < 0.05 & pcRes$Significant > pcRes$Expected),]
# 
# 	#elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
	resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
	elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
	elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
	elimRes[which(elimRes$q.value < 0.05 & elimRes$Significant > elimRes$Expected),]
	elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

	fishRes[,"cell_state"] <- tumor_tag[i]
	goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)

sigRes <- goRes %>%
		filter(q.value < 0.05 & Significant > Expected)

mult <- table(sigRes[,"Term"])
myterms <- mult[which(mult>1)]

plotRes <- goRes %>% filter(Term %in% names(myterms))


plotRes[,"logP"] <- -log10(as.numeric(plotRes$q.value))

plotRes <- plotRes %>% 
		   filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor"))  %>%
		   mutate(Term = recode(Term, "phospholipase C-activating G protein-coupled receptor signaling pathway" = "PLC-activating GPCR signaling pathway")) %>%
		   mutate(Term = fct_relevel(Term, rev(unique(Term[order(logP, decreasing = TRUE)]))))

 
# Plot the results
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_rec_sig_functional_enrichment.pdf",width=3.3,height=2.2)
ggplot(data = plotRes, aes(x = logP, y = Term, fill = cell_state)) +
geom_bar(stat="identity",position=position_dodge()) +
geom_vline(xintercept=-log10(0.05),linetype=2) +
scale_fill_manual(values = c("stemcell_tumor" = "#fb6a4a", "differentiated_tumor" = "#fcbba1")) + 
labs(x="-log10(adj. P-value)") +
theme_classic() +
theme(axis.text = element_text(size=7),
axis.title.x= element_text(size=7),
axis.title.y= element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
legend.title = element_blank(),
legend.position = "none") 
dev.off()







# Select 3 key myeloid markers
sub_cpm <- log2cpm[c("HMMR","ITGAX","CD68"),]

# Keep only myeloid cells from  primary IDHwt samples
pat_keep <- c("SM006", "SM012", "SM017", "SM018")
cell_keep <- clust_annot %>% filter(cell_type == "myeloid" & sample_id %in% pat_keep)
sub_cpm <- sub_cpm[,which(colnames(sub_cpm) %in% rownames(cell_keep))]


cd68_dens <- density(as.numeric(sub_cpm["CD68",]))		# Local minimum for gating = 0.301770430
cd163_dens <- density(as.numeric(sub_cpm["HMMR",]))		# Local minimum for gating = 9.407207e-02
cd11c_dens <- density(as.numeric(sub_cpm["ITGAX",]))		# Local minimum for gating = 3.010956e-01

cd68_thr <- 0	#0.301770430
cd163_thr <- 0	#9.407207e-02
cd11c_thr <- 0	#3.010956e-01
sum((as.numeric(sub_cpm["CD68",]) > cd68_thr) & (as.numeric(sub_cpm["HMMR",]) > cd163_thr) & (as.numeric(sub_cpm["ITGAX",]) > cd11c_thr)) # 48

long_cpm <- sub_cpm %>% 
rownames_to_column(var = "gene_symbol") %>%
pivot_longer(-gene_symbol, names_to = "cell_id", values_to = "cpm") %>%
inner_join(clust_annot %>% rownames_to_column(var = "cell_id"), by = "cell_id") %>%
select(-c(V1,V2,V3))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/heimberger_marker_distribution.pdf")
ggplot(long_cpm, aes(x = cpm, fill=gene_symbol)) +
geom_density() +
scale_fill_manual(values=c("HMMR" = "#fb6a4a", "ITGAX" = "#fcbba1", "CD68" = "#a50f15")) +
facet_grid(sample_id ~ gene_symbol,scales="free_y",switch="y") +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_text(size=7,angle = 0, hjust = 1),
strip.text.y.left = element_text(size=7,angle = 0, hjust = 1),
strip.background = element_blank(),
strip.placement = "outside",
legend.position="none")  
dev.off()

wide_cpm <- long_cpm %>%
pivot_wider(c(cell_id, dbCluster, cell_type, sample_id), names_from = "gene_symbol", values_from = "cpm")

cd68_pos <- wide_cpm %>%
			filter(CD68 > cd68_thr)
			
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/heimberger_marker_distribution_cd68pos.pdf")
ggplot(cd68_pos, aes(x = ITGAX, y = HMMR)) +
geom_bin2d(bins=100) +
geom_hline(yintercept = cd163_thr) + 
geom_vline(xintercept = cd11c_thr) + 
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_blank(),
axis.ticks.y=element_blank(),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.text.x = element_text(size=7,angle = 0, hjust = 1),
strip.text.y.left = element_text(size=7,angle = 0, hjust = 1),
strip.background = element_blank(),
strip.placement = "outside",
legend.position="none")  
dev.off()