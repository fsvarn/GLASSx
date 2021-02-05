
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
annot = tsne.data %>%
	rownames_to_column('cell') %>%
	mutate(cell_type = recode(dbCluster, `1` = "differentiated_tumor",  `2` = "myeloid", `3` = "stemcell_tumor",
							  `4` = "oligodendrocyte", `5` = "prolif_stemcell_tumor", `6` = "granulocyte", `7` = "endothelial",
							  `8` = "t_cell", `9` = "pericyte", `10` = "fibroblast", `11` = "b_cell", `12` = "dendritic_cell")) %>%
	column_to_rownames('cell')
	

# Get sample names
sample_id <- sapply(strsplit(rownames(annot), "-"), "[[", 3)
sample_id <- recode(sample_id, "0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006",
					"5" = "SM011", "6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018")
annot[,"sample_id"] <- sample_id
names(sample_id) <- rownames(annot)

# Convert ensembl ID to gene symbol
featuredata <- featuredata[rownames(log2cpm),]

gene <- featuredata[,"Associated.Gene.Name"]
names(gene) <- rownames(featuredata)

rownames(log2cpm) <- gene

# Read in bulk subtype results
con <- DBI::dbConnect(odbc::odbc(), "scgp")

bulk_res <- dbReadTable(con, Id(schema = "analysis", table="transcriptional_subtype"))
bulk_mes <- bulk_res %>% 
			filter(signature_name == "Mesenchymal") %>%
			mutate(sample_id = sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[2],x[3],sep=""))) #%>%
			#filter(sample_id %in% c("SM006","SM012","SM017","SM018","SM011"))

# Test receptor-ligand associations
myeloid_genes <- c("COL14A1", "GAL", "ICAM4", "TNFSF14", "OSM", "CCL2", "IL7R", "IL7R")
tumor_genes <- c("CD44", "GALR2","ITGB1", "LTBR", "OSMR", "CCR10", "IL7", "TSLP")
tumor_class <- c("differentiated_tumor", "differentiated_tumor", "differentiated_tumor", "differentiated_tumor", "differentiated_tumor", "stemcell_tumor", "differentiated_tumor", "differentiated_tumor")
rec_lig <- data.frame(myeloid_genes, tumor_genes, tumor_class, stringsAsFactors=FALSE)

myeloid_cells <- annot %>% filter(cell_type == "myeloid")
myeloid_cor <- tumor_cor <- list()
for(i in 1:nrow(rec_lig))
{
	myeloid_cpm <- t(log2cpm[rec_lig[i,"myeloid_genes"],rownames(myeloid_cells)]) %>% data.frame()
	myeloid_cpm[,"sample_id"] <- myeloid_cells[,"sample_id"]
	colnames(myeloid_cpm) <- c("cpm","sample_id")
	
	myeloid_cor[[i]] <- myeloid_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(myeloid_cor = cor(mean_cpm, enrichment_score, method="p"),
			  myeloid_pval = cor.test(mean_cpm, enrichment_score, method="p")$p.value) 
	
	
	tumor_cells <- annot %>% filter(cell_type == rec_lig[i,"tumor_class"])

	tumor_cpm <- t(log2cpm[rec_lig[i,"tumor_genes"],rownames(tumor_cells)]) %>% data.frame()
	tumor_cpm[,"sample_id"] <- tumor_cells[,"sample_id"]
	colnames(tumor_cpm) <- c("cpm","sample_id")
	
	
	tumor_cor[[i]] <- tumor_cpm %>%
	group_by(sample_id) %>%
	summarise(mean_cpm = mean(cpm)) %>%
	inner_join(bulk_mes, c("sample_id" = "sample_id")) %>%	
	summarise(tumor_cor = cor(mean_cpm, enrichment_score, method="p"),
			  tumor_pval = cor.test(mean_cpm, enrichment_score, method="p")$p.value)
}

myeloid_cor <- do.call(rbind,myeloid_cor)
tumor_cor <- do.call(rbind,tumor_cor)
res <- data.frame(rec_lig, myeloid_cor, tumor_cor)

# Plot the result (OSM/OSMR, best performer)

myeloid_cpm <- t(log2cpm["OSM",rownames(myeloid_cells)]) %>% data.frame()
myeloid_cpm[,"sample_id"] <- myeloid_cells[,"sample_id"]
colnames(myeloid_cpm) <- c("cpm","sample_id")

plot_myel <- myeloid_cpm %>%
group_by(sample_id) %>%
summarise(mean_cpm = mean(cpm)) %>%
inner_join(bulk_mes, c("sample_id" = "sample_id"))

tumor_cells <- annot %>% filter(cell_type == "differentiated_tumor")
tumor_cpm <- t(log2cpm["OSMR",rownames(tumor_cells)]) %>% data.frame()
tumor_cpm[,"sample_id"] <- tumor_cells[,"sample_id"]
colnames(tumor_cpm) <- c("cpm","sample_id")

plot_tumor <- tumor_cpm %>%
group_by(sample_id) %>%
summarise(mean_cpm = mean(cpm)) %>%
inner_join(bulk_mes, c("sample_id" = "sample_id"))

plot_res <- rbind(plot_myel, plot_tumor) %>% 
			mutate(cell = rep(c("Myeloid", "Diff.-like tumor"), each  = nrow(plot_myel))) %>%
			mutate(cell = fct_relevel(cell, "Myeloid", "Diff.-like tumor"))

breaks_fun <- function(x) {
  if (max(x) > 0.4) {
    seq(0, 0.6, 0.2)
  } else {
    seq(0, 0.2, 0.1)
  }
}

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/scgp_osm_osmr_bulk_mesenchymal.pdf",width=1.25,height=2.5)
ggplot(plot_res, aes(mean_cpm, enrichment_score/1000)) + 
geom_point()  +
geom_smooth(method="lm", se = FALSE, fullrange=TRUE) +
facet_wrap(vars(cell),scales="free_x",nrow=2, labeller = labeller(cell = c("Myeloid (OSM)", "Diff.-like ") = dose.labs, supp = supp.labs)) +
theme_bw() +
labs(x="Mean expression", y = "Bulk tumor mes. score") +
scale_x_continuous(breaks = breaks_fun, limits = c(0, NA)) +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text = element_text(size=7),
axis.title = element_text(size=7),
strip.background = element_blank(),
strip.text = element_text(size=7),
legend.position="none") 
dev.off()

