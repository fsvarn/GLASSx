
library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(topGO)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("t_cell",myinf1)]
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

geps <- read.delim(myinf1, row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))


g1 <- as.numeric(geps["GZMA",dat[,"tumor_barcode_a"]])
g2 <- as.numeric(geps["GZMA",dat[,"tumor_barcode_b"]])

p.val <- apply(geps, 1, function(x) wilcox.test(x[dat[,"tumor_barcode_a"]],x[dat[,"tumor_barcode_b"]])$p.value)
eff <- apply(geps, 1, function(x) log2(median(x[dat[,"tumor_barcode_b"]])/median(x[dat[,"tumor_barcode_a"]])))
q.val <- p.adjust(p.val, "BH")
res <- data.frame(p.val, q.val, eff)
res %>% filter(p.val < 0.05, eff < 0)

mes_t <- rownames(res %>% filter(p.val < 0.05, eff > log2(1.1)))

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
	filter(dbCluster == 8) %>%
	column_to_rownames('cell')
	
log2cpm_tcell <- log2cpm[,rownames(clust_annot)]
	

# Get sample names
sample_id <- sapply(strsplit(rownames(clust_annot), "-"), "[[", 3)
sample_id <- recode(sample_id, "0" = "UC917", "1" = "SM001", "2" = "SM002", "3" = "SM004", "4" = "SM006",
					"5" = "SM011", "6" = "SM008", "7" = "SM012", "8" = "SM015", "9" = "SM017", "10" = "SM018")
clust_annot[,"sample_id"] <- sample_id
names(sample_id) <- rownames(clust_annot)

# Convert ensembl ID to gene symbol
featuredata <- featuredata[rownames(log2cpm_tcell),]

gene <- featuredata[,"Associated.Gene.Name"]
names(gene) <- rownames(featuredata)

rownames(log2cpm_tcell) <- gene

small_sig <- mes_t[which(mes_t %in% rownames(log2cpm_tcell))]
sc_mes_score <- apply(log2cpm_tcell[small_sig,],2,mean)
sample_id <- clust_annot[names(sc_mes_score),"sample_id"]
plot_res <- data.frame(names(sc_mes_score), sc_mes_score, sample_id)
colnames(plot_res) <- c("cell_id","sc_mes_score","sample_id")

# Plot density and arrange by patient
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/myeloid_scmes_score_density.pdf",width=1.5,height=1)  
ggplot(plot_res, aes(x = sc_mes_score)) +
geom_density() +
theme_classic() +
theme(axis.text.x = element_text(size=7), axis.text.y=element_text(size=7),
axis.title = element_blank(), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.background = element_blank(),
legend.position="none") #+
#coord_cartesian(xlim=c(-35, 35))
dev.off()

thr <- density(sc_mes_score)$x[which.max(density(sc_mes_score)$y)]

plot_res[,"mes_t"] <- as.numeric(sc_mes_score > thr)

plot_res %>% group_by(sample_id) %>% summarise(mean = mean(sc_mes_score))  %>% arrange(mean)
plot_res %>% group_by(sample_id) %>% summarise(prop = sum(mes_t)/n()) %>% arrange(prop)

# Read in bulk subtype results
con <- DBI::dbConnect(odbc::odbc(), "scgp")

bulk_res <- dbReadTable(con, Id(schema = "analysis", table="transcriptional_subtype"))
bulk_mes <- bulk_res %>% 
			filter(signature_name == "Mesenchymal") %>%
			mutate(sample_id = sapply(strsplit(aliquot_barcode,"-"),function(x)paste(x[2],x[3],sep=""))) #%>%
			#filter(sample_id %in% c("SM006","SM012","SM017","SM018","SM011"))

plot_res <- plot_res %>% 
group_by(sample_id) %>% 
summarise(mean = mean(sc_mes_score)) %>% 
ungroup() %>%
inner_join(bulk_mes, by = "sample_id") %>%
as.data.frame()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/myeloid_sc_bulk_mes_scatterplot.pdf",width=2,height=1.5)  
ggplot(plot_res, aes(x = mean, y = enrichment_score)) +
geom_point() +
geom_smooth(method = lm, se = FALSE) + 
theme_classic() +
labs(x = "Mean myeloid cell mes. score", y = "Bulk mes. score") +
theme(axis.text = element_text(size=7),
axis.title = element_text(size=7), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.background = element_blank(),
legend.position="none") #+
#coord_cartesian(xlim=c(-35, 35))
dev.off()


##################################################
# Step 4: Perform functional enrichment analysis on the signature
##################################################

all_genes <- rownames(geps)
bg_genes <- as.numeric(all_genes %in% mes_t)
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
elimRes[which(elimRes$raw.p.value < 0.05 & elimRes$Significant > elimRes$Expected),]
elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

term_level <- rev(as.character(elimRes$Term))
plotRes <- elimRes %>% 
		   mutate(Term = as_factor(Term)) %>%
		   mutate(Term = fct_relevel(Term, term_level)) %>%
		   filter(q.value < 0.05 & elimRes$Significant > elimRes$Expected)

