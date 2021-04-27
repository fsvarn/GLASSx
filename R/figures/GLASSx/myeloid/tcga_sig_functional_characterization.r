library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(topGO)

rm(list=ls())

# Check how TCGA signature changes:
tcga_dat <- read.delim("data/res/CIBERSORTx/analysis/TCGA_myeloid_ts_result_v3.txt")
tcga_sig <- rownames(tcga_dat %>% filter(q.val < 0.05 & eff > log2(1.1)))


##################################################
# Step 1: Perform functional enrichment analysis on the signature
##################################################

all_genes <- rownames(tcga_dat)
bg_genes <- as.numeric(all_genes %in% tcga_sig)
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

# Parentchild
# resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
# pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
# pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
# pcRes[which(pcRes$q.value < 0.05 & pcRes$Significant > pcRes$Expected),]
# 
# Use this one: elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
# resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
# elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
# elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
# elimRes[which(elimRes$q.value < 0.05 & elimRes$Significant > elimRes$Expected),]
# elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))


goRes <- fishRes

plotRes <- goRes %>%
		filter(q.value < 0.01 &  Significant/Expected > 1) %>%
		arrange(desc(Significant/Expected)) %>%
		mutate(Term = recode(Term, "positive regulation of acute inflammatory response" = "pos. reg. of acute inflammatory response"))

plotRes[,"logP"] <- -log10(as.numeric(plotRes$q.value))

plotRes <- plotRes[1:15,]
term_level <- plotRes$Term[order(as.numeric(plotRes$raw.p.value),decreasing=TRUE)]

plotRes <- plotRes %>% 
		   arrange(desc(Significant/Expected)) %>%
		   mutate(Term = as_factor(Term)) %>%
		   mutate(Term = fct_relevel(Term, term_level))
		   
		   
# Plot the results
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/tcga_mes_sig_functional_enrichment.pdf",width=2.7,height=1.8)
ggplot(data = plotRes, aes(x = logP, y = Term)) +
geom_bar(stat="identity") +
geom_vline(xintercept=-log10(0.05),linetype=2) +
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


##################################################
# Step 2: Look at overlap between known signatures (Xue macrophages and glioma-specific)
##################################################

module_inf <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"


mac_module_df <- read.delim(module_inf,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- gsub("X","module",colnames(mac_module_df))

xue_jac <- xue_p.val <- rep(NA, ncol(mac_module_df))
for(i in 1:ncol(mac_module_df))
{
	# Xue macrophage modules
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mymod <- intersect(mymod, all_genes)
	
	# Jaccard Index
	num <- length(intersect(tcga_sig, mymod))
	denom <- length(unique(c(tcga_sig, mymod)))
	xue_jac[i] <- num/denom
	
	#Fisher
	xue_p.val[i] <- fisher.test(matrix(c(length(all_genes) - length(union(tcga_sig, mymod)),
	 								 length(setdiff(tcga_sig, mymod)), 
	 								 length(setdiff(mymod, tcga_sig)), 
	 								 length(intersect(tcga_sig, mymod))), nrow=2),alternative="greater")$p.value
	 								 
}

modDir <- "/projects/verhaak-lab/varnf/R/pkgs/v3.6/ssgsea.GBM.classification/data"
modinf <- paste(modDir,"/",dir(modDir),sep="")

glio_jac <- glio_p.val <- rep(NA, length(modinf))
for(i in 1:length(modinf))
{
	glioma_subtype_mod <- readLines(modinf[i])[11]
	glioma_subtype_mod <- gsub("msig.up.genes\\tc\\(", "", glioma_subtype_mod)
	glioma_subtype_mod <- gsub("'", "", glioma_subtype_mod)
	glioma_subtype_mod <- gsub(")", "", glioma_subtype_mod)
	glioma_subtype_mod <- unlist(strsplit(glioma_subtype_mod, ","))

	# Jaccard Index
	num <- length(intersect(tcga_sig, glioma_subtype_mod))
	denom <- length(unique(c(tcga_sig, glioma_subtype_mod)))
	glio_jac[i] <- num/denom

	#Fisher
	glio_p.val[i] <- fisher.test(matrix(c(length(all_genes) - length(union(tcga_sig, glioma_subtype_mod)),
									 length(setdiff(tcga_sig, glioma_subtype_mod)), 
									 length(setdiff(glioma_subtype_mod, tcga_sig)), 
									 length(intersect(tcga_sig, glioma_subtype_mod))), nrow=2),alternative="greater")$p.value
}

# Macrophage/microglia result

# Establish connection to db and get macrophage signature
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")
sigs <- dbReadTable(con, Id(schema = "ref", table = "immune_signatures"))
mac_sig <- sigs %>%
		filter(signature_set == "Muller") %>%
		filter(signature_name == "Macrophages") %>%
		.$gene_symbol
mg_sig <- sigs %>%
		filter(signature_set == "Muller") %>%
		filter(signature_name == "Microglia") %>%
		.$gene_symbol
		
# Jaccard Index
num <- length(intersect(tcga_sig, mac_sig))
denom <- length(unique(c(tcga_sig, mac_sig)))
mac_jac <- num/denom

#Fisher
mac_p.val <- fisher.test(matrix(c(length(all_genes) - length(union(tcga_sig, mac_sig)),
								 length(setdiff(tcga_sig, mac_sig)), 
								 length(setdiff(mac_sig, tcga_sig)), 
								 length(intersect(tcga_sig, mac_sig))), nrow=2),alternative="greater")$p.value
								 

# Jaccard Index
num <- length(intersect(tcga_sig, mg_sig))
denom <- length(unique(c(tcga_sig, mg_sig)))
mg_jac <- num/denom

#Fisher
mg_p.val <- fisher.test(matrix(c(length(all_genes) - length(union(tcga_sig, mg_sig)),
								 length(setdiff(tcga_sig, mg_sig)), 
								 length(setdiff(mg_sig, tcga_sig)), 
								 length(intersect(tcga_sig, mg_sig))), nrow=2),alternative="greater")$p.value
								 
												 