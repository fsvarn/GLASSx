###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHmut tumors excluding hypermutants
# Note: This drops the number of samples from 21 to 12 and no genes are significant following this exclusion
# Updated: 2021.03.16
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(gridExtra)
library(survival)
library(topGO)

#######################################################
rm(list=ls())
myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")

cell_state <- gsub("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze//CIBERSORTxHiRes_GLASS_", "", myinf1)
cell_state <- gsub("_Window48.txt", "", cell_state)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT rc.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status, mf2.coverage_adj_mut_freq AS hm
FROM analysis.tumor_rna_clinical_comparison rc
LEFT JOIN analysis.platinum_set ps ON ps.rna_barcode_a = rc.tumor_barcode_a AND ps.rna_barcode_b = rc.tumor_barcode_b
LEFT JOIN analysis.mut_freq mf2 ON ps.dna_barcode_b = mf2.aliquot_barcode 
JOIN clinical.cases cc ON cc.case_barcode = rc.case_barcode
WHERE idh_codel_subtype LIKE 'IDHmut%' AND subtype_a = subtype_b AND received_treatment AND mf2.coverage_adj_mut_freq < 10"


dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))


p <- se <- sigs <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	geps <- log10(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]	
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	nrow(geps)

	g1 <- geps[,dat[,"tumor_barcode_a"]]
	g2 <- geps[,dat[,"tumor_barcode_b"]]

	p.val <- eff <- rep(0, nrow(geps))
	names(p.val) <- names(eff) <- rownames(geps)
	for(j in 1:nrow(geps))
	{
		group1 <- as.numeric(g1[j,])
		group2 <- as.numeric(g2[j,])
	
		p.val[j] <- wilcox.test(group1,group2,paired=TRUE)$p.value
		eff[j] <- log2(mean(group2)/mean(group1))
		#eff[i] <- median(group2) - median(group1)
	}
	q.val <- p.adjust(p.val,"BH")
	idhmut_res <- data.frame(p.val, q.val, eff)
	idhmut_res <- idhmut_res[order(eff),]

	idhmut_res[,"logp"] <- -log10(idhmut_res[,"p.val"])
	idhmut_res[,"sig"] <- idhmut_res[,"q.val"] < 0.1
	
	idhmut_res <- idhmut_res[order(idhmut_res$p.val),]

	myoutf <- paste("data/res/CIBERSORTx/analysis/GLASS_idhmut_", cell_state[i],"_postreatment_result.txt",sep="")

	#write.table(idhmut_res, myoutf, sep="\t",quote=FALSE,row.names=TRUE) 
	
	sigs[[i]] <- idhmut_res

	# Plot heatmaps
	mygenes <- rownames(idhmut_res[which(idhmut_res[,"sig"]),])
	sig_matrix <- geps[mygenes,c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"])]
	
	sig_matrix <- t(apply(sig_matrix, 1, function(x) x - median(x)))
	
	grp1 <- sig_matrix[,dat[,"tumor_barcode_a"]]
	grp1 <- apply(grp1, 1, mean)
	grp2 <- sig_matrix[,dat[,"tumor_barcode_b"]]
	grp2 <- apply(grp2, 1, mean)
	
	expr <- c(grp1, grp2)
	# Rescale for color
	#expr[which(expr > quantile(expr,.99))] <- quantile(expr,.99)
	expr[which(expr > 0.33)] <- 0.33 # Plug the number in from the whole dataset (obtained below) back in for scaling
	gene_symbol <- c(names(grp1), names(grp2))
	timepoint <- c(rep("Initial", length(grp1)), rep("Recurrent", length(grp2)))
	
	plot_hm <- data.frame(gene_symbol, expr, timepoint)
	lev <- names(grp1)[order(grp2 - grp1)]
	plot_hm$gene_symbol <- factor(plot_hm$gene_symbol, levels = lev)
	plot_hm[,"cell"] <- cell_state[i]
	
	p[[i]] <- plot_hm
	
	se[[i]] <- ggplot(data = plot_hm, aes(x = timepoint, y = gene_symbol)) +
	geom_tile(aes(fill=expr)) +
	scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
	space = "Lab", na.value = "grey50", guide = "colourbar",
	aesthetics = "fill", limits = c(-0.33, 0.33)) + 
	theme_void() +
	theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	legend.title = element_blank(),
	legend.position = "none")	
}

plot_res <- do.call(rbind, p)

# Identify the 99th percentile globally
#quantile(plot_res$expr,.99) #0.33

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_treatment_pre_post_csx_heatmaps_scale.pdf",width=7, height =3)
grid.arrange(se[[1]],se[[2]],se[[3]],nrow=1)
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_treatment_pre_post_csx_legends.pdf",width=7, height =3)
ggplot(data = plot_hm, aes(x = timepoint, y = gene_symbol)) +
	geom_tile(aes(fill=expr)) +
	scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
	space = "Lab", na.value = "grey50", guide = "colourbar",
	aesthetics = "fill", limits = c(-0.33, 0.33)) + 
	theme_void() +
	theme(axis.text = element_blank(),
	axis.title= element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank())	
dev.off()



# GO enrichment

tumor_sigs <- sigs
tumor_tag <- mytag

# Diff gene function for GO enrichment
diffGenes <- function(bg_genes) {
	return(bg_genes == 1)}

goList <- list()
for(i in 1:length(tumor_sigs))
{
	cat("\r", i)
	mysig <- tumor_sigs[[i]]
	up <- rownames(mysig %>% filter(sig, eff > 0))
	
	all_genes <- rownames(mysig)
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
	resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
	fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
	fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")

	# Other options
	#Parentchild
# 	resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
# 	pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
# 	pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
# 	pcRes[which(pcRes$q.value < 0.05 & pcRes$Significant > pcRes$Expected),]
# 
# 	#elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
# 	resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
# 	elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
# 	elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
# 	elimRes[which(elimRes$q.value < 0.05 & elimRes$Significant > elimRes$Expected),]
# 	elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

	fishRes[,"cell_state"] <- tumor_tag[i]
	goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)

# Take the top 8 significant in both
sigRes <- goRes %>%
		filter(q.value < 0.05 & Significant > Expected)

mult <- table(sigRes[,"Term"])
myterms <- mult[which(mult>1)]

upRes <- goRes %>% filter(Term %in% names(myterms))


upRes[,"logP"] <- -log10(as.numeric(upRes$q.value))


# Get the top 8
keep <- upRes %>%
		filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor")) %>%
		group_by(Term) %>%
		summarise(p_rtg = mean(as.numeric(raw.p.value))) %>%
		as.data.frame() %>%
		.[1:8,"Term"]
		
upRes <- upRes %>% 
		   filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor"))  %>%
		   filter(Term %in% keep)

upRes$direction <- "up"

# Check down-regulated
goList <- list()
for(i in 1:length(tumor_sigs))
{
	cat("\r", i)
	mysig <- tumor_sigs[[i]]
	up <- rownames(mysig %>% filter(sig, eff < 0))
	
	all_genes <- rownames(mysig)
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
	resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
	fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
	fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")

	# Other options
	#Parentchild
# 	resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
# 	pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
# 	pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
# 	pcRes[which(pcRes$q.value < 0.05 & pcRes$Significant > pcRes$Expected),]
# 
# 	#elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
# 	resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
# 	elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
# 	elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
# 	elimRes[which(elimRes$q.value < 0.05 & elimRes$Significant > elimRes$Expected),]
# 	elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

	fishRes[,"cell_state"] <- tumor_tag[i]
	goList[[i]] <- fishRes
}

goRes <- do.call(rbind, goList)

sigRes <- goRes %>%
		filter(q.value < 0.05 & Significant > Expected)

mult <- table(sigRes[,"Term"])
myterms <- mult[which(mult>1)]

dnRes <- goRes %>% filter(Term %in% names(myterms))


dnRes[,"logP"] <- -log10(as.numeric(dnRes$q.value))

# Get the top 8 by average p-value
keep <- dnRes %>%
		filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor")) %>%
		group_by(Term) %>%
		summarise(p_rtg = mean(as.numeric(raw.p.value))) %>%
		as.data.frame() %>%
		.[1:8,"Term"]

dnRes <- dnRes %>% 
		   filter(cell_state %in% c("differentiated_tumor", "stemcell_tumor"))  %>%
		   filter(Term %in% keep)


dnRes$direction <- "dn"
plotRes <- rbind(upRes, dnRes)
 
term_order <- plotRes %>% filter(cell_state == "stemcell_tumor")  %>% dplyr::select(Term,logP) %>% .$Term %>% as.character()
plotRes <- plotRes %>% 
		   mutate(Term = fct_relevel(Term, rev(term_order))) %>%
		   mutate(direction = recode(direction, "up" = "Up","dn" = "Down")) %>%
		   mutate(direction = fct_relevel(direction, "Up", "Down"))

# Plot the results
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhmut_rec_sig_functional_enrichment.pdf",width=2.9867,height=2.2)
ggplot(data = plotRes, aes(x = logP, y = Term, fill = cell_state)) +
geom_bar(stat="identity",position=position_dodge()) +
geom_vline(xintercept=-log10(0.05),linetype=2) +
scale_fill_manual(values = c("stemcell_tumor" = "#fb6a4a", "differentiated_tumor" = "#fcbba1")) + 
labs(x="-log10(adj. P-value)") +
facet_wrap(direction~.,scales="free_y", ncol = 1) +
theme_classic() +
theme(axis.text = element_text(size=7),
axis.title.x= element_text(size=7),
axis.title.y= element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_blank(),
legend.title = element_blank(),
legend.position = "none") 
dev.off()
