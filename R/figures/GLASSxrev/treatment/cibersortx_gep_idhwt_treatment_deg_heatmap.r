###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDH-wild-type tumors
# Author: Frederick Varn
# Date: 2022.01.08
# Figures 4A, 4B, S4E
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

q <- "SELECT ps.*, cc.case_age_diagnosis_years, cc.case_overall_survival_mo, cc.case_vital_status
FROM analysis.tumor_rna_clinical_comparison ps
JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
WHERE idh_codel_subtype = 'IDHwt' AND received_treatment --AND subtype_a = subtype_b"

dat <- dbGetQuery(con,q)
dat <- dat %>% mutate(case_vital_status = recode(case_vital_status, "alive" = 1, "dead" = 0)) %>%
			   mutate(case_vital_status = as.numeric(case_vital_status))


p <- se <- sigs <- list()
for(i in 1:length(myinf1))
{
	cat("\r", i)
	geps <- read.delim(myinf1[i], row.names=1)
	#geps <- log10(geps+1)		No need because we are taking log2 fold-change later
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
	idhwt_res <- data.frame(p.val, q.val, eff)
	idhwt_res <- idhwt_res[order(eff),]

	idhwt_res[,"logp"] <- -log10(idhwt_res[,"p.val"])
	idhwt_res[,"sig"] <- idhwt_res[,"q.val"] < 0.05
	
	idhwt_res <- idhwt_res[order(idhwt_res$p.val),]

	myoutf <- paste("data/res/CIBERSORTx/analysis/GLASS_idhwt_", cell_state[i],"_postreatment_result.txt",sep="")

	#write.table(idhwt_res, myoutf, sep="\t",quote=FALSE,row.names=TRUE) 
	
	sigs[[i]] <- idhwt_res

	# Plot heatmaps
	mygenes <- rownames(idhwt_res[which(idhwt_res[,"sig"]),])
	sig_matrix <- geps[mygenes,c(dat[,"tumor_barcode_a"], dat[,"tumor_barcode_b"])]
	
	sig_matrix <- t(apply(sig_matrix, 1, function(x) log10(x) - median(log10(x))))
	
	grp1 <- sig_matrix[,dat[,"tumor_barcode_a"]]
	grp1 <- apply(grp1, 1, mean)
	grp2 <- sig_matrix[,dat[,"tumor_barcode_b"]]
	grp2 <- apply(grp2, 1, mean)
	
	expr <- c(grp1, grp2)
	# Rescale for color
	#expr[which(expr > quantile(expr,.99))] <- quantile(expr,.99)
	expr[which(expr > 0.31)] <- 0.31 # Plug the number in from the whole dataset (obtained below) back in for scaling
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
	aesthetics = "fill", limits = c(-0.31, 0.31)) + 
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
quantile(plot_res$expr,0.99) #0.31

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/treatment_pre_post_csx_heatmaps_scale.pdf",width=7, height =3)
grid.arrange(se[[1]],se[[2]],se[[3]],nrow=1)
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/treatment_pre_post_csx_legends.pdf",width=7, height =15)
ggplot(data = plot_hm, aes(x = timepoint, y = gene_symbol)) +
	geom_tile(aes(fill=expr)) +
	scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3",
	space = "Lab", na.value = "grey50", guide = "colourbar",
	aesthetics = "fill", limits = c(-0.31, 0.31)) + 
	theme_void() +
	theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=2),
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
	up <- rownames(mysig %>% filter(q.val < 0.05, eff > 0))
	
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

# Test
stemRes <- sigRes %>%
		   filter(cell_state == "stemcell_tumor")
diffRes <- sigRes %>%
		   filter(cell_state == "differentiated_tumor")
prolifRes <- sigRes %>%
			filter(cell_state == "prolif_stemcell_tumor")

plotRes <- stemRes[1:15,]


plotRes[,"logP"] <- -log10(as.numeric(plotRes$q.value))

termord <- rev(as.character(plotRes[order(plotRes$logP, decreasing=TRUE),"Term"]))

plotRes <- plotRes %>% 
		   mutate(Term = fct_relevel(Term, termord))

# Plot the results
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_rec_sig_functional_enrichment_stem.pdf",width=3.3,height=2.2)
ggplot(data = plotRes, aes(x = logP, y = Term, fill = cell_state)) +
geom_bar(stat="identity",position=position_dodge()) +
geom_vline(xintercept=-log10(0.05),linetype=2) +
scale_fill_manual(values = c("stemcell_tumor" = "#fb6a4a", "differentiated_tumor" = "#fcbba1", "prolif_stemcell_tumor" = "#a50f15")) + 
labs(x="-log10(adj. P-value)") +
#facet_wrap(grp~.,scales="free_y", ncol = 1) +
theme_classic() +
theme(axis.text = element_text(size=7),
axis.title.x= element_text(size=7),
axis.title.y= element_blank(),
panel.grid.major=element_blank(),
panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_blank(),
legend.title = element_blank(),
legend.position = "none")# +
#coord_cartesian(ylim = c(0,6.5))
dev.off()


plotRes <- rbind(diffRes[1:5,], prolifRes[1:5,])

plotRes[,"logP"] <- -log10(as.numeric(plotRes$q.value))

termord <- rev(as.character(plotRes[order(plotRes$logP, decreasing=TRUE),"Term"]))

plotRes <- plotRes %>% 
		   mutate(Term = fct_relevel(Term, termord)) %>%
		   mutate(cell_state = recode(cell_state, "differentiated_tumor" = "Differentiated-like", "prolif_stemcell_tumor" = "Proliferating stem-like"))


# Plot the results
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/idhwt_rec_sig_functional_enrichment_nonstem.pdf",width=3,height=2.2)
ggplot(data = plotRes, aes(x = logP, y = Term, fill = cell_state)) +
geom_bar(stat="identity",position=position_dodge()) +
geom_vline(xintercept=-log10(0.05),linetype=2) +
scale_fill_manual(values = c("stemcell_tumor" = "#fb6a4a", "Differentiated-like" = "#fcbba1", "Proliferating stem-like" = "#a50f15")) + 
labs(x="-log10(adj. P-value)") +
facet_wrap(cell_state~.,scales="free_y", ncol = 1) +
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




