###################################################
# Compare initial and recurrent genes in each cell type from initial and recurrent IDHwt tumors
# Updated: 2020.06.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(topGO)
library(biomaRt)

rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_mes_062320/CIBERSORTxGEP_Job27_GEPs_Filtered.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_mes_062320/CIBERSORTxGEP_Job27_GEPs_StdErrs.txt"
myinf3 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_nonmes_062320/CIBERSORTxGEP_Job28_GEPs_Filtered.txt"
myinf4 <- "/projects/verhaak-lab/GLASS-III/data/res/CIBERSORTx/CIBERSORTx_gep_group/CIBERSORTx_idhwt_nonmes_062320/CIBERSORTxGEP_Job28_GEPs_StdErrs.txt"

# Mesenchymal gene expression profiles
geps1 <- read.delim(myinf1, row.names=1)
stderr1 <- read.delim(myinf2, row.names=1)
stderr1 <- stderr1[rownames(geps1),]


# Non-mesenchymal gene expression profiles
geps2 <- read.delim(myinf3, row.names=1)
stderr2 <- read.delim(myinf4, row.names=1)
stderr2 <- stderr2[rownames(geps2),]

gene_qvals <- gene_diffs <- list()
for (i in 1:ncol(geps1))
{
	cat("\r",i)
	vBetaZ <- sapply(1:nrow(geps1), function(j) (geps1[j,i]-geps2[j,i])/sqrt(stderr1[j,i]^2+stderr2[j,i]^2))
	fc <- sapply(1:nrow(geps1), function(j) log2(geps1[j,i] +1 ) - log2(geps2[j,i] + 1))
	names(vBetaZ) <- names(fc) <- rownames(geps1)
	ZPs <- 2*pnorm(-abs(vBetaZ))
	Zqvals <- p.adjust(ZPs, method="BH")
	fc[which(Zqvals < 0.05 & !is.na(Zqvals))]
	
	gene_qvals[[i]] <- Zqvals[which(!is.na(Zqvals))]
	gene_diffs[[i]] <- fc[which(!is.na(Zqvals))]
	
}
names(gene_qvals) <- names(gene_diffs) <- colnames(geps1)

sapply(gene_diffs, length)
# 
# names(gene_diffs$myeloid[which(gene_qvals$myeloid < 0.05 & gene_diffs$myeloid > 0)])
# names(gene_diffs$myeloid[which(gene_qvals$myeloid < 0.05 & gene_diffs$myeloid < 0)])
# 
# cat(names(gene_diffs$t_cell[which(gene_qvals$t_cell < 0.05 & gene_diffs$t_cell > 0)]))
# cat(names(gene_diffs$t_cell[which(gene_qvals$t_cell < 0.05 & gene_diffs$t_cell < 0)]))
# 
# names(gene_diffs$dendritic_cell[which(gene_qvals$dendritic_cell < 0.05 & gene_diffs$dendritic_cell > 0)])
# names(gene_diffs$dendritic_cell[which(gene_qvals$dendritic_cell < 0.05 & gene_diffs$dendritic_cell < 0)])

cell <- "T cells"
gene_symbol <- names(gene_qvals$t_cell)
qval <- gene_qvals$t_cell
logfc <- gene_diffs$t_cell
logqval <- -log10(qval)
colour <- qval < 0.05
plot_res <- data.frame(gene_symbol, logfc, qval, logqval, colour)
plot_res <- plot_res[order(plot_res[,"logqval"],decreasing=TRUE),]

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_mes_nonmes_t_cell_volcano.pdf",width=2.5,height=3)
ggplot(plot_res, aes(x=logfc,y=logqval)) +
geom_vline(xintercept = 0, colour = "lightgray") +
geom_hline(yintercept = -log10(0.05), colour = "lightgray", linetype = "dashed") + 
geom_point(aes(colour=colour),size=1) +
scale_colour_manual(values=c("#27408B","#CD4F39")) +
labs(title = cell, x = "log2 fold change", y = "-log10(FDR)") +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.title=element_text(size=7,hjust=0.5),
axis.title=element_text(size=7),
axis.text.x=element_text(size=7),
axis.text.y=element_text(size=7),
legend.position="none") +
coord_cartesian(xlim = c(-18, 18))
dev.off()

# Make a GO enrichment plot

up <- names(gene_diffs$t_cell[which(gene_qvals$t_cell < 0.05 & gene_diffs$t_cell > 0)])
dn <- names(gene_diffs$t_cell[which(gene_qvals$t_cell < 0.05 & gene_diffs$t_cell < 0)])
bg <- names(gene_diffs$t_cell)

db <- useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids= getBM(attributes=c('go_id', 'external_gene_name', 'namespace_1003'), filters='external_gene_name', values=bg, mart=db)

# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])
 
# build the gene 2 GO annotation list (needed to create topGO object)
gene_2_GO=unstack(go_ids[,c(1,2)])
 
# remove any candidate genes without GO annotation
keep  <- dn %in% go_ids[,2]
keep  <- which(keep==TRUE)
candidate_list <- dn[keep]
 
# make named factor showing which genes are of interest
geneList=factor(as.integer(bg %in% candidate_list))
names(geneList)= bg

GOdata=new('topGOdata', ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene_2_GO)
	
# define test using the weight01 algorithm (default) with fisher
weight_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher') 

# generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
allGO=usedGO(GOdata)
all_res=GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))

#performing BH correction on our p values
p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)
 
# create the file with all the statistics from GO analysis
all_res_final=cbind(all_res,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]
 
#get list of significant GO before multiple testing correction
results.table.p= all_res_final[which(all_res_final$weightFisher<=0.001),]
 
#get list of significant GO after multiple testing correction
results.table.bh <- data.frame(all_res_final[which(all_res_final$p.adj<=0.05),])
umbrella_term <- results.table.bh[,"Term"]
umbrella_term[grep("metabol",umbrella_term)] <- "Metabolism"
umbrella_term[grep("stress",umbrella_term)] <- "Stress response"
umbrella_term[grep("oxy",umbrella_term)] <- "Hypoxia"
umbrella_term[grep("oxy",umbrella_term)] <- "Hypoxia"


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_mes_nonmes_t_cell_go_bar.pdf",width=5.5,height=3)
ggplot(results.table.bh, aes(x=logfdr,y=GO.biological.process.complete)) +
geom_bar(stat = "identity", fill="lightgray") +
geom_vline(xintercept = -log10(0.05), colour = "black") + 
labs(title = cell, x = "log2 fold change") +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.title=element_blank(),
axis.title.x=element_text(size=7),
axis.title.y = element_blank(),
axis.text.x=element_text(size=7),
axis.text.y=element_text(size=7),
legend.position="none")
dev.off()

pdf(file='/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_mes_nonmes_t_cell_topGOPlot_fullnames.pdf', height=12, width=12, paper='special', pointsize=18)
showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
dev.off()


# Make a GO enrichment plot from PANTHER
go_enr <- read.delim("/projects/verhaak-lab/GLASS-III/data/res/GO/cibersortx_mes_nonmes_tcell_GO.txt",stringsAsFactor=FALSE)
go_enr <-which(go_enr[,"GO.biological.process.complete"] == "regulation of transcription from RNA polymerase II promoter in response to stress ")
plot_go <- go_enr[,c("GO.biological.process.complete","Fold.Enrichment","FDR")]
logfdr <- -log10(plot_go[,"FDR"])
plot_go <- cbind(plot_go, logfdr)
plot_go[,"GO.biological.process.complete"] <- factor(plot_go[,"GO.biological.process.complete"],levels =rev(plot_go[,"GO.biological.process.complete"]))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/CIBERSORTx_mes_nonmes_t_cell_go_bar.pdf",width=5.5,height=3)
ggplot(plot_go, aes(x=logfdr,y=GO.biological.process.complete)) +
geom_bar(stat = "identity", fill="lightgray") +
geom_vline(xintercept = -log10(0.05), colour = "black") + 
labs(title = cell, x = "log2 fold change") +
theme_bw() +
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.title=element_blank(),
axis.title.x=element_text(size=7),
axis.title.y = element_blank(),
axis.text.x=element_text(size=7),
axis.text.y=element_text(size=7),
legend.position="none")
dev.off()
