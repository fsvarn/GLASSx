library(tidyverse)
library(ssgsea.GBM.classification)
library(umap)

#######################################################
rm(list=ls())

##################################################
# Step 1: Subtype each TCGA sample
##################################################

mes_path <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/group/tcga/MES/CIBERSORTxGEP_TCGA_MES_GEPs_Filtered.txt"
oth_path <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/group/tcga/OTH/CIBERSORTxGEP_TCGA_OTH_GEPs_Filtered.txt"

mese_path <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/group/tcga/MES/CIBERSORTxGEP_TCGA_MES_GEPs_StdErrs.txt"
othe_path <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/group/tcga/OTH/CIBERSORTxGEP_TCGA_OTH_GEPs_StdErrs.txt"

mes <- read.delim(mes_path,row.names=1)
oth <- read.delim(oth_path,row.names=1)

mese <- read.delim(mese_path,row.names=1)
othe <- read.delim(othe_path,row.names=1)

# Running DEG identifier per CIBERSORTx instructions
p.val <- q.val <- eff <- cell_state <- gene_symbol <- c()
for(i in 1:ncol(mes))
{
	cat("\r",i)
	vBetaZ <- sapply(1:nrow(mes), function(j) (mes[j,i] - oth[j,i])/sqrt(mese[j,i]^2 + othe[j,i]^2))
	zps <- 2*pnorm(-abs(vBetaZ))
	zqs <- p.adjust(zps, method="BH")
	
	fc <- sapply(1:nrow(mes), function(j) (log2((mes[j,i]+1)/(oth[j,i]+1))))
	eff <- c(eff, fc)
	p.val <- c(p.val, zps)
	q.val <- c(q.val, zqs)
	cell_state <- c(cell_state, rep(colnames(mes)[i],length(vBetaZ)))
	gene_symbol <- c(gene_symbol, rownames(mes))
}
logpval <- -log10(q.val)
sig <- q.val < 0.05
res <- data.frame(gene_symbol, eff, p.val, q.val, logpval, cell_state, sig)

# Plot a volcano plot for myeloid cells

plot_res <- res %>% 
			filter(cell_state == "myeloid") %>%
			filter(!is.na(logpval))

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_group_myeloid_volcano_nonmes_mes.pdf",width=2.5,height=3)  
ggplot(plot_res, aes(x=eff, y=logpval)) + 
geom_point(size=0.5,aes(colour = sig)) +
labs(x = "log2(FC))", y = "-log10(adj p-val)") +
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



plot_res <- res %>% 
			filter(cell_state == "differentiated_tumor") %>%
			filter(!is.na(logpval))
			
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_group_difftumor_volcano_nonmes_mes.pdf",width=2.5,height=3)  
ggplot(plot_res, aes(x=eff, y=logpval)) + 
geom_point(size=0.5,aes(colour = sig)) +
labs(x = "log2(FC))", y = "-log10(adj p-val)") +
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
