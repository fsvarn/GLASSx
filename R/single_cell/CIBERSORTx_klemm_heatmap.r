###################################################
# Plot a heatmap for the CIBERSORTx results from the Klemm purified cell data
# Source: Interrogation of the Microenvironmental Landscape in Brain Tumors Reveals Disease-Specific Alterations of Immune Cells, Cell 2020
# Updated: 2020.06.12
# Author: Frederick Varn
##################################################

library(tidyverse)
library(RColorBrewer)

rm(list=ls())

myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_klemm_mix_cibersortx_prop_06112020.txt"

prop <-read.delim(myinf1,stringsAsFactors=FALSE)
prop <- prop[,1:13]
prop[,1] <- sapply(strsplit(prop[,1],"_"),function(x) paste(x[c(length(x),length(x)-1)],collapse="_"))
plot_res <- prop %>% pivot_longer(-Mixture, values_to = "fraction", names_to = "sample")
plot_res <- data.frame(plot_res)

ord <- rev(c(grep("cd45n_",plot_res[,"Mixture"]), grep("mg_",plot_res[,"Mixture"]), grep("mdm_",plot_res[,"Mixture"]),
	grep("neutrophils_",plot_res[,"Mixture"]), grep("cd4_",plot_res[,"Mixture"]), grep("cd8_",plot_res[,"Mixture"])))
plot_res[,"Mixture"] <- factor(plot_res[,"Mixture"],levels=unique(plot_res[ord,"Mixture"]))

##################################################
# Plot the results
##################################################

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/klemm_cibersort_heatmap.pdf",width=5,height=6)
ggplot(plot_res, aes(x=sample,y=Mixture,fill=fraction)) +
geom_tile() +
scale_fill_gradient2(low="royalblue4", mid="#ffffff", high="tomato3",midpoint=0.5) +
theme_bw() +
theme(axis.title=element_blank(),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=5),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "right")
dev.off()