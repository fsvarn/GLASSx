###################################################
# Validate the Neftel synthetic mixture proportions with CIBERSORTx proportions
# Reference for dataset: An Integrative Model of Cellular States, Plasticity, and Genetics for Glioblastoma: Cell
# Updated: 2020.06.11
# Author: Frederick Varn
##################################################

library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(scales)

myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/SS2_neftel_synthetic_mix_cibersortx_prop_06092020.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/SS2_neftel_synthetic_mix_prop_06092020.txt"

prop <-read.delim(myinf1)
truth <- read.delim(myinf2)
truth <- t(truth)

comxx <- intersect(colnames(prop), colnames(truth))
prop <- prop[,comxx]
truth <- truth[,comxx]

# Get correlation coefficients
mycor <- rep(0, ncol(prop))
names(mycor) <- colnames(prop)
for(i in 1:ncol(prop))
{
	mycor[i] <- round(cor(prop[,i], truth[,i], method="p"),2)
}

colors <- brewer.pal(12,"Paired")

# Create scatterplots
se <- list()
for(i in 1:ncol(prop))
{
	cibersortx_prop <- prop[,i] * 100
	true_prop <- truth[,i] * 100
	plot_cor <- data.frame(cibersortx_prop, true_prop)
	
	mycoef <- format(mycor[i], nsmall = 2)
	names(mycoef) <- NULL
	myr <- deparse((bquote(italic("R") ~" = " ~ .(mycoef))))
	annotation_text <- data.frame(true_prop = 0.75 * max(true_prop), 
	cibersortx_prop = 0.02 * max(cibersortx_prop),
	myr)
	
	cell <- colnames(prop)[i]
	cell <- gsub("_"," ",cell)
	substr(cell,1,1) <- toupper(substr(cell,1,2))
	
	se[[i]] <- ggplot(plot_cor, aes(x=true_prop,y=cibersortx_prop)) +
	geom_point(colour="black",size=1) +
	geom_smooth(colour=colors[i],method='lm') +
	geom_text(data=annotation_text,label=myr, size=2.5, parse=TRUE) +
	labs(title = cell, x = "True fraction (%)", y = "CIBERSORTx fraction (%)") +
	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
	theme_bw() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.title=element_text(size=7,hjust=0.5),
	axis.title=element_text(size=7),
	axis.text.x=element_text(size=7),
	axis.text.y=element_text(size=7),
	legend.position="none") 
}

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/SS2_neftel_cibersortx_synth_mix_validation.pdf",width=4.9,height=3.5)
grid.arrange(se[[1]],se[[2]],se[[3]],se[[4]],se[[5]],se[[6]],nrow=2)
dev.off()