###################################################
# Validate the synthetic mixture proportions with CIBERSORTx proportions
# Updated: 2020.06.10
# Author: Frederick Varn
##################################################


library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(scales)

myinf1 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_synthetic_mix_cibersortx_prop_06092020.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/scgp/10x_scgp_synthetic_mix_prop_06092020.txt"

prop <-read.delim(myinf1)
truth <- read.delim(myinf2)
truth <- t(truth)

plot_data <- truth %>%
as.data.frame() %>%
rownames_to_column(var = "Mixture") %>%
pivot_longer(-Mixture, names_to = "cell_state", values_to = "true_proportion") %>%
inner_join(prop %>%
		   pivot_longer(-Mixture, names_to = "cell_state", values_to = "csx_proportion"),
		   by = c("Mixture","cell_state")) %>%
mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
					"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
					"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
					"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
					"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
					"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
 mutate(cell_state = as_factor(cell_state)) %>%
 mutate(cell_state = fct_relevel(cell_state, rev(c("B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
								"Oligodendrocyte", 
								"Endothelial", "Pericyte",
								"Fibroblast", 
								"Diff.-like", "Stem-like", "Prolif. stem-like"))))
		   
res <- plot_data %>%
		 group_by(cell_state) %>%
		 summarise(cor = cor(true_proportion, csx_proportion, method="p")) %>%
		 as.data.frame()

mycor <- round(as.numeric(res[,"cor"]),2)
names(mycor) <- as.character(res[,"cell_state"])

colors <- c("Stem-like" = "#fb6a4a", "Myeloid" = "#08519c", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15",
			"Oligodendrocyte" = "#2ca25f", "T cell" = "#6baed6", "Granulocyte" = "#bdd7e7", "Pericyte" = "#fee391", "Endothelial" = "#ffffd4", 
			"Dendritic cell" = "#3182bd", "Fibroblast" = "#feb24c", "B cell" = "#eff3ff")

cells <-  as.character(res[,"cell_state"])

# Create annotation text
myr <- myx <- myy <- rep(0,length(cells))
for(i in 1:length(cells))
{
	corval <- mycor[cells[i]]
	names(corval) <- NULL
	corval <- format(corval, nsmall = 2)
	myr[i] <- deparse((bquote(italic("R") ~" = " ~ .(corval))))

	sub_dat <- plot_data %>% filter(cell_state == cells[i])
	myx[i] <- 0.6 * max(sub_dat$true_proportion)
	myy[i] <- 0.05 * max(sub_dat$csx_proportion)
}

annotation_text <- data.frame(cell_state = cells,
							  true_proportion = myx,
							  csx_proportion = myy,
							  myr = myr)
							  
se <- ggplot(plot_data, aes(x=true_proportion*100,y=csx_proportion*100)) +
	geom_point(colour="black",size=1) +
	geom_smooth(method='lm',aes(colour = cell_state), se = FALSE) +
	facet_wrap(cell_state~.,scales="free", nrow=3) +
	scale_colour_manual(values = colors) + 
	geom_text(data = annotation_text, label=myr, size=2.5, parse=TRUE) +
	labs(x = "True fraction (%)", y = "CIBERSORTx fraction (%)") +
	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
	theme_classic() +
	theme(panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	plot.title=element_text(size=7,hjust=0.5),
	axis.title=element_text(size=7),
	axis.text.x=element_text(size=7),
	axis.text.y=element_text(size=7),
	strip.text = element_text(size=7),
	strip.background = element_blank(),
	legend.position="none") 

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cibersortx_synth_mix_validation.pdf",width=4,height=3.2)
se
dev.off()
# 
# Create scatterplots
# se <- list()
# for(i in 1:length(cells))
# {
# 	sub_dat <- plot_data %>% filter(cell_state == cells[i])
# 	
# 	mycoef <- format(mycor[cells[i]], nsmall = 2)
# 	names(mycoef) <- NULL
# 	myr <- deparse((bquote(italic("R") ~" = " ~ .(mycoef))))
# 	annotation_text <- data.frame(true_proportion = 0.75 * max(sub_dat$true_proportion), 
# 	csx_proportion = 0.02 * max(dat$csx_proportion),
# 	myr)
# 	
# 	se[[i]] <- ggplot(sub_dat, aes(x=true_proportion*100,y=csx_proportion*100)) +
# 	geom_point(colour="black",size=1) +
# 	geom_smooth(colour=colors[cells[i]],method='lm') +
# 	geom_text(data=annotation_text,label=myr, size=2.5, parse=TRUE) +
# 	labs(title = cells[i], x = "True fraction (%)", y = "CIBERSORTx fraction (%)") +
# 	scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) + 
# 	theme_bw() +
# 	theme(panel.grid.major = element_blank(),
# 	panel.grid.minor = element_blank(),
# 	plot.title=element_text(size=7,hjust=0.5),
# 	axis.title=element_blank(),
# 	axis.text.x=element_text(size=7),
# 	axis.text.y=element_text(size=7),
# 	legend.position="none") 
# }
# 
# pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cibersortx_synth_mix_validation.pdf",width=4,height=3.2)
# grid.arrange(se[[1]],se[[2]],se[[3]],se[[4]],se[[5]],se[[6]],se[[7]],se[[8]],se[[9]],se[[10]],se[[11]],se[[12]],nrow=3)
# dev.off()