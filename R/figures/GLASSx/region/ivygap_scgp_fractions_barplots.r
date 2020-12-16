
library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

#######################################################

rm(list=ls())

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/ivygap_scgp/CIBERSORTxGEP_ivygap_Fractions-Adjusted.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/columns-samples.csv"


dat <- read.delim(myinf1,header=TRUE)
samps <- read.csv(myinf2,stringsAsFactor=FALSE)


#Simplify analysis for now: Reference histology
samp_struct <- samps$structure_name
names(samp_struct) <- samps$rna_well_id

tumor_id <- samps$tumor_id
names(tumor_id) <- samps$rna_well_id

dat[,"structure"] <- samp_struct[as.character(dat[,"Mixture"])]
dat[,"tumor_id"] <- tumor_id[as.character(dat[,"Mixture"])]
dat <- dat[grep("reference histology", dat[,"structure"]),]

dat[,"structure"] <- gsub(" sampled by reference histology","",dat[,"structure"])

dat <- dat %>% 
		select(-c("P.value", "Correlation", "RMSE")) %>%
		pivot_longer(cols=-c("Mixture","structure","tumor_id"), names_to = "cell_state", values_to = "fraction") %>%
		mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		mutate(structure = recode(structure, 
		"Cellular Tumor" = "Cellular tumor", "Infiltrating Tumor" = "Infiltrating tumor", "Leading Edge" = "Leading edge")) %>%
		as.data.frame()
		
res <- dat %>%
	   group_by(structure, cell_state) %>%
	   summarise(avg = mean(fraction)*100) %>%
	   mutate(cell_state = as_factor(cell_state)) %>%
	   mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
									"Oligodendrocyte", 
									"Endothelial", "Pericyte",
									"Fibroblast", 
									"Diff.-like", "Stem-like", "Prolif. stem-like"))

   
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_scgp_fractions.pdf", width=2, height = 3) #,width=2.7,height=3)
ggplot(res, aes(x=structure, y = avg, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
labs(y = "Proportion (%)") +
theme(axis.text.x = element_text(size=7,angle=45,hjust=1),
	axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_blank(),
	strip.background = element_blank(),
	legend.position = "none") +
    coord_cartesian(ylim=c(0,100))
dev.off()

# Group by sample

samp_res <- dat %>%
	   mutate(fraction = fraction*100) %>%
	   mutate(cell_state = as_factor(cell_state)) %>%
	   mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
									"Oligodendrocyte", 
									"Endothelial", "Pericyte",
									"Fibroblast", 
									"Diff.-like", "Stem-like", "Prolif. stem-like")) %>%
	   mutate(Mixture = as_factor(Mixture)) %>%
	   mutate(tumor_id = as_factor(tumor_id)) %>%
	   mutate(structure = as_factor(structure)) %>%
	   mutate(structure = fct_relevel(structure, "Leading edge", "Infiltrating tumor", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation")) %>%
	   arrange(tumor_id, structure)

mixture_order <- as.character(unique(samp_res$Mixture))
samp_res <- samp_res %>%
	   mutate(Mixture = fct_relevel(Mixture, mixture_order)) 	   
	   

p1 <- ggplot(samp_res, aes(x=Mixture, y = fraction, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
labs(y = "Proportion (%)") +
facet_grid(.~tumor_id, scales="free_x",space = "free_x") +
theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),
	axis.title.y= element_text(size=7),
	axis.ticks.x = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_blank(),
	legend.position = "none") +
    coord_cartesian(ylim=c(0,100))

samp_res[,"placeholder"] <- "100"
p2 <- ggplot(samp_res, aes(x=Mixture, y = placeholder, fill = structure)) +
geom_tile() +
theme_classic() +
scale_fill_manual(values=c("Leading edge" = "#009999", "Infiltrating tumor"= "#ad4597","Cellular tumor" = "#01b050", "Pseudopalisading cells around necrosis"="#02ffcc", "Microvascular proliferation"="#c00000")) +
labs(y = "%") +
facet_grid(.~tumor_id, scales="free_x",space = "free_x") +
theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),
	axis.title.y= element_text(size=7),
	axis.ticks.x = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_blank(),
	strip.background = element_blank(),
	legend.position = "none")

p2_legend <- ggplot(samp_res, aes(x=Mixture, y = placeholder, fill = structure)) +
geom_tile() +
theme_classic() +
scale_fill_manual(values=c("Leading edge" = "#009999", "Infiltrating tumor"= "#ad4597","Cellular tumor" = "#01b050", "Pseudopalisading cells around necrosis"="#02ffcc", "Microvascular proliferation"="#c00000")) +
labs(y = "%") +
facet_grid(.~tumor_id, scales="free_x",space = "free_x") 


########################
## View test plot
########################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  else
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

nolabels <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + plot_grid + null_x + null_y + null_facet + null_legend
  else
    gg + plot_theme + null_x + null_facet + null_legend
} 

gg_cbind <- function(..., widths = NULL) {
  if(length(match.call()) - 2 != length(widths))
    message("Number of widths does not match number of columns")
  gg <- gtable_cbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$widths[panels] <- unit(widths, "null")
  return(gg)
}


gg_rbind <- function(..., heights = NULL, ncol = 2) {
  if(length(match.call()) - 3 != length(heights))
    message("Number of heights does not match number of rows")
  gg <- gtable_rbind(...)
  panels = gg$layout$t[grep("panel", gg$layout$name)]
  gg$heights[panels] <- unit(rep(heights,each = ncol), "null")
  return(gg)
}

## Extract legend
gg_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#Align figures for printing
gb1 <- ggplot_build(p1)
gb2 <- ggplot_build(p2)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)

g <- gtable:::rbind_gtable(gA, gB, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1*2]] <- unit(n1, "null")
g$heights[panels[2*2]] <- unit(n1*5, "null")

grid.newpage()

#Plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_scgp_fractions_by_samp.pdf", width=5.5, height = 2) #,width=2.7,height=3)
grid.draw(g)
dev.off()
					
					#Legends
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_scgp_fractions_by_samp_legends.pdf",width=7,height=7)
p2_legend
dev.off()

