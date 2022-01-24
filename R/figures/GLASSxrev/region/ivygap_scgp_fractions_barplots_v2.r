###################################################
# Compare CIBERSORTx profiles across Ivy GAP histological features
# Author: Frederick Varn
# Date: 2021.12.30
# Figures 2A, S2B
##################################################

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

tumor_name <- samps$tumor_name
names(tumor_name) <- samps$rna_well_id

dat[,"structure"] <- samp_struct[as.character(dat[,"Mixture"])]
dat[,"tumor_name"] <- tumor_name[as.character(dat[,"Mixture"])]
dat <- dat[grep("reference histology", dat[,"structure"]),]

dat[,"structure"] <- gsub(" sampled by reference histology","",dat[,"structure"])


dat <- dat %>% 
		select(-c("P.value", "Correlation", "RMSE")) %>%
		pivot_longer(cols=-c("Mixture","structure","tumor_name"), names_to = "cell_state", values_to = "fraction") %>%
		mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		mutate(structure = recode(structure, 
		"Cellular Tumor" = "Cellular tumor", "Infiltrating Tumor" = "Infiltrating tumor", "Leading Edge" = "Leading edge")) %>%
		mutate(Mixture = as.character(Mixture), tumor_name = as.character(tumor_name)) %>%
		mutate(cell_state = fct_relevel(cell_state, rev(c("B cell", "Granulocyte","T cell", "Dendritic cell","Myeloid",
												"Oligodendrocyte", "Endothelial", "Pericyte", "Fibroblast",
												"Diff.-like","Stem-like","Prolif. stem-like"))))

# Two way ANOVA
aov_res <- dat %>%
group_by(cell_state) %>%
summarise(struct_p = summary(aov(fraction~structure + tumor_name))[[1]]["structure",5],
		  patient_p = summary(aov(fraction~structure + tumor_name))[[1]]["tumor_name",5])
		  
plot_aov <- aov_res %>%
			pivot_longer(-cell_state, names_to = "source", values_to = "p_value") %>%
			mutate(source = recode(source, "struct_p" = "structure", "patient_p" = "patient")) %>%
			mutate(log_p = -log10(p_value))
			
		  
# Plot log10 p-values in faceted barplot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_anova_bar.pdf", width=3.25, height = 1.45)
ggplot(plot_aov, aes(x=cell_state, y=log_p, fill=source)) + 
geom_bar(stat="identity", position = position_dodge()) +
geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
scale_fill_manual(values = c("structure" = "#27408B", "patient" = "#CD4F39")) +
labs(y = "-log10(p-value)") +
theme_classic() +
theme(axis.text.x = element_text(size=7, angle=45, hjust = 1),
axis.text.y = element_text(size=7),
axis.title.x = element_blank(),
axis.title.y = element_text(size=7), 
panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
strip.background = element_blank(),
legend.position="right",
legend.text = element_text(size=7),
legend.title=element_blank())
dev.off()

# Create average profile for plotting
		
dat <- dat %>%
	   group_by(structure, cell_state) %>%
	   summarise(fraction = mean(fraction)) %>%
	   mutate(tumor_name = "Mean", Mixture = paste("Mean",structure,sep=" ")) %>%
	   ungroup() %>%
	   select(Mixture, structure, tumor_name, cell_state, fraction) %>%
	   bind_rows(dat) %>%
	   mutate(Mixture = as_factor(Mixture)) %>%
	   mutate(structure = fct_relevel(structure, "Leading edge", "Infiltrating tumor", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation")) %>%
	   mutate(cell_state = as_factor(cell_state)) %>%
	   mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
									"Oligodendrocyte", 
									"Endothelial", "Pericyte",
									"Fibroblast", 
									"Diff.-like", "Stem-like", "Prolif. stem-like"))

   
# pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_scgp_fractions.pdf", width=2, height = 3) #,width=2.7,height=3)
# ggplot(res, aes(x=structure, y = avg, fill = factor(cell_state))) +
# geom_bar(position="stack", stat="identity") +
# theme_classic() +
# scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
# 					 "Oligodendrocyte" = "#2ca25f",
# 					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
# 					 "Fibroblast" = "#feb24c",
# 					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
# labs(y = "Proportion (%)") +
# theme(axis.text.x = element_text(size=7,angle=45,hjust=1),
# 	axis.text.y = element_text(size=7),
# 	axis.title.x = element_blank(),
# 	axis.title.y= element_text(size=7),
# 	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
# 	strip.text = element_blank(),
# 	strip.background = element_blank(),
# 	legend.position = "none") +
#     coord_cartesian(ylim=c(0,100))
# dev.off()

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
	   mutate(tumor_name = as_factor(tumor_name)) %>%
	   mutate(structure = as_factor(structure)) %>%
	   mutate(structure = fct_relevel(structure, "Leading edge", "Infiltrating tumor", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation")) %>%
	   arrange(tumor_name, structure)

mixture_mut <- as.character(unique(samp_res$Mixture[which(samp_res$tumor_name == "W31-1-1")]))
mixture_mean <- as.character(unique(samp_res$Mixture[which(samp_res$tumor_name == "Mean")]))
mixture_other <- as.character(unique(samp_res$Mixture[which(!(samp_res$Mixture %in% c(mixture_mut, mixture_mean)))]))
mixture_order <- c(mixture_other, mixture_mut, mixture_mean)
tumor_order <- as.character(unique(samp_res$tumor_name))
tumor_order <- tumor_order[c(2:7, 9:length(tumor_order),8,1)]


# Remove the mean from plotting and indicate IDH status
samp_res <- samp_res %>%
	   mutate(Mixture = fct_relevel(Mixture, mixture_order), tumor_name = fct_relevel(tumor_name, tumor_order)) %>%
	   filter(tumor_name != "Mean") %>%
	   mutate(idh_status = as.character(tumor_name == "W31-1-1"))	   

p1 <- ggplot(samp_res, aes(x=Mixture, y = fraction, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
labs(y = "Proportion (%)") +
facet_grid(.~structure, scales="free_x",space = "free_x") +
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
p2 <- ggplot(samp_res, aes(x=Mixture, y = placeholder, fill =idh_status)) +
geom_tile() +
theme_classic() +
labs(y = "%") +
scale_fill_manual(values = c("white","black")) +
facet_grid(.~structure, scales="free_x",space = "free_x") +
theme(axis.text.x = element_blank(),
	axis.text.y = element_text(size=7),
	axis.title.x = element_blank(),
	axis.title.y= element_text(size=7),
	axis.ticks.x = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_blank(),
	strip.background = element_blank(),
	legend.position = "none")
	
p3 <- ggplot(samp_res, aes(x=Mixture, y = placeholder, fill =tumor_name)) +
geom_tile() +
theme_classic() +
labs(y = "%") +
scale_fill_manual(values = brewer.pal(10,"Spectral")) +
facet_grid(.~structure, scales="free_x",space = "free_x") +
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
facet_grid(.~tumor_name, scales="free_x",space = "free_x") 


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
gb3 <- ggplot_build(p3)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels]
g$heights[panels[1*5]] <- unit(n1*10, "null")
g$heights[panels[2*5]] <- unit(n1, "null")
g$heights[panels[3*5]] <- unit(n1, "null")

grid.newpage()

#Plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_scgp_fractions_by_region_v2.pdf", width=6.6, height = 1.88) #,width=2.7,height=3)
grid.draw(g)
dev.off()
					
#Legends
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/ivygap_scgp_fractions_by_region_legends.pdf",width=7,height=7)
p2_legend
dev.off()

