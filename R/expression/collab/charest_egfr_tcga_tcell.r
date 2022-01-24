library(tidyverse)
library(gridExtra)

rm(list=ls())
data <- read.delim("/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/tcga/CIBERSORTxGEP_TCGA_Fractions-Adjusted.txt",stringsAsFactor=FALSE)
egfr_status <-  read.delim("/projects/verhaak-lab/GLASS-III/data/dataset/TCGA/tcga_egfr_mutations_brennan_s5.txt",stringsAsFactor=FALSE)

data <- data[grep("\\.01A",data[,1]),]

egfr_id <- egfr_status[,"Sample"]
egfr_id <- sapply(strsplit(egfr_id, "\\."),function(x)x[2])
egfr_id <- str_pad(egfr_id, 4, pad = "0")
egfr_status[,"egfr_id"] <- egfr_id

tcga_id <- sapply(strsplit(data[,"Mixture"], "\\."),function(x)x[3])
data[,"tcga_id"] <- tcga_id
data <- data[which(tcga_id %in% egfr_id),]

overlap <- intersect(egfr_status$egfr_id, data$tcga_id)
data <- data[which(data$tcga_id %in% overlap),]
egfr_status <- egfr_status[which(egfr_status$egfr_id %in% overlap),]


# Definition of EGFR vIII: delta 2-7 transcript allele frequency > 10% (as used in Brennan et al Cell 2013)
# Definition of EGFR mutation: non-delta 2-7 mutation/amplification transcript/variant allele frequency > 10% (as used in Brennan et al Cell 2013)
# Definition of wild-type: not part of EGFRvIII or EGFR mutation groups and not focally amplified (as used in Brennan et al Cell 2013)

v3 <- egfr_status[which(egfr_status$delta..2.7 >= 0.1),"egfr_id"]

mut <- egfr_status[which(apply(egfr_status[,3:(ncol(egfr_status)-1)], 1, function(x)sum(x>=0.1))>=1),"egfr_id"]
mut <- mut[which(!mut %in% v3)]

wt <- egfr_status[which(egfr_status[,"EGFR.CNA"] != "Focal Amplification"),"egfr_id"]
wt <- wt[which(!wt %in% c(v3,mut))]

g1 <- data[which(data$tcga_id %in% v3),"t_cell"]
g2 <- data[which(data$tcga_id %in% wt),"t_cell"]
wilcox.test(g1,g2)

g1 <- data[which(data$tcga_id %in% v3),"t_cell"]
g2 <- data[which(data$tcga_id %in% mut),"t_cell"]
wilcox.test(g1,g2)

g1 <- data[which(data$tcga_id %in% mut),"t_cell"]
g2 <- data[which(data$tcga_id %in% wt),"t_cell"]
wilcox.test(g1,g2)

g1 <- data[which(data$tcga_id %in% c(v3,mut)),"t_cell"]
g2 <- data[which(data$tcga_id %in% wt),"t_cell"]
wilcox.test(g1,g2)

# Make plotting table
data <- data[,-which(colnames(data) %in% c("Mixture","P.value","Correlation","RMSE"))]

plot_dat <- data %>%
			pivot_longer(cols = -tcga_id, names_to = "cell_state", values_to = "fraction") %>%
			mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
			"differentiated_tumor" = "Differentiated tumor", "endothelial" = "Endothelial",
			"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
			"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
			"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Proliferating stem cell tumor",
			"stemcell_tumor" = "Stem cell tumor","t_cell" = "T cell")) %>%
			mutate(cell_state = as_factor(cell_state)) %>%
			mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
											"Oligodendrocyte", 
											"Endothelial", "Pericyte",
											"Fibroblast", 
											"Differentiated tumor", "Stem cell tumor", "Proliferating stem cell tumor"))

status <- rep("", nrow(plot_dat))
status[which(plot_dat$tcga_id %in% v3)] <- "EGFRvIII"
status[which(plot_dat$tcga_id %in% mut)] <- "EGFRmut"
status[which(plot_dat$tcga_id %in% wt)] <- "WT"
plot_dat$status <- status

# Boxplot
boxes <- ggplot(data = plot_dat %>% 
	   		  filter(cell_state == "T cell", status != ""),
	   aes(x = status, y = fraction*100)) +
geom_boxplot(outlier.shape=NA,colour="black") +
geom_jitter(size=0.5) +
labs(y = "T cell fraction (%)") +
theme_classic() +
theme(axis.text.x = element_text(size=7,angle=45, hjust =1),
	axis.text.y = element_text(size=7),
	axis.title.x= element_blank(),
	axis.title.y= element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 

plot_bars <- plot_dat %>%
			 filter(status != "") %>%
			 group_by(cell_state, status) %>%
			 summarise(mean = mean(fraction))

bars <- ggplot(plot_bars, aes(x=status, y = mean*100, fill = factor(cell_state))) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
labs(y = "Fraction (%)") +
theme(axis.text.x = element_text(size=7,angle=45, hjust =1),
	axis.text.y = element_text(size=7),	axis.title.x = element_blank(),axis.title.y = element_text(size=7),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 
	
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/TCGA_EGFR_comparison_tcell.pdf",width=2.2,height=1.5)
grid.arrange(boxes, bars, nrow = 1, ncol = 2, widths = c(1.3, 1))
dev.off()
