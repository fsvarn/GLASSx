##################################################
# Mixed effect models that link immune infiltrate to different driver alterations
# Also plots heatmap of these associations and a mutation co-occurrence matrix
# Updated: 2020.03.30
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(lme4)
library(car)
library(reshape)
library(grid)

##################################################
rm(list=ls())
# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
SELECT ds.*,im.signature_name, im.enrichment_score
FROM analysis.drivers_by_aliquot ds
JOIN analysis.analyte_sets an ON an.dna_barcode = ds.aliquot_barcode
JOIN analysis.davoli_immune_score im ON im.aliquot_barcode = an.rna_barcode
ORDER BY ds.aliquot_barcode, signature_name
"

dat <- dbGetQuery(con,q)

###################################################
# Step 1: Run mixed effect modeling for each mutation
##################################################

# Get list of unique driver alterations
drivers <- unique(unlist(strsplit(dat[,"snv_driver"],", ")))
drivers <- sapply(strsplit(drivers," "),function(x)x[1])
drivers <- unique(drivers)
drivers <- drivers[-which(is.na(drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
cells <- unique(dat[,"signature_name"])
genes <- drivers

# Create results matrix
snv_results <- matrix(NA, nrow = length(cells), ncol = length(genes))
rownames(snv_results) <- cells
colnames(snv_results) <- genes

snv_num <- rep(0,length(genes))
names(snv_num) <- genes

# nested for loop, one iteration for each cell, one for each mutation

for(i in 1:length(cells))
{
	# Subset data to only include that cell type's infiltration scores
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	for(j in 1:length(genes))
	{
		#Create a driver_status vector 
		driver_status <- rep(0, nrow(sub_dat))
		
		# Identify all samples that have a mutation in this gene
		indices <- grep(genes[j], sub_dat[,"snv_driver"])	
		
		# Change driver_status vector values to indicate driver is present at initial tumor
		driver_status[indices] <- 1
		if(sum(driver_status)==0){
			next}
		test_dat <- cbind(sub_dat, driver_status)
		
		lmm <- lmer(enrichment_score ~ driver_status + idh_codel_subtype + (1|case_barcode), data = test_dat)
		eff <- summary(lmm)[["coefficients"]][2]
		p.val <- Anova(lmm)[["Pr(>Chisq)"]][1]
		
		p.val <- ifelse(eff < 0, p.val*-1, p.val)
		snv_results[i,j] <- p.val
		snv_num[j] <- sum(driver_status)
	}
}

###################################################
# Step 2: Run mixed effect modeling for each copy number alteration
##################################################

# Get list of unique driver alterations
drivers <- unique(unlist(strsplit(dat[,"cnv_driver"],", ")))
drivers <- unique(drivers)
drivers <- drivers[-which(is.na(drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
cells <- unique(dat[,"signature_name"])
genes <- drivers

# Create results matrix
cnv_results <- matrix(NA, nrow = length(cells), ncol = length(genes))
rownames(cnv_results) <- cells
colnames(cnv_results) <- genes

cnv_num <- rep(0,length(genes))
names(cnv_num) <- genes

# nested for loop, one iteration for each cell, one for each mutation

for(i in 1:length(cells))
{
	# Subset data to only include that cell type's infiltration scores
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	for(j in 1:length(genes))
	{
		#Create a driver_status vector 
		driver_status <- rep(0, nrow(sub_dat))
		
		# Identify all samples that have a mutation in this gene
		indices <- grep(genes[j], sub_dat[,"cnv_driver"])	
		
		# Change driver_status vector values to indicate driver is present at initial tumor
		driver_status[indices] <- 1
		if(sum(driver_status)==0){
			next}
		test_dat <- cbind(sub_dat, driver_status)
		
		lmm <- lmer(enrichment_score ~ driver_status + idh_codel_subtype + (1|case_barcode), data = test_dat)
		eff <- summary(lmm)[["coefficients"]][2]
		p.val <- Anova(lmm)[["Pr(>Chisq)"]][1]
		
		p.val <- ifelse(eff < 0, p.val*-1, p.val)
		cnv_results[i,j] <- p.val
		cnv_num[j] <- sum(driver_status)
	}
}


###################################################
# Step 3: Plot the results and the mutual exclusivity of mutations
##################################################

#Figure 1: P-value heatmap

# Reshape data for ggplot2
snv_num <- snv_num[order(snv_num,decreasing=TRUE)]
cnv_num <- cnv_num[order(cnv_num,decreasing=TRUE)]
alt_order <- names(c(snv_num,cnv_num))

results <- cbind(snv_results,cnv_results)
plot_results <- melt(results)
colnames(plot_results) <- c("signature_name","alteration","p.value")
color <- rep("P > 0.05", nrow(plot_results)) 
color[which(plot_results[,"p.value"] > (-0.05) & plot_results[,"p.value"] < 0)] <- "Decrease (P < 0.05)"
color[which(plot_results[,"p.value"] < 0.05 & plot_results[,"p.value"] > 0)] <- "Increase (P < 0.05)"
plot_results <- cbind(plot_results,color)

# Filter out cell signatures that are too specific
plot_results <- plot_results[which(!plot_results[,"signature_name"] %in% c("Macrophages.M1","Macrophages.M2","CD8.effector.NK.cells")),] 

# Add column indicating whether the mutation is a mutation or CNA
alt_type <- rep("mut",nrow(plot_results))
alt_type[grep(" amp", plot_results[,"alteration"])] <- "cna"
alt_type[grep(" del", plot_results[,"alteration"])] <- "cna"
plot_results <- cbind(plot_results, alt_type)

# Set the orders
plot_results[,"alteration"] <- factor(plot_results[,"alteration"], levels = alt_order)
plot_results[,"alt_type"] <- factor(plot_results[,"alt_type"], levels = c("mut","cna"))

# Rename cells for plotting
plot_results <- plot_results %>% mutate(signature_name = recode_factor(signature_name,
															   "B.cells" = "B cells",
															   "CD4.mature" = "CD4+ T cells",
															   "CD8.effector" = "CD8+ T cells",
															   "Dendritic" = "Dendritic cells",
															   "Macrophages" = "Macrophages",
															   "NK.cells" = "NK cells",
															   "T.reg" = "T-regs"))

# Draw the plot

p1 <- ggplot(plot_results, aes(x=alteration,y=signature_name,fill=color)) +
geom_tile() +
facet_grid(.~alt_type, scales = "free_x") +
scale_fill_manual(values=c("royalblue4","tomato3", "#FFFFFF")) +
theme_bw() +
theme(axis.title=element_blank(),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_blank(),
strip.background = element_blank(),
legend.position = "none")


# convert ggplot object to grob object
gp <- ggplotGrob(p1)


# get gtable columns corresponding to the facets (5 & 7, in this case)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

# get the number of unique x-axis values per facet (10 & 12, in this case)
x.var <- sapply(ggplot_build(p1)$layout$panel_scales_x,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var

# plot result
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/mixed_effects_immune_mut_cna_heatmap.pdf",width=4,height=2)
grid.draw(gp)
dev.off()


#Figure 2: Mutation co-occurrence matrix ordered by CD8+ T cell infiltrate

uni_dat <- dat[which(dat[,"signature_name"]=="Macrophages"),]
alteration <- unique(plot_results[,"alteration"])

aliquot_barcode <- rep(uni_dat[,"aliquot_barcode"],each=length(alteration))
alteration <- rep(alteration, nrow(uni_dat))
present <- rep(0, length(aliquot_barcode))

# Create table for plotting that indicates whether a given alteration is in a sample
co_res <- data.frame(aliquot_barcode, alteration, present)
for(i in 1:nrow(co_res))
{
	myalt <- co_res[i,"alteration"]
	myali <- co_res[i,"aliquot_barcode"]
	
	sub_dat <-  uni_dat[which(uni_dat[,"aliquot_barcode"]==myali),]
	sub_alt <- paste(sub_dat[,"snv_driver"], sub_dat[,"cnv_driver"], sep = ", ")
	isthere <- grep(myalt, sub_alt)
	co_res[i,"present"] <- ifelse(length(isthere) > 0, 1, 0)
}

# Order the cases by macrophage infiltration level
mac_order <- uni_dat[order(uni_dat[,"enrichment_score"], decreasing=TRUE),"aliquot_barcode"]
co_res[,"aliquot_barcode"] <- factor(co_res[,"aliquot_barcode"],levels=mac_order)
co_res[,"present"] <- factor(co_res[,"present"],levels=c("0","1"))

all_alt_order <- c(snv_num, cnv_num)
all_alt_order <- all_alt_order[order(all_alt_order)]
co_res[,"alteration"] <- factor(co_res[,"alteration"],levels = names(all_alt_order))

# Draw the plot
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/mut_co_occurrence_macrophage.pdf",width=7,height=7)
ggplot(co_res, aes(x=aliquot_barcode,y=alteration,fill=present)) +
geom_tile() +
scale_fill_manual(values=c("white","black")) +
theme_bw() +
theme(axis.title=element_blank(),
axis.text.x=element_blank(),
axis.text.y=element_text(size=7),
strip.text = element_blank(),
strip.background = element_blank(),
legend.position = "none")
dev.off()
