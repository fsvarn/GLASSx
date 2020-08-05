###################################################
# Build a macrophage expression heatmap from Xue macrophage profiles split up by IDHmut subtype
# Updated: 2020.07.17
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(GSVA)
library(ggdendro)
library(grid)

#######################################################
rm(list=ls())
myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/CIBERSORTxHiRes_GLASS_myeloid_Window48.txt"

geps <- read.delim(myinf1, row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]

myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/xue_macrophage_2014/xue_macrophage_modules.txt"
mac_module_df <- read.delim(myinf2,stringsAsFactor=FALSE,row.names=1)
colnames(mac_module_df) <- paste("module_",1:49,sep="")

mac_modules <- list()
for(i in 1:ncol(mac_module_df))
{
	mymod <- mac_module_df[,i]
	mymod <- mymod[which(mymod != "")]
	mac_modules[[i]] <- mymod
}
names(mac_modules) <- colnames(mac_module_df)

res <- gsva(data.matrix(geps), mac_modules, method="ssgsea",parallel.sz=1)

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
WITH long_pairs AS
(
	SELECT ps.rna_barcode_a AS aliquot_barcode, signature_name, 'Initial' AS timepoint,  cc.case_sex, cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
	JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_a
	UNION
	SELECT ps.rna_barcode_b AS aliquot_barcode, signature_name, 'Recurrent' AS timepoint, cc.case_sex, cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	JOIN clinical.cases cc ON cc.case_barcode = ps.case_barcode
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_b
	JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_b
)
SELECT * 
FROM long_pairs
ORDER BY aliquot_batch
"
dat <- dbGetQuery(con,q)

module <- rownames(res)
res <- data.frame(module, res,stringsAsFactors=FALSE)
colnames(res) <- gsub("\\.","-",colnames(res))
plot_res <- res %>%
 			pivot_longer(-module, names_to = "aliquot_barcode", values_to = "es") %>%
 			inner_join(dat, "aliquot_barcode") %>%
 			as.data.frame()
 			
write.table(plot_res, "/projects/verhaak-lab/GLASS-III/data/res/macrophage_modules/CIBERSORTxHiRes_GLASS_myeloid_xue_modules.txt",sep="\t",quote=FALSE,row.names=FALSE)

# Rescale each measurement (independently) to have a mean of 0 and variance of 1
plot_scaled <- plot_res
for(i in 1:length(module))
{
	plot_scaled[which(plot_res[,"module"]==module[i]),"es"] <- scale(plot_scaled[which(plot_res[,"module"]==module[i]),"es"])
}

plot_scaled <- plot_scaled %>%
			   mutate(idhmut = recode(idh_codel_subtype, 'IDHwt' = 'IDHwt', 'IDHmut-noncodel' = 'IDHmut', 'IDHmut-codel' = 'IDHmut'))
plot_scaled[,"module"]  <- gsub("module_","",plot_scaled[,"module"])

#Make a matrix for clustering pan-glioma modules
plot_mat <- matrix(0,nrow=length(unique(plot_scaled[,"module"])),ncol=length(unique(plot_scaled[,"aliquot_barcode"])))
rownames(plot_mat) <- unique(plot_scaled[,"module"])
colnames(plot_mat) <- unique(plot_scaled[,"aliquot_barcode"])
for(i in 1:nrow(plot_res))
{
	mysigname <- as.character(plot_scaled[i,"module"])
	myalicode <- as.character(plot_scaled[i,"aliquot_barcode"])
	plot_mat[mysigname,myalicode] <- plot_scaled[i,"es"]
}
plot_mat <- t(plot_mat)

# Run clustering
mod_dendro <- as.dendrogram(hclust(d = dist(x = t(plot_mat))))
mod.order <- order.dendrogram(mod_dendro)

module_levels <- unique(plot_scaled[,"module"])
module_levels <- module_levels[mod.order]
plot_scaled[,"module"] <- factor(plot_scaled[,"module"], levels = module_levels, ordered=TRUE)

original_scaled <- plot_scaled

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Make IDHwt plot
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

plot_scaled <- original_scaled[which(original_scaled[,"idhmut"]=="IDHwt"),]

#Make a matrix for clustering
plot_mat <- matrix(0,nrow=length(unique(plot_scaled[,"module"])),ncol=length(unique(plot_scaled[,"aliquot_barcode"])))
rownames(plot_mat) <- unique(plot_scaled[,"module"])
colnames(plot_mat) <- unique(plot_scaled[,"aliquot_barcode"])
for(i in 1:nrow(plot_scaled))
{
	mysigname <- as.character(plot_scaled[i,"module"])
	myalicode <- as.character(plot_scaled[i,"aliquot_barcode"])
	plot_mat[mysigname,myalicode] <- plot_scaled[i,"es"]
}
plot_mat <- t(plot_mat)

# Run clustering
samp_clusters <- hclust(d = dist(x = plot_mat))
dendro <- as.dendrogram(samp_clusters)
#Reverse to get high on the right and low on the left
dendro <- rev(dendro)

# Create dendrogram plot
dendro.plot <- ggdendrogram(data = dendro) + 
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
dendro.plot


# Heatmap
# Data wrangling
# Extract the order of the samples in the dendrogram
plot.order <- order.dendrogram(dendro)

# Split data into two clusters and extract them for later
idhwt_clust <- cutree(samp_clusters,2)

# Order the levels according to their position in the cluster
plot_scaled[,"aliquot_barcode"] <- factor(plot_scaled[,"aliquot_barcode"],
                               levels = rownames(plot_mat)[plot.order], 
                               ordered = TRUE)
                           
# Height for sidebar                          
plot_scaled[,"status"] <- 1

#Create sidebars with sample information

# Initial/Recurrence
gg_timepoint <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = timepoint)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("white", "black")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
  
# IDH/codel status
gg_idh <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = idh_codel_subtype)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("#619CFF")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
  
# Transcriptional subtype
gg_ts <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = signature_name)) +
	geom_bar(stat="identity") +
	theme_void() +
  	scale_fill_manual(values=c("#008A22", "#8A0000", "#00458A")) +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Batch
gg_batch <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = aliquot_batch)) +
	geom_bar(stat="identity") +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
        
# Create heatmap plot

heatmap.plot <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = module)) +
  geom_tile(aes(fill=es)) +
  scale_fill_gradient2(low ="royalblue4", mid="white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(gg_batch)
gb3 <- ggplot_build(gg_timepoint)
gb4 <- ggplot_build(gg_idh)
gb5 <- ggplot_build(gg_ts)
gb6 <- ggplot_build(heatmap.plot)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)
n4 <- length(gb4$layout$panel_scales_y[[1]]$range$range)
n5 <- length(gb5$layout$panel_scales_y[[1]]$range$range)
n6 <- length(gb6$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)
gF <- ggplot_gtable(gb6)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")
g <- gtable:::rbind_gtable(g, gF, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n6/4, "null")
g$heights[panels[2]] <- unit(n6/16, "null")
g$heights[panels[3]] <- unit(n6/16,"null")
g$heights[panels[4]] <- unit(n6/16,"null")
g$heights[panels[5]] <- unit(n6/16,"null")
g$heights[panels[6]] <- unit(n6,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/macrophage_heatmap_idhwt.pdf",width=7,height=7, useDingbats = FALSE)
grid.draw(g)
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Make IDHmut plot
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

plot_scaled <- original_scaled[which(original_scaled[,"idhmut"]=="IDHmut"),]

#Make a matrix for clustering
plot_mat <- matrix(0,nrow=length(unique(plot_scaled[,"module"])),ncol=length(unique(plot_scaled[,"aliquot_barcode"])))
rownames(plot_mat) <- unique(plot_scaled[,"module"])
colnames(plot_mat) <- unique(plot_scaled[,"aliquot_barcode"])
for(i in 1:nrow(plot_scaled))
{
	mysigname <- as.character(plot_scaled[i,"module"])
	myalicode <- as.character(plot_scaled[i,"aliquot_barcode"])
	plot_mat[mysigname,myalicode] <- plot_scaled[i,"es"]
}
plot_mat <- t(plot_mat)

# Run clustering
samp_clusters <- hclust(d = dist(x = plot_mat))
dendro <- as.dendrogram(samp_clusters)
#Reverse to get high on the right and low on the left
dendro <- rev(dendro)

# Create dendrogram plot
dendro.plot <- ggdendrogram(data = dendro) + 
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
dendro.plot


# Heatmap
# Data wrangling
# Extract the order of the samples in the dendrogram
plot.order <- order.dendrogram(dendro)

# Split data into two clusters and extract them for later
idhmut_clust <- cutree(samp_clusters,2)

# Heatmap
# Data wrangling
# Extract the order of the samples in the dendrogram
plot.order <- order.dendrogram(dendro)
# Order the levels according to their position in the cluster
plot_scaled[,"aliquot_barcode"] <- factor(plot_scaled[,"aliquot_barcode"],
                               levels = rownames(plot_mat)[plot.order], 
                               ordered = TRUE)
                               

# Height for sidebar                          
plot_scaled[,"status"] <- 1

#Create sidebars with sample information

# Initial/Recurrence
gg_timepoint <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = timepoint)) +
	geom_bar(stat="identity") +
	scale_fill_manual(values=c("white", "black")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
  
# IDH/codel status
gg_idh <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = idh_codel_subtype)) +
	geom_bar(stat="identity") +
	theme_void() +
	scale_fill_manual(values=c("#F8766D","#00BA38")) +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
  
# Transcriptional subtype
gg_ts <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = signature_name)) +
	geom_bar(stat="identity") +
	theme_void() +
  	scale_fill_manual(values=c("#008A22", "#8A0000", "#00458A")) +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

# Batch
gg_batch <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = status, fill = aliquot_batch)) +
	geom_bar(stat="identity") +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
        
# Create heatmap plot

heatmap.plot <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = module)) +
  geom_tile(aes(fill=es)) +
  scale_fill_gradient2(low ="royalblue4", mid="white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
heatmap.plot

label.plot <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = module)) +
  geom_tile(aes(fill=es)) +
  scale_fill_gradient2(low ="royalblue4", mid="white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  theme_void() +
  theme(axis.text.y = element_text(size=7),
  		axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")


#Scale factor for IDHmut (to make size comparable to IDHwt)
sf <- (sum(original_scaled[,"idh_codel_subtype"]=="IDHmut-noncodel" | original_scaled[,"idh_codel_subtype"]=="IDHmut-codel")/length(module))/(sum(original_scaled[,"idh_codel_subtype"]=="IDHwt")/length(module))

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(gg_batch)
gb3 <- ggplot_build(gg_timepoint)
gb4 <- ggplot_build(gg_idh)
gb5 <- ggplot_build(gg_ts)
gb6 <- ggplot_build(heatmap.plot)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)
n4 <- length(gb4$layout$panel_scales_y[[1]]$range$range)
n5 <- length(gb5$layout$panel_scales_y[[1]]$range$range)
n6 <- length(gb6$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)
gF <- ggplot_gtable(gb6)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")
g <- gtable:::rbind_gtable(g, gF, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n6/4, "null")
g$heights[panels[2]] <- unit(n6/16, "null")
g$heights[panels[3]] <- unit(n6/16,"null")
g$heights[panels[4]] <- unit(n6/16,"null")
g$heights[panels[5]] <- unit(n6/16,"null")
g$heights[panels[6]] <- unit(n6,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/macrophage_heatmap_idhmut.pdf",width=7*sf,height=7, useDingbats = FALSE)
grid.draw(g)
dev.off()

# Get x-axis labels

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(gg_batch)
gb3 <- ggplot_build(gg_timepoint)
gb4 <- ggplot_build(gg_idh)
gb5 <- ggplot_build(gg_ts)
gb6 <- ggplot_build(label.plot)

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)
n4 <- length(gb4$layout$panel_scales_y[[1]]$range$range)
n5 <- length(gb5$layout$panel_scales_y[[1]]$range$range)
n6 <- length(gb6$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)
gF <- ggplot_gtable(gb6)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")
g <- gtable:::rbind_gtable(g, gF, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n6/4, "null")
g$heights[panels[2]] <- unit(n6/16, "null")
g$heights[panels[3]] <- unit(n6/16,"null")
g$heights[panels[4]] <- unit(n6/16,"null")
g$heights[panels[5]] <- unit(n6/16,"null")
g$heights[panels[6]] <- unit(n6,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/macrophage_heatmap_labelled.pdf",width=7,height=7, useDingbats = FALSE)
grid.draw(g)
dev.off()

# Cluster analysis
idhwt_clust_df <- data.frame(names(idhwt_clust), idhwt_clust)
colnames(idhwt_clust_df) <- c("aliquot_barcode", "cluster")
idhwt_res <- dat %>%
	inner_join(idhwt_clust_df, by = "aliquot_barcode")

c1 <- sum(idhwt_res[,"signature_name"] == "Mesenchymal" & idhwt_res[,"cluster"] == 1)
c2 <- sum(idhwt_res[,"signature_name"] == "Mesenchymal" & idhwt_res[,"cluster"] == 2)
c3 <- sum(idhwt_res[,"signature_name"] != "Mesenchymal" & idhwt_res[,"cluster"] == 1)
c4 <- sum(idhwt_res[,"signature_name"] != "Mesenchymal" & idhwt_res[,"cluster"] == 2)
ct <- matrix(c(c1,c2,c3,c4),nrow=2,ncol=2)
fisher.test(ct)
# P = 0.06 for Mesenchymal
# P = 0.06 for Proneural
# P = 0.72 for Classical
# P = 0.04 for Mesenchymal vs Proneural

# Cluster analysis
idhmut_clust_df <- data.frame(names(idhmut_clust), idhmut_clust)
colnames(idhmut_clust_df) <- c("aliquot_barcode", "cluster")
idhmut_res <- dat %>%
	inner_join(idhmut_clust_df, by = "aliquot_barcode")

c1 <- sum(idhmut_res[,"signature_name"] == "Classical" & idhmut_res[,"cluster"] == 1)
c2 <- sum(idhmut_res[,"signature_name"] == "Classical" & idhmut_res[,"cluster"] == 2)
c3 <- sum(idhmut_res[,"signature_name"] != "Classical" & idhmut_res[,"cluster"] == 1)
c4 <- sum(idhmut_res[,"signature_name"] != "Classical" & idhmut_res[,"cluster"] == 2)
ct <- matrix(c(c1,c2,c3,c4),nrow=2,ncol=2)
fisher.test(ct)
# P = 0.01 for Mesenchymal
# P = 0.02 for Proneural
# P = 1 for Classical

# Full cluster analysis (the clusters appear to be identical in each heatmap)
all_res <- rbind(idhwt_res, idhmut_res)

c1 <- sum(all_res[,"signature_name"] == "Proneural" & all_res[,"cluster"] == 1)
c2 <- sum(all_res[,"signature_name"] == "Proneural" & all_res[,"cluster"] == 2)
c3 <- sum(all_res[,"signature_name"] != "Proneural" & all_res[,"cluster"] == 1)
c4 <- sum(all_res[,"signature_name"] != "Proneural" & all_res[,"cluster"] == 2)
ct <- matrix(c(c1,c2,c3,c4),nrow=2,ncol=2)
fisher.test(ct)
# P = 2e-5 for Mesenchymal
# P = 2e-9 for Proneural
# P = 0.11 for Classical