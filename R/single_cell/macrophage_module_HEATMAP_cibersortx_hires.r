###################################################
# Build a macrophage expression heatmap from Xue macrophage profiles
# Updated: 2020.07.16
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
	SELECT ps.rna_barcode_a AS aliquot_barcode, signature_name, 'Initial' AS timepoint,  cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ps.rna_barcode_a
	JOIN biospecimen.aliquots al ON al.aliquot_barcode = ps.rna_barcode_a
	UNION
	SELECT ps.rna_barcode_b AS aliquot_barcode, signature_name, 'Recurrent' AS timepoint, cs.idh_codel_subtype, al.aliquot_batch
	FROM analysis.platinum_set ps
	JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
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

# Rescale each measurement (independently) to have a mean of 0 and variance of 1
plot_scaled <- plot_res
for(i in 1:length(module))
{
	plot_scaled[which(plot_res[,"module"]==module[i]),"es"] <- scale(plot_scaled[which(plot_res[,"module"]==module[i]),"es"])
}

#Make a matrix for clustering
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
dendro <- as.dendrogram(hclust(d = dist(x = plot_mat)))
#Reverse to get high on the right and low on the left
dendro <- rev(dendro)

######################## 
## Common plotting elements
########################

plot_theme    <- theme_bw(base_size = 10) + theme(axis.title = element_text(size = 10),
                                                  axis.text = element_text(size=10),
                                                  panel.grid.major = element_blank(),
                                                  panel.grid.minor = element_blank(),
                                                  panel.background = element_rect(fill = "transparent"),
                                                  axis.line = element_blank())
null_legend   <- theme(legend.position = 'none')
null_x        <- theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
null_y        <- theme(axis.title.y=element_blank(),
                       axis.text.y=element_text(size=7),
                       axis.ticks.y=element_blank())
bottom_x      <- theme(axis.text.x=element_blank())
null_facet    <- theme(strip.background = element_blank(),
                       strip.text.x = element_blank())
top_margin    <- theme(plot.margin= unit(c(1, 1, 0.1, 1), "lines")) ## Top, Right, Bottom, Left
middle_margin <- theme(plot.margin= unit(c(0, 1, 0.1, 1), "lines"))
bottom_margin <- theme(plot.margin= unit(c(0, 1, 1, 1), "lines"))

########################
## View test plot
########################

testPlot <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  else
    gg + plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

nolabels <- function(gg, grid = TRUE) {
  if(grid)
    gg + plot_theme + null_x + null_y + null_facet + null_legend
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

########################
## Blank plot
########################

gg_blank <-
  ggplot(data.frame()) +
  geom_blank()
  
# Dendrogram
dendro.plot <- ggdendrogram(data = dendro) + 
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 

# Heatmap
# Data wrangling
# Extract the order of the samples in the dendrogram
plot.order <- order.dendrogram(dendro)
# Order the levels according to their position in the cluster
plot_scaled[,"aliquot_barcode"] <- factor(plot_scaled[,"aliquot_barcode"],
                               levels = rownames(plot_mat)[plot.order], 
                               ordered = TRUE)
plot_scaled[,"status"] <- 1
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

# Initial/Recurrence
gg_timepoint <-
  ggplot(plot_scaled, aes(x=aliquot_barcode, fill = timepoint)) +
  geom_tile(aes(y=status)) +
  scale_fill_manual(values=c("white","black")) +
  theme(axis.title.y=element_blank())
  
# IDH/codel status
gg_idh <-
  ggplot(plot_scaled, aes(x=aliquot_barcode, fill = idh_codel_subtype)) +
  geom_tile(aes(y=status)) +
  theme(axis.title.y=element_blank())
  
# Transcriptional subtype
gg_ts <-
  ggplot(plot_scaled, aes(x=aliquot_barcode, fill = signature_name)) +
  geom_tile(aes(y=status)) +
  scale_fill_manual(values=c("#008A22", "#8A0000", "#00458A"))
  theme(axis.title.y=element_blank())
  
# Batch
gg_batch <-
  ggplot(plot_scaled, aes(x=aliquot_barcode, fill = aliquot_batch)) +
  geom_tile(aes(y=status)) +
  theme(axis.title.y=element_blank())

#Reorder alignment plot by ordering case_barcode by where there aliquots first appear in dendrogram

#Align figures for printing
gb1 <- ggplot_build(nolabels(dendro.plot))
gb2 <- ggplot_build(nolabels(gg_batch))
gb3 <- ggplot_build(nolabels(gg_timepoint))
gb4 <- ggplot_build(nolabels(gg_idh))
gb5 <- ggplot_build(nolabels(gg_ts))
gb6 <- ggplot_build(nolabels(heatmap.plot))

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


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/macrophage_heatmap.pdf",width=7,height=7)
grid.draw(g)
dev.off()

