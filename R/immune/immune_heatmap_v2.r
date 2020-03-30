library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggdendro)
library(grid)

#######################################################

rm(list=ls()) 

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "
SELECT ps.case_barcode, 
ps.tumor_barcode_a, 
ps.tumor_barcode_b, 
bs1.sample_type AS sample_type_a,
bs2.sample_type AS sample_type_b,
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs1.grade AS grade_a,
cs2.grade AS grade_b, 
cs1.idh_codel_subtype AS subtype_a, 
cs2.idh_codel_subtype AS subtype_b
FROM analysis.rna_silver_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.tumor_barcode_b AND im2.signature_name = im1.signature_name
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
WHERE cs1.idh_codel_subtype IS NOT NULL
ORDER BY 1, 2, 6
"

dat <- dbGetQuery(con,q)

dat <- dat[which(dat[,"signature_name"] %in% c("B.cells","CD4.mature","CD8.effector","Dendritic","Macrophages","NK.cells","T.reg")),]
cells <- unique(dat[,"signature_name"])


#Figures

#Make a long format table for heatmap
es <- c(dat[,"es_a"],dat[,"es_b"])
signature_name <- rep(dat[,"signature_name"],2)
subtype <- c(dat[,"subtype_a"],dat[,"subtype_b"])
status <- c(rep("Initial",nrow(dat)), rep("Recurrent",nrow(dat)))
case_barcode <- rep(dat[,"case_barcode"],2)
aliquot_barcode <- c(dat[,"tumor_barcode_a"],dat[,"tumor_barcode_b"])
plot_res <- data.frame(aliquot_barcode, case_barcode, signature_name, es, subtype, status)

# Rescale each measurement (independently) to have a mean of 0 and variance of 1
plot_scaled <- plot_res
for(i in 1:length(cells))
{
	plot_scaled[which(plot_scaled[,"signature_name"]==cells[i]),"es"] <- scale(plot_scaled[which(plot_scaled[,"signature_name"]==cells[i]),"es"])
}

#Make a matrix for clustering
plot_mat <- matrix(0,nrow=length(unique(signature_name)),ncol=length(unique(aliquot_barcode)))
rownames(plot_mat) <- unique(signature_name)
colnames(plot_mat) <- unique(aliquot_barcode)
for(i in 1:nrow(plot_res))
{
	mysigname <- as.character(plot_scaled[i,"signature_name"])
	myalicode <- as.character(plot_scaled[i,"aliquot_barcode"])
	plot_mat[mysigname,myalicode] <- plot_scaled[i,"es"]
}
plot_mat <- t(plot_mat)

# Run clustering
dendro <- as.dendrogram(hclust(d = dist(x = plot_mat)))
#Reverse to get high on the right and low on the left
dendro <- rev(dendro)

# Create dendrogram plot
#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_dendrogram_idhmut.pdf",width=4,height=1)
dendro.plot <- ggdendrogram(data = dendro) + 
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") 
dendro.plot
#dev.off()
# Heatmap

# Data wrangling
# Extract the order of the samples in the dendrogram
plot.order <- order.dendrogram(dendro)
# Order the levels according to their position in the cluster
plot_scaled[,"aliquot_barcode"] <- factor(plot_scaled[,"aliquot_barcode"],
                               levels = rownames(plot_mat)[plot.order], 
                               ordered = TRUE)

# Create heatmap plot

#Old figure using heatmap.2
#heatmap.2(t(plot_mat[plot.order,]),trace="none",density.info="none", dendrogram="none", col=colorRampPalette(c("navy","white","red"),bias=1,space="rgb"))
#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_heatmap_idhmut.pdf",width=4,height=1)
heatmap.plot <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = signature_name)) +
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
#dev.off()

#print(heatmap.plot, vp = viewport(x = 0.5, y = 0.45, width = 1.0, height = 0.9))
#print(dendro.plot, vp = viewport(x = 0.525, y = 0.9, width = 1.0, height = 0.2))

#Reorder alignment plot by ordering case_barcode by where there aliquots first appear in dendrogram

#Define case order based on where its first aliquot appears in the dendrogram
aliquot_order <- rownames(plot_mat)[plot.order]
case_order <- rep("",length(aliquot_order))
for(i in 1:length(aliquot_order))
{
	sub_scaled <- plot_scaled[which(plot_scaled[,"aliquot_barcode"]==aliquot_order[i]),]
	case_order[i] <- as.character(sub_scaled[1,"case_barcode"])
}
case_order <- unique(case_order)
case_order <- rev(case_order)

plot_scaled[,"case_barcode"] <- factor(plot_scaled[,"case_barcode"],
                               levels = case_order,
                               ordered = TRUE)

#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_alignment_idhmut.pdf",width=4,height=4)
align.plot <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = case_barcode)) +
  geom_line(aes(group = case_barcode)) +
  geom_tile(aes(fill = status)) +
  scale_fill_manual(values=c("royalblue4","tomato3")) +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
align.plot 
#dev.off()

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(heatmap.plot)
gb3 <- ggplot_build(align.plot)

n1 <- length(gb1$layout$panel_params[[1]]$y.labels)
n2 <- length(gb2$layout$panel_params[[1]]$y.labels)
n3 <- length(gb3$layout$panel_params[[1]]$y.labels)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n3 * 0.5, "null")
g$heights[panels[2]] <- unit(n3 * 1,"null")
g$heights[panels[3]] <- unit(n3 * 3,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_all.pdf",width=7,height=7)
grid.draw(g)
dev.off()


#Draw scatterplot of simplicity scores
uni_dat <- dat[,c("tumor_pair_barcode","ss_a","ss_b","switch","subtype_a","subtype_b","idh_codel_subtype")]
uni_dat <- uni_dat[-which(duplicated(uni_dat)),]
cor1 <- round(cor(uni_dat[,"ss_a"],uni_dat[,"ss_b"]),2)
cor.test(uni_dat[,"ss_a"],uni_dat[,"ss_b"])

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/init_rec_simplicity_scatterplot.pdf",width=1.5,height=1.5)
p1 <- ggplot(uni_dat,aes(x = ss_a, y = ss_b)) +
geom_point(size=0.5, aes(colour = idh_codel_subtype)) +
geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
labs(x = "Initial", y = "Recurrent", title="Simplicity score") +
annotate(geom="text", size=2.5, 
	x=max(uni_dat[,"ss_a"]) - (0.23*(max(uni_dat[,"ss_a"]) - min(uni_dat[,"ss_a"]))), 
	y=max(uni_dat[,"ss_b"]) - (0.95*(max(uni_dat[,"ss_b"]) - min(uni_dat[,"ss_b"]))), 
	label=deparse(bquote(italic("R")~" = "~.(cor1))),parse=TRUE) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 
p1
dev.off()

