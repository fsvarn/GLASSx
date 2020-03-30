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

original_scaled <- plot_scaled

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Make IDHwt plot
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

plot_scaled <- original_scaled[which(original_scaled[,"subtype"]=="IDHwt"),]

#Make a matrix for clustering
plot_mat <- matrix(0,nrow=length(unique(plot_scaled[,"signature_name"])),ncol=length(unique(plot_scaled[,"aliquot_barcode"])))
rownames(plot_mat) <- unique(plot_scaled[,"signature_name"])
colnames(plot_mat) <- unique(plot_scaled[,"aliquot_barcode"])
for(i in 1:nrow(plot_scaled))
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
                               
#Extract cluster information from dendrogram
cut_dendro <- cut(dendro,11)
clusters <- rep(0,nrow(plot_scaled))
clusters[which(plot_scaled[,"aliquot_barcode"] %in% labels(cut_dendro[[2]][[2]]))] <- 1
clusters <- as.factor(clusters)

#Add cluster information to scaled table
plot_scaled[,"clusters"] <- clusters

#Create sidebar with cluster information
heights <- rep(1,length(clusters))
clusters <- data.frame(plot_scaled[,"aliquot_barcode"],plot_scaled[,"subtype"],clusters,heights)
colnames(clusters) <- c("aliquot_barcode","subtype","cluster","height")
clusters <- clusters[-which(duplicated(clusters)),]

#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_sidebars_idhwt.pdf",width=4,height=1)
sidebars <- ggplot(data = clusters, aes(x = aliquot_barcode, y = height)) +
	geom_bar(aes(fill=cluster),stat="identity") +
	scale_fill_manual(values=c("#26ABE2", "#BD1E2D")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
sidebars
#dev.off()

# Create heatmap plot

#Old figure using heatmap.2
#heatmap.2(t(plot_mat[plot.order,]),trace="none",density.info="none", dendrogram="none", col=colorRampPalette(c("navy","white","red"),bias=1,space="rgb"))
#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_heatmap_idhwt.pdf",width=4,height=1)
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

#Identify cases that switched clusters
uni.aliquots <- plot_scaled[,c("aliquot_barcode","case_barcode","status","subtype","clusters")]
uni.aliquots <- uni.aliquots[-which(duplicated(uni.aliquots)),]
uni.cases <- as.character(unique(uni.aliquots[,"case_barcode"]))
cluster_change <- rep("",length(uni.cases))
names(cluster_change) <- as.character(uni.cases)
for(i in 1:length(uni.cases))
{
	mycase <- uni.cases[i]
	sub.aliquots <-uni.aliquots[which(uni.aliquots[,"case_barcode"]==mycase),c("case_barcode","aliquot_barcode","status","clusters")]
	c1 <- as.numeric(sub.aliquots[which(sub.aliquots[,"status"]=="Initial"),"clusters"])
	c2 <- as.numeric(sub.aliquots[which(sub.aliquots[,"status"]=="Recurrent"),"clusters"])
	if (c1 == c2){
		cluster_change[mycase] <- "None"
	} else if (c2 < c1){
		cluster_change[mycase] <- "Decrease"
	} else if (c2 > c1){
		cluster_change[mycase] <- "Increase"}	
}
cluster_change <- rep(cluster_change,2) #Double vector so that initial and recurrent tumors are tied to it
uni.aliquots <- data.frame(uni.aliquots,cluster_change)
uni.aliquots[,"cluster_change"] <- factor(uni.aliquots[,"cluster_change"],levels=c("None","Decrease","Increase"))
uni.aliquots[,c("aliquot_barcode","case_barcode","clusters","cluster_change")]


#Plot change over time for below heatmap
#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_alignment_idhmut.pdf",width=4,height=4)
align.plot <- ggplot(data = uni.aliquots, aes(x = aliquot_barcode, y = case_barcode)) +
  geom_line(aes(group = case_barcode,colour=cluster_change),size=0.25) +
  geom_tile(aes(fill = status)) +
  scale_fill_manual(values=c("royalblue4","tomato3")) +
  scale_colour_manual(values=c("black","#26ABE2", "#BD1E2D")) +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
align.plot 
#dev.off()

#Check whether cluster change truly associates with immune change
#Increase
myinc <- unique(as.character(uni.aliquots[which(uni.aliquots[,"cluster_change"]=="Increase"),"case_barcode"]))
inc_dat <- dat[which(dat[,"case_barcode"] %in% myinc),]	
inc.pval <- inc.eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_inc <- inc_dat[which(inc_dat[,"signature_name"]==cells[i]),]
	inc.pval[i] <- wilcox.test(sub_inc[,"es_a"],sub_inc[,"es_b"],paired=TRUE)$p.value
	inc.eff[i] <- median(sub_inc[,"es_b"]) - median(sub_inc[,"es_a"])
}
inc.res <- data.frame(cells,inc.eff,inc.pval)

#Decrease
mydec <- unique(as.character(uni.aliquots[which(uni.aliquots[,"cluster_change"]=="Decrease"),"case_barcode"]))
dec_dat <- dat[which(dat[,"case_barcode"] %in% mydec),]	
dec.pval <- dec.eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dec <- dec_dat[which(dec_dat[,"signature_name"]==cells[i]),]
	dec.pval[i] <- wilcox.test(sub_dec[,"es_a"],sub_dec[,"es_b"],paired=TRUE)$p.value
	dec.eff[i] <- median(sub_dec[,"es_b"]) - median(sub_dec[,"es_a"])
}
dec.res <- data.frame(cells,dec.eff,dec.pval)

#No change
mync <- unique(as.character(uni.aliquots[which(uni.aliquots[,"cluster_change"]=="None"),"case_barcode"]))
nc_dat <- dat[which(dat[,"case_barcode"] %in% mync),]	
nc.pval <- nc.eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_nc <- nc_dat[which(nc_dat[,"signature_name"]==cells[i]),]
	nc.pval[i] <- wilcox.test(sub_nc[,"es_a"],sub_nc[,"es_b"],paired=TRUE)$p.value
	nc.eff[i] <- median(sub_nc[,"es_b"]) - median(sub_nc[,"es_a"])
}
nc.res <- data.frame(cells,nc.eff,nc.pval)

#Build cluster change table
change <- c(rep("Increase",length(myinc)),rep("Decrease",length(mydec)),rep("None",length(mync)))
change_case <- c(myinc,mydec,mync)
change_res <- data.frame(change_case,change)
colnames(change_res) <- c("case_barcode","change")

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(sidebars)
gb3 <- ggplot_build(heatmap.plot)
gb4 <- ggplot_build(align.plot)

n1 <- length(gb1$layout$panel_params[[1]]$y.labels)
n2 <- length(gb2$layout$panel_params[[1]]$y.labels)
n3 <- length(gb3$layout$panel_params[[1]]$y.labels)
n4 <- length(gb4$layout$panel_params[[1]]$y.labels)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n3/3, "null")
g$heights[panels[2]] <- unit(n3/10, "null")
g$heights[panels[3]] <- unit(n3,"null")
g$heights[panels[4]] <- unit(n3,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_idhwt.pdf",width=4,height=4)
grid.draw(g)
dev.off()

#Heatmap with row labels (cells)
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_idhwt_label.pdf",width=4,height=4)
heatmap.plot <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = signature_name)) +
  geom_tile(aes(fill=es)) +
  scale_fill_gradient2(low ="royalblue4", mid="white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  theme_void() +
  theme(axis.text.x = element_blank(),
		axis.text.y = element_text(size=7),
        axis.ticks = element_blank(),
        legend.position = "none")
heatmap.plot
dev.off()



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Make IDHmut plot
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

plot_scaled <- original_scaled[which(original_scaled[,"subtype"] %in% c("IDHmut-noncodel","IDHmut-codel")),]

#Make a matrix for clustering
plot_mat <- matrix(0,nrow=length(unique(plot_scaled[,"signature_name"])),ncol=length(unique(plot_scaled[,"aliquot_barcode"])))
rownames(plot_mat) <- unique(plot_scaled[,"signature_name"])
colnames(plot_mat) <- unique(plot_scaled[,"aliquot_barcode"])
for(i in 1:nrow(plot_scaled))
{
	mysigname <- as.character(plot_scaled[i,"signature_name"])
	myalicode <- as.character(plot_scaled[i,"aliquot_barcode"])
	plot_mat[mysigname,myalicode] <- plot_scaled[i,"es"]
}
plot_mat <- t(plot_mat)

# Run clustering
dendro <- as.dendrogram(hclust(d = dist(x = plot_mat)))
#Reverse to get high on the right and low on the left
#dendro <- rev(dendro)

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
                               
#Extract cluster information from dendrogram
cut_dendro <- cut(dendro,11)
clusters <- rep(0,nrow(plot_scaled))
clusters[which(plot_scaled[,"aliquot_barcode"] %in% labels(cut_dendro[[2]][[2]]))] <- 1
clusters <- as.factor(clusters)

#Add cluster information to scaled table
plot_scaled[,"clusters"] <- clusters

#Create sidebar with cluster information
heights <- rep(1,length(clusters))
clusters <- data.frame(plot_scaled[,"aliquot_barcode"],plot_scaled[,"subtype"],clusters,heights)
colnames(clusters) <- c("aliquot_barcode","subtype","cluster","height")
clusters <- clusters[-which(duplicated(clusters)),]

#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_sidebars_idhwt.pdf",width=4,height=1)
sidebars <- ggplot(data = clusters, aes(x = aliquot_barcode, y = height)) +
	geom_bar(aes(fill=cluster),stat="identity") +
	scale_fill_manual(values=c("#26ABE2", "#BD1E2D")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
sidebars
#dev.off()

#Create sidebar with codel information
subtype_bar <- ggplot(data = clusters, aes(x = aliquot_barcode, y = height)) +
	geom_bar(aes(fill=subtype),stat="identity") +
	scale_fill_manual(values=c("#F8766D", "#00BA38")) +
	theme_void() +
	theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
subtype_bar

# Create heatmap plot

#Old figure using heatmap.2
#heatmap.2(t(plot_mat[plot.order,]),trace="none",density.info="none", dendrogram="none", col=colorRampPalette(c("navy","white","red"),bias=1,space="rgb"))
#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_heatmap_idhwt.pdf",width=4,height=1)
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

#Identify cases that switched clusters
uni.aliquots <- plot_scaled[,c("aliquot_barcode","case_barcode","status","subtype","clusters")]
uni.aliquots <- uni.aliquots[-which(duplicated(uni.aliquots)),]
uni.cases <- as.character(unique(uni.aliquots[,"case_barcode"]))
cluster_change <- rep("",length(uni.cases))
names(cluster_change) <- as.character(uni.cases)
for(i in 1:length(uni.cases))
{
	mycase <- uni.cases[i]
	sub.aliquots <-uni.aliquots[which(uni.aliquots[,"case_barcode"]==mycase),c("case_barcode","aliquot_barcode","status","clusters")]
	c1 <- as.numeric(sub.aliquots[which(sub.aliquots[,"status"]=="Initial"),"clusters"])
	c2 <- as.numeric(sub.aliquots[which(sub.aliquots[,"status"]=="Recurrent"),"clusters"])
	if (c1 == c2){
		cluster_change[mycase] <- "None"
	} else if (c2 < c1){
		cluster_change[mycase] <- "Decrease"
	} else if (c2 > c1){
		cluster_change[mycase] <- "Increase"}	
}
cluster_change <- rep(cluster_change,2) #Double vector so that initial and recurrent tumors are tied to it
uni.aliquots <- data.frame(uni.aliquots,cluster_change)
uni.aliquots[,"cluster_change"] <- factor(uni.aliquots[,"cluster_change"],levels=c("None","Decrease","Increase"))
uni.aliquots[,c("aliquot_barcode","case_barcode","clusters","cluster_change")]


#Plot change over time for below heatmap
#pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_alignment_idhmut.pdf",width=4,height=4)
align.plot <- ggplot(data = uni.aliquots, aes(x = aliquot_barcode, y = case_barcode)) +
  geom_line(aes(group = case_barcode,colour=cluster_change),size=0.25) +
  geom_tile(aes(fill = status)) +
  scale_fill_manual(values=c("royalblue4","tomato3")) +
  scale_colour_manual(values=c("black","#26ABE2", "#BD1E2D")) +
  theme_void() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")
align.plot 
#dev.off()

#Check whether cluster change truly associates with immune change
#Increase
myinc <- unique(as.character(uni.aliquots[which(uni.aliquots[,"cluster_change"]=="Increase"),"case_barcode"]))
inc_dat <- dat[which(dat[,"case_barcode"] %in% myinc),]	
inc.pval <- inc.eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_inc <- inc_dat[which(inc_dat[,"signature_name"]==cells[i]),]
	inc.pval[i] <- wilcox.test(sub_inc[,"es_a"],sub_inc[,"es_b"],paired=TRUE)$p.value
	inc.eff[i] <- median(sub_inc[,"es_b"]) - median(sub_inc[,"es_a"])
}
inc.res <- data.frame(cells,inc.eff,inc.pval)

#Decrease
mydec <- unique(as.character(uni.aliquots[which(uni.aliquots[,"cluster_change"]=="Decrease"),"case_barcode"]))
dec_dat <- dat[which(dat[,"case_barcode"] %in% mydec),]	
dec.pval <- dec.eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dec <- dec_dat[which(dec_dat[,"signature_name"]==cells[i]),]
	dec.pval[i] <- wilcox.test(sub_dec[,"es_a"],sub_dec[,"es_b"],paired=TRUE)$p.value
	dec.eff[i] <- median(sub_dec[,"es_b"]) - median(sub_dec[,"es_a"])
}
dec.res <- data.frame(cells,dec.eff,dec.pval)

#No change
mync <- unique(as.character(uni.aliquots[which(uni.aliquots[,"cluster_change"]=="None"),"case_barcode"]))
nc_dat <- dat[which(dat[,"case_barcode"] %in% mync),]	
nc.pval <- nc.eff <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_nc <- nc_dat[which(nc_dat[,"signature_name"]==cells[i]),]
	nc.pval[i] <- wilcox.test(sub_nc[,"es_a"],sub_nc[,"es_b"],paired=TRUE)$p.value
	nc.eff[i] <- median(sub_nc[,"es_b"]) - median(sub_nc[,"es_a"])
}
nc.res <- data.frame(cells,nc.eff,nc.pval)

#Add to cluster change table and write to db
change <- c(rep("Increase",length(myinc)),rep("Decrease",length(mydec)),rep("None",length(mync)))
change_case <- c(myinc,mydec,mync)
change_add <- data.frame(change_case,change)
colnames(change_add) <- c("case_barcode","change")
change_res <- rbind(change_res,change_add)
dbWriteTable(con, Id(schema="analysis", table="immune_cluster_change"), change_res, overwrite=TRUE)


#Scale factor for IDHmut (to make size comparable to IDHwt)
sf <- sum(dat[,"subtype_a"]=="IDHmut-noncodel" | dat[,"subtype_a"]=="IDHmut-codel")/length(cells)/(sum(dat[,"subtype_a"]=="IDHwt")/length(cells))

#Align figures for printing
gb1 <- ggplot_build(dendro.plot)
gb2 <- ggplot_build(sidebars)
gb3 <- ggplot_build(heatmap.plot)
gb4 <- ggplot_build(align.plot)

n1 <- length(gb1$layout$panel_params[[1]]$y.labels)
n2 <- length(gb2$layout$panel_params[[1]]$y.labels)
n3 <- length(gb3$layout$panel_params[[1]]$y.labels)
n4 <- length(gb4$layout$panel_params[[1]]$y.labels)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1]] <- unit(n3/3, "null")
g$heights[panels[2]] <- unit(n3/10, "null")
g$heights[panels[3]] <- unit(n3,"null")
g$heights[panels[4]] <- unit(n3,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_idhmut.pdf",width=4*sf,height=4)
grid.draw(g)
dev.off()

#Heatmap with row labels (cells)
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_heatmap_idhmut_label.pdf",width=4,height=4)
heatmap.plot <- ggplot(data = plot_scaled, aes(x = aliquot_barcode, y = signature_name)) +
  geom_tile(aes(fill=es)) +
  scale_fill_gradient2(low ="royalblue4", mid="white", high = "tomato3",
  space = "Lab", na.value = "grey50", guide = "colourbar",
  aesthetics = "fill") +
  theme_void() +
  theme(axis.text.x = element_blank(),
		axis.text.y = element_text(size=7),
        axis.ticks = element_blank())
heatmap.plot
dev.off()
