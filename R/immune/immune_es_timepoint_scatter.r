library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(survival)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "SELECT ps.*, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b,
cs.idh_codel_subtype
FROM analysis.rna_silver_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.tumor_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.tumor_barcode_b AND im2.signature_name = im1.signature_name
JOIN clinical.subtypes cs ON cs.case_barcode = ps.case_barcode
ORDER BY 2, 5
"

dat <- dbGetQuery(con, q)

dat[,"recur_status"] <- rep(1,nrow(dat))

cells <- unique(dat[,"signature_name"])

mycor <- mypval <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	mycor[i] <- cor(sub_dat[,"es_a"],sub_dat[,"es_b"],method="p")
	mypval[i] <- cor.test(sub_dat[,"es_a"],sub_dat[,"es_b"],method="p")$p.value
}

res <- data.frame(cells,mycor,mypval)

#Scatterplot and waterfall plot of correlations

#Waterfall plot
wf_res <- res[which(res[,"cells"] %in% c("B.cells","CD4.mature","CD8.effector","Dendritic","Macrophages","NK.cells","T.reg")),]
wf_res[,"cells"] <- c("B cells", "CD4+ T cells", "CD8+ T cells", "Dendritic cells", "Macrophages","NK cells","Tregs")
wf_res[,"cells"] <- factor(wf_res[,"cells"],levels=wf_res[order(wf_res[,"mycor"],decreasing=TRUE),"cells"])


pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/init_rec_immune_waterfall.pdf",width=3,height=1.75)
ggplot(wf_res,aes(x = cells, y = mycor)) +
geom_bar(stat="identity") +
labs(y = "Pearson coefficient") +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0,0.6))
dev.off()

#Scatterplots (2 examples)
mdat <- dat[which(dat[,"signature_name"] == "Macrophages"),] 
tdat <- dat[which(dat[,"signature_name"] == "CD8.effector"),] 

cor1 <- round(wf_res[which(wf_res[,"cells"]=="Macrophages"),"mycor"],2)
cor2 <- round(wf_res[which(wf_res[,"cells"]=="CD8+ T cells"),"mycor"],2)

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/init_rec_immune_scatterplot.pdf",width=3,height=1.5)

p1 <- ggplot(mdat,aes(x = es_a, y = es_b)) +
geom_point(size=0.5, aes(colour = idh_codel_subtype)) +
geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
labs(x = "Initial", y = "Recurrent", title="Macrophages") +
annotate(geom="text", size=2.5, 
	x=max(mdat[,"es_a"]) - (0.23*(max(mdat[,"es_a"]) - min(mdat[,"es_a"]))), 
	y=max(mdat[,"es_b"]) - (0.95*(max(mdat[,"es_b"]) - min(mdat[,"es_b"]))), 
	label=deparse(bquote(italic("R")~" = "~.(cor1))),parse=TRUE) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 

p2 <- ggplot(tdat,aes(x = es_a, y = es_b)) +
geom_point(size=0.5, aes(colour = idh_codel_subtype)) +
geom_smooth(method='lm',formula= y~x, colour="black",size=0.65) +
labs(x = "Initial", y = "Recurrent", title="CD8+ T cells") +
annotate(geom="text", size=2.5, 
	x=max(tdat[,"es_a"]) - (0.23*(max(tdat[,"es_a"]) - min(tdat[,"es_a"]))), 
	y=max(tdat[,"es_b"]) - (0.95*(max(tdat[,"es_b"]) - min(tdat[,"es_b"]))), 
	label=deparse(bquote(italic("R")~" = "~.(cor2))), parse=TRUE) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 

grid.arrange(p1,p2,ncol=2)
dev.off()
