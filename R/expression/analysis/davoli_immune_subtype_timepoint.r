library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(ggbeeswarm)


#######################################################

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
ORDER BY 1, 2, 6
"

dat <- dbGetQuery(con,q)

#See if there are significant changes across samples over time

cells <- unique(dat[,"signature_name"])

eff <- p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	eff[i] <- median(sub_dat[,"es_b"]) - median(sub_dat[,"es_a"])
	p.value[i] <- wilcox.test(sub_dat[,"es_a"], sub_dat[,"es_b"],paired=TRUE)$p.value
}
time_res <- data.frame(cells,eff,p.value)

#See if there are significant changes  over time when stratifying by subtype 

cells <- unique(dat[,"signature_name"])
subtypes <- unique(dat[,"subtype_a"])

eff <- p.value <- rep(0,length(cells)*length(subtypes))
ct <- 1
for(i in 1:length(subtypes))
{
	sub_dat <- dat[which(dat[,"subtype_a"]==subtypes[i]),]

	for(j in 1:length(cells))
	{
		cell_dat <- sub_dat[which(sub_dat[,"signature_name"]==cells[j]),]
		eff[ct] <- median(cell_dat[,"es_b"]) - median(cell_dat[,"es_a"])
		p.value[ct] <- wilcox.test(cell_dat[,"es_a"], cell_dat[,"es_b"],paired=TRUE)$p.value
		ct <- ct + 1
	}
}
cell <- rep(cells,length(subtypes))
subtype <- rep(subtypes,each=length(cells))
time_subtype_res <- data.frame(subtype,cells,eff,p.value)


#See if there are significant differences between subtypes from initial/recurrent samples

init.p.value <- rec.p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	sub_dat[,"subtype_a"] <- as.factor(sub_dat[,"subtype_a"])
	sub_dat[,"subtype_b"] <- as.factor(sub_dat[,"subtype_b"])
	
	init.p.value[i] <- kruskal.test(es_a ~ subtype_a, data = sub_dat)$p.value
	rec.p.value[i] <- kruskal.test(es_b ~ subtype_b, data = sub_dat) $p.value
}
subtype_kw <- data.frame(cells,init.p.value,rec.p.value)


#See whether IDHwt or IDHmut are higher

init.eff <- init.p.value <- rec.eff <- rec.p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	#Subtypes are the same in both initial and recurrent tumors
	g1 <- sub_dat[which(sub_dat[,"subtype_a"]=="IDHwt"),]
	g2 <- sub_dat[which(sub_dat[,"subtype_a"]=="IDHmut-noncodel"),]

	init.eff[i] <- median(g1[,"es_a"]) - median(g2[,"es_a"])
	init.p.value[i] <- wilcox.test(g1[,"es_a"],g2[,"es_a"])$p.value
	rec.eff[i] <- median(g1[,"es_b"]) - median(g2[,"es_b"])
	rec.p.value[i] <- wilcox.test(g1[,"es_b"],g2[,"es_b"])$p.value
}
idh_diffs1 <- data.frame(cells,init.eff,init.p.value,rec.eff,rec.p.value)

init.eff <- init.p.value <- rec.eff <- rec.p.value <- rep(0,length(cells))
for(i in 1:length(cells))
{
	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	#Subtypes are the same in both initial and recurrent tumors
	g1 <- sub_dat[which(sub_dat[,"subtype_a"]=="IDHwt"),]
	g2 <- sub_dat[which(sub_dat[,"subtype_a"]!="IDHwt"),]

	init.eff[i] <- median(g1[,"es_a"]) - median(g2[,"es_a"])
	init.p.value[i] <- wilcox.test(g1[,"es_a"],g2[,"es_a"])$p.value
	rec.eff[i] <- median(g1[,"es_b"]) - median(g2[,"es_b"])
	rec.p.value[i] <- wilcox.test(g1[,"es_b"],g2[,"es_b"])$p.value
}
idh_diffs2 <- data.frame(cells,init.eff,init.p.value,rec.eff,rec.p.value)

#Figures
es <- c(dat[,"es_a"],dat[,"es_b"])
signature_name <- rep(dat[,"signature_name"],2)
subtype <- c(dat[,"subtype_a"],dat[,"subtype_b"])
status <- c(rep("Initial",nrow(dat)), rep("Recurrent",nrow(dat)))
case_barcode <- rep(dat[,"case_barcode"],2)
plot_res <- data.frame(case_barcode, signature_name, es, subtype, status)

es <- c(dat[,"es_a"],dat[,"es_b"])
signature_name <- rep(dat[,"signature_name"],2)
subtype <- c(dat[,"subtype_a"],dat[,"subtype_b"])
status <- c(rep("Initial",nrow(dat)), rep("Recurrent",nrow(dat)))
case_barcode <- rep(dat[,"case_barcode"],2)
plot_res <- data.frame(case_barcode, signature_name, es, subtype, status)

plot_res <- plot_res[which(plot_res[,"signature_name"] %in% c("B.cells","CD4.mature","CD8.effector","Dendritic","Macrophages","NK.cells","T.reg")),]

#Plot showing all data
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/immune_cell_subtype_timepoint.pdf",width=7,height=2)
ggplot(plot_res,aes(x = status, y = es, fill=subtype)) +
geom_boxplot(outlier.size=0) +
#geom_beeswarm(size=0.3) +
#stat_summary(fun.y=median,geom="line",colour="red") +
labs(y = "Enrichment score") +
facet_grid(.~signature_name) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none")
dev.off()

#Plot showing ladder plots between initial and recurrence for all subtypes
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/all_subtypes_ladders.pdf",width=7,height=2)
ggplot(plot_res,aes(x = status, y = es, group= case_barcode, colour=subtype)) +
geom_line(size=0.45,alpha=0.4) +
geom_point(size=1,colour="black") +
#stat_summary(fun.y=median,geom="line",colour="red") +
labs(y = "Enrichment score") +
facet_grid(.~signature_name) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none")
dev.off()

#Plot showing the one significant difference (CD4+ T cells in IDHwt)
idhwt_cd4 <- plot_res[which(plot_res[,"subtype"] == "IDHwt" & plot_res[,"signature_name"] == "CD4.mature"),]

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/idhwt_cd4_ladderplot.pdf",width=1.5,height=2)
ggplot(idhwt_cd4,aes(x = status, y = es, colour=subtype)) +
geom_boxplot(outlier.size=0,colour="black") +
geom_line(size=0.8,alpha=0.4,aes(group=case_barcode)) +
geom_point(size=1,colour="black") +
#stat_summary(fun.y=median,geom="line",colour="red") +
scale_colour_manual(values="#619CFF") +
labs(y = "Enrichment score") +
facet_grid(.~signature_name) +
theme_bw() +
theme(axis.text.x= element_text(size=7,angle=45,hjust=1),axis.text.y= element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") +
coord_cartesian(ylim=c(-0.25,0.2))
dev.off()

