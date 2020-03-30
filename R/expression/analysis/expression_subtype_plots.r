library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(ggbeeswarm)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
WITH subtype_label AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	ts1.signature_name, 
	ts1.p_value AS pval_a,
	ts2.p_value AS pval_b,
	CASE WHEN ts1.p_value < 0.05 THEN ts1.signature_name ELSE NULL END AS subtype_a,
	CASE WHEN ts2.p_value < 0.05 THEN ts2.signature_name ELSE NULL END AS subtype_b
	FROM analysis.rna_silver_set ss
	JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
	WHERE ts1.p_value < 0.05 OR ts2.p_value < 0.05
),
subtype_switch AS
(
	SELECT sl.tumor_pair_barcode,
	sl.case_barcode,
	string_agg(subtype_a,',') AS subtype_a,
	string_agg(subtype_b,',') AS subtype_b
	FROM subtype_label sl
	GROUP BY tumor_pair_barcode, sl.case_barcode
	ORDER BY 1
)
SELECT ss.*, 
ts1.signature_name, 
ts1.enrichment_score AS es_a,
ts2.enrichment_score AS es_b,
ts1.p_value AS pval_a,
ts2.p_value AS pval_b,
si1.simplicity_score AS ss_a,
si2.simplicity_score AS ss_b,
CASE WHEN sw.subtype_a IS NULL THEN 'Mixed' ELSE sw.subtype_a END AS subtype_a,
CASE WHEN sw.subtype_b IS NULL THEN 'Mixed' ELSE sw.subtype_b END AS subtype_b,
CASE WHEN sw.subtype_a = sw.subtype_b THEN 'None' ELSE 'Switch' END AS switch,
su.treatment_tmz,
su.treatment_radiotherapy,
su.idh_codel_subtype
FROM analysis.rna_silver_set ss
JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.surgeries su ON su.sample_barcode = al1.sample_barcode
JOIN analysis.simplicity_score si1 ON ss.tumor_barcode_a = si1.aliquot_barcode
JOIN analysis.simplicity_score si2 ON ss.tumor_barcode_b = si2.aliquot_barcode
JOIN subtype_switch sw ON ss.tumor_pair_barcode = sw.tumor_pair_barcode
WHERE idh_codel_subtype IS NOT NULL -- ts1.p_value < 0.05 OR ts2.p_value < 0.05
ORDER BY 1, 2 
"

dat <- dbGetQuery(con, q)

plot_res1 <- dat[,c("case_barcode","tumor_barcode_a","signature_name","pval_a","ss_a","switch","idh_codel_subtype")]
plot_res2 <- dat[,c("case_barcode","tumor_barcode_b","signature_name","pval_b","ss_b","switch","idh_codel_subtype")]
status <- c(rep("Initial",nrow(plot_res1)),rep("Recurrent",nrow(plot_res2)))
colnames(plot_res1) <- colnames(plot_res2) <- c("case_barcode","aliquot_barcode","signature_name","pval","ss","switch","idh_codel_subtype")
plot_res <- rbind(plot_res1,plot_res2)
plot_res <- cbind(plot_res,status)

aliquots <- unique(plot_res[,"aliquot_barcode"])
fraction <- c()
for(i in 1:length(aliquots))
{
	sub_dat <- plot_res[which(plot_res[,"aliquot_barcode"]==aliquots[i]),]
	tmp_fraction <- as.numeric(sub_dat[,"pval"] < 0.05)/sum(as.numeric(sub_dat[,"pval"] < 0.05))
	if(is.na(sum(tmp_fraction)))
	{
		tmp_fraction[which(sub_dat[,"pval"]==min(sub_dat[,"pval"]))] <- 1
		tmp_fraction[which(sub_dat[,"pval"]!=min(sub_dat[,"pval"]))] <- 0
		tmp_fraction <- tmp_fraction/sum(tmp_fraction)
	}
	fraction <- c(fraction,tmp_fraction)
}

plot_res <- cbind(plot_res,fraction)
plot_res <- plot_res[which(!is.na(plot_res[,"idh_codel_subtype"])),]
plot_res[,"pval"] <- -log(plot_res[,"pval"])

#Set order for x-axis of plot (will be by initial tumor simplicity score high to low)
myorder <- dat[,c("case_barcode","ss_a")]
myorder <- myorder[-which(duplicated(myorder)),]
myorder <- myorder[order(myorder[,"ss_a"],decreasing=TRUE),"case_barcode"]
plot_res[,"case_barcode"] <- factor(plot_res[,"case_barcode"], levels = myorder)

#Set order for transcript subtypes (should match IDH-codel subtypes)
plot_res[,"signature_name"] <- factor(plot_res[,"signature_name"], levels = c("Mesenchymal","Classical","Proneural"))

######################## 
## Common plotting elements
########################

plot_grid     <- facet_grid(. ~ idh_codel_subtype, scales = "free_x", space = "free")
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

########################
## Blank plot
########################

gg_blank <-
  ggplot(data.frame()) +
  geom_blank()

gg_simplicity_score_ini <-
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Mesenchymal"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, alpha= -1 * ss)) +
  scale_fill_manual(values=c("black")) +
  theme(axis.title.y=element_blank())
  
gg_transcript_subtype_ini <-
  ggplot(plot_res %>% filter(status=="Initial"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#00458a","#008a22","#8a0000")) +
  theme(axis.title.y=element_blank())
  
gg_simplicity_score_rec <-
  ggplot(plot_res %>% filter(status=="Recurrent",signature_name=="Mesenchymal"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, alpha= -1 * ss)) +
  scale_fill_manual(values=c("black")) +
  theme(axis.title.y=element_blank())

gg_transcript_subtype_rec <-
  ggplot(plot_res %>% filter(status=="Recurrent"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#00458a","#008a22","#8a0000")) +
  theme(axis.title.y=element_blank())
  
gg_subtype_switch <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Mesenchymal"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=switch)) +
  scale_fill_manual(values=c("white","black")) +
  theme(axis.title.y=element_blank())
  
#Align figures for printing
gb1 <- ggplot_build(nolabels(gg_simplicity_score_ini))
gb2 <- ggplot_build(nolabels(gg_transcript_subtype_ini))
gb3 <- ggplot_build(nolabels(gg_simplicity_score_rec))
gb4 <- ggplot_build(nolabels(gg_transcript_subtype_rec))
gb5 <- ggplot_build(nolabels(gg_subtype_switch))

n1 <- length(gb1$layout$panel_params[[1]]$y.labels)
n2 <- length(gb2$layout$panel_params[[1]]$y.labels)
n3 <- length(gb3$layout$panel_params[[1]]$y.labels)
n4 <- length(gb4$layout$panel_params[[1]]$y.labels)
n5 <- length(gb5$layout$panel_params[[1]]$y.labels)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1*3]] <- unit(n2/3, "null")
g$heights[panels[2*3]] <- unit(n2, "null")
g$heights[panels[3*3]] <- unit(n2/3,"null")
g$heights[panels[4*3]] <- unit(n2,"null")
g$heights[panels[5*3]] <- unit(n2/3,"null")

grid.newpage()

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/transcriptional_subtype_change.pdf",width=7,height=3)
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


#Boxplot of subtype switching
sub_dat <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-codel"),]
wilcox.test(sub_dat[which(sub_dat[,"switch"]=="Switch"),"ss_a"],sub_dat[which(sub_dat[,"switch"]=="None"),"ss_a"])

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/subtype_simplicity_boxplot.pdf",width=2,height=2)
p1 <- ggplot(dat, aes(x = idh_codel_subtype, y = ss_a, fill=switch)) +
geom_boxplot() +
#geom_jitter(position=position_jitterdodge(jitter.width=0.4),size=0.4) +
#geom_beeswarm(size=0.2,dodge.width=1) +
labs(y = "Simplicity score") +
stat_summary(fun.y=median,geom="line",colour="black") +
scale_fill_manual(values=c("lightblue","pink")) +
#facet_grid(.~idh_codel_subtype) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0.,1.1))
p1
dev.off()


#Old figures
# pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/transcriptional_subtype_change.pdf",width=7,height=3)
# grid.arrange(nolabels(gg_simplicity_score_ini),nolabels(gg_transcript_subtype_ini),nolabels(gg_simplicity_score_rec),nolabels(gg_transcript_subtype_rec),nrow=4)
# #grid.arrange(nolabels(gg_transcript_subtype_ini),nolabels(gg_transcript_subtype_rec),nrow=2)
# dev.off()

# 
# pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/transcriptional_subtype_change.pdf",width=2,height=2)
# ggplot(plot_res,aes(x = case_barcode, y = fraction, fill=signature_name)) +
# geom_bar(stat="identity",position=position_stack()) +
# #scale_colour_manual(values=c("#00BA38","#619CFF")) +
# labs(y = "Transcriptional subtype") +
# facet_grid(idh_codel_subtype ~ status) +
# theme_bw() +
# theme(axis.text= element_text(size=7),
# axis.title = element_text(size=7),
# plot.title = element_text(size=7,hjust=0.5),
# panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
# legend.position="none")  +
# coord_cartesian(ylim=c(0,1))
# dev.off()


