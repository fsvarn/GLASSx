library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

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
),
roman_to_int(grade, grade_int) AS 
(
	VALUES ('I'::text,1), ('II'::text,2), ('III'::text,3), ('IV'::text,4)
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
CASE
	WHEN r2.grade_int > r1.grade_int THEN 'Grade up'::text
	WHEN r2.grade_int = r1.grade_int THEN 'Grade stable'::text
	WHEN r2.grade_int < r1.grade_int THEN 'Grade down'::text
	ELSE NULL::text
END AS grade_change,
su1.treatment_tmz,
su1.treatment_radiotherapy,
CASE WHEN su1.treatment_chemotherapy_other LIKE '%Nivolumab%' OR su1.treatment_chemotherapy_other LIKE '%Pembrolizumab%' THEN true ELSE false END AS treatment_pd1 ,
su1.idh_codel_subtype
FROM analysis.rna_silver_set ss
JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = ss.tumor_barcode_a
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = ss.tumor_barcode_b
JOIN clinical.surgeries su1 ON su1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries su2 ON su2.sample_barcode = al2.sample_barcode
JOIN roman_to_int r1 ON r1.grade = su1.grade::text
JOIN roman_to_int r2 ON r2.grade = su2.grade::text
JOIN analysis.simplicity_score si1 ON ss.tumor_barcode_a = si1.aliquot_barcode
JOIN analysis.simplicity_score si2 ON ss.tumor_barcode_b = si2.aliquot_barcode
JOIN subtype_switch sw ON ss.tumor_pair_barcode = sw.tumor_pair_barcode
WHERE su1.idh_codel_subtype IS NOT NULL -- ts1.p_value < 0.05 OR ts2.p_value < 0.05
ORDER BY 1, 2 
"

dat <- dbGetQuery(con, q)

plot_res1 <- dat[,c("case_barcode","tumor_barcode_a","signature_name","pval_a","ss_a","switch","grade_change","treatment_tmz","treatment_radiotherapy","treatment_pd1","idh_codel_subtype")]
plot_res2 <- dat[,c("case_barcode","tumor_barcode_b","signature_name","pval_b","ss_b","switch","grade_change","treatment_tmz","treatment_radiotherapy","treatment_pd1","idh_codel_subtype")]
status <- c(rep("Initial",nrow(plot_res1)),rep("Recurrent",nrow(plot_res2)))
colnames(plot_res1) <- colnames(plot_res2) <- c("case_barcode","aliquot_barcode","signature_name","pval","ss","switch","grade_change","tmz","rt","aPD1","idh_codel_subtype")
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

#Set order for status
plot_res[,"status"] <- factor(plot_res[,"status"], levels = c("Recurrent","Initial"))

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

gg_simplicity_score <-
  ggplot(plot_res %>% filter(signature_name=="Mesenchymal"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, alpha= -1 * ss)) +
  scale_fill_manual(values=c("black")) +
  theme(axis.title.y=element_blank())
  
gg_transcript_subtype_pro <-
  ggplot(plot_res %>% filter(signature_name=="Proneural"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#00458a")) +
  theme(axis.title.y=element_blank())
  
gg_transcript_subtype_cla <-
  ggplot(plot_res %>% filter(signature_name=="Classical"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#008a22")) +
  theme(axis.title.y=element_blank())
  
gg_transcript_subtype_mes <-
  ggplot(plot_res %>% filter(signature_name=="Mesenchymal"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#8a0000")) +
  theme(axis.title.y=element_blank())
  
gg_subtype_switch <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=switch)) +
  scale_fill_manual(values=c("white","black")) +
  theme(axis.title.y=element_blank())
  
gg_grade_change <- 
  ggplot(plot_res %>% filter(status=="Recurrent",signature_name=="Proneural"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=grade_change)) +
  scale_fill_manual(values=c("white","black")) +
  theme(axis.title.y=element_blank())
  
gg_tmz <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=tmz)) +
  scale_fill_manual(values=c("white","black"),na.value="grey50") +
  theme(axis.title.y=element_blank())
  
gg_rt <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=rt)) +
  scale_fill_manual(values=c("white","black"),na.value="grey50") +
  theme(axis.title.y=element_blank())
  
gg_pd1 <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=aPD1)) +
  scale_fill_manual(values=c("white","black"),na.value="grey50") +
  theme(axis.title.y=element_blank())
  
#Align figures for printing
gb1 <- ggplot_build(nolabels(gg_simplicity_score))
gb2 <- ggplot_build(nolabels(gg_transcript_subtype_pro))
gb3 <- ggplot_build(nolabels(gg_transcript_subtype_cla))
gb4 <- ggplot_build(nolabels(gg_transcript_subtype_mes))
gb5 <- ggplot_build(nolabels(gg_subtype_switch))
gb6 <- ggplot_build(nolabels(gg_grade_change))
gb7 <- ggplot_build(nolabels(gg_tmz))
gb8 <- ggplot_build(nolabels(gg_rt))
gb9 <- ggplot_build(nolabels(gg_pd1))

n1 <- length(gb1$layout$panel_params[[1]]$y.labels)
n2 <- length(gb2$layout$panel_params[[1]]$y.labels)
n3 <- length(gb3$layout$panel_params[[1]]$y.labels)
n4 <- length(gb4$layout$panel_params[[1]]$y.labels)
n5 <- length(gb5$layout$panel_params[[1]]$y.labels)
n6 <- length(gb6$layout$panel_params[[1]]$y.labels)
n7 <- length(gb6$layout$panel_params[[1]]$y.labels)
n8 <- length(gb7$layout$panel_params[[1]]$y.labels)
n9 <- length(gb8$layout$panel_params[[1]]$y.labels)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)
gF <- ggplot_gtable(gb6)
gG <- ggplot_gtable(gb7)
gH <- ggplot_gtable(gb8)
gI <- ggplot_gtable(gb9)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")
g <- gtable:::rbind_gtable(g, gF, "last")
g <- gtable:::rbind_gtable(g, gG, "last")
g <- gtable:::rbind_gtable(g, gH, "last")
g <- gtable:::rbind_gtable(g, gI, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1*3]] <- unit(n1, "null")
g$heights[panels[2*3]] <- unit(n1, "null")
g$heights[panels[3*3]] <- unit(n1,"null")
g$heights[panels[4*3]] <- unit(n1,"null")
g$heights[panels[5*3]] <- unit(n1/2,"null")
g$heights[panels[6*3]] <- unit(n1/2,"null")
g$heights[panels[7*3]] <- unit(n1/2,"null")
g$heights[panels[8*3]] <- unit(n1/2,"null")
g$heights[panels[9*3]] <- unit(n1/2,"null")

grid.newpage()

#Plot
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/transcriptional_subtype_change_v2.pdf",width=7,height=3)
grid.draw(g)
dev.off()

#Legends
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/transcriptional_subtype_change_legends.pdf",width=7,height=7)
grid.arrange(gg_legend(gg_simplicity_score),
	gg_legend(gg_transcript_subtype_pro),
	gg_legend(gg_transcript_subtype_cla),
	gg_legend(gg_transcript_subtype_mes),
	gg_legend(gg_tmz),ncol=5)
dev.off()

#P-values vs subtype switch
test_dat <- plot_res %>% filter(status=="Recurrent",signature_name=="Proneural")
myvars <- c("grade_change","tmz","rt","aPD1")
p.val <- rep(0,length(myvars))
for(i in 1:length(myvars))
{
	g1 <- nrow(test_dat[which(test_dat[,"switch"]=="Switch" & test_dat[,myvars[i]]==1),])
	g2 <- nrow(test_dat[which(test_dat[,"switch"]=="Switch" & test_dat[,myvars[i]]==0),])
	g3 <- nrow(test_dat[which(test_dat[,"switch"]=="None" & test_dat[,myvars[i]]==1),])
	g4 <- nrow(test_dat[which(test_dat[,"switch"]=="None" & test_dat[,myvars[i]]==0),])
	
	ct <- matrix(c(g1,g2,g3,g4),nrow=2)
	p.val[i] <- fisher.test(ct)$p.value
}

#Test for mesenchymal gain after radiation
test_dat <- dat %>% filter(signature_name=="Mesenchymal")
mes_b <- grep("Mesenchymal",test_dat[,"subtype_b"])
mes_a <- grep("Mesenchymal",test_dat[,"subtype_a"])
myrows <- mes_b[which(!mes_b %in% mes_a)]
mes_gain <- rep(0,nrow(test_dat))
mes_gain[myrows] <- 1
test_dat[,"mes_gain"] <- mes_gain

#test_dat <- test_dat[grep("Proneural",test_dat[,"subtype_a"])]
myvars <- c("grade_change","treatment_tmz","treatment_radiotherapy","treatment_pd1")
p.val <- rep(0,length(myvars))
for(i in 1:length(myvars))
{
	g1 <- nrow(test_dat[which(test_dat[,"mes_gain"]==1 & test_dat[,myvars[i]]==1),])
	g2 <- nrow(test_dat[which(test_dat[,"mes_gain"]==1 & test_dat[,myvars[i]]==0),])
	g3 <- nrow(test_dat[which(test_dat[,"mes_gain"]==0 & test_dat[,myvars[i]]==1),])
	g4 <- nrow(test_dat[which(test_dat[,"mes_gain"]==0 & test_dat[,myvars[i]]==0),])
	
	ct <- matrix(c(g1,g2,g3,g4),nrow=2)
	p.val[i] <- fisher.test(ct)$p.value
}


#Boxplot of subtype switching
sub_dat <- dat[which(dat[,"idh_codel_subtype"]=="IDHmut-codel"),]
wilcox.test(sub_dat[which(sub_dat[,"switch"]=="Switch"),"ss_a"],sub_dat[which(sub_dat[,"switch"]=="None"),"ss_a"])

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/subtype_simplicity_boxplot.pdf",width=2.6,height=2)
p1 <- ggplot(dat, aes(x = idh_codel_subtype, y = ss_a, colour=switch)) +
geom_pointrange(position=position_dodge(width=0.5),stat="summary",
	fun.ymin=function(z) { quantile(z,0.25) },
	fun.ymax=function(z) { quantile(z,0.75) },
	fun.y=median) +
labs(y = "Simplicity score") +
scale_colour_manual(values=c("royalblue4","tomato3")) +
theme_bw() +
theme(axis.text= element_text(size=7),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") +
coord_cartesian(ylim=c(0.,1.1))
p1
dev.off()

