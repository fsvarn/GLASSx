###################################################
# Create stacked barplots (transcriptional classifier/simplicity score/CIBERSORTx) of each GLASS sample
# Updated: 2020.07.06
# Author: Frederick Varn
##################################################
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
	--WHERE ts1.p_value < 0.05 OR ts2.p_value < 0.05
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
ci1.cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
su1.idh_codel_subtype
FROM analysis.rna_silver_set ss
JOIN analysis.transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b AND ts1.signature_name = ts2.signature_name
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
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

# Correlate CIBERSORTx scores with similarity score (for visual purposes)

cell_states <- unique(dat[,"cell_state"])
init_cscor <- rec_cscor <- rep(0, length(cell_states))
names(init_cscor) <- names(rec_cscor) <- cell_states
for(i in 1:length(cell_states))
{
	mystate <- cell_states[i]
	substate <- dat[which(dat[,"cell_state"] == mystate),c("tumor_pair_barcode","ss_a","ss_b","cell_state","fraction_a","fraction_b")]
	substate <- distinct(substate)
		
	init_cscor[i] <- cor(substate[,"fraction_a"], substate[,"ss_a"],method="s")
	rec_cscor[i] <- cor(substate[,"fraction_b"], substate[,"ss_b"],method="s")

}

# Build CIBERSORTx scores with similarity score (for visual purposes)

plot_res1 <- dat[,c("case_barcode","tumor_barcode_a","signature_name","pval_a","ss_a","switch","grade_change","treatment_tmz","treatment_radiotherapy","treatment_pd1","idh_codel_subtype","cell_state","fraction_a")]
plot_res2 <- dat[,c("case_barcode","tumor_barcode_b","signature_name","pval_b","ss_b","switch","grade_change","treatment_tmz","treatment_radiotherapy","treatment_pd1","idh_codel_subtype","cell_state","fraction_b")]
status <- c(rep("Initial",nrow(plot_res1)),rep("Recurrent",nrow(plot_res2)))
colnames(plot_res1) <- colnames(plot_res2) <- c("case_barcode","aliquot_barcode","signature_name","pval","ss","switch","grade_change","tmz","rt","aPD1","idh_codel_subtype","cell_state", "proportion")
plot_res <- rbind(plot_res1,plot_res2)
plot_res <- cbind(plot_res,status)

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

#Rename the cell states
plot_res <- plot_res %>%
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
  ggplot(plot_res %>% filter(signature_name=="Mesenchymal", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, alpha= -1 * ss)) +
  scale_fill_manual(values=c("black")) +
  theme(axis.title.y=element_blank())
  
gg_transcript_subtype_pro <-
  ggplot(plot_res %>% filter(signature_name=="Proneural", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#00458a")) +
  theme(axis.title.y=element_blank())
  
gg_transcript_subtype_cla <-
  ggplot(plot_res %>% filter(signature_name=="Classical", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#008a22")) +
  theme(axis.title.y=element_blank())
  
gg_transcript_subtype_mes <-
  ggplot(plot_res %>% filter(signature_name=="Mesenchymal", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=status, fill=signature_name, alpha= pval)) +
  scale_fill_manual(values=c("#8a0000")) +
  theme(axis.title.y=element_blank())
  
gg_cell_state_init <- 
  ggplot(plot_res %>% filter(status=="Initial", signature_name=="Mesenchymal"), aes(x=case_barcode, y = proportion, fill = cell_state)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
						 "Fibroblast" = "#feb24c",
						 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
  theme(axis.title.y=element_blank())

gg_cell_state_rec <- 
  ggplot(plot_res %>% filter(status=="Recurrent", signature_name=="Mesenchymal"), aes(x=case_barcode, y = proportion, fill = cell_state)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
						 "Fibroblast" = "#feb24c",
						 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
  theme(axis.title.y=element_blank())
   
gg_subtype_switch <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=switch)) +
  scale_fill_manual(values=c("white","black")) +
  theme(axis.title.y=element_blank())
  
gg_grade_change <- 
  ggplot(plot_res %>% filter(status=="Recurrent",signature_name=="Proneural", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=grade_change)) +
  scale_fill_manual(values=c("white","black")) +
  theme(axis.title.y=element_blank())
  
gg_tmz <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=tmz)) +
  scale_fill_manual(values=c("white","black"),na.value="grey50") +
  theme(axis.title.y=element_blank())
  
gg_rt <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=rt)) +
  scale_fill_manual(values=c("white","black"),na.value="grey50") +
  theme(axis.title.y=element_blank())
  
gg_pd1 <- 
  ggplot(plot_res %>% filter(status=="Initial",signature_name=="Proneural", cell_state == "Stem cell tumor"), aes(x=case_barcode)) +
  geom_tile(aes(y=signature_name, fill=aPD1)) +
  scale_fill_manual(values=c("white","black"),na.value="grey50") +
  theme(axis.title.y=element_blank())
  
#Align figures for printing
gb1 <- ggplot_build(nolabels(gg_simplicity_score))
gb2 <- ggplot_build(nolabels(gg_transcript_subtype_pro))
gb3 <- ggplot_build(nolabels(gg_transcript_subtype_cla))
gb4 <- ggplot_build(nolabels(gg_transcript_subtype_mes))
gb5 <- ggplot_build(nolabels(gg_cell_state_init))
gb6 <- ggplot_build(nolabels(gg_cell_state_rec))
gb7 <- ggplot_build(nolabels(gg_subtype_switch))
gb8 <- ggplot_build(nolabels(gg_grade_change))
gb9 <- ggplot_build(nolabels(gg_tmz))
gb10 <- ggplot_build(nolabels(gg_rt))
gb11 <- ggplot_build(nolabels(gg_pd1))

n1 <- length(gb1$layout$panel_scales_y[[1]]$range$range)
n2 <- length(gb2$layout$panel_scales_y[[1]]$range$range)
n3 <- length(gb3$layout$panel_scales_y[[1]]$range$range)
n4 <- length(gb4$layout$panel_scales_y[[1]]$range$range)
n5 <- length(gb5$layout$panel_scales_y[[1]]$range$range)
n6 <- length(gb6$layout$panel_scales_y[[1]]$range$range)
n7 <- length(gb7$layout$panel_scales_y[[1]]$range$range)
n8 <- length(gb8$layout$panel_scales_y[[1]]$range$range)
n9 <- length(gb9$layout$panel_scales_y[[1]]$range$range)
n10 <- length(gb10$layout$panel_scales_y[[1]]$range$range)
n11 <- length(gb11$layout$panel_scales_y[[1]]$range$range)

gA <- ggplot_gtable(gb1)
gB <- ggplot_gtable(gb2)
gC <- ggplot_gtable(gb3)
gD <- ggplot_gtable(gb4)
gE <- ggplot_gtable(gb5)
gF <- ggplot_gtable(gb6)
gG <- ggplot_gtable(gb7)
gH <- ggplot_gtable(gb8)
gI <- ggplot_gtable(gb9)
gJ <- ggplot_gtable(gb10)
gK <- ggplot_gtable(gb11)

g <- gtable:::rbind_gtable(gA, gB, "last")
g <- gtable:::rbind_gtable(g, gC, "last")
g <- gtable:::rbind_gtable(g, gD, "last")
g <- gtable:::rbind_gtable(g, gE, "last")
g <- gtable:::rbind_gtable(g, gF, "last")
g <- gtable:::rbind_gtable(g, gG, "last")
g <- gtable:::rbind_gtable(g, gH, "last")
g <- gtable:::rbind_gtable(g, gI, "last")
g <- gtable:::rbind_gtable(g, gJ, "last")
g <- gtable:::rbind_gtable(g, gK, "last")

panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels[1*3]] <- unit(n1, "null")
g$heights[panels[2*3]] <- unit(n1, "null")
g$heights[panels[3*3]] <- unit(n1,"null")
g$heights[panels[4*3]] <- unit(n1,"null")
g$heights[panels[5*3]] <- unit(n1*3.2,"null")
g$heights[panels[6*3]] <- unit(n1*3.2,"null")
g$heights[panels[7*3]] <- unit(n1/2,"null")
g$heights[panels[8*3]] <- unit(n1/2,"null")
g$heights[panels[9*3]] <- unit(n1/2,"null")
g$heights[panels[10*3]] <- unit(n1/2,"null")
g$heights[panels[11*3]] <- unit(n1/2,"null")

grid.newpage()

#Plot
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/transcriptional_subtype_change_v4.pdf",width=7,height=5)
grid.draw(g)
dev.off()

#Legends
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/transcriptional_subtype_change_legends_v4.pdf",width=7,height=7)
grid.arrange(gg_legend(gg_simplicity_score),
	gg_legend(gg_transcript_subtype_pro),
	gg_legend(gg_transcript_subtype_cla),
	gg_legend(gg_transcript_subtype_mes),
	gg_legend(gg_cell_state_rec),
	gg_legend(gg_tmz),ncol=5)
dev.off()


# Initial versus recurrence

barplot_res <- plot_res %>% 
			   group_by(cell_state, status, idh_codel_subtype) %>% 
			   summarise(proportion = mean(proportion)) %>%
			   mutate(status = as_factor(status)) %>%
		       mutate(status = fct_relevel(status, "Initial","Recurrent"))
			   

mysubtype <- unique(plot_res[,"idh_codel_subtype"])   
mycellstate <- unique(plot_res[,"cell_state"])
p.val <- matrix(0, nrow=length(mycellstate), ncol=length(mysubtype))
rownames(p.val) <- mycellstate
colnames(p.val) <- mysubtype
for(i in 1:length(mycellstate))
{
	for(j in 1:length(mysubtype))
	{
		dat1 <- plot_res %>% filter(cell_state == mycellstate[i] & idh_codel_subtype == mysubtype[j] & status == "Initial")
		dat2 <- plot_res %>% filter(cell_state == mycellstate[i] & idh_codel_subtype == mysubtype[j] & status == "Recurrent")
		p.val[i,j] <- wilcox.test(dat1[,"proportion"], dat2[,"proportion"], paired=TRUE)$p.value
	}
}
			   
pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/cibersortx_stacked_barplot_timepoint_subtype.pdf",width=5,height=3.5)  
ggplot(barplot_res, aes(fill=cell_state, y=proportion, x=status)) + 
geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
						 "Oligodendrocyte" = "#2ca25f",
						 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
						 "Fibroblast" = "#feb24c",
						 "Stem cell tumor" = "#fb6a4a", "Differentiated tumor" = "#fcbba1", "Proliferating stem cell tumor" = "#a50f15")) +
labs(y = "Proportion") +
facet_grid(~idh_codel_subtype) +
theme_bw() +
theme(axis.title.x=element_blank(),
axis.title.y=element_text(size=7),
axis.text.x=element_text(size=7,angle=45,hjust=1),
axis.text.y=element_text(size=7),
strip.text = element_text(size=7, hjust=0.5),
strip.background = element_blank(),
legend.title = element_blank(),
legend.text = element_text(size=7),
legend.position = "right")
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

