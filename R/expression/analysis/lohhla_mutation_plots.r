library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(ggbeeswarm)


#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
WITH collapse AS
(
	SELECT *,
	CASE WHEN pval < 0.05 THEN 'loss' ELSE 'none' END AS loss_status
	FROM variants.lohhla
	WHERE coverage_filter=10
),
agg AS
(
	SELECT pair_barcode, string_agg(loss_status,'; ') AS any_loss
	FROM collapse
	GROUP BY pair_barcode
)
SELECT ag.*, 
CASE WHEN any_loss LIKE '%loss%' THEN 1 ELSE 0 END AS loss_status,
mf.mt_count, mf.neoag_count, cs.idh_codel_subtype
FROM agg ag
JOIN analysis.pairs pa ON ag.pair_barcode = pa.pair_barcode
JOIN analysis.neoag_freq mf ON mf.aliquot_barcode = pa.tumor_barcode
JOIN analysis.diamond_set gs ON mf.aliquot_barcode = gs.tumor_barcode_a OR mf.aliquot_barcode = gs.tumor_barcode_b
JOIN clinical.subtypes cs ON cs.case_barcode = gs.case_barcode
ORDER BY 4 DESC
"

dat <- dbGetQuery(con, q)
dat[,"loss_status"] <- factor(dat[,"loss_status"],levels=c("0","1"))
dat[,"idh_codel_subtype"] <- as.factor(dat[,"idh_codel_subtype"])

pdf("/projects/varnf/GLASS-III/GLASS-III/figures/analysis/lohhla_mutation_beeswarm.pdf",width=3,height=4)
p1 <- ggplot(dat,aes(x = idh_codel_subtype, y = mt_count, colour=loss_status)) +
geom_beeswarm(size=0.3,dodge.width=0.5) +
stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", position=position_dodge(width=0.5),
	width=0.3, size = 0.25,colour="black",aes(group=loss_status)) +
scale_colour_manual(values=c("royalblue4","tomato3")) +
labs(y = "Total NS mutations") +
scale_y_continuous(trans='log10') +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")

p2 <- ggplot(dat,aes(x = idh_codel_subtype, y = neoag_count, colour=loss_status)) +
geom_beeswarm(size=0.3,dodge.width=0.5) +
stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", position=position_dodge(width=0.5),
	width=0.3, size = 0.25,colour="black",aes(group=loss_status)) +
scale_colour_manual(values=c("royalblue4","tomato3")) +
labs(y = "Total neoantigen-coding mutations") +
scale_y_continuous(trans='log10') +
theme_bw() +
theme(axis.text.x=element_text(size=7),axis.text.y = element_text(size=7),
axis.title.x = element_blank(),axis.title.y = element_text(size=7),
plot.title = element_blank(),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none")

grid.arrange(p1,p2,nrow=2)
dev.off()