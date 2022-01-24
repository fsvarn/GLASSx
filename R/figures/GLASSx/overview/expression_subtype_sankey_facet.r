###################################################
# Subtype switching Sankey plot, faceted by IDH subtype (used in Fig. 1)
# Updated: 2020.08.24
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)
library(ggalluvial)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
WITH full_data AS
(
	SELECT ss.tumor_pair_barcode, 
	ss.case_barcode,
	ts1.signature_name AS subtype_a, 
	ts2.signature_name AS subtype_b,
	CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
	FROM analysis.rna_silver_set ss
	JOIN analysis.top_transcriptional_subtype ts1 ON ts1.aliquot_barcode = ss.tumor_barcode_a
	JOIN analysis.top_transcriptional_subtype ts2 ON ts2.aliquot_barcode = ss.tumor_barcode_b 
	JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
)
SELECT subtype_a, subtype_b, idh_status, COUNT(*)
FROM full_data
GROUP BY subtype_a, subtype_b, idh_status
"

dat <- dbGetQuery(con,q)
dat[,"count"] <- as.numeric(dat[,"count"])

# Fisher's exact tests

c1 <- sum(dat %>% filter(subtype_b == "Classical") %>% .$count)
c2 <- sum(dat %>% filter(subtype_b != "Classical") %>% .$count)
c3 <- sum(dat %>% filter(subtype_a == "Classical") %>% .$count)
c4 <- sum(dat %>% filter(subtype_a != "Classical") %>% .$count)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))

c1 <- sum(dat %>% filter(subtype_a != "Classical" & subtype_a == subtype_b) %>% .$count)
c2 <- sum(dat %>% filter(subtype_a != "Classical" & subtype_a != subtype_b) %>% .$count)
c3 <- sum(dat %>% filter(subtype_a == "Classical" & subtype_a == subtype_b) %>% .$count)
c4 <- sum(dat %>% filter(subtype_a == "Classical" & subtype_a != subtype_b) %>% .$count)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))


wtdat <- dat %>% filter(idh_status == "IDHwt")
mutdat <- dat %>% filter(idh_status == "IDHmut")

c1 <- sum(wtdat %>% filter(subtype_b == "Mesenchymal") %>% .$count)
c2 <- sum(wtdat %>% filter(subtype_b != "Mesenchymal") %>% .$count)
c3 <- sum(wtdat %>% filter(subtype_a == "Mesenchymal") %>% .$count)
c4 <- sum(wtdat %>% filter(subtype_a != "Mesenchymal") %>% .$count)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))

c1 <- sum(wtdat %>% filter(subtype_b == "Classical") %>% .$count)
c2 <- sum(wtdat %>% filter(subtype_b != "Classical") %>% .$count)
c3 <- sum(wtdat %>% filter(subtype_a == "Classical") %>% .$count)
c4 <- sum(wtdat %>% filter(subtype_a != "Classical") %>% .$count)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))

c1 <- sum(wtdat %>% filter(subtype_a != "Mesenchymal" & subtype_a == subtype_b) %>% .$count)
c2 <- sum(wtdat %>% filter(subtype_a != "Mesenchymal" & subtype_a != subtype_b) %>% .$count)
c3 <- sum(wtdat %>% filter(subtype_a == "Mesenchymal" & subtype_a == subtype_b) %>% .$count)
c4 <- sum(wtdat %>% filter(subtype_a == "Mesenchymal" & subtype_a != subtype_b) %>% .$count)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))

c1 <- sum(wtdat %>% filter(subtype_a != "Classical" & subtype_a == subtype_b) %>% .$count)
c2 <- sum(wtdat %>% filter(subtype_a != "Classical" & subtype_a != subtype_b) %>% .$count)
c3 <- sum(wtdat %>% filter(subtype_a == "Classical" & subtype_a == subtype_b) %>% .$count)
c4 <- sum(wtdat %>% filter(subtype_a == "Classical" & subtype_a != subtype_b) %>% .$count)

fisher.test(matrix(c(c1,c2,c3,c4),nrow=2))

# Check number switching
sum(wtdat %>% filter(subtype_a != subtype_b) %>% .$count)/sum(wtdat$count)
sum(mutdat %>% filter(subtype_a == "Proneural", subtype_b == "Proneural") %>% .$count)/sum(mutdat$count)
sum(dat %>% filter(subtype_a != subtype_b) %>% .$count)/sum(dat$count)


dat <- dat %>% 
	   mutate(subtype_a = recode(subtype_a, "Mesenchymal" = "Mes.", "Proneural" = "Pro.", "Classical" = "Class.")) %>%
	   mutate(subtype_b = recode(subtype_b, "Mesenchymal" = "Mes.", "Proneural" = "Pro.", "Classical" = "Class."))

p1 <- ggplot(dat %>% filter(idh_status == 'IDHwt'), aes(axis1 = subtype_a, axis2 = subtype_b, y=count)) +
geom_alluvium(aes(fill = subtype_a), width = 1/12, knot.pos = 0.4) +
geom_stratum(width = 1/4, fill = "white", color = "grey") +
scale_x_discrete(limits = c("Initial", "Recurrent"),expand = c(0,0)) +
scale_y_continuous(expand = c(0,0)) +
scale_fill_manual(values=c("#008A22","#8A0000","#00458A")) +
ggtitle("IDHwt") +
geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
theme_bw() +
theme(axis.text.x= element_text(size=7), axis.text.y=element_blank(),
axis.ticks = element_blank(),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust = 0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 

p2 <- ggplot(dat %>% filter(idh_status == 'IDHmut'), aes(axis1 = subtype_a, axis2 = subtype_b, y=count)) +
geom_alluvium(aes(fill = subtype_a), width = 1/12, knot.pos = 0.4) +
geom_stratum(width = 1/4, fill = "white", color = "grey") +
scale_x_discrete(limits = c("Initial", "Recurrent"),expand = c(0,0)) +
scale_y_continuous(expand = c(0,0)) +
scale_fill_manual(values=c("#008A22","#8A0000","#00458A")) +
ggtitle("IDHmut") +
geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.8) +
theme_bw() +
theme(axis.text.x= element_text(size=7), axis.text.y=element_blank(),
axis.ticks = element_blank(),
axis.title = element_blank(),
plot.title = element_text(size=7, hjust = 0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
legend.position="none") 

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/subtype_sankey_plot_facet.pdf",width=3.5,height=1.6)
grid.arrange(p1, p2, ncol=2)
dev.off()
