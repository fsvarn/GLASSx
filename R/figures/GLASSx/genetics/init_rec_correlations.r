###################################################
# Compare cell state fractions between initial and recurrent tumors
# Updated: 2020.11.19
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(survival)
library(DescTools)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
SELECT ss.*, ci1.cell_state,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
"

dat <- dbGetQuery(con, q)

dat <- dat %>% mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		mutate(cell_state = fct_relevel(cell_state, "Prolif. stem-like", "Stem-like","Diff.-like",
										"Fibroblast", "Pericyte","Endothelial", "Oligodendrocyte",
										"Myeloid", "Dendritic cell", "T cell", "Granulocyte","B cell"))
		
tumor_dat <- dat %>% 
		filter(grepl("-like", cell_state)) %>%
		group_by(case_barcode) %>%
		summarise(cell_state, fraction_a = (fraction_a/sum(fraction_a))*100, fraction_b = (fraction_b/sum(fraction_b))*100, residual = fraction_b - fraction_a, idh_status) %>%
	   	mutate(cell_state = fct_relevel(cell_state, c("Prolif. stem-like", "Diff.-like" ,"Stem-like"))) %>%
		as.data.frame()
	
		
# Get correlation coefficients

# res <- dat %>% group_by(cell_state, idh_status) %>% 
# 		summarise(cor = cor(fraction_a,fraction_b,method="p"), p.value=cor.test(fraction_a, fraction_b)$p.value, n = n()) %>%
# 		as.data.frame()

res <- dat %>% group_by(cell_state, idh_status) %>% 
		summarise(cor = CCC(fraction_a,fraction_b)$rho.c$est, n = n()) %>%
		as.data.frame()

res2 <- dat %>% group_by(cell_state, idh_status) %>% 
		summarise(pval = cor.test(fraction_a,fraction_b)$p.value, n = n()) %>%
		as.data.frame()
		
avg <- res %>% 
	   group_by(cell_state) %>%
	   summarise(avg = mean(cor)) %>%
	   arrange(desc(avg))

# res <- res %>%
# 	   mutate(cell_state = fct_relevel(cell_state, avg$cell_state))

# tumor_res <- tumor_dat %>% group_by(cell_state, idh_status) %>% 
# 		summarise(cor = cor(fraction_a,fraction_b,method="p"), p.value=cor.test(fraction_a, fraction_b)$p.value, n = n()) %>%
# 		as.data.frame()

tumor_res <- tumor_dat %>% group_by(cell_state, idh_status) %>% 
		summarise(cor = CCC(fraction_a,fraction_b)$rho.c$est, n = n()) %>%
		as.data.frame()

# Residuals
#test <- tumor_dat %>% filter(cell_state == "Prolif. stem-like", idh_status == "IDHmut")
# test <- tumor_dat %>% filter(cell_state == "Prolif. stem-like", idh_status == "IDHmut")
# X <- cbind(rep(1, nrow(test)), test$fraction_a)
# H = X %*% solve(t(X) %*% X) %*% t(X)
# lev <- diag(H)
# MSE <- mean(test$residual^2)
# 
# test$residual/(sqrt(MSE * 1 - lev))
# 
# 
# ks.test(x=test$residual, y="pnorm", mean = mean(test$residual), sd = sd(test$residual))
# 
# lm_obj <- lm(fraction_b~0+fraction_a, test)
# 		as.data.frame()
# 		
# Heatmap of coefficients

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_cor_hm.pdf", width=2.75, height = 1.15) #,width=2.7,height=3)
ggplot(res, aes(x = cell_state, y = idh_status, fill=cor)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") 
dev.off()
		
		
# Legends
dummy_res <- res
dummy_res[1,3] <- 1
dummy_res[2,3] <- -1
	
pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_cor_hm_legend.pdf", width=2.75, height = 1.15) #,width=2.7,height=3)
ggplot(dummy_res, aes(x = cell_state, y = idh_status, fill=cor)) +
geom_tile() +
scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title = element_blank(),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white")) 
dev.off()
		
# Barplot of coefficients

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_cor_bar.pdf", width=5.6, height = 2.4) #,width=2.7,height=3)
ggplot(res, aes(x = cell_state, y = cor, fill=cell_state)) +
geom_bar(stat="identity") +
labs(y = "Lin's concordance") +
scale_fill_manual(values=c("B cell" = "#eff3ff", "Granulocyte" = "#bdd7e7", "T cell" = "#6baed6", "Dendritic cell" = "#3182bd", "Myeloid" = "#08519c",
					 "Oligodendrocyte" = "#2ca25f",
					 "Endothelial" = "#ffffd4", "Pericyte" = "#fee391",
					 "Fibroblast" = "#feb24c",
					 "Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
facet_wrap(vars(idh_status), nrow=1, strip.position="top") +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") 
dev.off()

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_cor_bar_grouped.pdf",width=2.7,height=1.6)
ggplot(res, aes(x = cell_state, y = cor, fill=idh_status)) +
geom_bar(position="dodge",stat="identity") +
labs(y = "Pearson coefficient (R)") +
scale_fill_manual(values=c("IDHwt" = "#619CFF", "IDHmut" = "#00BA38")) +
#facet_wrap(vars(idh_status), nrow=2, strip.position="top") +
theme_bw() +
theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
axis.title.x = element_blank(), axis.title.y = element_text(size=7),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") 
dev.off()

# Scatterplots of tumor cell fraction

idhmut_cor <- tumor_res %>% filter(idh_status == "IDHmut")
idhmut_cor$fraction_a = 70
idhmut_cor$fraction_b = 5
idhmut_cor$label = sapply(round(idhmut_cor$cor,2), function(x)deparse(bquote(italic("R")~" = "~.(x))))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_cor_idhmut_scatterplot.pdf",width=3,height=1.4)
ggplot(tumor_dat %>% filter(idh_status == "IDHmut"),aes(x = fraction_a, y = fraction_b)) +
geom_point(size=0.5) +
geom_abline(aes(slope = 1, intercept = 0, colour=cell_state)) +
facet_wrap(vars(cell_state), nrow=1, strip.position="top") +
labs(x = "Initial tumor proportion (%)", y = "Recurrent tumor proportion (%)") +
geom_text(data = idhmut_cor, aes(label = label),size=2.5, parse=TRUE) +
# annotate(geom="text", size=2.5, 
#  	x=0.05, y = 0.95, label=deparse(bquote(italic("R")~" = "~.(cor))),parse=TRUE) +
scale_colour_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
theme_bw() +
theme(
axis.text.x= element_text(size=7),
axis.text.y= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") + 
coord_cartesian(xlim=c(0,100), ylim=c(0,100))
dev.off()



idhwt_cor <- tumor_res %>% filter(idh_status == "IDHwt")
idhwt_cor$fraction_a = 70
idhwt_cor$fraction_b = 5
idhwt_cor$label = sapply(round(idhwt_cor$cor,2), function(x)deparse(bquote(italic("R")~" = "~.(x))))


pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/cell_state_cor_idhwt_scatterplot.pdf",width=3,height=1.4)
ggplot(tumor_dat %>% filter(idh_status == "IDHwt"),aes(x = fraction_a, y = fraction_b)) +
geom_point(size=0.5) +
geom_abline(aes(slope = 1, intercept = 0, colour=cell_state)) +
facet_wrap(vars(cell_state), nrow=1, strip.position="top") +
labs(x = "Initial tumor proportion (%)", y = "Recurrent tumor proportion (%)") +
geom_text(data = idhwt_cor, aes(label = label),size=2.5, parse=TRUE) +
# annotate(geom="text", size=2.5, 
#  	x=0.05, y = 0.95, label=deparse(bquote(italic("R")~" = "~.(cor))),parse=TRUE) +
scale_colour_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
theme_bw() +
theme(
axis.text.x= element_text(size=7),
axis.text.y= element_text(size=7),
axis.title = element_text(size=7),
plot.title = element_text(size=7,hjust=0.5),
panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
strip.text = element_text(size=7),
strip.background = element_rect(colour="white",fill="white"),
legend.position="none") + 
coord_cartesian(xlim=c(0,100), ylim=c(0,100))
dev.off()


			
			