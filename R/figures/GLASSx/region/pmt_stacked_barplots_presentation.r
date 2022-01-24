
library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(grid)
library(gridExtra)

#######################################################

rm(list=ls())

myinf1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/ivygap_scgp/CIBERSORTxGEP_ivygap_Fractions-Adjusted.txt"
myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/ivygap/columns-samples.csv"


dat <- read.delim(myinf1,header=TRUE)
samps <- read.csv(myinf2,stringsAsFactor=FALSE)


#Simplify analysis for now: Reference histology
samp_struct <- samps$structure_name
names(samp_struct) <- samps$rna_well_id

tumor_name <- samps$tumor_name
names(tumor_name) <- samps$rna_well_id

dat[,"structure"] <- samp_struct[as.character(dat[,"Mixture"])]
dat[,"tumor_name"] <- tumor_name[as.character(dat[,"Mixture"])]
dat <- dat[grep("reference histology", dat[,"structure"]),]

dat[,"structure"] <- gsub(" sampled by reference histology","",dat[,"structure"])


dat <- dat %>% 
		select(-c("P.value", "Correlation", "RMSE")) %>%
		pivot_longer(cols=-c("Mixture","structure","tumor_name"), names_to = "cell_state", values_to = "fraction") %>%
		mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell")) %>%
		mutate(structure = recode(structure, 
		"Cellular Tumor" = "Cellular tumor", "Infiltrating Tumor" = "Infiltrating tumor", "Leading Edge" = "Leading edge")) %>%
		mutate(Mixture = as.character(Mixture), tumor_name = as.character(tumor_name))
# Create average profile for plotting
		
dat <- dat %>%
	   group_by(structure, cell_state) %>%
	   summarise(fraction = mean(fraction)) %>%
	   mutate(tumor_name = "Mean", Mixture = paste("Mean",structure,sep=" ")) %>%
	   ungroup() %>%
	   select(Mixture, structure, tumor_name, cell_state, fraction) %>%
	   bind_rows(dat) %>%
	   mutate(Mixture = as_factor(Mixture)) %>%
	   mutate(structure = fct_relevel(structure, "Leading edge", "Infiltrating tumor", "Cellular tumor", "Pseudopalisading cells around necrosis", "Microvascular proliferation")) %>%
	   mutate(cell_state = as_factor(cell_state)) %>%
	   mutate(cell_state = fct_relevel(cell_state, "B cell", "Granulocyte", "T cell", "Dendritic cell", "Myeloid", 
									"Oligodendrocyte", 
									"Endothelial", "Pericyte",
									"Fibroblast", 
									"Diff.-like", "Stem-like", "Prolif. stem-like"))

mean_values <- dat %>% filter(tumor_name == "Mean")

#---------------------------------------------------------------------------


#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ag1.signature_name AS subtype_a, 
ag2.signature_name AS subtype_b,
ci1.cell_state AS cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_scgp ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
"

dat1 <- dbGetQuery(con,q)

dat1 <- dat1 %>%
mutate(cell_state = recode(cell_state, "b_cell" = "B cell", "dendritic_cell" = "Dendritic cell",
		"differentiated_tumor" = "Diff.-like", "endothelial" = "Endothelial",
		"fibroblast" = "Fibroblast", "granulocyte" = "Granulocyte",
		"myeloid" = "Myeloid", "oligodendrocyte" = "Oligodendrocyte",
		"pericyte" = "Pericyte", "prolif_stemcell_tumor" = "Prolif. stem-like",
		"stemcell_tumor" = "Stem-like","t_cell" = "T cell"))


# Figure 1: No adjustment:

unadj_res <- dat1 %>%
filter(idh_status == "IDHwt" & !(subtype_a == "Proneural" & subtype_b =="Classical")) %>%
filter(cell_state == "Stem-like" | cell_state =="Prolif. stem-like" | cell_state == "Diff.-like") %>%
group_by(tumor_pair_barcode, idh_status) %>%
summarise(cell_state, subtype_a, subtype_b, fraction_a = fraction_a/sum(fraction_a), fraction_b = fraction_b/sum(fraction_b)) %>%
ungroup() %>%
group_by(cell_state, subtype_a, subtype_b) %>%
summarise(p.val=t.test(fraction_a,fraction_b, paired=TRUE)$p.value, fraction_a = mean(fraction_a) , fraction_b = mean(fraction_b), n = n()) %>%
mutate(trans = paste(subtype_a,subtype_b,sep="_")) %>%
filter(subtype_a != subtype_b, n > 5)

fraction <- c(unadj_res$fraction_a, unadj_res$fraction_b)
trans <- rep(unadj_res$trans,2)
cell_state <- rep(unadj_res$cell_state, 2)
timepoint <- rep(c("Initial","Recurrent"), each = nrow(unadj_res))
plot_res <- data.frame(cell_state, fraction, trans, timepoint)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/malig_subtype_switch_bars_unadj.pdf",width=2.5,height=1.5)
ggplot(plot_res, aes(x=timepoint, y = fraction*100, fill = cell_state)) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
facet_grid(.~trans) +
labs(y = "Fraction (%)") +
scale_fill_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
theme(axis.text.x = element_text(size=7,angle=45,hjust=1),
	axis.text.y = element_text(size=7),
	axis.title.y = element_text(size=7),
	axis.title.x = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 
dev.off()		  

#---------------------------------------------------------------------------


#Read in data for analysis
q <- "
SELECT ss.tumor_pair_barcode, 
ss.case_barcode,
ss.tumor_barcode_a,
ss.tumor_barcode_b,
ag1.signature_name AS subtype_a, 
ag2.signature_name AS subtype_b,
ci1.cell_state AS cell_state,
ci1.fraction AS fraction_a,
ci2.fraction AS fraction_b,
cs.idh_codel_subtype,
CASE WHEN cs.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.rna_silver_set ss
JOIN analysis.cibersortx_ivygap ci1 ON ci1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_ivygap ci2 ON ci2.aliquot_barcode = ss.tumor_barcode_b AND ci2.cell_state = ci1.cell_state
JOIN analysis.top_transcriptional_subtype ag1 ON ag1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.top_transcriptional_subtype ag2 ON ag2.aliquot_barcode = ss.tumor_barcode_b 
JOIN clinical.subtypes cs ON cs.case_barcode = ss.case_barcode
"

dat2 <- dbGetQuery(con,q)


# Adjust for all non-CT  features:
# Leading edge
le_fractions <- mean_values %>% filter(structure == "Leading edge") %>% data.frame()

new_fraction_a <- new_fraction_b <- rep(0, nrow(dat1))
for(i in 1:nrow(dat1))
{
	cell_state <- dat1[i,"cell_state"]
	
	# Representation of a cell state in the leading edge
	region_factor <- le_fractions[which(le_fractions[,"cell_state"] == cell_state), "fraction"]
	
	# Representation of the leading edge in a sample
	region_fraction_a <- dat2[which(dat2$tumor_pair_barcode == dat1[i,"tumor_pair_barcode"] & dat2$cell_state == "LE"), "fraction_a"]
	region_fraction_b <- dat2[which(dat2$tumor_pair_barcode == dat1[i,"tumor_pair_barcode"] & dat2$cell_state == "LE"), "fraction_b"]
	
	# Cell state fraction normalized for no leading edge
	new_fraction_a[i] = dat1[i,"fraction_a"] - region_factor * region_fraction_a
	new_fraction_b[i] = dat1[i,"fraction_b"] - region_factor * region_fraction_b
	
}

new_dat1 <- data.frame(dat1, new_fraction_a, new_fraction_b)

# Microvascular proliferation

le_fractions <- mean_values %>% filter(structure == "Microvascular proliferation") %>% data.frame()

new_fraction_a <- new_fraction_b <- rep(0, nrow(dat1))
for(i in 1:nrow(dat1))
{
	cell_state <- dat1[i,"cell_state"]
	
	# Representation of a cell state in the leading edge
	region_factor <- le_fractions[which(le_fractions[,"cell_state"] == cell_state), "fraction"]
	
	# Representation of the leading edge in a sample
	region_fraction_a <- dat2[which(dat2$tumor_pair_barcode == dat1[i,"tumor_pair_barcode"] & dat2$cell_state == "CTmvp"), "fraction_a"]
	region_fraction_b <- dat2[which(dat2$tumor_pair_barcode == dat1[i,"tumor_pair_barcode"] & dat2$cell_state == "CTmvp"), "fraction_b"]
	
	# Cell state fraction normalized for no leading edge
	new_fraction_a[i] = new_dat1[i,"new_fraction_a"] - region_factor * region_fraction_a
	new_fraction_b[i] = new_dat1[i,"new_fraction_b"] - region_factor * region_fraction_b
	
}
new_dat1 <- data.frame(dat1, new_fraction_a, new_fraction_b)

# Pseudopalisading cells around necrosis

le_fractions <- mean_values %>% filter(structure == "Pseudopalisading cells around necrosis") %>% data.frame()

new_fraction_a <- new_fraction_b <- rep(0, nrow(dat1))
for(i in 1:nrow(dat1))
{
	cell_state <- dat1[i,"cell_state"]
	
	# Representation of a cell state in the leading edge
	region_factor <- le_fractions[which(le_fractions[,"cell_state"] == cell_state), "fraction"]
	
	# Representation of the leading edge in a sample
	region_fraction_a <- dat2[which(dat2$tumor_pair_barcode == dat1[i,"tumor_pair_barcode"] & dat2$cell_state == "CTpan"), "fraction_a"]
	region_fraction_b <- dat2[which(dat2$tumor_pair_barcode == dat1[i,"tumor_pair_barcode"] & dat2$cell_state == "CTpan"), "fraction_b"]
	
	# Cell state fraction normalized for no leading edge
	new_fraction_a[i] = new_dat1[i,"new_fraction_a"] - region_factor * region_fraction_a
	new_fraction_b[i] = new_dat1[i,"new_fraction_b"] - region_factor * region_fraction_b
	
}


new_dat1 <- data.frame(dat1, new_fraction_a, new_fraction_b)
new_dat1[which(new_dat1$new_fraction_a < 0), "new_fraction_a"] <- 0
new_dat1[which(new_dat1$new_fraction_b < 0), "new_fraction_b"] <- 0


# Figure 4: Post-adjustment all non-CT features:

adj_res <- new_dat1 %>%
filter(idh_status == "IDHwt" & !(subtype_a == "Proneural" & subtype_b =="Classical")) %>%
filter(cell_state == "Stem-like" | cell_state =="Prolif. stem-like" | cell_state == "Diff.-like") %>%
group_by(tumor_pair_barcode, idh_status) %>%
summarise(cell_state, subtype_a, subtype_b, new_fraction_a = new_fraction_a/sum(new_fraction_a), new_fraction_b = new_fraction_b/sum(new_fraction_b)) %>%
ungroup() %>%
group_by(cell_state, subtype_a, subtype_b) %>%
summarise(fraction_a = mean(new_fraction_a,na.rm=TRUE) , fraction_b = mean(new_fraction_b,na.rm=TRUE), n = n()) %>%
mutate(trans = paste(subtype_a,subtype_b,sep="_")) %>%
filter(subtype_a != subtype_b, n > 5)

fraction <- c(adj_res$fraction_a, adj_res$fraction_b)
trans <- rep(adj_res$trans,2)
cell_state <- rep(adj_res$cell_state, 2)
timepoint <- rep(c("Initial","Recurrent"), each = nrow(adj_res))
plot_res <- data.frame(cell_state, fraction, trans, timepoint)

pdf("/projects/verhaak-lab/GLASS-III/figures/analysis/malig_subtype_switch_bars_adj.pdf",width=2.5,height=1.5)
ggplot(plot_res, aes(x=timepoint, y = fraction*100, fill = cell_state)) +
geom_bar(position="stack", stat="identity") +
theme_classic() +
facet_grid(.~trans) +
labs(y = "Fraction (%)") +
scale_fill_manual(values=c("Stem-like" = "#fb6a4a", "Diff.-like" = "#fcbba1", "Prolif. stem-like" = "#a50f15")) +
theme(axis.text.x = element_text(size=7,angle=45,hjust=1),
	axis.text.y = element_text(size=7),
	axis.title.y = element_text(size=7),
	axis.title.x = element_blank(),
	panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	strip.text = element_text(size=7),
	strip.background = element_rect(colour="white",fill="white"),
	legend.position = "none") 
dev.off()		  
