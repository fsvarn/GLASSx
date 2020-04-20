###################################################
# How does each transcriptional simplicity associate with cancer cell fraction of drivers?
# Updated: 2020.03.31
# Author: Frederick Varn
##################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(reshape)
library(grid)


##################################################
rm(list=ls())

# Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in data
q <- 
"
WITH subtype_rank AS
(
	SELECT *,
	RANK() OVER (PARTITION BY aliquot_barcode ORDER BY p_value ASC) AS p_rank
	FROM analysis.transcriptional_subtype
),
top_rank AS
(
	SELECT *
	FROM subtype_rank
	WHERE p_rank = 1
),
agg AS
(
	SELECT aliquot_barcode, 
	string_agg(signature_name,',') AS subtype 
	FROM top_rank
	GROUP BY aliquot_barcode
)
SELECT da.*, ag.subtype, ss.simplicity_score, es.purity
FROM analysis.drivers_by_aliquot da
JOIN analysis.analyte_sets an ON an.dna_barcode = da.aliquot_barcode
JOIN analysis.platinum_set ps ON ps.dna_barcode_a = da.aliquot_barcode
JOIN agg ag ON ag.aliquot_barcode = an.rna_barcode 
JOIN analysis.estimate es ON es.aliquot_barcode = an.rna_barcode
JOIN analysis.simplicity_score ss ON ss.aliquot_barcode = an.rna_barcode
WHERE idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con,q)

subtypes <- c("Mesenchymal","Classical","Proneural")

###################################################
# Step 1: Examine how simplicity associates with subtype
##################################################

p.val <- eff <- rep(0, length(subtypes))
names(p.val) <- names(eff) <- subtypes
for(i in 1:length(subtypes))
{
	g1 <- dat[grep(subtypes[i], dat[,"subtype"]),"simplicity_score"]
	g2 <- dat[grep(subtypes[i], dat[,"subtype"], invert=TRUE),"simplicity_score"]
	
	p.val[i] <- wilcox.test(g1, g2)$p.value
	eff[i] <- median(g1)
}

###################################################
# Step 2: Examine how simplicity associates with subtype
##################################################

purity_cor <- cor.test(dat[,"purity"], dat[,"simplicity_score"], method="s")

###################################################
# Step 3: Identify the driver alterations that associate with subtype
##################################################

# Start with copy number:
#---------------------------

# Get list of unique driver alterations
cnv_drivers <- unique(unlist(strsplit(dat[,"cnv_driver"],", ")))
cnv_drivers <- unique(cnv_drivers)
cnv_drivers <- cnv_drivers[-which(is.na(cnv_drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
cnv_genes <- cnv_drivers

# Create results matrix
cnv_results <- matrix(NA, nrow=length(cnv_genes), ncol=length(subtypes))
rownames(cnv_results) <- cnv_genes
colnames(cnv_results) <- subtypes

sig_cnv <- c()
sig_cnv_subtype <- c()
for(i in 1:length(cnv_genes))
{
	#Create a driver_status vector 
	driver_status <- rep(0, nrow(dat))
	
	# Identify all samples that have a mutation in this gene
	indices <- grep(cnv_genes[i], dat[,"cnv_driver"])	
	
	# Change driver_status vector values to indicate driver is present at initial tumor
	driver_status[indices] <- 1
	if(sum(driver_status)==0){
		next}
	test_dat <- cbind(dat, driver_status)
	
	for(j in 1:length(subtypes))
	{
		g1 <- nrow(test_dat[which(test_dat[,"driver_status"]==1 & test_dat[,"subtype"]==subtypes[j]),])
		g2 <- nrow(test_dat[which(test_dat[,"driver_status"]==0 & test_dat[,"subtype"]==subtypes[j]),])
		g3 <- nrow(test_dat[which(test_dat[,"driver_status"]==1 & test_dat[,"subtype"]!=subtypes[j]),])
		g4 <- nrow(test_dat[which(test_dat[,"driver_status"]==0 & test_dat[,"subtype"]!=subtypes[j]),])
		
		p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value
		eff <- g1/(g1+g2) > g3/(g3+g4)
		p.val <- ifelse(eff, p.val, p.val*-1)
		cnv_results[i,j] <- p.val
		
		# Identify all CNVs that are associated with a specific subtype
		# Ignoring CNVs that are significantly depleted, as this is likely due to enrichment in other subtypes
		if(p.val < 0.1 & p.val > 0 & sum(driver_status) > 5){
		sig_cnv <- c(sig_cnv, cnv_genes[i])
		sig_cnv_subtype <- c(sig_cnv_subtype, subtypes[j])
		}
	}
	
}

sig_cnv_drivers <- data.frame(sig_cnv, sig_cnv_subtype)
colnames(sig_cnv_drivers) <- c("alteration","subtype")

#---------------------------------------

# Now look at single nucleotide variants:
#---------------------------

# Get list of unique driver alterations
snv_drivers <- unique(unlist(strsplit(dat[,"snv_driver"],", ")))
snv_drivers <- sapply(strsplit(snv_drivers," "),function(x)x[1])
snv_drivers <- unique(snv_drivers)
snv_drivers <- snv_drivers[-which(is.na(snv_drivers))]

# Select variables to test:
# Testing how specific cell infiltration changes in response to the presence/absence of different genes
snv_genes <- snv_drivers

# Create results matrix
snv_results <- matrix(NA, nrow=length(snv_genes), ncol=length(subtypes))
rownames(snv_results) <- snv_genes
colnames(snv_results) <- subtypes

sig_snv <- c()
sig_snv_subtype <- c()
for(i in 1:length(snv_genes))
{
	#Create a driver_status vector 
	driver_status <- rep(0, nrow(dat))
	
	# Identify all samples that have a mutation in this gene
	indices <- grep(snv_genes[i], dat[,"snv_driver"])	
	
	# Change driver_status vector values to indicate driver is present at initial tumor
	driver_status[indices] <- 1
	if(sum(driver_status)==0){
		next}
	test_dat <- cbind(dat, driver_status)
	
	for(j in 1:length(subtypes))
	{
		g1 <- nrow(test_dat[which(test_dat[,"driver_status"]==1 & test_dat[,"subtype"]==subtypes[j]),])
		g2 <- nrow(test_dat[which(test_dat[,"driver_status"]==0 & test_dat[,"subtype"]==subtypes[j]),])
		g3 <- nrow(test_dat[which(test_dat[,"driver_status"]==1 & test_dat[,"subtype"]!=subtypes[j]),])
		g4 <- nrow(test_dat[which(test_dat[,"driver_status"]==0 & test_dat[,"subtype"]!=subtypes[j]),])
		
		p.val <- fisher.test(matrix(c(g1,g2,g3,g4),nrow=2,ncol=2))$p.value
		eff <- g1/(g1+g2) > g3/(g3+g4)
		p.val <- ifelse(eff, p.val, p.val*-1)
		snv_results[i,j] <- p.val
		
		# Identify all SNVs that are associated with a specific subtype
		# Ignoring SNVs that are significantly depleted, as this is likely due to enrichment in other subtypes
		if(p.val < 0.1 & p.val > 0 & sum(driver_status) > 5){
		sig_snv <- c(sig_snv, snv_genes[i])
		sig_snv_subtype <- c(sig_snv_subtype, subtypes[j])
		}
	}
	
}

sig_snv_drivers <- data.frame(sig_snv, sig_snv_subtype)
colnames(sig_snv_drivers) <- c("alteration","subtype")

# Combine the two lists
subtype_drivers <- rbind(sig_cnv_drivers, sig_snv_drivers)

###################################################
# Step 4: Test how the cancer cell fraction of each alteration associates with subtype and simplicity
##################################################

# Identify samples that are clonal for ANY of the subtype drivers

ccf_matrix <- matrix(0, nrow = nrow(dat), ncol = nrow(subtype_drivers))
colnames(ccf_matrix) <- subtype_drivers[,"alteration"]

for(i in 1:nrow(sig_cnv_drivers))
{
	mycnv <- as.character(sig_cnv_drivers[i,"alteration"])
	mysubtype <- sig_cnv_drivers[i,"subtype"]
	
	# Get the CCF for the alteration specifically
	cnv_split <- strsplit(dat[,"cnv_driver"],",")
	cnv_ccf_split <- strsplit(dat[,"cnv_driver_ccf"],",")
	cnv_ind <- lapply(cnv_split,function(x)grep(mycnv,x))

	cnv_ccf <- rep(0, length(cnv_ind))
	for(i in 1:length(cnv_ind))
	{
		if(length(cnv_ind[[i]]) == 0){
			next}
		
		cnv_ccf[i] <- as.numeric(cnv_ccf_split[[i]][cnv_ind[[i]]])
	}
	ccf_matrix[,mycnv] <- cnv_ccf
}

# Now test the single nucleotide variants

for(i in 1:nrow(sig_snv_drivers))
{
	mysnv <- as.character(sig_snv_drivers[i,"alteration"])
	mysubtype <- sig_snv_drivers[i,"subtype"]
	
	# Get the CCF for the alteration specifically
	snv_split <- strsplit(dat[,"snv_driver"],",")
	snv_ccf_split <- strsplit(dat[,"snv_driver_ccf"],",")
	snv_ind <- lapply(snv_split,function(x)grep(mysnv,x))

	snv_ccf <- rep(0, length(snv_ind))
	for(i in 1:length(snv_ind))
	{
		if(length(snv_ind[[i]]) == 0){
			next}
		
		snv_ccf[i] <- as.numeric(snv_ccf_split[[i]][snv_ind[[i]]])
	}
	ccf_matrix[,mysnv] <- snv_ccf
}

# Test how the presence of each clonal alteration associates with simplicity

ccf_p.val <- ccf_median <- clonal_num <- rep(0, ncol(ccf_matrix))
names(ccf_p.val) <- names(ccf_median) <- names(clonal_num) <- colnames(ccf_matrix)
for(i in 1:ncol(ccf_matrix))
{
	g1 <- dat[which(ccf_matrix[,i] >= 0.75),"simplicity_score"]
	g2 <- dat[which(ccf_matrix[,i] >= 0 & ccf_matrix[,i] < 0.75),"simplicity_score"]
	ccf_p.val[i] <- wilcox.test(g1,g2)$p.value
	ccf_median[i] <- median(g1) - median(g2)
	clonal_num[i] <- length(g1)
} 

any_clonal <- apply(ccf_matrix, 1, function(x)sum(x > 0.75, na.rm=TRUE))
clonal_p.val <- wilcox.test(dat[which(any_clonal %in% c(1,2)),"simplicity_score"], dat[which(any_clonal == 0),"simplicity_score"])


aliquot_barcode <- dat[,"aliquot_barcode"]
subtype <- dat[,"subtype"]

ccf_matrix <- data.frame(aliquot_barcode, ccf_matrix, subtype)
plot_ccf <- melt(ccf_matrix)
simplicity_score <- rep(dat[,"simplicity_score"],6)
plot_ccf[,"simplicity_score"] <- simplicity_score

# Part 1: Heatmap of CCFs in the driver alterations
# Part 2: Sidebar with dominant subtype
# Part 3: Side bar with IDH codel subtype
# Part 4: Bargraph of simplicity score
# Whole plot is sorted by simplicity score