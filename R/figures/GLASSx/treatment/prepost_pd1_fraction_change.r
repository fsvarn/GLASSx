library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(topGO)


rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "
SELECT tc.*, CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_ivygap sc1 ON sc1.aliquot_barcode = tc.tumor_barcode_a
JOIN analysis.cibersortx_ivygap sc2 ON sc2.aliquot_barcode = tc.tumor_barcode_b AND sc1.cell_state = sc2.cell_state
WHERE received_pd1 AND idh_codel_subtype = 'IDHwt'
ORDER BY 1
"
dat <- dbGetQuery(con, q)


dat %>% 
group_by(cell_state) %>%
summarise(p.val = t.test(fraction_a, fraction_b, paired=TRUE)$p.value,
		  eff = median(fraction_b - fraction_a))


dat %>%
filter(base::grepl("_tumor", cell_state)) %>%
group_by(tumor_pair_barcode) %>%
summarise(cell_state,
fraction_a = fraction_a/sum(fraction_a), 
fraction_b = fraction_b/sum(fraction_b)) %>%
group_by(cell_state) %>%
summarise(p.val = t.test(fraction_a, fraction_b, paired=TRUE)$p.value,
		  eff = median(fraction_b - fraction_a))

dat %>%
filter(!base::grepl("_tumor", cell_state)) %>%
group_by(tumor_pair_barcode) %>%
summarise(cell_state,
fraction_a = fraction_a/sum(fraction_a), 
fraction_b = fraction_b/sum(fraction_b)) %>%
group_by(cell_state) %>%
summarise(p.val = t.test(fraction_a, fraction_b, paired=TRUE)$p.value,
		  eff = median(fraction_b - fraction_a))


q <- "
SELECT tc.*, CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
WHERE received_pd1 
ORDER BY 1
"

# Get treatment info
dat <- dbGetQuery(con,q)

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
#myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")


geps <- read.delim(myinf1[11], row.names=1)
geps <- log10(geps+1)
rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
colnames(geps) <- gsub("\\.","-",colnames(geps))

eff <- apply(geps, 1, function(x) log2(mean(x[dat$tumor_barcode_b])/mean(x[dat$tumor_barcode_a])))
p.val <- apply(geps, 1, function(x) wilcox.test(x[dat$tumor_barcode_b], x[dat$tumor_barcode_a], paired = TRUE)$p.value)
q.val <- p.adjust(p.val,"BH")

res <- data.frame(eff, p.val, q.val)
res <- res[order(p.val),]

