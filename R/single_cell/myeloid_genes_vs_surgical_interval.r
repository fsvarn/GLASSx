library(tidyverse)
library(odbc)
library(DBI)
library(grid)
library(survival)
library(topGO)

#######################################################
rm(list=ls())

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

q <- "SELECT tc.*, cc.case_age_diagnosis_years, 
CASE WHEN tc.idh_codel_subtype = 'IDHwt' THEN 'IDHwt' ELSE 'IDHmut' END AS idh_status,
sc1.cell_state, sc1.fraction AS fraction_a, sc2.fraction AS fraction_b, sc2.fraction - sc1.fraction AS diff
FROM analysis.tumor_rna_clinical_comparison tc
JOIN analysis.rna_silver_set ss ON tc.tumor_pair_barcode = ss.tumor_pair_barcode
JOIN analysis.cibersortx_scgp sc1 ON sc1.aliquot_barcode = ss.tumor_barcode_a
JOIN analysis.cibersortx_scgp sc2 ON sc2.aliquot_barcode = ss.tumor_barcode_b AND sc2.cell_state = sc1.cell_state
JOIN clinical.cases cc ON cc.case_barcode = tc.case_barcode
WHERE sc1.cell_state = 'myeloid' AND idh_codel_subtype = 'IDHwt'
"

dat <- dbGetQuery(con, q)

myDir1 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/GLASS-freeze/"
mytag <- dir(myDir1)
myinf1 <- mytag[grep("_Window48.txt", mytag)]
myinf1 <- myinf1[grep("myeloid|stemcell|differentiated",myinf1)]		# These are the signatures we have some confidence in
mytag <- gsub("CIBERSORTxHiRes_GLASS_", "", myinf1)
mytag <- gsub("_Window48.txt", "", mytag)
myinf1 <- paste(myDir1, myinf1, sep = "/")


gep_list <- good_list <- bad_list <- list()

geps <- read.delim(myinf1[2], row.names=1)
colnames(geps) <- gsub("\\.","-",colnames(geps))
geps <- log10(geps+1)

rem <- apply(geps,1,function(x)sum(is.na(x)))
geps <- geps[-which(rem==ncol(geps)),]
vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
geps <- geps[which(vars > 0),]
all_genes <- rownames(geps)

geps <- t(geps)
geps <- geps[dat$tumor_barcode_b,]

cohort <- substr(rownames(geps), 1, 7)

eff <- apply(geps, 2, function(x) summary(lm(dat$fraction_b ~ x + cohort))$coefficients[2,1])
p.val <- apply(geps, 2, function(x) summary(lm(dat$fraction_b ~ x + cohort))$coefficients[2,4])
q.val <- p.adjust(p.val, "BH")

frac_res <- data.frame(eff, p.val,q.val) %>% rownames_to_column("gene") %>% arrange(p.val)

eff <- apply(geps, 2, function(x) summary(lm(dat$surgical_interval ~ x + cohort))$coefficients[2,1])
p.val <- apply(geps, 2, function(x) summary(lm(dat$surgical_interval ~ x + cohort))$coefficients[2,4])
q.val <- p.adjust(p.val, "BH")

surg_res <- data.frame(eff, p.val,q.val) %>% rownames_to_column("gene") %>% arrange(p.val)

surg_genes <- surg_res %>% filter(q.val < 0.2 & eff < 0) %>% .$gene

#######################################################

# Test in CGGA693

myinf2 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/cgga_693_scgp/CIBERSORTxHiRes_cgga_myeloid_Window48.txt"
myinf3 <- "/projects/verhaak-lab/GLASS-III/results/cibersortx/hires/cgga_325_scgp/CIBERSORTxHiRes_cgga_myeloid_Window48.txt"
#myinf2 <- "/projects/verhaak-lab/GLASS-III/data/dataset/cgga/CGGA.mRNAseq_693.RSEM-genes.20200506.txt"

myinf4 <- "/projects/verhaak-lab/GLASS-III/data/dataset/cgga/CGGA.mRNAseq_693_clinical.20200506.txt"
myinf5 <- "/projects/verhaak-lab/GLASS-III/data/dataset/cgga/CGGA.mRNAseq_325_clinical.20200506.txt"

expr1 <- read.delim(myinf2, row.names=1)
expr1 <- log10(expr1+1)
rem <- apply(expr1,1,function(x)sum(is.na(x)))
expr1 <- expr1[-which(rem==ncol(expr1)),]
vars <- apply(expr1,1,function(x)var(x,na.rm=TRUE))
expr1 <- expr1[which(vars > 0),]

expr2 <- read.delim(myinf3, row.names=1)
expr2 <- log10(expr2+1)
rem <- apply(expr2,1,function(x)sum(is.na(x)))
expr2 <- expr2[-which(rem==ncol(expr2)),]
vars <- apply(expr2,1,function(x)var(x,na.rm=TRUE))
expr2 <- expr2[which(vars > 0),]

sub_genes1 <- intersect(surg_genes, rownames(expr1))
rec_score1 <- apply(expr1, 2, function(x) mean(x[sub_genes1]))

sub_genes2 <- intersect(surg_genes, rownames(expr2))
rec_score2 <- apply(expr2, 2, function(x) mean(x[sub_genes2]))

rec_score <- c(rec_score1, rec_score2)

info1 <- read.delim(myinf4)
info1 <- info1 %>% filter(PRS_type != "Primary")
info1[,"cohort"] <- "693"

info2 <- read.delim(myinf5)
info2 <- info2 %>% filter(PRS_type != "Primary")
info2[,"cohort"] <- "325"

info <- rbind(info1, info2)

sub_score <- rec_score[info$CGGA_ID]
treatment <- info$Chemo_status..TMZ.treated.1.un.treated.0. | info$Radio_status..treated.1.un.treated.0.
surv <- data.frame(info, treatment, sub_score)

summary(coxph(Surv(OS, Censor..alive.0..dead.1.) ~ sub_score + IDH_mutation_status + Age + treatment + cohort, data = surv))

bg_genes <- as.numeric(all_genes %in% surg_genes)
names(bg_genes) <- all_genes

diffGenes <- function(bg_genes) {
	return(bg_genes == 1)}

# Functional enrichment of signature
sampleGOdata <- new("topGOdata",
				description = "Simple session", 
				ontology = "BP",
			    allGenes = bg_genes,
			    geneSelectionFun = diffGenes, 
			    nodeSize = 10,
			    annot = annFUN.org, mapping="org.Hs.eg.db", ID="symbol")

# Fishers test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher",)
fishRes <- GenTable(sampleGOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score), numChar=120)
fishRes[,"q.value"] <- p.adjust(fishRes[,"raw.p.value"],"BH")

# #Parentchild
# resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
# pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
# pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
# pcRes[which(pcRes$q.value < 0.05 & pcRes$Significant > pcRes$Expected),]

# Use this one: elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
elimRes[which(elimRes$raw.p.value < 0.05 & elimRes$Significant > elimRes$Expected),]
elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

term_level <- rev(as.character(fishRes$Term))
plotRes <- fishRes %>% 
		   mutate(Term = as_factor(Term)) %>%
		   mutate(Term = fct_relevel(Term, term_level)) %>%
		   filter(q.value < 0.1 & fishRes$Significant > fishRes$Expected)


