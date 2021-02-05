
# Functional annotation in each tumor cell state

for(i in 1:length(sigs))
{
	mysig <- sigs[[1]]
	mysig %>%filter(sig & eff > 0) %>% arrange(q.val)
	
	tmp <- read.csv("/projects/verhaak-lab/GLASS-III/data/dataset/scgp/neuron_cell_states_comp_kcj_20210122.csv")
	tmp <- tmp[1:100,2]
	test <- rownames(mysig %>% filter(sig))
	
 	fisher.test(matrix(c(length(all_genes) - length(union(test, tmp)),
	 								 length(setdiff(test, tmp)), 
	 								 length(setdiff(tmp, test)), 
	 								 length(intersect(test, tmp))), nrow=2),alternative="greater")$p.value
	
	all_genes <- rownames(mysig)
	bg_genes <- as.numeric(all_genes %in% rownames(mysig %>% filter(sig)))
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
	fishRes[which(fishRes$q.value < 0.05 & fishRes$Significant > fishRes$Expected),]

	#Parentchild
	resultPC <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
	pcRes <- GenTable(sampleGOdata, raw.p.value = resultPC, topNodes = length(resultPC@score), numChar = 120)
	pcRes[,"q.value"] <- p.adjust(pcRes[,"raw.p.value"],"BH")
	pcRes[which(pcRes$q.value < 0.1 & pcRes$Significant > pcRes$Expected),]

	# Use this one: elim algorithm with fisher statistic. elim algorithm removes false positives by getting adjusting for co-dependent GO terms 
	resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
	elimRes <- GenTable(sampleGOdata, raw.p.value = resultElim, topNodes = length(resultElim@score), numChar = 120)
	elimRes[,"q.value"] <- p.adjust(elimRes[,"raw.p.value"],"BH")
	elimRes[which(elimRes$q.value < 0.05 & elimRes$Significant > elimRes$Expected),]
	elimRes[,"logp"] <- -log10(as.numeric(elimRes$raw.p.value))

	term_level <- rev(as.character(elimRes$Term))
	plotRes <- elimRes %>% 
			   mutate(Term = as_factor(Term)) %>%
			   mutate(Term = fct_relevel(Term, term_level)) %>%
			   filter(q.value < 0.05 & elimRes$Significant > elimRes$Expected)
}




# Survival associations

survdat <- dat

for(i in 1:length(myinf1))
{
	cat("\r", i)
	tag <- cell_state[i]
	
	geps <- read.delim(myinf1[i], row.names=1)
	geps <- log10(geps+1)
	rem <- apply(geps,1,function(x)sum(is.na(x)))
	geps <- geps[-which(rem==ncol(geps)),]	
	vars <- apply(geps,1,function(x)var(x,na.rm=TRUE))
	geps <- geps[which(vars > 0),]
	colnames(geps) <- gsub("\\.","-",colnames(geps))
	
	sig <- sigs[[i]]
	
	mysig <- rownames(sig %>% filter(q.val < 0.1))
	
	coll <- apply(geps[mysig,], 2, mean)
	
	score_a <- coll[dat$tumor_barcode_a]
	score_b <- coll[dat$tumor_barcode_b]
	score_diff <-  coll[dat$tumor_barcode_b] - coll[dat$tumor_barcode_a]
	survdat <- data.frame(survdat, score_a, score_b, score_diff)

	cols <- c("score_a", "score_b", "score_diff")
	cols <- setNames(cols, paste0(tag, ".", cols, sep=""))

	survdat <- survdat %>% rename_(.dots = cols)
}

survdat[,"event"] <- 1

summary(coxph(Surv(surgical_interval, event) ~ differentiated_tumor.score_a + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ stemcell_tumor.score_a + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ prolif_stemcell_tumor.score_a + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ myeloid.score_a + case_age_diagnosis_years, data =survdat))

summary(coxph(Surv(surgical_interval, event) ~ differentiated_tumor.score_b + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ stemcell_tumor.score_b + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ prolif_stemcell_tumor.score_b + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ myeloid.score_b + case_age_diagnosis_years, data =survdat))

summary(coxph(Surv(surgical_interval, event) ~ differentiated_tumor.score_diff + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ stemcell_tumor.score_diff + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ prolif_stemcell_tumor.score_diff + case_age_diagnosis_years, data =survdat))
summary(coxph(Surv(surgical_interval, event) ~ myeloid.score_diff + case_age_diagnosis_years, data =survdat))

summary(coxph(Surv(case_overall_survival_mo, case_vital_status) ~ differentiated_tumor.score_a + case_age_diagnosis_years, data =dat))


thr <- ifelse(survdat$stemcell_tumor.score_a > mean(survdat$stemcell_tumor.score_a), 1, 0)
survdat$thr <- thr
summary(coxph(Surv(surgical_interval, event) ~ thr + case_age_diagnosis_years, data =survdat))
