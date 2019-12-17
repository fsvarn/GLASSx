#######################################################

library(tidyverse)
library(odbc)
library(DBI)
library(RColorBrewer)
library(gridExtra)

#######################################################

#Establish connection
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

#Read in data
q <- "SELECT ps.case_barcode, 
im1.signature_name, 
im1.enrichment_score AS es_a,
im2.enrichment_score AS es_b, 
an1.prop_aneuploidy AS prop_aneuploidy_a, 
an2.prop_aneuploidy AS prop_aneuploidy_b,
an1.aneuploidy_score AS aneuploidy_score_a, 
an2.aneuploidy_score AS aneuploidy_score_b,
cs1.idh_codel_subtype
FROM analysis.platinum_set ps
JOIN analysis.davoli_immune_score im1 ON im1.aliquot_barcode = ps.rna_barcode_a
JOIN analysis.davoli_immune_score im2 ON im2.aliquot_barcode = ps.rna_barcode_b AND im2.signature_name = im1.signature_name
JOIN analysis.gatk_aneuploidy an1 ON an1.aliquot_barcode = ps.dna_barcode_a
JOIN analysis.gatk_aneuploidy an2 ON an2.aliquot_barcode = ps.dna_barcode_b
JOIN biospecimen.aliquots al1 ON al1.aliquot_barcode = im1.aliquot_barcode
JOIN biospecimen.aliquots al2 ON al2.aliquot_barcode = im2.aliquot_barcode
JOIN clinical.surgeries cs1 ON cs1.sample_barcode = al1.sample_barcode
JOIN clinical.surgeries cs2 ON cs2.sample_barcode = al2.sample_barcode
JOIN biospecimen.samples bs1 ON bs1.sample_barcode = al1.sample_barcode
JOIN biospecimen.samples bs2 ON bs2.sample_barcode = al2.sample_barcode
--WHERE cs1.idh_codel_subtype LIKE 'IDHwt'
ORDER BY case_barcode
"

dat <- dbGetQuery(con,q)

cells <- unique(dat[,"signature_name"])

score_tp_cor <- score_tp_pval <- score_rec_cor <- score_rec_pval <- prop_tp_cor <- prop_tp_pval <- prop_rec_cor <- prop_rec_pval <- rep(0,length(cells))
for(i in 1:length(cells))
{

	sub_dat <- dat[which(dat[,"signature_name"]==cells[i]),]
	
	score_tp_cor[i] <- cor(sub_dat[,"aneuploidy_score_a"],sub_dat[,"es_a"],method="s")
	score_tp_pval[i] <- cor.test(sub_dat[,"aneuploidy_score_a"],sub_dat[,"es_a"],method="s")$p.value
	
	score_rec_cor[i] <- cor(sub_dat[,"aneuploidy_score_b"],sub_dat[,"es_b"],method="s")
	score_rec_pval[i] <- cor.test(sub_dat[,"aneuploidy_score_b"],sub_dat[,"es_b"],method="s")$p.value
	
		
	prop_tp_cor[i] <- cor(sub_dat[,"prop_aneuploidy_a"],sub_dat[,"es_a"],method="s")
	prop_tp_pval[i] <- cor.test(sub_dat[,"prop_aneuploidy_a"],sub_dat[,"es_a"],method="s")$p.value
	
	prop_rec_cor[i] <- cor(sub_dat[,"prop_aneuploidy_b"],sub_dat[,"es_b"],method="s")
	prop_rec_pval[i] <- cor.test(sub_dat[,"prop_aneuploidy_b"],sub_dat[,"es_b"],method="s")$p.value
	
	
}

res <- data.frame(cells, score_tp_cor, score_tp_pval, score_rec_cor, score_rec_pval, prop_tp_cor, prop_tp_pval, prop_rec_cor, prop_rec_pval)

# IDHwt (n = 63); Proportion of aneuploidy is significantly associated with lower DCs
#                    cells score_tp_cor score_tp_pval score_rec_cor
# 1         Macrophages.M1  -0.07391952     0.5647790  -0.060249496
# 2         Macrophages.M2  -0.02823619     0.8261254  -0.108739452
# 3            Macrophages  -0.03800464     0.7674449  -0.143001817
# 4  CD8.effector.NK.cells   0.08558335     0.5048238  -0.039464630
# 5              Dendritic  -0.01336480     0.9172009  -0.182103498
# 6                  T.reg   0.03307182     0.7969370   0.008662377
# 7                B.cells   0.06560904     0.6094361   0.183990831
# 8               NK.cells  -0.17486024     0.1704673   0.007839693
# 9           CD8.effector  -0.02804179     0.8273040  -0.125677062
# 10            CD4.mature   0.15252887     0.2327069  -0.009194702
#    score_rec_pval  prop_tp_cor prop_tp_pval prop_rec_cor prop_rec_pval
# 1       0.6390226 -0.028657834   0.82320791 -0.087630660    0.49464852
# 2       0.3962549 -0.032402074   0.80054110 -0.186038139    0.14433387
# 3       0.2635369 -0.027841782   0.82816809 -0.185534101    0.14544372
# 4       0.7587765 -0.002976190   0.98162513 -0.146531136    0.25181044
# 5       0.1531675 -0.211501536   0.09610912 -0.250315024    0.04785529
# 6       0.9462793  0.118423579   0.35437547 -0.007008533    0.95652471
# 7       0.1488815 -0.030769969   0.81040243  0.093247087    0.46729286
# 8       0.9513748 -0.083909370   0.51227714 -0.096823359    0.45030995
# 9       0.3263691 -0.007632488   0.95266680 -0.213568231    0.09283060
# 10      0.9429833  0.144537250   0.25768571 -0.006000456    0.96277320

#IDHmut-noncodel (n = 18); Proportion of aneuploidy is significantly associated with higher CD8.effector.NK.cells
#                    cells score_tp_cor score_tp_pval score_rec_cor
# 1             CD4.mature   0.30548834     0.2176635  -0.366953582
# 2           CD8.effector   0.29179403     0.2400391   0.112269084
# 3               NK.cells   0.25387134     0.3093726   0.168403627
# 4                B.cells  -0.12219534     0.6290751   0.006237171
# 5                  T.reg  -0.06952493     0.7839890   0.560305894
# 6         Macrophages.M1  -0.04108291     0.8714206  -0.074846056
# 7              Dendritic  -0.22332251     0.3730458  -0.108110970
# 8  CD8.effector.NK.cells   0.30548834     0.2176635   0.287949411
# 9            Macrophages  -0.17591915     0.4850176  -0.218300998
# 10        Macrophages.M2  -0.20120094     0.4233758  -0.329530554
#    score_rec_pval prop_tp_cor prop_tp_pval prop_rec_cor prop_rec_pval
# 1      0.13415525  0.15583075   0.53564005   0.01341589     0.9606477
# 2      0.65738328  0.43034056   0.07612890   0.19298246     0.4413564
# 3      0.50414694  0.31888545   0.19681508   0.39525284     0.1055265
# 4      0.98040403  0.19298246   0.44135639   0.06707946     0.7922816
# 5      0.01558387 -0.10010320   0.69253961   0.28998968     0.2422781
# 6      0.76787016  0.04231166   0.86932594   0.12899897     0.6091089
# 7      0.66937753 -0.09391125   0.71092178   0.15376677     0.5411403
# 8      0.24657423  0.43240454   0.07461608   0.48194014     0.0446757
# 9      0.38416539 -0.13312693   0.59753834  -0.13519092     0.5917882
# 10     0.18175464 -0.17234262   0.49260688  -0.18885449     0.4513725

#IDHmut-codel (n = 5)
#                    cells score_tp_cor score_tp_pval score_rec_cor
# 1         Macrophages.M1    0.9746794    0.00481823           0.4
# 2         Macrophages.M2    0.8720816    0.05385422           0.3
# 3            Macrophages    0.9746794    0.00481823           0.3
# 4  CD8.effector.NK.cells    0.4616903    0.43376616          -0.1
# 5              Dendritic    0.5642881    0.32172334          -0.3
# 6                  T.reg    0.3590924    0.55281475          -0.1
# 7                B.cells    0.3590924    0.55281475           0.6
# 8               NK.cells   -0.6668859    0.21889398          -0.2
# 9           CD8.effector    0.6155870    0.26899777           0.4
# 10            CD4.mature    0.6668859    0.21889398           0.3
#    score_rec_pval prop_tp_cor prop_tp_pval prop_rec_cor prop_rec_pval
# 1       0.5166667         0.1    0.9500000          0.3     0.6833333
# 2       0.6833333        -0.3    0.6833333          0.1     0.9500000
# 3       0.6833333         0.1    0.9500000          0.1     0.9500000
# 4       0.9500000         0.4    0.5166667          0.3     0.6833333
# 5       0.6833333        -0.5    0.4500000          0.0     1.0000000
# 6       0.9500000         0.6    0.3500000          0.3     0.6833333
# 7       0.3500000         0.6    0.3500000          0.7     0.2333333
# 8       0.7833333         0.6    0.3500000          0.1     0.9500000
# 9       0.5166667         0.3    0.6833333          0.7     0.2333333
# 10      0.6833333        -0.2    0.7833333         -0.1     0.9500000

#IDHmut together (n = 23)
#                    cells score_tp_cor score_tp_pval score_rec_cor
# 1         Macrophages.M1   0.18903321    0.38767468    0.09052679
# 2         Macrophages.M2   0.03639142    0.86906326   -0.11788378
# 3            Macrophages   0.08238613    0.70861826   -0.09699298
# 4  CD8.effector.NK.cells   0.34824568    0.10343329    0.22283517
# 5              Dendritic  -0.09451661    0.66794496   -0.19398597
# 6                  T.reg   0.10058184    0.64793390    0.39841734
# 7                B.cells  -0.07278284    0.74137972    0.18005877
# 8               NK.cells   0.09603292    0.66292100    0.10246438
# 9           CD8.effector   0.36441964    0.08734119    0.15717838
# 10            CD4.mature   0.37048488    0.08181653   -0.16364457
#    score_rec_pval prop_tp_cor prop_tp_pval prop_rec_cor prop_rec_pval
# 1      0.68122992  0.04249012   0.84780855  0.268774704    0.21412138
# 2      0.59216794 -0.19169960   0.37917508  0.025691700    0.90834558
# 3      0.65974724 -0.13438735   0.53942370  0.021739130    0.92267601
# 4      0.30677760  0.44466403   0.03474094  0.370553360    0.08252017
# 5      0.37512700 -0.14426877   0.50969711  0.055335968    0.80206503
# 6      0.05969909 -0.01185771   0.95859672  0.174901186    0.42299715
# 7      0.41100472  0.28656126   0.18442942  0.198616601    0.36193023
# 8      0.64176965  0.42984190   0.04184738  0.224308300    0.30205433
# 9      0.47384621  0.33596838   0.11733098  0.153162055    0.48365654
# 10     0.45560858  0.04841897   0.82662634  0.005928854    0.98019193

