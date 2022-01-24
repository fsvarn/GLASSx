library(odbc)
library(DBI)
library(tidyverse)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in aggregate pathologist annotations
myinf1 <- "data/pathology_annot/Histo_review_aggregate_v2.txt"

dat <- read.delim(myinf1)

dat_long <- dat %>%
  select(-Notes) %>%
  pivot_longer(-c(CASES,Pathologist),names_to ="feature", values_to = "metric")

feat_class <- rep("Primary",nrow(dat_long))
feat_class[grepl("NecTreatEff|Fibrosis|DepopTum", dat_long$feature)] <- "Recurrent"
dat_long$feat_class <- feat_class

metric_class <- rep("Percent",nrow(dat_long))
metric_class[grepl("_score", dat_long$feature)] <- "Score"
dat_long$metric_class <- metric_class

dat_long <- dat_long %>%
  mutate(feature = fct_relevel(feature, unique(dat_long$feature))) %>%
  filter(Pathologist=="JH")

dat_comp <- dat_long %>%
  select(-c(Pathologist, feat_class, metric_class)) %>%
  pivot_wider(names_from = feature, values_from = metric)

# Transcriptomic analysis

# Read in IvyGAP table
ivy <- dbReadTable(con, Id(schema="analysis",table="cibersortx_ivygap"))


ivy$aliquot_barcode <- sapply(strsplit(ivy$aliquot_barcode,"-"),function(x)paste(x[1:4],collapse="-"))
ivy <- ivy %>% filter(aliquot_barcode %in% unique(dat_long$CASES))

ivy$fraction <- ivy$fraction * 100

ivy_comp <- ivy %>%
  pivot_wider(names_from = cell_state, values_from = fraction)

# Compare both metrics
full_comp <- dat_comp %>%
  inner_join(ivy_comp, by = c("CASES" = "aliquot_barcode")) %>%
  column_to_rownames("CASES")

le_comp <- full_comp[,c("LE","LE_percent")]
ct_comp <- full_comp[,c("CT","CT_percent")]
pan_comp <- full_comp[,c("CTpan","Necrosis_percent")]
colnames(le_comp) <- colnames(ct_comp) <- colnames(pan_comp) <- c("transcriptomic", "pathology")

plot_comp <- rbind(le_comp, ct_comp, pan_comp)
cor(plot_comp, method="p")
cor(plot_comp, method="s")

plot_comp[,"sample_barcode"] <- rep(rownames(le_comp), 3)
rownames(plot_comp) <- NULL
status <- rep("Recurrent", nrow(plot_comp))
status[which(grepl("-TP", plot_comp$sample_barcode))] <- "Initial"
plot_comp[,"status"] <- status
plot_comp[,"cohort"] <- substring(plot_comp$sample_barcode,6,7)
plot_comp[,"feature"] <- rep(c("LE","CT","Necrosis"), each = nrow(le_comp))
plot_comp <- plot_comp[,c(3, 1, 2, 6, 4, 5)]

plot_comp %>% 
group_by(feature, status, cohort) %>%
summarise(cor(transcriptomic, pathology,method="s"))

pdf("figures/pathology/feature_cor_HF_JH.pdf", height = 6, width = 4)
ggplot(data = plot_comp %>% filter(cohort == "HF"), aes(x = pathology, y = transcriptomic)) +
  geom_point() +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Pathology", y = "Transcriptomic") +
  facet_grid(feature~status, scales= "free", space = "free") +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
dev.off()

pdf("figures/pathology/feature_cor_LU_JH.pdf", height = 6, width = 4)
ggplot(data = plot_comp %>% filter(cohort == "LU"), aes(x = pathology, y = transcriptomic)) +
  geom_point() +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Pathology", y = "Transcriptomic") +
  facet_grid(feature~status) +
  theme_classic() +
  theme(axis.text = element_text(size=12),
        axis.title= element_text(size=12),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(colour="white",fill="white")) 
dev.off()

plot_comp %>% 
  group_by(feature, status) %>%
  summarise(r=cor(transcriptomic, pathology), s = cor(transcriptomic, pathology, method="s"), p = cor.test(transcriptomic, pathology,method="s")$p.value)

plot_comp[,"case_barcode"] <- substring(plot_comp$sample_barcode, 1,10)
test <- plot_comp %>% 
        filter(status == "Initial", feature == "LE") %>%
        inner_join(plot_comp %>% filter(status == "Recurrent", feature == "LE"), by = "case_barcode")

up <- test %>% filter(transcriptomic.y < transcriptomic.x)
wilcox.test(up$pathology.x, up$pathology.y,paired=TRUE)

q <- "SELECT ts.*
FROM analysis.top_transcriptional_subtype ts
JOIN analysis.rna_silver_set ss1 ON ts.aliquot_barcode = ss1.tumor_barcode_a

UNION

SELECT ts.*
FROM analysis.top_transcriptional_subtype ts
JOIN analysis.rna_silver_set ss1 ON ts.aliquot_barcode = ss1.tumor_barcode_b"
top_ts <- dbGetQuery(con, q)
top_ts <- top_ts %>% 
filter(signature_name == "Mesenchymal") %>%
mutate(sample_barcode = substring(aliquot_barcode, 1, 15))

mes_status <- rep("No", nrow(plot_comp))
mes_status[which(plot_comp$sample_barcode %in% top_ts$sample_barcode)] <- "Mes"
plot_comp[,"mes_status"] <- mes_status

g1 <- plot_comp[which(plot_comp$mes_status == "Mes" & plot_comp$feature == "Necrosis" & plot_comp$cohort == "HF"), "pathology"]
g2 <- plot_comp[which(plot_comp$mes_status != "Mes" & plot_comp$feature == "Necrosis" & plot_comp$cohort == "HF"), "pathology"]

g1 <- dat_long[which(dat_long$CASES %in% top_ts$sample_barcode & dat_long$feature == "PN_percent"), ]$metric
g2 <- dat_long[which(!(dat_long$CASES %in% top_ts$sample_barcode) & dat_long$feature == "PN_percent"),]$metric
t.test(g1,g2)

pdf("figures/pathology/feature_cor_JH.pdf", height = 3.5, width = 4)
ggplot(data = plot_comp, aes(x = pathology, y = transcriptomic)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "Pathology", y = "Transcriptomic") +
  facet_wrap(status ~ feature, scales = "free") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white")) 
dev.off()

