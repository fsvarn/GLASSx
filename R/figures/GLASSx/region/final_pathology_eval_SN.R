library(odbc)
library(DBI)
library(tidyverse)
library(corrr)

rm(list=ls())
setwd("/Users/varnf/Documents/Projects/GLASS/GLASS-III/GLASS-III/")

# Connect to DB
con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

# Read in aggregate pathologist annotations
myinf1 <- "data/pathology_annot/neuropath_template_final_aggregate.txt"

dat <- read.delim(myinf1)

dat <- dat %>%
       mutate(sample_barcode = recode(sample_barcode,
         "Pt10_primary" = "GLSS-SN-0010-TP", "Pt10_recurrence" = "GLSS-SN-0010-R1",
         "Pt16_primary" = "GLSS-SN-0016-TP", "Pt16_recurrence" = "GLSS-SN-0016-R1",
         "Pt17_primary" = "GLSS-SN-0017-TP", "Pt17_recurrence" = "GLSS-SN-0017-R1",
         "Pt6_primary" = "GLSS-SN-0015-TP", "Pt6_recurrence" = "GLSS-SN-0015-R2",
         "Pt8_primary" = "GLSS-SN-0009-TP", "Pt8_recurrence" = "GLSS-SN-0009-R1"))

# Pathologist agreement

path_dat <- dat %>% 
  filter(folder == "SN") %>%
  select(-c(folder, annotations, priority, total)) %>%
  pivot_longer(-c(sample_barcode, timepoint, pathologist), names_to = "feature", values_to = "fraction") %>%
  filter(!(timepoint == "Primary" & feature %in% c("depopulated_tumor", "fibrosis", "treatment_effec_other", "treatment_necrosis")))
path_coll <- path_dat %>%
  group_by(sample_barcode, feature, timepoint) %>%
  summarise(mean = mean(fraction), sd = sd(fraction), cov = sd(fraction)/mean(fraction)) %>%
  mutate(feature = fct_relevel(feature, "histologically_normal","leading_edge","cellular_tumor", "necrosis", 
                               "depopulated_tumor", "fibrosis", "treatment_necrosis", "treatment_effec_other"))

path_cov <- path_coll %>%
  group_by(feature) %>%
  summarise(mean_cov = mean(cov,na.rm=TRUE))

pdf("figures/pathology/SN_pathology_concordance.pdf",height = 4, width = 6.5)
ggplot(path_dat, aes(x = feature, y = fraction)) +
  geom_boxplot() +
  facet_grid(. ~ sample_barcode, scale = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title = element_blank(),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position="none") 
dev.off()

pdf("figures/pathology/SN_pathology_cov.pdf",height = 4, width = 4)
ggplot(path_coll, aes(x = feature, y = cov)) +
  geom_hline(yintercept = 1, colour = "gray") +
  geom_boxplot() +
  theme_bw() +
  facet_grid(.~timepoint, scales = "free_x", space = "free_x") +
  labs(y = "Coefficient of variation") + 
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position="none") 
dev.off()

dat_comb <- dat
dat_comb$leading_edge = dat_comb$histologically_normal + dat_comb$leading_edge
dat_comb <- dat_comb %>% select(-histologically_normal)
path_dat <- dat_comb %>% 
  filter(folder == "SN") %>%
  select(-c(folder, annotations, priority, total)) %>%
  pivot_longer(-c(sample_barcode, timepoint, pathologist), names_to = "feature", values_to = "fraction") %>%
  filter(!(timepoint == "Primary" & feature %in% c("depopulated_tumor", "fibrosis", "treatment_effec_other", "treatment_necrosis")))
path_coll <- path_dat %>%
  group_by(sample_barcode, feature, timepoint) %>%
  summarise(mean = mean(fraction), sd = sd(fraction), cov = sd(fraction)/mean(fraction)) %>%
  mutate(feature = fct_relevel(feature, "histologically_normal","leading_edge","cellular_tumor", "necrosis", 
                               "depopulated_tumor", "fibrosis", "treatment_necrosis", "treatment_effec_other"))

pdf("figures/pathology/SN_pathology_lecomb_cov.pdf",height = 4, width = 4)
ggplot(path_coll, aes(x = feature, y = cov)) +
  geom_hline(yintercept = 1, colour = "gray") +
  geom_boxplot() +
  theme_bw() +
  facet_grid(.~timepoint, scales = "free_x", space = "free_x") +
  labs(y = "Coefficient of variation") + 
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position="none") 
dev.off()

# Transcriptomic analysis

# Read in IvyGAP table
ivy <- dbReadTable(con, Id(schema="analysis",table="cibersortx_ivygap"))
ivy$sample_barcode <- sapply(strsplit(ivy$aliquot_barcode,"-"),function(x)paste(x[1:4],collapse="-"))
ivy$fraction <- ivy$fraction 
ivy <- ivy %>%
  pivot_wider(names_from = "cell_state", values_from = "fraction")

full_comp <- dat %>% 
  inner_join(ivy, by = "sample_barcode") %>%
  select(-c(aliquot_barcode, annotations, priority, total))

#--------------------------
# General analyses (focused more on primary-only)
# Working terms:
# Leading edge per pathology = Histologic normal + leading edge
# Cellular tumor IvyGAP % combines the three parts of the cellular tumor region (CT, CTmvp, CTpan)
# Additionally correlate CTpan with necrosis to show that they are found together

# Removing non-cellular content
full_comp_nonec <- full_comp[,-which(colnames(full_comp)=="necrosis" | colnames(full_comp)=="treatment_necrosis" | colnames(full_comp)=="fibrosis")]
# Renormalize area
full_comp_nonec[,c("histologically_normal", "leading_edge", "cellular_tumor", "depopulated_tumor","treatment_effec_other")] <- t(apply(full_comp_nonec[,c("histologically_normal", "leading_edge", "cellular_tumor", "depopulated_tumor","treatment_effec_other")], 1, function(x)x/sum(x)))
# Add necrosis back in as its own "score"
full_comp_nonec[,c("necrosis","treatment_necrosis")] <- full_comp[,c("necrosis","treatment_necrosis")]/100
full_comp_nonec <- full_comp_nonec[,c("sample_barcode", "timepoint", "folder", "pathologist", 
                                      "histologically_normal", "leading_edge", "cellular_tumor", "depopulated_tumor","treatment_effec_other",
                                      "necrosis", "treatment_necrosis", "LE", "CT", "CTpan", "CTmvp")]

plot_comp_nonec <- full_comp_nonec %>%
  filter(folder == "SN") %>%
  mutate(leading_edge = leading_edge + histologically_normal) %>%
  group_by(sample_barcode) %>%
  select(-c(histologically_normal, timepoint, folder, pathologist)) %>%
  summarise_if(is_numeric, mean) %>%
  pivot_longer(-sample_barcode, names_to = "feature", values_to = "fraction") %>%
  filter(!feature %in% c("necrosis", "treatment_necrosis"))
plot_comp_nonec$timepoint <- ifelse(grepl("-TP", plot_comp_nonec$sample_barcode),"Primary", "Recurrent")
plot_comp_nonec$source <- ifelse(plot_comp_nonec$feature %in% c("LE","CT","CTpan","CTmvp"), "Transcriptome","Pathology")
plot_comp_nonec$case_barcode <- substring(plot_comp_nonec$sample_barcode, 1, 12)
plot_comp_nonec <- plot_comp_nonec %>%
  mutate(feature = recode(feature, "leading_edge" = "Leading edge", "cellular_tumor" = "Cellular tumor", 
                          "depopulated_tumor" = "Depopulated tumor", "treatment_effec_other" = "Treatment effect other",
                          "LE" = "Leading edge", "CT" = "Cellular tumor", "CTpan" = "Cellular tumor", "CTmvp" = "Cellular tumor"))

ggplot(plot_comp_nonec, aes(fill=feature, y=fraction, x=case_barcode)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_manual(values=c("Total" = "#8B8B8B", "Subtotal" = "#ECECEC")) +
  facet_grid(source~timepoint) +
  labs(y = "Feature proportion (%)") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7,angle=45,hjust=1),
        axis.text.y=element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_blank(),
        legend.text = element_text(size=7))




# Goal distinguish CT from LE and show that PAN is found in more necrotic tumors
avg_cors <- full_comp_nonec %>%
  filter(folder == "SN") %>%
  group_by(pathologist,timepoint) %>%
  summarise(le_cor = cor(LE, histologically_normal + leading_edge), 
            ct_cor = cor(cellular_tumor, CT+CTpan+CTmvp), 
            nec_cor = cor(necrosis, CTpan)) %>%
  group_by(timepoint) %>%
  summarise(le_cor = mean(le_cor), ct_cor = mean(ct_cor), nec_cor = mean(nec_cor))

# Plot the correlation heatmap (no recurrent features)
plot_cors <- full_comp_nonec %>%
  filter(folder == "SN") %>%
  group_by(pathologist,timepoint) %>%
  summarise(leading_edge = cor(LE, histologically_normal + leading_edge), 
            cellular_tumor = cor(cellular_tumor, CT+CTpan+CTmvp), 
            necrosis = cor(necrosis, CTpan)) %>%
  pivot_longer(-c(pathologist, timepoint), values_to = "cor", names_to = "feature")

ggplot(plot_cors, aes(x = feature, y = pathologist, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
  facet_grid(. ~ timepoint) +
  theme_bw() +
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title = element_blank(),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position="none") 

# Mean pathologist correlation plot
mean_pathologist <- full_comp_nonec %>%
  filter(folder == "SN") %>%
  transmute(sample_barcode = sample_barcode,
            timepoint = timepoint,
            leading_edge = histologically_normal + leading_edge, 
            cellular_tumor = cellular_tumor, 
            depopulated_tumor = depopulated_tumor,
            treatment_effec_other = treatment_effec_other,
            necrosis = necrosis,
            treatment_necrosis = treatment_necrosis,
            LE = LE,
            CT = CT + CTpan + CTmvp, 
            CTpan = CTpan) %>%
  group_by(sample_barcode,) %>%
  summarise_if(is.numeric,mean) 
mean_results <- mean_pathologist %>%
  mutate(timepoint = ifelse(grepl("TP",sample_barcode),"Primary", "Recurrent")) %>%
  group_by(timepoint) %>%
  select_if(is.numeric) %>%
  group_map(~ correlate(.x)) 

pri_mean <- mean_results[[1]] %>%
  select(term, LE, CT, CTpan) %>%
  pivot_longer(-term, names_to = "cibersort", values_to = "cor") %>%
  mutate(timepoint = "Primary") %>%
  filter(!(term %in% c("depopulated_tumor", "treatment_effec_other", "treatment_necrosis", "necrosis", "LE", "CT", "CTpan"))) %>% 
  filter(!(cibersort =="CTpan"))

rec_mean <- mean_results[[2]] %>%
  select(term, LE, CT, CTpan) %>%
  pivot_longer(-term, names_to = "cibersort", values_to = "cor") %>%
  mutate(timepoint = "Recurrent") %>%
  filter(!(term %in% c("LE", "CT", "CTpan", "necrosis", "treatment_necrosis")))  %>%
  filter(!(cibersort=="CTpan"))

plot_mean <- rbind(pri_mean, rec_mean) %>%
  mutate(term = fct_relevel(term, "leading_edge", "cellular_tumor", "depopulated_tumor", "treatment_effec_other")) %>%
  mutate(cibersort = fct_relevel(cibersort, rev(c("LE", "CT")))) %>%
  mutate(cibersort = recode(cibersort, "LE" = "Leading edge", "CT" = "Cellular tumor")) %>%
  mutate(term = recode(term, "leading_edge" = "Leading edge", "cellular_tumor" = "Cellular tumor", 
                       "depopulated_tumor" = "Depopulated tumor", "treatment_effec_other" = "Treatment effect other"))


pdf("figures/pathology/SN_mean_correlations.pdf",height = 2, width = 3.75)
ggplot(plot_mean, aes(x = term, y = cibersort, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2("R",low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
  facet_grid(. ~ timepoint, scale = "free_x", space= "free_x") +
  labs(x = "Pathology (% area on slide)", y = "CIBERSORT") +
  theme_bw() +
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title = element_text(size=7),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.15, 'inch')) 
dev.off()

# Necrosis features
plot_nec <- mean_pathologist[,c("sample_barcode", "necrosis", "treatment_necrosis", "CTpan")]
plot_nec$timepoint <- ifelse(grepl("-TP", plot_nec$sample_barcode), "Primary", "Recurrent")
plot_nec <- plot_nec %>% mutate(timepoint = fct_relevel(timepoint, "Primary", "Recurrent"))
plot_nec$necrosis <- plot_nec$necrosis
plot_nec$treatment_necrosis <- plot_nec$treatment_necrosis
plot_nec$total_necrosis <- plot_nec$necrosis + plot_nec$treatment_necrosis

pdf("figures/pathology/total_necrosis_cor_glass_legend.pdf", height = 2, width = 2)
ggplot(data = plot_nec, aes(x = total_necrosis, y = CTpan*100)) +
  geom_point(aes(shape = timepoint)) +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Necrosis % (Pathology)", y = "PAN % (CIBERSORT)") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.text=element_text(size=7),
        legend.title = element_blank())
dev.off()

pdf("figures/pathology/total_necrosis_cor_glass.pdf", height = 2, width = 2)
ggplot(data = plot_nec, aes(x = total_necrosis, y = CTpan*100)) +
  geom_point(aes(shape = timepoint)) +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Necrosis % (Pathology)", y = "PAN % (CIBERSORT)") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position="none")
dev.off()

pdf("figures/pathology/necrosis_cor_glass.pdf", height = 2, width = 2)
ggplot(data = plot_nec, aes(x = necrosis, y = CTpan*100)) +
  geom_point(aes(shape = timepoint)) +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Necrosis % (Pathology)", y = "PAN % (CIBERSORT)") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position="none")
dev.off()

cor(plot_nec$CTpan, plot_nec$necrosis)
cor(plot_nec$CTpan, plot_nec$total_necrosis)

# Attempt at combining recurrent features
# Working terms:
# Leading edge per pathology = Histologic normal + leading edge
# Cellular tumor per pathology = Cellular tumor + depopulated tumor
# Necrosis per pathology = Necrosis + treatment necrosis
# Cellular tumor IvyGAP % combines the three parts of the cellular tumor region (CT, CTmvp, CTpan)
# Additionally correlate CTpan with necrosis to show that they are found together

plot_cors <- full_comp_nonec %>%
  filter(folder == "SN") %>%
  group_by(pathologist,timepoint) %>%
  summarise(leading_edge = cor(LE, histologically_normal + leading_edge), 
            cellular_tumor = cor(cellular_tumor+depopulated_tumor, CT+CTpan+CTmvp), 
            necrosis = cor(necrosis + treatment_necrosis, CTpan)) %>%
  pivot_longer(-c(pathologist, timepoint), values_to = "cor", names_to = "feature")

ggplot(plot_cors, aes(x = feature, y = pathologist, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2(low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
  facet_grid(. ~ timepoint) +
  theme_bw() +
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title = element_blank(),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position="none") 

# Done

# Remaining questions
# 1. Are my groupings reasonable?
# 2. Should PAN be correlated separately with necrosis and treatment necrosis?
# 3. Should I repeat recurrent analysis by removing the recurrence-only features and renormalizing?

# Repeat analysis by removing recurrent features and renormalizing:

# Removing non-cellular content
full_comp_norec <- full_comp[,-which(colnames(full_comp)=="necrosis" | colnames(full_comp)=="treatment_necrosis" | colnames(full_comp)=="fibrosis" | 
                                       colnames(full_comp) == "depopulated_tumor" | colnames(full_comp) == "treatment_effec_other")]
# Renormalize area
full_comp_norec[,c("histologically_normal", "leading_edge", "cellular_tumor")] <- t(apply(full_comp_norec[,c("histologically_normal", "leading_edge", "cellular_tumor")], 1, function(x)x/sum(x)))
# Add necrosis back in as its own "score"
full_comp_norec[,c("necrosis","treatment_necrosis")] <- full_comp[,c("necrosis","treatment_necrosis")]
full_comp_norec <- full_comp_norec[,c("sample_barcode", "timepoint", "folder", "pathologist", 
                                      "histologically_normal", "leading_edge", "cellular_tumor",
                                      "necrosis", "treatment_necrosis", "LE", "CT", "CTpan", "CTmvp")]

# Goal distinguish CT from LE and show that PAN is found in more necrotic tumors
avg_cors_norec <- full_comp_norec %>%
  filter(folder == "SN") %>%
  group_by(pathologist,timepoint) %>%
  summarise(le_cor = cor(LE, histologically_normal + leading_edge), 
            ct_cor = cor(cellular_tumor, CT+CTpan+CTmvp), 
            nec_cor = cor(necrosis+treatment_necrosis, CTpan)) %>%
  group_by(timepoint) %>%
  summarise(le_cor = mean(le_cor), ct_cor = mean(ct_cor), nec_cor = mean(nec_cor))


# Transcriptional subtype test

ts <- dbReadTable(con, Id(schema="analysis",table="top_transcriptional_subtype"))
ts$sample_barcode <- sapply(strsplit(ts$aliquot_barcode,"-"),function(x)paste(x[1:4],collapse="-"))

ts_comp <- dat %>%
  group_by(sample_barcode) %>%
  summarise_if(is_numeric, mean) %>%
  select(-c(folder, timepoint, annotations, priority)) %>%
  inner_join(ts, by = "sample_barcode") %>%
  select(-c(aliquot_barcode, total, enrichment_score, p_value))

ts_comp %>%
  filter(grepl("-LU-",sample_barcode), signature_name == "Mesenchymal")
ts_comp %>%
  filter(grepl("-LU-",sample_barcode), signature_name == "Classical")
ts_comp %>%
  filter(grepl("-LU-",sample_barcode), signature_name == "Proneural")


ts_comp %>%
  filter(grepl("-HF-",sample_barcode), signature_name == "Mesenchymal")
ts_comp %>%
  filter(grepl("-HF-",sample_barcode), signature_name == "Classical")
ts_comp %>%
  filter(grepl("-HF-",sample_barcode), signature_name == "Proneural")

ts_plot <- ts_comp %>%
  mutate(cohort = substring(sample_barcode, 6,7),
         timepoint = substring(sample_barcode, 14, 15)) %>%
  mutate(timepoint = recode(timepoint, "R1" = "Recurrent", "R2" = "Recurrent", "TP" = "Primary", "R3" = "Recurrent")) %>%
  group_by(cohort, signature_name, timepoint) %>%
  summarise_if(is_numeric, mean) %>%
  pivot_longer(-c(cohort:timepoint), values_to = "percent", names_to = "feature")

ggplot(ts_plot, aes(fill=feature, y=percent, x=timepoint)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_manual(values=c("Total" = "#8B8B8B", "Subtotal" = "#ECECEC")) +
  facet_grid(signature_name~cohort) +
  labs(y = "Feature proportion (%)") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=7),
        axis.text.x=element_text(size=7,angle=45,hjust=1),
        axis.text.y=element_text(size=7),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_blank(),
        legend.text = element_text(size=7))

# Visualize the data with scatterplots
leading_edge <- full_comp_nonec$leading_edge + full_comp_nonec$histologically_normal
cellular_tumor <- full_comp_nonec$cellular_tumor
pathology <- c(leading_edge, cellular_tumor)

LE<- full_comp_nonec$LE
CT <- full_comp_nonec$CT + full_comp_nonec$CTpan + full_comp_nonec$CTmvp
cibersort <- c(LE, CT)

feature <- c(rep("Leading edge", nrow(full_comp_nonec)), rep("Cellular tumor", nrow(full_comp_nonec)))
plot_nonec <- data.frame(full_comp_nonec[,c("sample_barcode","timepoint","folder","pathologist")], pathology, cibersort, feature)

ggplot(plot_nonec %>% filter(folder=="SN",timepoint=="Primary"), aes(x=pathology,y=cibersort, colour=feature)) +
  geom_point(size=1) +
  geom_smooth(method='lm',se=FALSE) +
  facet_wrap(.~pathologist, scale="free") +
  labs(x = "Pathology fraction (%)", y = "CIBERSORT fraction (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=7,hjust=0.5),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=7),
        legend.position="none") 

# Create a simple correlation heatmap






#--------------------------
# Analyze which recurrent features are being detected

rec_comp <- full_comp %>% 
  filter(folder == "SN") %>%
  filter(timepoint == "Recurrent") %>%
  group_by(pathologist) %>%
  summarise(tn_le = cor(treatment_necrosis, LE), 
            dp_le = cor(depopulated_tumor, LE), 
            to_le = cor(treatment_necrosis, LE),
            tn_ct = cor(treatment_effec_other, CT), 
            dp_ct = cor(depopulated_tumor, CT), 
            to_ct = cor(treatment_necrosis, CT),
            tn_pan = cor(treatment_effec_other, CTpan), 
            dp_pan = cor(depopulated_tumor, CTpan), 
            to_pan = cor(treatment_effec_other, CTpan),
            tn_mvp = cor(treatment_necrosis, CTmvp), 
            dp_mvp = cor(depopulated_tumor, CTmvp), 
            to_mvp = cor(treatment_effec_other, CTmvp))


full_comp_nonec <- full_comp[,-which(colnames(full_comp)=="necrosis" | colnames(full_comp)=="treatment_necrosis" | colnames(full_comp)=="fibrosis")]


# Combine relevant recurrent features and get correlations:

lu_comp <- full_comp %>% 
  filter(folder == "SN") %>%
  #select(-c(fibrosis, depopulated_tumor, treatment_necrosis, treatment_effec_other)) %>%
  group_by(timepoint, pathologist) %>%
  summarise(le_cor = cor(LE, histologically_normal + leading_edge), nec_cor = cor(necrosis+depopulated_tumor, CTpan), ct_cor = cor(cellular_tumor+leading_edge+histologically_normal, CT+LE))

# Try removing non-cellular content
full_comp_nonec <- full_comp[,-which(colnames(full_comp)=="necrosis" | colnames(full_comp)=="treatment_necrosis" | colnames(full_comp)=="fibrosis")]
full_comp_nonec[,4:8] <- t(apply(full_comp_nonec[,4:8], 1, function(x)x/sum(x)))

full_comp_nonec %>%
  filter(folder == "SN") %>%
  group_by(pathologist,timepoint) %>%
  summarise(le_cor = cor(LE, histologically_normal + leading_edge,method="p"), ct_cor = cor(cellular_tumor, CT+CTpan+CTmvp,method="p")) %>%
  group_by(timepoint) %>%
  summarise(median(le_cor), median(ct_cor))

# Try removing non-cellular content and recurrence-specific features
full_comp_nonec <- full_comp %>%
  select(-c(necrosis, treatment_necrosis, fibrosis, depopulated_tumor, treatment_effec_other))
full_comp_nonec[,4:8] <- t(apply(full_comp_nonec[,4:6], 1, function(x)x/sum(x)))

full_comp_nonec %>%
  filter(folder == "SN") %>%
  group_by(pathologist,timepoint) %>%
  summarise(le_cor = cor(LE, histologically_normal + leading_edge,method="p"), ct_cor = cor(cellular_tumor, CT+CTpan+CTmvp,method="p")) %>%
  group_by(timepoint) %>%
  summarise(median(le_cor), median(ct_cor))

full_comp_nonec %>%
  filter(folder == "HF") %>%
  group_by(pathologist,timepoint) %>%
  summarise(le_cor = cor(LE, histologically_normal + leading_edge,method="s"), ct_cor = cor(cellular_tumor, CT+CTpan+CTmvp,method="s")) %>%
  ungroup(pathologist) %>%
  summarise(median(le_cor), median(ct_cor))



full_comp_nonec %>%
  filter(folder == "SN", timepoint == "Primary") %>%
  group_by(pathologist) %>%
  select_if(is.numeric) %>%
  group_map(~correlate(.x)) 


hf_comp <- full_comp %>% 
  filter(folder == "HF") %>%
  select(-c(fibrosis, depopulated_tumor, treatment_necrosis, treatment_effec_other)) %>%
  group_by(timepoint, pathologist) %>%
  summarise(le_cor = cor(LE, histologically_normal + leading_edge), nec_cor = cor(necrosis, CTpan), ct_cor = cor(cellular_tumor, CT))


