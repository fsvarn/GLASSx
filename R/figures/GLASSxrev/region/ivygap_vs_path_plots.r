#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Final pathologist vs CIBERSORT plots
# Includes Leeds samples with features grouped based on periphery/core and recurrence
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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
dat <- dat %>% filter(folder == "LU")

# Transcriptomic analysis

# Read in IvyGAP table
ivy <- dbReadTable(con, Id(schema="analysis",table="cibersortx_ivygap"))
ivy$sample_barcode <- sapply(strsplit(ivy$aliquot_barcode,"-"),function(x)paste(x[1:4],collapse="-"))
ivy$fraction <- ivy$fraction 

ivy <- ivy %>%
  pivot_wider(names_from = "cell_state", values_from = "fraction") %>%
  group_by(aliquot_barcode) %>%
  mutate(CTtotal = sum(CT, CTpan, CTmvp)) %>%
  pivot_longer(-c(aliquot_barcode, sample_barcode), values_to = "cibersort_percent", names_to = "cibersort_feature")


#--------------------------

# Estimate percent histologically normal for the response letter
histnorm_estimates <- dat %>%
  select(-c(annotations, priority)) %>%
  pivot_longer(-c(sample_barcode, timepoint, folder, pathologist), names_to = "path_feature", values_to = "path_percent") %>%
  filter(!(path_feature %in% c("necrosis", "fibrosis", "treatment_necrosis", "total"))) %>%
  group_by(sample_barcode, timepoint, folder, pathologist) %>%
  mutate(total = sum(path_percent)) %>%
  ungroup() %>%
  mutate(path_percent = path_percent/total) %>%
  select(-total) %>%
  group_by(folder, path_feature, sample_barcode) %>% 
  summarise(mean = mean(path_percent)) %>%
  filter(folder == "LU", path_feature == "histologically_normal")

# Subset to cellular regions, combine recurrence-specific features
long_dat <- dat %>%
  mutate(periphery = leading_edge + histologically_normal) %>%
  select(-c(annotations, priority, histologically_normal, leading_edge)) %>%
  pivot_longer(-c(sample_barcode, timepoint, folder, pathologist), names_to = "path_feature", values_to = "path_percent") %>%
  filter(!(path_feature %in% c("necrosis", "fibrosis", "treatment_necrosis", "total"))) %>%
  group_by(sample_barcode, timepoint, folder, pathologist) %>%
  mutate(total = sum(path_percent)) %>%
  ungroup() %>%
  mutate(path_percent = path_percent/total) %>%
  select(-total)

# Add summary rows
rec_feat <- long_dat %>%
  filter(path_feature %in% c("depopulated_tumor", "treatment_effect_other"), timepoint == "Recurrent") %>%
  group_by(sample_barcode, timepoint, folder, pathologist) %>%
  summarise(path_percent = sum(path_percent)) %>%
  mutate(path_feature = "total_rec_feat") %>%
  relocate(path_feature, .before = path_percent)

rec_ct <- long_dat %>%
  filter(path_feature %in% c("depopulated_tumor", "treatment_effect_other", "cellular_tumor"), timepoint == "Recurrent") %>%
  group_by(sample_barcode, timepoint, folder, pathologist) %>%
  summarise(path_percent = sum(path_percent)) %>%
  mutate(path_feature = "ct_and_rec") %>%
  relocate(path_feature, .before = path_percent)

long_dat <- long_dat %>%
  bind_rows(rec_feat) %>%
  bind_rows(rec_ct)

# Create a mean pathologist value and correlate everything against everything (use this, easier to understand)

mean_pathologist <- long_dat %>%
  inner_join(ivy, by = "sample_barcode") %>%
  select(-aliquot_barcode) %>%
  group_by(sample_barcode, folder, timepoint, path_feature, cibersort_feature) %>%
  summarise(path_percent = mean(path_percent), cibersort_percent = mean(cibersort_percent)) %>%
  ungroup() %>%
  group_by(folder, timepoint, path_feature, cibersort_feature) %>%
  summarise(cor = cor(path_percent, cibersort_percent)) %>%
  filter(!is.na(cor)) 

# Select relevant values to plot
plot_cor <- mean_pathologist %>%
  filter(cibersort_feature %in% c("LE", "CTtotal")) %>%
  filter(path_feature %in% c("periphery","cellular_tumor","total_rec_feat", "ct_and_rec")) %>%
  mutate(cibersort_feature = recode(cibersort_feature, "CTtotal" = "Tumor core", "LE" = "Periphery")) %>%
  mutate(path_feature = recode(path_feature, "cellular_tumor" = "Cellular tumor", "periphery" = "Periphery",
                               "total_rec_feat" = "Treatment effect", "ct_and_rec" = "Tumor core")) %>%
  mutate(path_feature = as_factor(path_feature)) %>%
  mutate(path_feature = fct_relevel(path_feature, "Periphery", "Tumor core", "Treatment effect"))

# Simplify for the paper
plot_cor_simple <- plot_cor %>%
  filter(!(timepoint=="Recurrent" & (path_feature %in% c("Cellular tumor", "Treatment effect")))) %>%
  mutate(path_feature = recode(path_feature, "Cellular tumor" = "Tumor core"))

pdf("figures/pathology/mean_pathologist_cor_LU_v2.pdf",height = 1.5, width = 2.7)
ggplot(plot_cor_simple %>% filter(folder=="LU"), aes(x = path_feature, y = cibersort_feature, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2("R",low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
  facet_grid(. ~ timepoint, scale = "free_x", space= "free_x") +
  labs(x = "Pathology feature", y = "IvyGAP feature") +
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Breakdown of recurrence-specific features vs IvyGAP

rec_cor <- mean_pathologist %>%
  filter(timepoint == "Recurrent") %>%
  filter(cibersort_feature != "CTtotal") %>%
  filter(!path_feature %in% c("total_rec_feat", "ct_and_rec")) %>%
  #  mutate(cibersort_feature = recode(cibersort_feature, "CT" = "Cellular tumor", "LE" = "Leading edge")) %>%
  mutate(path_feature = recode(path_feature, "cellular_tumor" = "Cellular tumor", "periphery" = "Periphery",
                               "depopulated_tumor" = "Depopulated tumor", "treatment_effect_other" = "Treatment effect (other)")) %>%
  mutate(path_feature = as_factor(path_feature)) %>%
  mutate(path_feature = fct_relevel(path_feature, "Periphery", "Cellular tumor", "Depopulated tumor", "Treatment effect (other)")) %>%
  mutate(cibersort_feature = fct_relevel(cibersort_feature, "CTmvp", "CTpan", "CT", "LE"))


pdf("figures/pathology/mean_pathologist_rec_correlation_subset_LU.pdf",height = 2, width = 1.8)
ggplot(rec_cor %>% filter(folder=="LU", path_feature %in% c("Depopulated tumor", "Treatment effect (other)")), aes(x = path_feature, y = cibersort_feature, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2("R",low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
  facet_grid(. ~ timepoint, scale = "free_x", space= "free_x") +
  labs(x = "Pathology feature", y = "CIBERSORT feature") +
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Add the necrosis plots

nec_dat <- dat %>%
  select(c(sample_barcode, timepoint, folder, necrosis, fibrosis, treatment_necrosis, pathologist)) %>%
  pivot_longer(-c(sample_barcode, timepoint, folder, pathologist), names_to = "path_feature", values_to = "path_percent")

# Correlate necrosis against features
mean_pathologist_cor_nec <- nec_dat %>%
  inner_join(ivy, by = "sample_barcode") %>%
  select(-aliquot_barcode) %>%
  group_by(sample_barcode, folder, timepoint, path_feature, cibersort_feature) %>%
  summarise(path_percent = mean(path_percent), cibersort_percent = mean(cibersort_percent)) %>%
  ungroup() %>%
  group_by(folder, timepoint, path_feature, cibersort_feature) %>%
  summarise(cor = cor(path_percent, cibersort_percent)) %>%
  filter(!is.na(cor)) 

plot_nec <- nec_dat %>%
  filter(folder == "LU") %>%
  inner_join(ivy, by = "sample_barcode") %>%
  filter(path_feature == "necrosis", cibersort_feature == "CTpan") %>%
  select(sample_barcode, timepoint, folder, pathologist, path_feature, path_percent, cibersort_feature, cibersort_percent) %>%
  group_by(sample_barcode, timepoint, folder, path_feature, cibersort_feature) %>%
  summarise(path_percent = mean(path_percent), cibersort_percent = mean(cibersort_percent))

pdf("figures/pathology/mean_pathologist_necrosis_cor_legend.pdf",height = 1, width = 1)
ggplot(data = plot_nec, aes(x = path_percent, y = cibersort_percent*100)) +
  geom_point(aes(shape = timepoint), size = 1) +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Necrosis % (Pathology)", y = "CTpan % (CIBERSORT)") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.text=element_text(size=7),
        legend.title = element_blank())
dev.off()

pdf("figures/pathology/mean_pathologist_necrosis_cor.pdf",height = 1, width = 1)
ggplot(data = plot_nec, aes(x = path_percent, y = cibersort_percent*100)) +
  geom_point(aes(shape = timepoint), size = 1) +
  geom_smooth(method="lm",se = FALSE) +
  labs(x = "Necrosis % (Pathology)", y = "CTpan % (CIBERSORT)") +
  theme_classic() +
  theme(axis.text = element_text(size=7),
        axis.title= element_blank(),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.position = "none") +
  scale_y_continuous(labels = scales::number_format(accuracy = 1), limits = c(0, NA))
dev.off()

