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
mean_estimates <- dat %>%
  select(-c(annotations, priority)) %>%
  pivot_longer(-c(sample_barcode, timepoint, folder, pathologist), names_to = "path_feature", values_to = "path_percent") %>%
  filter(!(path_feature %in% c("necrosis", "fibrosis", "treatment_necrosis", "total"))) %>%
  filter(!(path_feature=="treatment_effect_other" & path_percent > 0 & timepoint == "Primary")) %>% # Effects two HF samples where small levels of recurrence-specific feature was annotated in primary
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
            mutate(leading_edge = leading_edge + histologically_normal) %>%
            select(-c(annotations, priority, histologically_normal)) %>%
            pivot_longer(-c(sample_barcode, timepoint, folder, pathologist), names_to = "path_feature", values_to = "path_percent") %>%
            filter(!(path_feature %in% c("necrosis", "fibrosis", "treatment_necrosis", "total"))) %>%
            filter(!(path_feature=="treatment_effect_other" & path_percent > 0 & timepoint == "Primary")) %>% # Effects two HF samples where small levels of recurrence-specific feature was annotated in primary
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

# Correlate everything against everything
cor_dat <- long_dat %>%
           inner_join(ivy, by = "sample_barcode") %>%
           select(-aliquot_barcode) %>%
           group_by(pathologist, folder, timepoint, path_feature, cibersort_feature) %>%
           summarise(cor = cor(path_percent, cibersort_percent)) %>%
           filter(!is.na(cor))

mean_cor <- cor_dat %>%
            group_by(folder, timepoint, path_feature, cibersort_feature) %>%
            summarise(cor = mean(cor)) 
median_cor <- cor_dat %>%
  group_by(folder, timepoint, path_feature, cibersort_feature) %>%
  summarise(cor = median(cor)) 

# Select relevant values to plot
plot_cor <- median_cor %>%
            filter(cibersort_feature %in% c("LE", "CTtotal")) %>%
            filter(path_feature %in% c("leading_edge","cellular_tumor","total_rec_feat", "ct_and_rec")) %>%
            mutate(cibersort_feature = recode(cibersort_feature, "CTtotal" = "Cellular tumor (all)", "LE" = "Leading edge")) %>%
            mutate(path_feature = recode(path_feature, "cellular_tumor" = "Cellular tumor", "leading_edge" = "Leading edge",
                                                       "total_rec_feat" = "Treatment effect", "ct_and_rec" = "Cellular tumor + treatment effect")) %>%
            mutate(path_feature = as_factor(path_feature)) %>%
            mutate(path_feature = fct_relevel(path_feature, "Leading edge", "Cellular tumor", "Treatment effect"))

# all cohorts
ggplot(plot_cor, aes(x = path_feature, y = cibersort_feature, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2("R",low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
  facet_grid(folder ~ timepoint, scale = "free_x", space= "free_x") +
  labs(x = "Pathology features (% area on slide)", y = "CIBERSORT features (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title = element_text(size=7),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.15, 'inch')) 

# Show full breakdown to the reviewer
pdf("figures/pathology/median_correlation_rev_LU.pdf",height = 1.7, width = 3.25)
ggplot(plot_cor %>% filter(folder=="LU", path_feature != "Cellular tumor + treatment effect"), aes(x = path_feature, y = cibersort_feature, fill=cor)) +
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

# Simplify for the paper
plot_cor_simple <- plot_cor %>%
                   filter(!(timepoint=="Recurrent" & (path_feature %in% c("Cellular tumor", "Treatment effect"))))

pdf("figures/pathology/median_correlation_LU.pdf",height = 2.2, width = 3)
ggplot(plot_cor_simple %>% filter(folder=="LU"), aes(x = path_feature, y = cibersort_feature, fill=cor)) +
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

# Breakdown of recurrence-specific features vs IvyGAP

rec_cor <- median_cor %>%
  filter(timepoint == "Recurrent") %>%
  filter(cibersort_feature != "CTtotal") %>%
  filter(!path_feature %in% c("total_rec_feat", "ct_and_rec")) %>%
#  mutate(cibersort_feature = recode(cibersort_feature, "CT" = "Cellular tumor", "LE" = "Leading edge")) %>%
  mutate(path_feature = recode(path_feature, "cellular_tumor" = "Cellular tumor", "leading_edge" = "Leading edge",
                               "depopulated_tumor" = "Depopulated tumor", "treatment_effect_other" = "Treatment effect (other)")) %>%
  mutate(path_feature = as_factor(path_feature)) %>%
  mutate(path_feature = fct_relevel(path_feature, "Leading edge", "Cellular tumor", "Depopulated tumor", "Treatment effect (other)")) %>%
  mutate(cibersort_feature = fct_relevel(cibersort_feature, "CTmvp", "CTpan", "CT", "LE"))

ggplot(rec_cor , aes(x = path_feature , y = cibersort_feature, fill=cor)) +
  geom_tile() +
  scale_fill_gradient2("R",low ="royalblue4", mid = "white", high = "tomato3", midpoint = 0, limits = c(-1,1)) +
  facet_grid(folder ~ ., scale = "free_x", space= "free_x") +
  labs(x = "Pathology features (% area on slide)", y = "CIBERSORT features (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(size=7), axis.text.x= element_text(size=7,angle=45,hjust=1),
        axis.title = element_text(size=7),
        panel.grid.major=element_blank(),#panel.grid.minor=element_blank(),
        strip.text = element_text(size=7),
        strip.background = element_rect(colour="white",fill="white"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=7), 
        legend.key.size = unit(0.15, 'inch')) 

pdf("figures/pathology/median_rec_correlation_LU.pdf",height = 2.2, width = 2.3)
ggplot(rec_cor %>% filter(folder=="LU"), aes(x = path_feature, y = cibersort_feature, fill=cor)) +
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

pdf("figures/pathology/median_rec_correlation_subset_LU.pdf",height = 2, width = 2)
ggplot(rec_cor %>% filter(folder=="LU", cibersort_feature != "LE", path_feature != "Leading edge"), aes(x = path_feature, y = cibersort_feature, fill=cor)) +
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



