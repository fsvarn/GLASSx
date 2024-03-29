# Used this script to look at the distribution of mean coverage across all WXS samples when adding the H2 cohort
library(tidyverse)

## Load WGS metrics
wgsmetricsf = list.files("results/align/wgsmetrics", recursive = T, pattern = "WgsMetrics.txt", full.names = T)
wgsmetrics = wgsmetricsf %>%
  map(function(f) data.table::fread(sprintf("sed -n '/METRICS CLASS/,/HISTOGRAM/p' %s | sed '/^#/d'", f)))
names(wgsmetrics) = gsub(".WgsMetrics.txt","",basename(wgsmetricsf))
wgsmetrics = wgsmetrics %>% 
  map2_df(names(wgsmetrics), ~mutate(.x, aliquot_id = .y)) 

## Load duplication metrics
dupmetricsf = list.files("results/align/markduplicates", recursive = T, pattern = "metrics.txt", full.names = T)
dupmetrics = dupmetricsf %>%
  map(function(f) data.table::fread(sprintf("sed -n '/METRICS CLASS/,/HISTOGRAM/p' %s | sed '/^#/d'", f), colClasses = c("ESTIMATED_LIBRARY_SIZE"="int")))
names(dupmetrics) = gsub(".metrics.txt","",basename(dupmetricsf))
dupmetrics = dupmetrics %>% 
  map2_df(names(dupmetrics), ~mutate(.x, aliquot_id = .y)) 

## Calculate summaries
metrics = dupmetrics %>% 
  full_join(wgsmetrics) %>%
  mutate(SAMPLE_TYPE = substr(aliquot_id,14,15),
         TOTAL_READS = 2*READ_PAIRS_EXAMINED + UNPAIRED_READS_EXAMINED + UNMAPPED_READS + SECONDARY_OR_SUPPLEMENTARY_RDS,
         PERCENT_MAPPED = (UNPAIRED_READS_EXAMINED + 2*READ_PAIRS_EXAMINED) / (2*READ_PAIRS_EXAMINED + UNPAIRED_READS_EXAMINED + UNMAPPED_READS),
         GB_PER_1X_COVERAGE_FFPE = ((TOTAL_READS / MEDIAN_COVERAGE) * 150)/(1E9),
         COVERAGE_REMAIN_30X = 30 - MEDIAN_COVERAGE,
         COVERAGE_REMAIN_40X = 40 - MEDIAN_COVERAGE,
         COVERAGE_REMAIN_50X = 50 - MEDIAN_COVERAGE,
         COVERAGE_REMAIN_60X = 60 - MEDIAN_COVERAGE,
         GB_TO_FFPE_30X = COVERAGE_REMAIN_30X * GB_PER_1X_COVERAGE_FFPE,
         GB_TO_FFPE_40X = COVERAGE_REMAIN_40X * GB_PER_1X_COVERAGE_FFPE,
         GB_TO_FFPE_50X = COVERAGE_REMAIN_50X * GB_PER_1X_COVERAGE_FFPE,
         GB_TO_FFPE_60X = COVERAGE_REMAIN_60X * GB_PER_1X_COVERAGE_FFPE,
         GB_TO_FROZEN_30X = (3.3) * COVERAGE_REMAIN_30X / 0.85,
         GB_TO_FROZEN_40X = (3.3) * COVERAGE_REMAIN_40X / 0.85,
         GB_TO_FROZEN_50X = (3.3) * COVERAGE_REMAIN_50X / 0.85,
         GB_TO_FROZEN_60X = (3.3) * COVERAGE_REMAIN_60X / 0.85) %>%
  select(aliquot_id, SAMPLE_TYPE,
         LIBRARY, TOTAL_READS, MEAN_COVERAGE, PERCENT_DUPLICATION, PERCENT_MAPPED,
         COVERAGE_REMAIN_30X, COVERAGE_REMAIN_40X, COVERAGE_REMAIN_50X, COVERAGE_REMAIN_60X,
         GB_TO_FFPE_30X, GB_TO_FFPE_40X, GB_TO_FFPE_50X, GB_TO_FFPE_60X,
         GB_TO_FROZEN_30X, GB_TO_FROZEN_40X, GB_TO_FROZEN_50X, GB_TO_FROZEN_60X)

frame_metrics <- data.frame(metrics)

#Look at WXS mean_coverage metrics
wxs_metrics <- frame_metrics[grep("-WXS-",frame_metrics[,"aliquot_id"]),]
# Remove WXS samples that are probably WGS samples (would not expect 30x coverage from DF cohort)
wxs_metrics <- wxs_metrics[-which(wxs_metrics[,"MEAN_COVERAGE"] > 30),]
# Remove NA coverage WXS
wxs_metrics <- wxs_metrics[-which(is.na(wxs_metrics[,"MEAN_COVERAGE"])),]

# Set thresholds
mean_coverage <- mean(wxs_metrics[,"MEAN_COVERAGE"])
sd_coverage <- sd(wxs_metrics[,"MEAN_COVERAGE"])
min_threshold <- mean_coverage - 2*sd_coverage

# Identify failed coverage
coverage_fail <- wxs_metrics %>%
				 filter(MEAN_COVERAGE < min_threshold)

round(sum(metrics$GB_TO_FFPE_30X[metrics$GB_TO_FFPE_30X > 0]))
round(sum(metrics$GB_TO_FFPE_40X[metrics$GB_TO_FFPE_40X > 0]))
round(sum(metrics$GB_TO_FFPE_50X[metrics$GB_TO_FFPE_50X > 0]))
round(sum(metrics$GB_TO_FFPE_60X[metrics$GB_TO_FFPE_60X > 0]))

round(sum(metrics$GB_TO_FROZEN_30X[metrics$GB_TO_FROZEN_30X > 0]))
round(sum(metrics$GB_TO_FROZEN_40X[metrics$GB_TO_FROZEN_40X > 0]))
round(sum(metrics$GB_TO_FROZEN_50X[metrics$GB_TO_FROZEN_50X > 0]))
round(sum(metrics$GB_TO_FROZEN_60X[metrics$GB_TO_FROZEN_60X > 0]))
