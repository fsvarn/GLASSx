library(DBI)
library(tidyverse)

con <- DBI::dbConnect(odbc::odbc(), "GLASSv3")

segfiles <- list.files("results/cnv/callsegments", full.names = TRUE)
segfiles <- segfiles[grep(".called.seg",segfiles)]


segs <- parallel::mclapply(segfiles, function(f) {
  dat <- read.delim(f, comment.char = "@", as.is= TRUE)
  dat <- dat %>%
    mutate(aliquot_barcode = substr(basename(f),1,30), pos = sprintf("[%s,%s]", START, END)) %>%
    select(aliquot_barcode, chrom = CONTIG, pos, num_points = NUM_POINTS_COPY_RATIO, log2_copy_ratio = MEAN_LOG2_COPY_RATIO, cnv_call = CALL)
  return(dat)
}, mc.cores = 8)
segs <- data.table::rbindlist(segs) %>% as.data.frame()

# Add the new HF cohort
# hf_segs <- segs[grep("GLSS-HF-",segs[,"aliquot_barcode"]),]
# hf_segs <- hf_segs[-grep("GLSS-HF-2",hf_segs[,"aliquot_barcode"]),]
# hf_segs <- hf_segs[-grep("GLSS-HF-3",hf_segs[,"aliquot_barcode"]),]
# 
# hf_segs[which(hf_segs[,"chrom"] == "X"),"chrom"] <- 23
# hf_segs[which(hf_segs[,"chrom"] == "Y"),"chrom"] <- 24

# Add the new SN cohort
sn_segs <- segs[grep("GLSS-SN-",segs[,"aliquot_barcode"]),]

sn_segs[which(sn_segs[,"chrom"] == "X"),"chrom"] <- 23
sn_segs[which(sn_segs[,"chrom"] == "Y"),"chrom"] <- 24

# Check for conflicts with aliquots table
q <-"SELECT * FROM biospecimen.aliquots WHERE aliquot_batch = 'GLSS-SN-WGS'"
aliquots <- dbGetQuery(con, q)

# Two aliquots were removed early in the process due to high contamination issues
sn_segs_final <- sn_segs[which(sn_segs[,"aliquot_barcode"] %in% aliquots[,"aliquot_barcode"]),]

dbWriteTable(con, Id(schema="variants",table="gatk_seg"), sn_segs_final, append=TRUE)
