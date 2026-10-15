#finds the percent of cpgs per samplE(out of all the valid samples in the matrix) that have methylation above a certain level
library(bsseq)
library(GenomicRanges)
library(DSS)
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(dplyr)

rd= readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds")
str(rd)
dat <- getData(rd)
dim(dat)
colnames(dat)[1:20]
sample_ids <- rd@sample.ids
cov_cols <- paste0("coverage_", sample_ids)
cs_cols  <- paste0("numCs_", sample_ids)
ts_cols  <- paste0("numTs_", sample_ids)
str(cov_cols)
str(cs_cols)
# Genome-wide methylation for each sample
genome_meth <- sapply(seq_along(sample_ids), function(i) {
  sum(dat[[cs_cols[i]]], na.rm = TRUE) /
    sum(dat[[cov_cols[i]]], na.rm = TRUE)
})
# Put into a data frame
genome_meth_df <- data.frame(
  sample = sample_ids,
  methylation = genome_meth
)
genome_meth_df
dim(dat)
head(dat[, 1:15])
colnames(dat)[1:20]
cov_cols[1:5]
cs_cols[1:5]
cov_cols[1] %in% colnames(dat)
cs_cols[1] %in% colnames(dat)
sample_ids <- rd@sample.ids
genome_meth <- sapply(seq_along(sample_ids), function(i) {
  cov_col <- paste0("coverage", i)
  cs_col  <- paste0("numCs", i)
  sum(dat[[cs_col]], na.rm = TRUE) /
    sum(dat[[cov_col]], na.rm = TRUE)
})
genome_meth_df <- data.frame(
  sample = sample_ids,
  methylation = genome_meth
)
genome_meth_df$methylation_percent <-
  genome_meth_df$methylation * 100
genome_meth_df

fai <- fread(
  "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.fai",
  header = FALSE
)
colnames(fai)[1:2] <- c("chr", "length")
head(fai)
genome_length <- sum(fai$length)

sample_ids <- rd@sample.ids
thresholds <- c(0.75, 0.85, 0.90)
results <- lapply(seq_along(sample_ids), function(i) {
  cov <- dat[[paste0("coverage", i)]]
  cs  <- dat[[paste0("numCs", i)]]
  # Methylation at each CpG
  meth <- cs / cov
  # Only CpGs with valid coverage
  valid <- !is.na(meth) & cov > 0
  data.frame(
    sample = sample_ids[i],
    n_CpGs = sum(valid),
    n_75 = sum(meth[valid] >= 0.75),
    n_85 = sum(meth[valid] >= 0.85),
    n_90 = sum(meth[valid] >= 0.90),
    pct_75 = sum(meth[valid] >= 0.75) / sum(valid) * 100,
    pct_85 = sum(meth[valid] >= 0.85) / sum(valid) * 100,
    pct_90 = sum(meth[valid] >= 0.90) / sum(valid) * 100
  )
})
genome_meth_thresholds <- rbindlist(results)
genome_meth_thresholds
write.csv(
  genome_meth_thresholds,
  file = "/scratch/project_2010912/ezel/salmon/data_files/genome_meth_thresholds.csv",
  row.names = FALSE
)