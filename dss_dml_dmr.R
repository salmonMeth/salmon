###
library(bsseq)
library(GenomicRanges)
library(DSS)
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(dplyr)


path_db_normalized="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_meth_normalized.rds"
db_norm=readRDS(path_db_normalized)

df <- getData(db_norm)
Cov <- as.matrix(df[, db_norm@coverage.index])
M <- as.matrix(df[, db_norm@numCs.index])
treatment <- db_norm@treatment


rownames(M) <- NULL
rownames(Cov) <- NULL
colnames(M) <- db_norm@sample.ids
colnames(Cov) <- db_norm@sample.ids

keep <- rowSums(is.na(M)) == 0 &
  rowSums(is.na(Cov)) == 0
M <- M[keep, , drop = FALSE]
Cov <- Cov[keep, , drop = FALSE]


gr <- GRanges(
  seqnames = df$chr,
  ranges = IRanges(df$start, df$end)
)


BSobj <- BSseq(
  gr = gr,
  M = M,
  Cov = Cov,
  sampleNames = db_norm@sample.ids
)

cat("Starting DMLtest 1\n")
flush.console()
dmlTest_no_smoothing <- DMLtest(
  BSobj,
  group1 = which(treatment == 0),
  group2 = which(treatment == 1)
)

cat("Starting DMLtest 2\n")
flush.console()
dmlTest_smoothing <- DMLtest(
  BSobj,
  group1 = which(treatment == 0),
  group2 = which(treatment == 1),smoothing=TRUE)

# Save the DML test object
saveRDS(dmlTest_no_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/DMLtest_no_smoothing.rds")
saveRDS(dmlTest_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/DMLtest_smoothing.rds")


cat("Starting callDMR 1\n")
flush.console()
dmrs_no_smoothing <- callDMR(
  dmlTest_no_smoothing, delta=0.1
)


cat("Starting callDMR 2\n")
flush.console() 

dmrs_smoothing <- callDMR(
  dmlTest_smoothing, delta=0.1
)

saveRDS(dmrs_no_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/DSS_DMRs_no_smoothing.rds")
write.csv(dmrs_no_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/DSS_DMRs_no_smoothing.csv", row.names = FALSE)
#####
saveRDS(dmrs_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/DSS_DMRs_smoothing.rds")
write.csv(dmrs_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/DSS_DMRs_smoothing.csv", row.names = FALSE)

#######################################
#######################################
#compare the results for cpgs

library(data.table)
test_smooth=readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/DMLtest_smoothing.rds")
#get the dss results
test_smooth_fdr=test_smooth[test_smooth$fdr<.01,]
test_diff_10=test_smooth_fdr[abs(test_smooth_fdr$diff) >.1,]
#methlkit results
methylKit_normalized_overdisper_diff_10=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_normalized_with_overdispersion_difference_10.csv")
####
mk <- as.data.table(methylKit_normalized_overdisper_diff_10)
dss <- as.data.table(test_diff_10)

library(data.table)

mk <- as.data.table(methylKit_normalized_overdisper_diff_10)
dss <- as.data.table(test_diff_10)

# Create a common coordinate column
mk[, coord := paste(chr, start, sep = ":")]
dss[, coord := paste(chr, pos, sep = ":")]

# Shared CpGs
common <- intersect(mk$coord, dss$coord)

# Unique to each method
mk_only <- setdiff(mk$coord, dss$coord)
dss_only <- setdiff(dss$coord, mk$coord)

cat("methylKit:", nrow(mk), "\n")
cat("DSS:", nrow(dss), "\n")
cat("Common:", length(common), "\n")
cat("methylKit only:", length(mk_only), "\n")
cat("DSS only:", length(dss_only), "\n")

common_chr <- mk_common[, .N, by = chr][order(chr)]

mk_unique_chr <- mk_unique[, .N, by = chr][order(chr)]

dss_unique_chr <- dss_unique[, .N, by = chr][order(chr)]

mk_common  <- mk[coord %in% common]
dss_common <- dss[coord %in% common]

mk_unique  <- mk[coord %in% mk_only]
dss_unique <- dss[coord %in% dss_only]

# Create a common coordinate column
mk[, coord := paste(chr, start, sep = ":")]
dss[, coord := paste(chr, pos, sep = ":")]

# Shared CpGs
common <- intersect(mk$coord, dss$coord)

# Unique to each method
mk_only <- setdiff(mk$coord, dss$coord)
dss_only <- setdiff(dss$coord, mk$coord)

cat("methylKit:", nrow(mk), "\n")
cat("DSS:", nrow(dss), "\n")
cat("Common:", length(common), "\n")
cat("methylKit only:", length(mk_only), "\n")
cat("DSS only:", length(dss_only), "\n")

library(data.table)

common_chr <- mk_common[, .(common = .N), by = chr]

mk_unique_chr <- mk_unique[, .(methylKit_only = .N), by = chr]

dss_unique_chr <- dss_unique[, .(DSS_only = .N), by = chr]

chr_distribution <- Reduce(
  function(x, y) merge(x, y, by = "chr", all = TRUE),
  list(common_chr, mk_unique_chr, dss_unique_chr)
)

chr_distribution[is.na(chr_distribution)] <- 0

chr_distribution <- chr_distribution[order(chr)]

chr_distribution

##
library(data.table)

dss_common[, diff_percent := diff * 100]

comparison <- merge(
  mk_common,
  dss_common[, .(
    coord,
    diff_percent,
    phi1,
    phi2,
    pval,
    fdr
  )],
  by = "coord"
)

comparison[, diff_error := meth.diff - diff_percent]
comparison[, abs_diff_error := abs(diff_error)]

library(data.table)
dat <- getData(db_norm)
Cov <- as.matrix(dat[, db_norm@coverage.index])
coords <- paste(dat$chr, dat$start, sep = ":")
idx <- match(comparison$coord, coords)
Cov_common <- Cov[idx, ]
comparison[, meanCov := rowMeans(Cov_common, na.rm = TRUE)]
comparison[, abs_disagreement :=
             abs((-meth.diff) - diff_percent)]
plot(
  log10(comparison$meanCov),
  comparison$abs_disagreement,
  pch = 16,
  cex = 0.7,
  xlab = "log10(mean coverage)",
  ylab = "Absolute disagreement (%)"
)
#############segs
diff_segs=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_segments_with_overdispersion_difference_10.csv")

library(GenomicRanges)

## Segment ranges
seg_gr <- GRanges(
  seqnames = diff_segs$chr,
  ranges = IRanges(diff_segs$start, diff_segs$end)
)

## methylKit CpGs
mk_gr <- GRanges(
  seqnames = methylKit_normalized_overdisper_diff_10$chr,
  ranges = IRanges(
    methylKit_normalized_overdisper_diff_10$start,
    methylKit_normalized_overdisper_diff_10$start
  )
)

## DSS CpGs
dss_gr <- GRanges(
  seqnames = test_diff_10$chr,
  ranges = IRanges(
    test_diff_10$pos,
    test_diff_10$pos
  )
)

## Common CpGs
common_gr <- GRanges(
  seqnames = comparison$chr,
  ranges = IRanges(
    comparison$start,
    comparison$start
  )
)

mk_hits <- findOverlaps(mk_gr, seg_gr)
dss_hits <- findOverlaps(dss_gr, seg_gr)
common_hits <- findOverlaps(common_gr, seg_gr)

n_mk_in_seg <- length(unique(queryHits(mk_hits)))
n_dss_in_seg <- length(unique(queryHits(dss_hits)))
n_common_in_seg <- length(unique(queryHits(common_hits)))

n_mk_in_seg
n_dss_in_seg
n_common_in_seg
100 * n_mk_in_seg / length(mk_gr)

100 * n_dss_in_seg / length(dss_gr)

100 * n_common_in_seg / length(common_gr)
mk_per_seg <- tabulate(subjectHits(mk_hits), nbins = length(seg_gr))
dss_per_seg <- tabulate(subjectHits(dss_hits), nbins = length(seg_gr))
common_per_seg <- tabulate(subjectHits(common_hits), nbins = length(seg_gr))

summary(mk_per_seg)
summary(dss_per_seg)
summary(common_per_seg)


######getting the DSS DMR
diff_seg_dss=read.csv("/scratch/project_2010912/ezel/salmon/data_files/DSS_DMRs_smoothing.csv")
library(GenomicRanges)

## methylKit differential segments
mk_seg_gr <- GRanges(
  seqnames = diff_segs$chr,
  ranges = IRanges(diff_segs$start, diff_segs$end)
)

## DSS DMRs
dss_seg_gr <- GRanges(
  seqnames = diff_seg_dss$chr,
  ranges = IRanges(diff_seg_dss$start, diff_seg_dss$end)
)
hits <- findOverlaps(mk_seg_gr, dss_seg_gr)

library(GenomicRanges)

mk_idx  <- queryHits(hits)
dss_idx <- subjectHits(hits)

# Overlap intervals
ov <- pintersect(
  mk_seg_gr[mk_idx],
  dss_seg_gr[dss_idx]
)

# Lengths
overlap_bp <- width(ov)
mk_len  <- width(mk_seg_gr[mk_idx])
dss_len <- width(dss_seg_gr[dss_idx])

# Union length
union_bp <- mk_len + dss_len - overlap_bp

# Summary table
overlap_summary <- data.frame(
  mk_chr   = as.character(seqnames(mk_seg_gr[mk_idx])),
  mk_start = start(mk_seg_gr[mk_idx]),
  mk_end   = end(mk_seg_gr[mk_idx]),
  mk_len   = mk_len,
  
  dss_start = start(dss_seg_gr[dss_idx]),
  dss_end   = end(dss_seg_gr[dss_idx]),
  dss_len   = dss_len,
  
  overlap_bp = overlap_bp,
  
  pct_mk  = 100 * overlap_bp / mk_len,
  pct_dss = 100 * overlap_bp / dss_len,
  
  jaccard = overlap_bp / union_bp
)

overlap_summary

##compare the segments found
diff_segs$length <- diff_segs$end - diff_segs$start + 1

###################
#############
path_segs= "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_segments.rds"
all_segments= readRDS(path_segs)
#check the distribution of NA values(from the OG data) that fall into the methylkit segments
df=readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds")
dat <- getData(df)
library(GenomicRanges)
cpg_gr <- GRanges(
  seqnames = dat$chr,
  ranges = IRanges(dat$start, dat$start)
)

####find the cpgs with NA values
cov_cols <- db@coverage.index
cs_cols  <- db@numCs.index
ts_cols  <- db@numTs.index

has_na <-
  rowSums(is.na(dat[, cov_cols])) +
  rowSums(is.na(dat[, cs_cols])) +
  rowSums(is.na(dat[, ts_cols])) > 0

mcols(cpg_gr)$has_na <- has_na


all_seg_gr <- GRanges(
  seqnames = all_segments$chr,
  ranges = IRanges(
    start = all_segments$start,
    end = all_segments$end
  )
)

na_samples_per_cpg <- rowSums(is.na(dat[, cov_cols]))
mcols(cpg_gr)$na_samples <- na_samples_per_cpg

#we want to get more detailed stats regarding the NA values in the segments
library(GenomicRanges)
library(data.table)
#change accordingly
na_prop_threshold <- 0.10
n_samples <- length(db@sample.ids)
na_sample_threshold <- ceiling(na_prop_threshold * n_samples)

cov_cols <- db@coverage.index

na_samples_per_cpg <- rowSums(is.na(dat[, cov_cols]))

mcols(cpg_gr)$na_samples <- na_samples_per_cpg
mcols(cpg_gr)$affected <-
  na_samples_per_cpg >= na_sample_threshold

hits <- findOverlaps(all_seg_gr, cpg_gr)

seg_list <- split(subjectHits(hits), queryHits(hits))

##

segment_summary <- data.table(
  chr = as.character(seqnames(all_seg_gr)),
  start = start(all_seg_gr),
  end = end(all_seg_gr)
)

segment_summary[, segment_bp := width(all_seg_gr)]

segment_summary[, `:=`(
  n_cpg = 0L,
  affected_cpg = 0L,
  prop_affected_cpg = NA_real_,
  missing_sample_values = 0L,
  normalized_missing = NA_real_,
  max_missing_samples = 0L,
  median_missing_samples = NA_real_
)]

for(i in seq_along(seg_list)){
  
  idx <- seg_list[[i]]
  
  n_cpg <- length(idx)
  
  miss <- mcols(cpg_gr)$na_samples[idx]
  
  segment_summary$n_cpg[i] <- n_cpg
  
  segment_summary$affected_cpg[i] <-
    sum(miss >= na_sample_threshold)
  
  segment_summary$prop_affected_cpg[i] <-
    sum(miss >= na_sample_threshold) / n_cpg
  
  segment_summary$missing_sample_values[i] <-
    sum(miss)
  
  segment_summary$normalized_missing[i] <-
    sum(miss) / (n_cpg * n_samples)
  
  segment_summary$max_missing_samples[i] <-
    max(miss)
  
  segment_summary$median_missing_samples[i] <-
    median(miss)
}

chrom_summary <-
  segment_summary[, .(
    
    n_segments = .N,
    
    covered_bp = sum(segment_bp),
    
    total_cpg = sum(n_cpg),
    
    affected_cpg = sum(affected_cpg),
    
    prop_affected_cpg =
      sum(affected_cpg) / sum(n_cpg),
    
    mean_seg_prop_affected =
      mean(prop_affected_cpg),
    
    median_seg_prop =
      median(prop_affected_cpg),
    
    p90_seg_prop =
      quantile(prop_affected_cpg,0.9),
    
    max_seg_prop =
      max(prop_affected_cpg),
    #(num_tot_NA/(num_cpg sample_size)) 
    mean_norm_missing=
      mean(normalized_missing),
    
    median_norm_missing =
      median(normalized_missing),
    
    pct_segments_with_affected =
      mean(affected_cpg>0)*100,
    
    median_missing_samples =
      median(median_missing_samples),
    
    max_missing_samples =
      max(max_missing_samples)
    
  ), by=chr]

## ---------- CpG-level distribution within each chromosome ----------

cpg_dt <- data.table(
  chr = as.character(seqnames(cpg_gr)),
  na_samples = mcols(cpg_gr)$na_samples
)

## Keep only CpGs that fall inside a segment
hits <- findOverlaps(cpg_gr, all_seg_gr)

cpg_in_seg <- cpg_dt[queryHits(hits)]

## Add threshold indicators
cpg_in_seg[, ge10 := na_samples >= ceiling(0.10 * n_samples)]
cpg_in_seg[, ge20 := na_samples >= ceiling(0.20 * n_samples)]
cpg_in_seg[, ge50 := na_samples >= ceiling(0.50 * n_samples)]
cpg_in_seg[, ge90 := na_samples >= ceiling(0.90 * n_samples)]

## Summarize by chromosome
cpg_distribution_chr <-
  cpg_in_seg[, .(
    
    n_cpg = .N,
    
    mean_missing_samples = mean(na_samples),
    
    median_missing_samples = median(na_samples),
    
    p90_missing_samples = quantile(na_samples, 0.90),
    
    p99_missing_samples = quantile(na_samples, 0.99),
    
    max_missing_samples = max(na_samples),
    
    pct_cpg_ge10 = 100 * mean(ge10),
    
    pct_cpg_ge20 = 100 * mean(ge20),
    
    pct_cpg_ge50 = 100 * mean(ge50),
    
    pct_cpg_ge90 = 100 * mean(ge90)
    
  ), by = chr]

setnames(
  chrom_summary,
  "median_missing_samples",
  "median_segment_missing_samples"
)

setnames(
  cpg_distribution_chr,
  "median_missing_samples",
  "median_cpg_missing_samples"
)

chrom_summary <- merge(
  chrom_summary,
  cpg_distribution_chr,
  by = "chr",
  all.x = TRUE
)
fwrite(
  chrom_summary,
  "chromosome_missingness_summary.csv"
)