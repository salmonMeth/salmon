library(methylKit)

base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"
path_diff_meth =  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_segments"
# Find files
cov_files <- list.files(base_dir, pattern = "_segment_coverage\\.csv$", full.names = TRUE)
cs_files  <- list.files(base_dir, pattern = "_segment_cs\\.csv$", full.names = TRUE)

# Sort to make sure chromosomes match
cov_files <- sort(cov_files)
cs_files  <- sort(cs_files)

# Read and combine all chromosomes
cov_list <- lapply(cov_files, read.csv)
cs_list  <- lapply(cs_files, read.csv)

coverage_df <- do.call(rbind, cov_list)
cs_df       <- do.call(rbind, cs_list)

names(coverage_df) <- sub("^X", "", names(coverage_df))
names(cs_df) <- sub("^X", "", names(cs_df))
# Check columns
head(coverage_df)
head(cs_df)

sample_ids <- setdiff(names(coverage_df), c("chr", "start", "end"))
new_treatment <- ifelse(grepl("T1", sample_ids), 0, 1)


df <- merge(
  coverage_df,
  cs_df,
  by = c("chr", "start", "end"),
  suffixes = c("_cov", "_cs")
)

# Get sample IDs from coverage columns
sample_ids <- sub("_cov$", "", grep("_cov$", names(df), value = TRUE))

# Rename coverage columns
names(df)[grep("_cov$", names(df))] <- paste0("coverage_", sample_ids)

# Rename Cs columns
names(df)[grep("_cs$", names(df))] <- paste0("numCs_", sample_ids)

# Create numTs columns
for (s in sample_ids) {
  df[[paste0("numTs_", s)]] <-
    df[[paste0("coverage_", s)]] -
    df[[paste0("numCs_", s)]]
}

# Add strand column (required by methylKit)
df$strand <- "*"

# Reorder columns for methylKit
df <- df[, c(
  "chr",
  "start",
  "end",
  "strand",
  unlist(lapply(sample_ids, function(s)
    c(
      paste0("coverage_", s),
      paste0("numCs_", s),
      paste0("numTs_", s)
    )
  ))
)]

mb <- new(
  "methylBase",
  df,
  sample.ids = sample_ids,
  assembly = "-",
  context = "CpG",
  treatment = new_treatment,
  resolution = "region"
)

##############
#possibly unnecessary bits
n <- length(sample_ids)
coverage_idx <- seq(5, by = 3, length.out = n)
numCs_idx    <- seq(6, by = 3, length.out = n)
numTs_idx    <- seq(7, by = 3, length.out = n)
mb@coverage.index <- coverage_idx
mb@numCs.index <- numCs_idx
mb@numTs.index <- numTs_idx
#############

mbdb <- makeMethylDB(mb, "methylDB")


diff <- calculateDiffMeth(mbdb)

differ <- 15

dmr <- getMethylDiff(
  diff,
  differ,
  qvalue = 0.01
)

write.csv(
  getData(dmr),
  file = paste0(path_diff_meth, "_difference_", differ, ".csv"),
  row.names = FALSE
)