library(data.table)

data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/no_extras"

# Methylation range threshold
threshold <- 0.1  

cs_files <- list.files(
  data_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)

all_segment_ranges <- list()

for (cs_file in cs_files) {
  
  # Find matching coverage file
  prefix <- sub("_segment_cs\\.csv$", "", basename(cs_file))
  cov_file <- file.path(
    data_dir,
    paste0(prefix, "_segment_coverage.csv")
  )
  
  if (!file.exists(cov_file)) {
    warning("Missing coverage file for: ", basename(cs_file))
    next
  }
  
  cs  <- fread(cs_file)
  cov <- fread(cov_file)
  
  coord_cols <- intersect(
    c("chr", "start", "end", "segment", "segment_id"),
    intersect(names(cs), names(cov))
  )
  
  sample_cols <- intersect(
    setdiff(names(cs), coord_cols),
    setdiff(names(cov), coord_cols)
  )
  
  if (length(sample_cols) == 0) {
    warning("No shared sample columns in: ", basename(cs_file))
    next
  }
  
  if (nrow(cs) != nrow(cov)) {
    warning("Different row counts in: ", basename(cs_file))
    next
  }
  
  if (length(coord_cols) > 0 &&
      !identical(cs[, ..coord_cols], cov[, ..coord_cols])) {
    warning("Coordinates/order differ in: ", basename(cs_file))
    next
  }
  
#get methyl percent
  cs_mat  <- as.matrix(cs[, ..sample_cols])
  cov_mat <- as.matrix(cov[, ..sample_cols])
  
  # Calculate methylation proportions
  meth <- cs_mat / cov_mat
  meth[is.na(cov_mat) | cov_mat <= 0] <- NA_real_
  meth[!is.finite(meth)] <- NA_real_
  
#get rage
    ranges <- apply(meth, 1, function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    }
  })
  
  result <- copy(cs[, ..coord_cols])
  result[, chromosome_file := prefix]
  result[, methylation_range := ranges]
  result[, passes_threshold := !is.na(ranges) & ranges > threshold]
  
  all_segment_ranges[[prefix]] <- result
  
  message(
    prefix, ": ",
    sum(result$passes_threshold, na.rm = TRUE),
    " segments pass threshold"
  )
}

segment_ranges <- rbindlist(all_segment_ranges, fill = TRUE)
selected_segments <- segment_ranges[passes_threshold == TRUE]


cat("Total segments:", nrow(segment_ranges), "\n")
cat("Segments above threshold:", nrow(selected_segments), "\n")
cat("Threshold:", threshold, "\n")


#for the loci

library(data.table)
meth=readRDS("/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/simulated_T1_parental_methylation_1000per_chr.rds")

threshold <- 0.1

meth_mat <- meth$methylation_probability

ranges <- apply(meth_mat, 1, function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  }
})

segment_ranges <- copy(meth$loci)

segment_ranges[, methylation_range := ranges]

segment_ranges[, passes_threshold :=
                 !is.na(methylation_range) &
                 methylation_range > threshold]

selected_segments <- segment_ranges[passes_threshold == TRUE]

cat("Total loci:", nrow(segment_ranges), "\n")
cat("Loci above threshold:", nrow(selected_segments), "\n")
cat("Threshold:", threshold, "\n")