library(data.table)

# Directory containing all chromosome-specific CSV files
data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/cs_cov"

# Methylation range threshold
threshold <- 0.1  # 20 percentage points

# Find all CS files
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
  
  # Identify shared sample columns, excluding coordinates
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
  
  # Check matching coordinates if available
  if (length(coord_cols) > 0 &&
      !identical(cs[, ..coord_cols], cov[, ..coord_cols])) {
    warning("Coordinates/order differ in: ", basename(cs_file))
    next
  }
  
  # Convert counts and coverage to matrices
  cs_mat  <- as.matrix(cs[, ..sample_cols])
  cov_mat <- as.matrix(cov[, ..sample_cols])
  
  # Calculate methylation proportions
  meth <- cs_mat / cov_mat
  meth[is.na(cov_mat) | cov_mat <= 0] <- NA_real_
  meth[!is.finite(meth)] <- NA_real_
  
  # Calculate range for each segment across samples
  ranges <- apply(meth, 1, function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    }
  })
  
  # Store results for this chromosome
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

# Combine all chromosomes
segment_ranges <- rbindlist(all_segment_ranges, fill = TRUE)

# Extract qualifying segments
selected_segments <- segment_ranges[passes_threshold == TRUE]

# Save results
fwrite(
  segment_ranges,
  file.path(data_dir, "all_segment_methylation_ranges.csv")
)

fwrite(
  selected_segments,
  file.path(data_dir, "selected_segments_range_above_threshold.csv")
)

# Summary
cat("Total segments:", nrow(segment_ranges), "\n")
cat("Segments above threshold:", nrow(selected_segments), "\n")
cat("Threshold:", threshold, "\n")