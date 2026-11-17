
library(data.table)

dharma_dir <- file.path(
  "/scratch/project_2010912/ezel/salmon/data_files",
  "simulated_data/sim_segs_complex/DHARMa_T1_sample50"
)

output_dir <- file.path(
  dharma_dir,
  "randomly_assembled_datasets_100.3"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

n_datasets <- 500L

#base_seed <- 20261115L
#base_seed <- 34209126L
base_seed <- 98210004L


alpha <- 0.05
under_threshold <- 0.8
over_threshold <- 1.2

result_files <- list.files(
  dharma_dir,
  pattern = "^DHARMa_results_.*_T1\\.csv$",
  full.names = TRUE
)

if (length(result_files) == 0L) {
  stop("No DHARMa result files found.")
}

cat("DHARMa result files:", length(result_files), "\n")

dharma <- rbindlist(
  lapply(result_files, fread),
  fill = TRUE
)

cat("Total DHARMa rows:", nrow(dharma), "\n")

# Check required columns
required_columns <- c(
  "key",
  "segment_id",
  "combination",
  "replicate",
  "status",
  "dispersion_statistic",
  "dispersion_p_value"
)

missing_columns <- setdiff(
  required_columns,
  names(dharma)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

eligible <- dharma[
  status == "success" &
    !is.na(dispersion_statistic) &
    !is.na(dispersion_p_value)
]

cat(
  "Successfully diagnosed models:",
  nrow(eligible),
  "\n"
)

# Check unique model identity
if (anyDuplicated(eligible$key)) {
  stop("Duplicate model keys found in successful results.")
}

# Check unique segment-combination-replicate identities
if (anyDuplicated(
  eligible[, .(segment_id, combination, replicate)]
)) {
  stop(
    "Duplicate segment-combination-replicate combinations found."
  )
}
segments <- sort(unique(eligible$segment_id))

n_segments <- length(segments)

cat("Eligible segments:", n_segments, "\n")

# One row per available segment-parameter combination
combination_pool <- unique(
  eligible[, .(segment_id, combination)]
)

combination_counts <- combination_pool[
  ,
  .(n_combinations = .N),
  by = segment_id
]

cat("\nParameter combinations per segment:\n")
print(summary(combination_counts$n_combinations))

# Every segment must have at least one available combination
if (!setequal(segments, combination_counts$segment_id)) {
  stop("Some eligible segments have no available combinations.")
}
classify_dispersion <- function(statistic, p_value) {
  
  classification <- rep(
    NA_character_,
    length(statistic)
  )
  
  classification[p_value >= alpha] <- "unflagged"
  
  classification[
    p_value < alpha &
      statistic < under_threshold
  ] <- "underdispersed"
  
  classification[
    p_value < alpha &
      statistic > over_threshold
  ] <- "overdispersed"
  
  classification[
    p_value < alpha &
      statistic >= under_threshold &
      statistic <= over_threshold
  ] <- "other_significant"
  
  classification
}

set.seed(base_seed)

all_selections <- vector("list", n_datasets)

for (dataset_index in seq_len(n_datasets)) {
  
  cat(
    "\nGenerating dataset",
    dataset_index,
    "of",
    n_datasets,
    "\n"
  )
 
  selected_combinations <- combination_pool[
    ,
    .SD[sample.int(.N, 1L)],
    by = segment_id
  ]
  
  candidates <- eligible[
    selected_combinations,
    on = .(segment_id, combination),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]
  
  selected <- candidates[
    ,
    .SD[sample.int(.N, 1L)],
    by = segment_id
  ]

  
  if (nrow(selected) != n_segments) {
    stop(
      "Incorrect number of segments in dataset ",
      dataset_index
    )
  }
  
  if (uniqueN(selected$segment_id) != n_segments) {
    stop(
      "Duplicate or missing segments in dataset ",
      dataset_index
    )
  }
  
  
  selected[
    ,
    classification := classify_dispersion(
      dispersion_statistic,
      dispersion_p_value
    )
  ]
  
  if (anyNA(selected$classification)) {
    stop(
      "Some segments could not be classified in dataset ",
      dataset_index
    )
  }
  
  selected[
    ,
    dataset := dataset_index
  ]
  
  # Keep useful columns
  selected <- selected[
    ,
    .(
      dataset,
      segment_id,
      key,
      combination,
      replicate,
      dispersion_statistic,
      dispersion_p_value,
      classification
    )
  ]
  
  all_selections[[dataset_index]] <- selected
  
}

selections <- rbindlist(
  all_selections,
  use.names = TRUE
)

setorder(
  selections,
  dataset,
  segment_id
)

cat("\nAll datasets generated.\n")

cat("Total selected rows:", nrow(selections), "\n")

cat(
  "Expected rows:",
  n_datasets * n_segments,
  "\n"
)

fwrite(
  selections,
  file.path(
    output_dir,
    "selected_models_all_100_datasets.csv"
  )
)

classification_levels <- c(
  "unflagged",
  "underdispersed",
  "overdispersed",
  "other_significant"
)

counts_long <- selections[
  ,
  .N,
  by = .(dataset, classification)
]

count_grid <- CJ(
  dataset = seq_len(n_datasets),
  classification = classification_levels
)

counts_long <- merge(
  count_grid,
  counts_long,
  by = c("dataset", "classification"),
  all.x = TRUE
)

counts_long[is.na(N), N := 0L]

counts <- dcast(
  counts_long,
  dataset ~ classification,
  value.var = "N"
)

# ------------------------------------------------------------
# 9. Add percentages and total
# ------------------------------------------------------------

counts[
  ,
  total_segments := rowSums(
    .SD
  ),
  .SDcols = classification_levels
]

for (category in classification_levels) {
  
  percentage_name <- paste0(
    category,
    "_percent"
  )
  
  counts[
    ,
    (percentage_name) :=
      100 * get(category) / total_segments
  ]
}

# Verify totals
if (any(counts$total_segments != n_segments)) {
  stop("Some assembled datasets have incorrect totals.")
}

fwrite(
  counts,
  file.path(
    output_dir,
    "classification_counts_100_datasets.csv"
  )
)

summary_results <- rbindlist(
  lapply(classification_levels, function(category) {
    
    values <- counts[[category]]
    
    percentages <- counts[[
      paste0(category, "_percent")
    ]]
    
    data.table(
      classification = category,
      
      mean_count = mean(values),
      sd_count = sd(values),
      min_count = min(values),
      median_count = median(values),
      max_count = max(values),
      
      mean_percent = mean(percentages),
      sd_percent = sd(percentages),
      min_percent = min(percentages),
      median_percent = median(percentages),
      max_percent = max(percentages)
    )
    
  })
)

fwrite(
  summary_results,
  file.path(
    output_dir,
    "classification_summary_across_100_datasets.csv"
  )
)


cat("\n========================================\n")
cat("RESAMPLING COMPLETED\n")
cat("========================================\n")

cat("Number of assembled datasets:", n_datasets, "\n")
cat("Segments per dataset:", n_segments, "\n")
cat("Total selected models:", nrow(selections), "\n")

cat("\nClassification counts for first 10 datasets:\n")
print(head(counts, 10))

cat("\nSummary across all datasets:\n")
print(summary_results)

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nFinished.\n")
