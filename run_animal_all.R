library(data.table)
library(brms)


data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/animal_model/"

output_T1 <- file.path(data_dir, "final_models/animal_model_T1_all")
output_T2 <- file.path(data_dir, "final_models/animal_model_T2_all")

dir.create(output_T1, recursive = TRUE, showWarnings = FALSE)
dir.create(output_T2, recursive = TRUE, showWarnings = FALSE)


task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID is not set.")
}

chunk_size <- 10L

dat_T1 <- readRDS(
  file.path(data_dir, "keep_segments_model_data_T1.rds")
)

dat_T2 <- readRDS(
  file.path(data_dir, "keep_segments_model_data_T2.rds")
)

A_T1 <- readRDS(
  file.path(data_dir, "A_T1.rds")
)

A_T2 <- readRDS(
  file.path(data_dir, "A_T2.rds")
)



test_segments <- intersect(
  unique(dat_T1$segment_id),
  unique(dat_T2$segment_id)
)

test_segments <- sort(test_segments)

n_segments <- length(test_segments)

cat("Total shared segments:", n_segments, "\n")

# Save segment list
writeLines(
  test_segments,
  file.path(data_dir, "all_shared_test_segments.txt")
)


start_index <- task_id * chunk_size + 1L

end_index <- min(
  start_index + chunk_size - 1L,
  n_segments
)

if (start_index > n_segments) {
  stop("Array task exceeds the number of segments.")
}

segments_this_task <- test_segments[start_index:end_index]

cat("\n====================================\n")
cat("ARRAY TASK:", task_id, "\n")
cat("Segment indices:", start_index, "-", end_index, "\n")
cat("Number of segments:", length(segments_this_task), "\n")
cat("====================================\n")

formula <- brms::bf(
  meth | trials(coverage) ~
    1 + Weight_sc +
    (1 | gr(animal, cov = A))
)



fit_one_segment <- function(seg_id, timepoint, dat, A, output_dir, seed) {
  
  cat("\n------------------------------------\n")
  cat("Timepoint:", timepoint, "\n")
  cat("Segment:", seg_id, "\n")
  cat("------------------------------------\n")
  
  file_prefix <- paste0(timepoint, "_", seg_id)
  
  model_file <- file.path(
    output_dir,
    paste0("model_", file_prefix, ".rds")
  )
  
  summary_file <- file.path(
    output_dir,
    paste0("summary_", file_prefix, ".txt")
  )
  
  diagnostics_file <- file.path(
    output_dir,
    paste0("diagnostics_", file_prefix, ".csv")
  )
  
  divergences_file <- file.path(
    output_dir,
    paste0("divergences_", file_prefix, ".txt")
  )
  
  # Skip if all outputs already exist
  if (
    file.exists(model_file) &&
    file.exists(summary_file) &&
    file.exists(diagnostics_file) &&
    file.exists(divergences_file)
  ) {
    cat("Already completed. Skipping:", seg_id, "\n")
    return(invisible(NULL))
  }
  
  # Extract segment data
  dat_seg <- as.data.table(
    dat[dat$segment_id == seg_id, ]
  )
  
  if (nrow(dat_seg) == 0L) {
    stop("No observations found for ", timepoint, ": ", seg_id)
  }
  
  # Match relationship matrix to animals
  ids <- sort(unique(as.character(dat_seg$animal)))
  
  missing_ids <- setdiff(ids, rownames(A))
  
  if (length(missing_ids) > 0L) {
    stop(
      "Animals missing from relationship matrix for ",
      timepoint, ": ", seg_id
    )
  }
  
  A_seg <- A[ids, ids, drop = FALSE]
  
  dat_seg[, animal := factor(
    as.character(animal),
    levels = ids
  )]
  
  stopifnot(
    identical(levels(dat_seg$animal), rownames(A_seg)),
    identical(rownames(A_seg), colnames(A_seg))
  )
  
  # Fit model
  model_test <- brms::brm(
    formula = formula,
    data = dat_seg,
    data2 = list(A = A_seg),
    family = brms::beta_binomial(link = "logit"),
    chains = 4,
    cores = 1,
    iter = 10000,
    warmup = 3000,
    thin = 5,
    seed = seed,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    refresh = 100
  )
  
  # Save model
  saveRDS(
    model_test,
    model_file
  )
  
  # Save summary
  capture.output(
    summary(model_test),
    file = summary_file
  )
  
  # Diagnostics
  rhat_values <- brms::rhat(model_test)
  ess_values <- brms::neff_ratio(model_test)
  
  diagnostics <- data.table(
    parameter = names(rhat_values),
    Rhat = as.numeric(rhat_values),
    ESS_ratio = as.numeric(ess_values)
  )
  
  fwrite(
    diagnostics,
    diagnostics_file
  )
  
  # Divergences
  np <- as.data.frame(brms::nuts_params(model_test))
  
  divergences <- sum(
    np$Parameter == "divergent__" & np$Value == 1
  )
  
  writeLines(
    paste("Divergent transitions:", divergences),
    divergences_file
  )
  
  cat("Completed:", timepoint, seg_id, "\n")
  cat("Divergences:", divergences, "\n")
  
  # Release memory
  rm(
    model_test,
    dat_seg,
    A_seg,
    np,
    rhat_values,
    ess_values,
    diagnostics
  )
  
  gc()
  
  invisible(NULL)
}


for (i in seq_along(segments_this_task)) {
  
  seg_id <- segments_this_task[i]
  
  segment_index <- start_index + i - 1L
  
  cat("\n\n####################################\n")
  cat("SEGMENT", segment_index, "OF", n_segments, "\n")
  cat("####################################\n")
  
  # T1
  fit_one_segment(
    seg_id = seg_id,
    timepoint = "T1",
    dat = dat_T1,
    A = A_T1,
    output_dir = output_T1,
    seed = 123 + segment_index
  )
  
  gc()
  
  # T2
  fit_one_segment(
    seg_id = seg_id,
    timepoint = "T2",
    dat = dat_T2,
    A = A_T2,
    output_dir = output_T2,
    seed = 10000 + segment_index
  )
  
  gc()
}

cat("\n====================================\n")
cat("ARRAY TASK", task_id, "COMPLETED\n")
cat("====================================\n")
