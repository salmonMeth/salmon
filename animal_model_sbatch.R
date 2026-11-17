
library(data.table)
library(brms)

# Directories
data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/animal_model"

output_dir <- file.path(
  data_dir,
  "animal_model_151-251"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load saved data
dat_T1 <- readRDS(
  file.path(data_dir, "keep_segments_model_data_T1.rds")
)

A_T1 <- readRDS(
  file.path(data_dir, "A_T1.rds")
)

dat_T2 <- readRDS(
  file.path(data_dir, "keep_segments_model_data_T2.rds")
)

# Select first 20 shared segments
test_segments <- intersect(
  unique(dat_T1$segment_id),
  unique(dat_T2$segment_id)
)

test_segments <- test_segments[151:251]

cat("Number of segments:", length(test_segments), "\n")
writeLines(
  test_segments,
  file.path(output_dir, "test_segments.txt")
)

# Model formula
formula <- brms::bf(
  meth | trials(coverage) ~
    1 + Weight_sc +
    (1 | gr(animal, cov = A))
)

# Fit each segment separately
for (i in seq_along(test_segments)) {
  
  seg_id <- test_segments[i]
  
  cat("\n----------------------------------\n")
  cat("Fitting segment", i, "of", length(test_segments), "\n")
  cat("Segment:", seg_id, "\n")
  cat("----------------------------------\n")
  
  dat_seg <- as.data.table(
    dat_T1[dat_T1$segment_id == seg_id, ]
  )
  
  # Match relationship matrix to animals
  ids <- sort(unique(as.character(dat_seg$animal)))
  
  A_seg <- A_T1[ids, ids, drop = FALSE]
  
  dat_seg[, animal := factor(
    as.character(animal),
    levels = ids
  )]
  
  stopifnot(
    identical(levels(dat_seg$animal), rownames(A_seg)),
    identical(rownames(A_seg), colnames(A_seg))
  )
  
  # Fit model
  model_test_T1 <- brms::brm(
    formula = formula,
    data = dat_seg,
    data2 = list(A = A_seg),
    family = brms::beta_binomial(link = "logit"),
    chains = 4,
    cores = 1,
    iter = 6000,
    warmup = 3000,
    thin = 5,
    seed = 123 + i,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    refresh = 100
  )
  
  file_prefix <- paste0("T1_", seg_id)
  
  # Save model
  saveRDS(
    model_test_T1,
    file.path(output_dir, paste0("model_", file_prefix, ".rds"))
  )
    capture.output(
    summary(model_test_T1),
    file = file.path(output_dir, paste0("summary_", file_prefix, ".txt"))
  )
  
  rhat_values <- brms::rhat(model_test_T1)
  ess_values <- brms::neff_ratio(model_test_T1)
  
  diagnostics <- data.table(
    parameter = names(rhat_values),
    Rhat = as.numeric(rhat_values),
    ESS_ratio = as.numeric(ess_values)
  )
  
  fwrite(
    diagnostics,
    file.path(output_dir, paste0("diagnostics_", file_prefix, ".csv"))
  )
  np <- as.data.frame(brms::nuts_params(model_test_T1))
  
  divergences <- sum(
    np$Parameter == "divergent__" & np$Value == 1
  )
  
  writeLines(
    paste("Divergent transitions:", divergences),
    file.path(output_dir, paste0("divergences_", file_prefix, ".txt"))
  )
  
  cat("Completed:", seg_id, "\n")
  cat("Divergent transitions:", divergences, "\n")
  
  # Release memory
  rm(model_test_T1, dat_seg, A_seg, np,
     rhat_values, ess_values, diagnostics)
  gc()
}

cat("\nAll segment models completed.\n")
