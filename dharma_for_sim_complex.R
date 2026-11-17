
# ============================================================
# DHARMa diagnostics for simulated T1 beta-binomial models
#
# Runs:
#   1. testUniformity
#   2. testDispersion
#   3. testOutliers
#
# One chromosome per Slurm array task.
# Checkpointing: completed model keys are not repeated.
# ============================================================

library(data.table)
library(glmmTMB)
library(DHARMa)


simulation_dir <- file.path(
  "/scratch/project_2010912/ezel/salmon/data_files",
  "simulated_data/sim_segs_complex"
)

model_output_dir <- file.path(
  simulation_dir,
  "model_fits_T1_sample50"
)

dharma_output_dir <- file.path(
  simulation_dir,
  "DHARMa_T1_sample50"
)

dir.create(
  dharma_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

n_sim <- 250L
batch_size <- 25L
base_seed <- 98765L


chromosomes <- list.dirs(
  model_output_dir,
  recursive = FALSE,
  full.names = FALSE
)

chromosomes <- sort(
  chromosomes[grepl("^NC_[0-9]+\\.[0-9]+$", chromosomes)]
)

if (length(chromosomes) == 0L) {
  stop("No chromosome directories found.")
}

task_id_text <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if (!nzchar(task_id_text)) {
  stop("SLURM_ARRAY_TASK_ID is not set.")
}

task_id <- as.integer(task_id_text)

if (is.na(task_id) ||
    task_id < 1L ||
    task_id > length(chromosomes)) {
  stop("Invalid array task ID: ", task_id_text)
}

chr <- chromosomes[task_id]

cat("========================================\n")
cat("DHARMa diagnostics\n")
cat("Chromosome:", chr, "\n")
cat("Array task:", task_id, "of", length(chromosomes), "\n")
cat("DHARMa simulations per model:", n_sim, "\n")
cat("========================================\n")

# ============================================================
# 3. Input and output files
# ============================================================

chr_model_dir <- file.path(
  model_output_dir,
  chr,
  "models"
)

model_results_file <- file.path(
  model_output_dir,
  chr,
  paste0("results_", chr, "_T1.csv")
)

dharma_results_file <- file.path(
  dharma_output_dir,
  paste0("DHARMa_results_", chr, "_T1.csv")
)

dharma_errors_file <- file.path(
  dharma_output_dir,
  paste0("DHARMa_errors_", chr, "_T1.csv")
)

if (!dir.exists(chr_model_dir)) {
  stop("Model directory does not exist: ", chr_model_dir)
}

if (!file.exists(model_results_file)) {
  stop("Model results CSV not found: ", model_results_file)
}


model_metadata <- fread(model_results_file)

required_metadata <- c(
  "key",
  "segment_id",
  "combination",
  "replicate",
  "model_file"
)

missing_metadata <- setdiff(
  required_metadata,
  names(model_metadata)
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing columns in model results CSV: ",
    paste(missing_metadata, collapse = ", ")
  )
}

model_metadata[, model_basename := basename(model_file)]

if (anyDuplicated(model_metadata$model_basename)) {
  stop("Duplicate model filenames in model results CSV.")
}

model_files <- list.files(
  chr_model_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

cat("Saved models found:", length(model_files), "\n")

if (length(model_files) == 0L) {
  cat("No models found for this chromosome.\n")
  quit(save = "no", status = 0)
}

model_info <- data.table(
  model_file = model_files,
  model_basename = basename(model_files)
)

# Match each RDS to its original model results row
model_info[, metadata_index := match(
  model_basename,
  model_metadata$model_basename
)]

unmatched_models <- model_info[
  is.na(metadata_index)
]

if (nrow(unmatched_models) > 0L) {
  stop(
    nrow(unmatched_models),
    " model files could not be matched to the model results CSV."
  )
}

model_info[, `:=`(
  key = model_metadata$key[metadata_index],
  segment_id = model_metadata$segment_id[metadata_index],
  combination = model_metadata$combination[metadata_index],
  replicate = model_metadata$replicate[metadata_index]
)]


if (file.exists(dharma_results_file)) {
  
  old_results <- fread(dharma_results_file)
  
  if (!"key" %in% names(old_results)) {
    stop("Existing DHARMa results CSV has no key column.")
  }
  
  completed_keys <- unique(old_results$key)
  
} else {
  
  completed_keys <- character(0)
  
}

cat("Already completed:", length(completed_keys), "\n")
#HELPERS

model_seed <- function(key) {
  
  chars <- as.integer(charToRaw(key))
  
  as.integer(
    (sum(as.double(chars) * seq_along(chars)) +
       base_seed) %% 2147483646 + 1
  )
}

make_key_dt <- function(key_value, ...) {
  
  dt <- data.table(...)
  
  dt[, key := key_value]
  
  setcolorder(
    dt,
    c("key", setdiff(names(dt), "key"))
  )
  
  dt
}

run_dharma_test <- function(test_function, residuals) {
  
  warning_messages <- character(0)
  
  test_result <- tryCatch(
    
    withCallingHandlers(
      
      test_function(residuals),
      
      warning = function(w) {
        
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        
        invokeRestart("muffleWarning")
      }
    ),
    
    error = function(e) e
  )
  
  if (inherits(test_result, "error")) {
    
    return(list(
      statistic = NA_real_,
      p_value = NA_real_,
      error = conditionMessage(test_result),
      warnings = paste(
        unique(warning_messages),
        collapse = " | "
      )
    ))
  }
  
  statistic <- NA_real_
  
  if (!is.null(test_result$statistic)) {
    statistic <- as.numeric(test_result$statistic[1])
  }
  
  p_value <- NA_real_
  
  if (!is.null(test_result$p.value)) {
    p_value <- as.numeric(test_result$p.value[1])
  }
  
  list(
    statistic = statistic,
    p_value = p_value,
    error = NA_character_,
    warnings = paste(
      unique(warning_messages),
      collapse = " | "
    )
  )
}


# Append a batch to CSV
append_batch <- function(x, file) {
  
  if (length(x) == 0L) {
    return(invisible(NULL))
  }
  
  dt <- rbindlist(
    x,
    fill = TRUE
  )
  
  if (nrow(dt) > 0L) {
    
    fwrite(
      dt,
      file,
      append = file.exists(file),
      col.names = !file.exists(file)
    )
  }
  
  invisible(NULL)
}


pending_results <- list()
pending_errors <- list()

n_processed <- 0L
n_success <- 0L
n_failed <- 0L
n_skipped <- 0L

for (i in seq_len(nrow(model_info))) {
  
  model_file <- model_info$model_file[i]
  key <- model_info$key[i]
  
  # Skip models already recorded
  if (key %in% completed_keys) {
    n_skipped <- n_skipped + 1L
    next
  }
  
  segment_id <- model_info$segment_id[i]
  combination <- model_info$combination[i]
  replicate <- model_info$replicate[i]
  
  cat(
    "\nModel", i, "of", nrow(model_info),
    "|", key, "\n"
  )
  
  n_processed <- n_processed + 1L

  
  result <- make_key_dt(
    key,
    segment_id = segment_id,
    combination = combination,
    replicate = replicate,
    model_file = model_file,
    n_sim = n_sim,
    seed = model_seed(key),
    status = NA_character_,
    uniformity_statistic = NA_real_,
    uniformity_p_value = NA_real_,
    dispersion_statistic = NA_real_,
    dispersion_p_value = NA_real_,
    outlier_statistic = NA_real_,
    outlier_p_value = NA_real_,
    uniformity_error = NA_character_,
    dispersion_error = NA_character_,
    outlier_error = NA_character_,
    warnings = NA_character_,
    simulation_error = NA_character_
  )
  error_rows <- list()
  

  
  fit <- tryCatch(
    readRDS(model_file),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    
    result[, `:=`(
      status = "model_load_error",
      simulation_error = conditionMessage(fit)
    )]
    
    error_rows[[1]] <- make_key_dt(
      key,
      model_file = model_file,
      stage = "model_load",
      error = conditionMessage(fit)
    )
    
  } else {
#dhMA
    
    warning_messages <- character(0)
    
    set.seed(model_seed(key))
    
    simulated_residuals <- tryCatch(
      
      withCallingHandlers(
        
        DHARMa::simulateResiduals(
          fittedModel = fit,
          n = n_sim,
          refit = FALSE,
          plot = FALSE
        ),
        
        warning = function(w) {
          
          warning_messages <<- c(
            warning_messages,
            conditionMessage(w)
          )
          
          invokeRestart("muffleWarning")
        }
      ),
      
      error = function(e) e
    )
    
    result[, warnings := paste(
      unique(warning_messages),
      collapse = " | "
    )]
    
    if (inherits(simulated_residuals, "error")) {
      
      result[, `:=`(
        status = "simulation_error",
        simulation_error = conditionMessage(simulated_residuals)
      )]
      
      error_rows[[length(error_rows) + 1L]] <- make_key_dt(
        key,
        model_file = model_file,
        stage = "simulateResiduals",
        error = conditionMessage(simulated_residuals)
      )
      
    } else {
      
      # ------------------------------------------------------
      # Uniformity test
      # ------------------------------------------------------
      
      uniformity <- run_dharma_test(
        DHARMa::testUniformity,
        simulated_residuals
      )
      
      result[, `:=`(
        uniformity_statistic = uniformity$statistic,
        uniformity_p_value = uniformity$p_value,
        uniformity_error = uniformity$error
      )]
 
      
      dispersion <- run_dharma_test(
        DHARMa::testDispersion,
        simulated_residuals
      )
      
      result[, `:=`(
        dispersion_statistic = dispersion$statistic,
        dispersion_p_value = dispersion$p_value,
        dispersion_error = dispersion$error
      )]

      
      outliers <- run_dharma_test(
        DHARMa::testOutliers,
        simulated_residuals
      )
      
      result[, `:=`(
        outlier_statistic = outliers$statistic,
        outlier_p_value = outliers$p_value,
        outlier_error = outliers$error
      )]

      
      all_warnings <- c(
        warning_messages,
        uniformity$warnings,
        dispersion$warnings,
        outliers$warnings
      )
      
      result[, warnings := paste(
        unique(all_warnings[nzchar(all_warnings)]),
        collapse = " | "
      )]

      
      test_errors <- c(
        uniformity$error,
        dispersion$error,
        outliers$error
      )
      
      if (all(is.na(test_errors))) {
        
        result[, status := "success"]
        
      } else {
        
        result[, status := "partial"]
        
      }
      
      # Save individual test errors
      test_names <- c(
        "uniformity",
        "dispersion",
        "outliers"
      )
      
      for (test_index in seq_along(test_names)) {
        
        if (!is.na(test_errors[test_index])) {
          
          error_rows[[length(error_rows) + 1L]] <- make_key_dt(
            key,
            model_file = model_file,
            stage = test_names[test_index],
            error = test_errors[test_index]
          )
        }
      }
      
      rm(simulated_residuals)
    }
    
    rm(fit)
  }

  
  pending_results[[length(pending_results) + 1L]] <- result
  
  if (length(error_rows) > 0L) {
    pending_errors <- c(
      pending_errors,
      error_rows
    )
  }
  
  if (result$status == "success") {
    n_success <- n_success + 1L
  } else {
    n_failed <- n_failed + 1L
  }
  

  
  if (length(pending_results) >= batch_size) {
    
    append_batch(
      pending_results,
      dharma_results_file
    )
    
    append_batch(
      pending_errors,
      dharma_errors_file
    )
    
    completed_keys <- c(
      completed_keys,
      vapply(
        pending_results,
        function(x) as.character(x$key[1]),
        character(1)
      )
    )
    
    pending_results <- list()
    pending_errors <- list()
    
    cat(
      "\nCheckpoint saved.",
      "Processed:", n_processed,
      "| Success:", n_success,
      "| Failed/partial:", n_failed,
      "| Skipped:", n_skipped,
      "\n"
    )
    
    gc(verbose = FALSE)
  }
  
  if (n_processed %% 10L == 0L) {
    cat(
      "Progress:",
      n_processed,
      "new models processed\n"
    )
  }
  
} # End model loop

# ============================================================
# 9. Save remaining results
# ============================================================

append_batch(
  pending_results,
  dharma_results_file
)

append_batch(
  pending_errors,
  dharma_errors_file
)

# ============================================================
# 10. Final summary
# ============================================================

cat("\n========================================\n")
cat("DHARMa chromosome task completed\n")
cat("Chromosome:", chr, "\n")
cat("Models found:", nrow(model_info), "\n")
cat("Newly processed:", n_processed, "\n")
cat("Successful:", n_success, "\n")
cat("Failed or partial:", n_failed, "\n")
cat("Skipped as already completed:", n_skipped, "\n")
cat("========================================\n")

cat("Results:\n", dharma_results_file, "\n")

if (file.exists(dharma_errors_file)) {
  cat("Errors:\n", dharma_errors_file, "\n")
}

cat("\nFinished.\n")
