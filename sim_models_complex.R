
# ============================================================
# Fit simulated T1 beta-binomial models
# Save one glmmTMB RDS per successful dataset
# Process one chromosome per Slurm array task
# Restart-safe; maximum 50 datasets per segment
# ============================================================

library(data.table)
library(glmmTMB)

# ============================================================
# 1. Paths and settings
# ============================================================

simulation_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/sim_segs_complex"

metadata_file <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"

output_dir <- file.path(
  simulation_dir,
  "model_fits_T1_sample50"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

max_datasets_per_segment <- 50L
minimum_samples <- 42L
base_seed <- 123L

# ============================================================
# 2. Load metadata
# ============================================================

meta <- fread(metadata_file)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta_T1 <- meta[TimePoint == "T1"]

if (nrow(meta_T1) == 0L) {
  stop("No T1 samples found in metadata.")
}

if (anyDuplicated(meta_T1$sample)) {
  stop("Duplicate T1 sample IDs found in metadata.")
}

meta_T1[, Weight.sc := as.numeric(scale(Weight.g))]

sample_metadata <- meta_T1[
  ,
  .(sample, Weight.sc)
]

# ============================================================
# 3. Find simulation files
# ============================================================

simulation_files <- list.files(
  simulation_dir,
  pattern = "_simulations\\.rds$",
  full.names = TRUE
)

if (length(simulation_files) == 0L) {
  stop("No simulation files found.")
}

file_info <- data.table(
  file = simulation_files,
  segment_id = sub(
    "_simulations\\.rds$",
    "",
    basename(simulation_files)
  )
)

file_info[, chr := sub(
  "^(NC_[0-9]+\\.[0-9]+)_.*$",
  "\\1",
  segment_id
)]

if (anyDuplicated(file_info$segment_id)) {
  stop("Duplicate segment IDs found in simulation files.")
}

all_chromosomes <- sort(unique(file_info[["chr"]]))

cat("Simulation files:", nrow(file_info), "\n")
cat("Chromosomes:", length(all_chromosomes), "\n")

# ============================================================
# 4. Select chromosome for this Slurm array task
# ============================================================

task_id_text <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if (nzchar(task_id_text)) {
  
  task_id <- as.integer(task_id_text)
  
  if (is.na(task_id) ||
      task_id < 1L ||
      task_id > length(all_chromosomes)) {
    stop("Invalid SLURM_ARRAY_TASK_ID: ", task_id_text)
  }
  
  chromosomes_to_process <- all_chromosomes[task_id]
  
} else {
  
  # If running interactively, process all chromosomes
  chromosomes_to_process <- all_chromosomes
  
}

cat("\nAssigned chromosome(s):\n")
print(chromosomes_to_process)

# ============================================================
# 5. Helper functions
# ============================================================

get_random_sd <- function(fit, group) {
  
  vc <- VarCorr(fit)$cond
  
  if (!group %in% names(vc)) {
    return(NA_real_)
  }
  
  as.numeric(attr(vc[[group]], "stddev")[1])
}


segment_seed <- function(segment_id) {
  
  chars <- utf8ToInt(segment_id)
  
  as.integer(
    (sum(as.double(chars) * seq_along(chars)) +
       base_seed) %% 2147483646 + 1
  )
}


model_filename <- function(segment_id, combo, rep_id) {
  
  paste0(
    segment_id,
    "__comb", combo,
    "__rep", rep_id,
    ".rds"
  )
}


append_csv <- function(x, file) {
  
  if (nrow(x) == 0L) {
    return(invisible(NULL))
  }
  
  fwrite(
    x,
    file,
    append = file.exists(file),
    col.names = !file.exists(file)
  )
}


fit_one_dataset <- function(dat) {
  
  warning_messages <- character(0)
  
  fit <- tryCatch(
    
    withCallingHandlers(
      
      glmmTMB(
        cbind(meth, unmeth) ~
          Weight.sc +
          (1 | ID_ma) +
          (1 | ID_pa),
        data = dat,
        family = betabinomial(link = "logit")
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
  
  if (inherits(fit, "error")) {
    
    return(list(
      fit = NULL,
      error = conditionMessage(fit),
      warnings = paste(
        unique(warning_messages),
        collapse = " | "
      )
    ))
  }
  
  list(
    fit = fit,
    error = NULL,
    warnings = paste(
      unique(warning_messages),
      collapse = " | "
    )
  )
}


make_result_row <- function(
    fit,
    dat,
    key,
    segment_id,
    combo,
    rep_id,
    model_file,
    warning_text = ""
) {
  
  result <- data.table(
    estimated_sd_ma = get_random_sd(fit, "ID_ma"),
    estimated_sd_pa = get_random_sd(fit, "ID_pa"),
    estimated_phi = as.numeric(sigma(fit)),
    weight_estimate = unname(
      fixef(fit)$cond["Weight.sc"]
    ),
    AIC = AIC(fit),
    BIC = BIC(fit),
    logLik = as.numeric(logLik(fit)),
    convergence = fit$fit$convergence,
    pdHess = fit$sdr$pdHess,
    warnings = warning_text
  )
  
  result[, `:=`(
    segment_id = segment_id,
    combination = combo,
    replicate = rep_id,
    n_samples = uniqueN(dat$sample),
    n_mothers = uniqueN(dat$ID_ma),
    n_fathers = uniqueN(dat$ID_pa),
    true_sd_ma = unique(dat$sd_ma)[1],
    true_sd_pa = unique(dat$sd_pa)[1],
    true_sd_env = unique(dat$sd_env)[1],
    true_phi = unique(dat$phi)[1],
    model_file = model_file
  )]
  
  # IMPORTANT:
  # Do not pass key = key to data.table().
  # "key" is a reserved argument in data.table().
  # Add the column afterwards using :=.
  
  result[, key := key]
  
  setcolorder(
    result,
    c(
      "key",
      "segment_id",
      "combination",
      "replicate",
      "model_file",
      setdiff(
        names(result),
        c(
          "key",
          "segment_id",
          "combination",
          "replicate",
          "model_file"
        )
      )
    )
  )
  
  result
}

# ============================================================
# 6. Process assigned chromosomes
# ============================================================

for (current_chr in chromosomes_to_process) {
  
  cat("\n========================================\n")
  cat("Starting chromosome:", current_chr, "\n")
  cat("========================================\n")
  
  chr_dir <- file.path(
    output_dir,
    current_chr
  )
  
  model_dir <- file.path(
    chr_dir,
    "models"
  )
  
  dir.create(
    model_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  results_file <- file.path(
    chr_dir,
    paste0("results_", current_chr, "_T1.csv")
  )
  
  errors_file <- file.path(
    chr_dir,
    paste0("errors_", current_chr, "_T1.csv")
  )
  
  # ----------------------------------------------------------
  # Get only the files belonging to this chromosome
  # ----------------------------------------------------------
  
  chr_files <- file_info[
    file_info[["chr"]] == current_chr
  ]
  
  cat("Segments assigned:", nrow(chr_files), "\n")
  
  if (nrow(chr_files) == 0L) {
    cat("No simulation files for this chromosome.\n")
    next
  }
  
  # ----------------------------------------------------------
  # Repair old result CSV:
  # Keep only rows whose corresponding model RDS exists.
  #
  # Old results without model RDS must be refitted so that
  # we have actual models available for DHARMa.
  # ----------------------------------------------------------
  
  if (file.exists(results_file)) {
    
    old_results <- fread(results_file)
    
    if ("model_file" %in% names(old_results)) {
      
      valid_old_results <- old_results[
        file.exists(model_file)
      ]
      
      if (nrow(valid_old_results) < nrow(old_results)) {
        
        cat(
          "Removing",
          nrow(old_results) - nrow(valid_old_results),
          "result rows without saved models.\n"
        )
      }
      
      old_results <- valid_old_results
      
    } else {
      
      cat(
        "Previous results have no model_file column;",
        "they will be regenerated from saved models.\n"
      )
      
      old_results <- old_results[0]
    }
    
    # Keep only one row per key
    if (nrow(old_results) > 0L) {
      setorder(old_results, key)
      old_results <- unique(old_results, by = "key")
      fwrite(old_results, results_file)
    } else {
      file.remove(results_file)
    }
  }
  
  # Read keys with both a result row and saved model
  if (file.exists(results_file)) {
    
    existing_results <- fread(results_file)
    
    completed_keys <- unique(existing_results[["key"]])
    
  } else {
    
    completed_keys <- character(0)
  }
  
  # ----------------------------------------------------------
  # Process each segment
  # ----------------------------------------------------------
  
  for (file_index in seq_len(nrow(chr_files))) {
    
    simulation_file <- chr_files[["file"]][file_index]
    segment_id <- chr_files[["segment_id"]][file_index]
    
    cat(
      "\n----------------------------------------\n"
    )
    cat(
      "Segment", file_index, "of", nrow(chr_files),
      ":", segment_id, "\n"
    )
    
    sim <- as.data.table(readRDS(simulation_file))
    
    required_columns <- c(
      "sample",
      "ID_ma",
      "ID_pa",
      "meth",
      "unmeth",
      "coverage",
      "sd_ma",
      "sd_pa",
      "sd_env",
      "phi",
      "combination",
      "replicate"
    )
    
    missing_columns <- setdiff(
      required_columns,
      names(sim)
    )
    
    if (length(missing_columns) > 0L) {
      
      cat(
        "Skipping segment; missing columns:",
        paste(missing_columns, collapse = ", "),
        "\n"
      )
      
      next
    }
    
    if (nrow(sim) == 0L) {
      cat("Empty simulation file. Skipping.\n")
      next
    }
    
    # --------------------------------------------------------
    # Add sample weights
    # --------------------------------------------------------
    
    sim[, Weight.sc := sample_metadata$Weight.sc[
      match(sample, sample_metadata$sample)
    ]]
    
    # --------------------------------------------------------
    # Filter invalid observations
    # --------------------------------------------------------
    
    sim <- sim[
      is.finite(meth) &
        is.finite(unmeth) &
        is.finite(coverage) &
        coverage > 0 &
        meth >= 0 &
        unmeth >= 0 &
        meth + unmeth == coverage &
        is.finite(Weight.sc) &
        !is.na(ID_ma) &
        !is.na(ID_pa)
    ]
    
    if (nrow(sim) == 0L) {
      cat("No valid observations. Skipping.\n")
      next
    }
    
    sim[, ID_ma := factor(ID_ma)]
    sim[, ID_pa := factor(ID_pa)]
    
    # Keep grouping columns as character
    sim[, combination := as.character(combination)]
    sim[, replicate := as.character(replicate)]
    
    # --------------------------------------------------------
    # Identify eligible datasets
    # --------------------------------------------------------
    
    dataset_qc <- sim[
      ,
      .(
        n_samples = uniqueN(sample),
        duplicate_samples = anyDuplicated(sample) > 0L
      ),
      by = .(combination, replicate)
    ]
    
    eligible <- dataset_qc[
      n_samples >= minimum_samples &
        !duplicate_samples
    ]
    
    cat("Eligible datasets:", nrow(eligible), "\n")
    
    if (nrow(eligible) == 0L) {
      cat("No eligible datasets. Skipping segment.\n")
      next
    }
    
    # --------------------------------------------------------
    # Select up to 50 datasets
    # --------------------------------------------------------
    
    selected_file <- file.path(
      chr_dir,
      paste0(
        "selected_datasets_",
        segment_id,
        ".csv"
      )
    )
    
    if (file.exists(selected_file)) {
      
      selected <- fread(selected_file)
      
      selected[, combination := as.character(combination)]
      selected[, replicate := as.character(replicate)]
      
      cat(
        "Loaded existing selection:",
        nrow(selected),
        "datasets\n"
      )
      
      # Verify that the saved selection is still eligible
      eligible_keys <- paste(
        eligible[["combination"]],
        eligible[["replicate"]],
        sep = "__"
      )
      
      selected_keys <- paste(
        selected[["combination"]],
        selected[["replicate"]],
        sep = "__"
      )
      
      if (any(!selected_keys %in% eligible_keys)) {
        stop(
          "Saved selection contains ineligible datasets: ",
          segment_id
        )
      }
      
    } else {
      
      set.seed(segment_seed(segment_id))
      
      n_selected <- min(
        max_datasets_per_segment,
        nrow(eligible)
      )
      
      selected <- eligible[
        sample.int(
          nrow(eligible),
          size = n_selected,
          replace = FALSE
        )
      ]
      
      setorder(
        selected,
        combination,
        replicate
      )
      
      fwrite(selected, selected_file)
      
      cat(
        "Selected",
        nrow(selected),
        "datasets\n"
      )
    }
    
    if (nrow(selected) == 0L) {
      cat("No selected datasets. Skipping.\n")
      next
    }
    
    # --------------------------------------------------------
    # Fit selected datasets
    # --------------------------------------------------------
    
    for (j in seq_len(nrow(selected))) {
      
      combo <- selected[["combination"]][j]
      rep_id <- selected[["replicate"]][j]
      
      key <- paste(
        segment_id,
        combo,
        rep_id,
        sep = "__"
      )
      
      model_file <- file.path(
        model_dir,
        model_filename(
          segment_id,
          combo,
          rep_id
        )
      )
      
      # ------------------------------------------------------
      # If model exists but result is missing, recover result
      # ------------------------------------------------------
      
      if (file.exists(model_file)) {
        
        if (key %in% completed_keys) {
          next
        }
        
        cat("Recovering result from saved model:", key, "\n")
        
        fit <- readRDS(model_file)
        
        dat <- sim[
          combination == combo &
            replicate == rep_id
        ]
        
        recovered_result <- make_result_row(
          fit = fit,
          dat = dat,
          key = key,
          segment_id = segment_id,
          combo = combo,
          rep_id = rep_id,
          model_file = model_file,
          warning_text = ""
        )
        
        append_csv(
          recovered_result,
          results_file
        )
        
        completed_keys <- c(completed_keys, key)
        
        rm(fit, dat, recovered_result)
        next
      }
      
      # ------------------------------------------------------
      # Fit only if model RDS does not already exist
      # ------------------------------------------------------
      
      dat <- sim[
        combination == combo &
          replicate == rep_id
      ]
      
      if (uniqueN(dat$sample) < minimum_samples ||
          anyDuplicated(dat$sample) > 0L) {
        
        cat("Dataset failed QC at fitting stage. Skipping.\n")
        next
      }
      
      cat(
        "\nFitting", j, "of", nrow(selected),
        ":", key, "\n"
      )
      
      fit_output <- fit_one_dataset(dat)
      
      if (is.null(fit_output$fit)) {
        
        # IMPORTANT:
        # Do not use key = key inside data.table().
        # "key" is a reserved data.table argument.
        
        error_row <- data.table(
          segment_id = segment_id,
          combination = combo,
          replicate = rep_id,
          error = fit_output$error,
          warnings = fit_output$warnings
        )
        
        error_row[, key := key]
        
        setcolorder(
          error_row,
          c(
            "key",
            setdiff(names(error_row), "key")
          )
        )
        
        append_csv(
          error_row,
          errors_file
        )
        
        cat(
          "Fit failed:",
          fit_output$error,
          "\n"
        )
        
        rm(fit_output, dat)
        next
      }
      
      fit <- fit_output$fit
      
      # ------------------------------------------------------
      # Save model FIRST
      # ------------------------------------------------------
      
      saveRDS(
        fit,
        model_file
      )
      
      cat(
        "Saved model:",
        basename(model_file),
        "\n"
      )
      
      # ------------------------------------------------------
      # Extract and save results
      # ------------------------------------------------------
      
      result <- make_result_row(
        fit = fit,
        dat = dat,
        key = key,
        segment_id = segment_id,
        combo = combo,
        rep_id = rep_id,
        model_file = model_file,
        warning_text = fit_output$warnings
      )
      
      append_csv(
        result,
        results_file
      )
      
      completed_keys <- c(
        completed_keys,
        key
      )
      
      cat(
        "Completed:", key,
        "| Model exists:", file.exists(model_file),
        "| Successful fits:", length(completed_keys),
        "\n"
      )
      
      rm(
        fit,
        fit_output,
        result,
        dat
      )
      
    } # End selected dataset loop
    
    # --------------------------------------------------------
    # Segment summary
    # --------------------------------------------------------
    
    segment_model_files <- list.files(
      model_dir,
      pattern = paste0(
        "^",
        gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", segment_id),
        "__comb.*\\.rds$"
      ),
      full.names = TRUE
    )
    
    cat(
      "\nFinished segment:", segment_id,
      "| Selected:", nrow(selected),
      "| Saved models:", length(segment_model_files),
      "\n"
    )
    
    rm(
      sim,
      selected,
      eligible,
      dataset_qc
    )
    
    gc(verbose = FALSE)
    
  } # End segment loop
  
  # ----------------------------------------------------------
  # Chromosome summary
  # ----------------------------------------------------------
  
  chromosome_models <- list.files(
    model_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )
  
  chromosome_selections <- list.files(
    chr_dir,
    pattern = "^selected_datasets_.*\\.csv$"
  )
  
  cat("\n========================================\n")
  cat("Finished chromosome:", current_chr, "\n")
  cat("Segments with selection files:",
      length(chromosome_selections), "\n")
  cat("Saved model RDS files:",
      length(chromosome_models), "\n")
  cat("========================================\n")
  
} # End chromosome loop

cat("\nAll assigned chromosomes finished.\n")
