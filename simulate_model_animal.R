
library(data.table)
library(glmmTMB)


simulation_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/sim_segs_complex"

metadata_file <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"

output_dir <- file.path(simulation_dir, "model_fits_T1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


meta <- fread(metadata_file)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta_T1 <- meta[TimePoint == "T1"]

# If TimePoint is stored as something other than "T1",
# this also handles values such as "1_T1" only if the
# metadata itself uses "T1" in TimePoint.
if (nrow(meta_T1) == 0) {
  stop("No T1 rows found in metadata. Check the TimePoint column.")
}

if (anyDuplicated(meta_T1$sample)) {
  stop("Duplicate T1 sample IDs found in metadata.")
}

meta_T1[, Weight.sc := as.numeric(scale(Weight.g))]

sample_metadata <- meta_T1[
  ,
  .(sample, ID_ma, ID_pa, Weight.sc)
]


simulation_files <- list.files(
  simulation_dir,
  pattern = "_simulations\\.rds$",
  full.names = TRUE
)

if (length(simulation_files) == 0) {
  stop("No simulation RDS files found in simulation_dir.")
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

if (any(grepl("^NC_", file_info$chr) == FALSE)) {
  stop("Could not extract chromosome names from some simulation filenames.")
}

chromosomes <- sort(unique(file_info$chr))

cat("Simulation files:", nrow(file_info), "\n")
cat("Chromosomes:", length(chromosomes), "\n")
print(chromosomes)

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if (nzchar(task_id)) {
  task_id <- as.integer(task_id)
  
  if (is.na(task_id) || task_id < 1) {
    stop("Invalid SLURM_ARRAY_TASK_ID.")
  }
  
  if (task_id > length(chromosomes)) {
    cat("Array task", task_id, "has no chromosome assigned. Exiting.\n")
    quit(save = "no", status = 0)
  }
  
  chromosomes <- chromosomes[task_id]
}

cat("Chromosome(s) assigned to this task:\n")
print(chromosomes)

# ============================================================
# Helper: extract random-effect SDs
# ============================================================

get_random_sd <- function(fit, group) {
  
  vc <- VarCorr(fit)$cond
  
  if (!group %in% names(vc)) {
    return(NA_real_)
  }
  
  as.numeric(attr(vc[[group]], "stddev")[1])
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
      warnings = paste(unique(warning_messages), collapse = " | ")
    ))
  }
  
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
    warnings = paste(unique(warning_messages), collapse = " | ")
  )
  
  list(
    fit = fit,
    result = result,
    error = NULL,
    warnings = paste(unique(warning_messages), collapse = " | ")
  )
}

# ============================================================
# Process each chromosome
# ============================================================

for (chr in chromosomes) {
  
  cat("\n========================================\n")
  cat("Starting chromosome:", chr, "\n")
  cat("========================================\n")
  
  chr_dir <- file.path(output_dir, chr)
  dir.create(chr_dir, recursive = TRUE, showWarnings = FALSE)
  
  results_file <- file.path(
    chr_dir,
    paste0("results_", chr, "_T1.csv")
  )
  
  errors_file <- file.path(
    chr_dir,
    paste0("errors_", chr, "_T1.csv")
  )
  
  models_file <- file.path(
    chr_dir,
    paste0("models_", chr, "_T1.rds")
  )
  
  # ----------------------------------------------------------
  # Load previous results and models for resuming
  # ----------------------------------------------------------
  
  if (file.exists(results_file)) {
    chromosome_results <- fread(results_file)
  } else {
    chromosome_results <- data.table()
  }
  
  if (file.exists(errors_file)) {
    chromosome_errors <- fread(errors_file)
  } else {
    chromosome_errors <- data.table()
  }
  
  if (file.exists(models_file)) {
    chromosome_models <- readRDS(models_file)
  } else {
    chromosome_models <- list()
  }
  
  if (nrow(chromosome_results) > 0) {
    completed <- unique(chromosome_results$key)
  } else {
    completed <- character(0)
  }
  
  # ----------------------------------------------------------
  # Simulation files for this chromosome
  # ----------------------------------------------------------
  
  chr_files <- file_info[file_info$chr == chr]
  
  cat("Segments in chromosome:", nrow(chr_files), "\n")
  
  for (file_index in seq_len(nrow(chr_files))) {
    
    simulation_file <- chr_files$file[file_index]
    segment_id <- chr_files$segment_id[file_index]
    
    cat(
      "\nSegment", file_index, "of", nrow(chr_files),
      ":", segment_id, "\n"
    )
    
    sim <- as.data.table(readRDS(simulation_file))
    
    required_columns <- c(
      "sample", "ID_ma", "ID_pa",
      "meth", "unmeth", "coverage",
      "sd_ma", "sd_pa", "sd_env",
      "phi", "combination", "replicate"
    )
    
    missing_columns <- setdiff(
      required_columns,
      names(sim)
    )
    
    if (length(missing_columns) > 0) {
      warning(
        "Skipping ", segment_id,
        "; missing columns: ",
        paste(missing_columns, collapse = ", ")
      )
      next
    }
    
    # --------------------------------------------------------
    # Add sample metadata and scaled weight
    # --------------------------------------------------------
    
    sim[, Weight.sc := sample_metadata$Weight.sc[
      match(sample, sample_metadata$sample)
    ]]
    
    # --------------------------------------------------------
    # Keep valid observations
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
    
    if (nrow(sim) == 0) {
      cat("No valid observations; skipping segment.\n")
      next
    }
    
    sim[, ID_ma := factor(ID_ma)]
    sim[, ID_pa := factor(ID_pa)]
    
    # --------------------------------------------------------
    # Identify combination-replicate datasets
    # --------------------------------------------------------
    
    datasets <- unique(
      sim[, .(combination, replicate)]
    )
    
    cat("Datasets in segment:", nrow(datasets), "\n")
    
    for (j in seq_len(nrow(datasets))) {
      
      combo <- datasets$combination[j]
      rep_id <- datasets$replicate[j]
      
      key <- paste(
        segment_id,
        combo,
        rep_id,
        sep = "__"
      )
      
      # Skip successfully completed fits
      if (key %in% completed) {
        next
      }
      
      dat <- sim[
        sim$combination == combo &
          sim$replicate == rep_id
      ]
      
      # Require enough samples
      if (uniqueN(dat$sample) < 42) {
        next
      }
      
      # Ensure one observation per sample
      if (anyDuplicated(dat$sample)) {
        chromosome_errors <- rbind(
          chromosome_errors,
          data.table(
            segment_id = segment_id,
            combination = combo,
            replicate = rep_id,
            error = "Duplicate sample observations",
            warnings = ""
          ),
          fill = TRUE
        )
        
        fwrite(chromosome_errors, errors_file)
        next
      }
      
      # ------------------------------------------------------
      # Fit model
      # ------------------------------------------------------
      
      fit_output <- fit_one_dataset(dat)
      
      if (is.null(fit_output$fit)) {
        
        chromosome_errors <- rbind(
          chromosome_errors,
          data.table(
            segment_id = segment_id,
            combination = combo,
            replicate = rep_id,
            error = fit_output$error,
            warnings = fit_output$warnings
          ),
          fill = TRUE
        )
        
        fwrite(chromosome_errors, errors_file)
        
        next
      }
      
      fit <- fit_output$fit
      
      # ------------------------------------------------------
      # Save model
      # ------------------------------------------------------
      
      chromosome_models[[key]] <- fit
      
      saveRDS(
        chromosome_models,
        models_file
      )
      
      # ------------------------------------------------------
      # Save result row
      # ------------------------------------------------------
      
      result <- copy(fit_output$result)
      
      result[, `:=`(
        key = key,
        segment_id = segment_id,
        combination = combo,
        replicate = rep_id,
        n_samples = uniqueN(dat$sample),
        n_mothers = uniqueN(dat$ID_ma),
        n_fathers = uniqueN(dat$ID_pa),
        true_sd_ma = unique(dat$sd_ma)[1],
        true_sd_pa = unique(dat$sd_pa)[1],
        true_sd_env = unique(dat$sd_env)[1],
        true_phi = unique(dat$phi)[1]
      )]
      
      # Put identifying columns first
      setcolorder(
        result,
        c(
          "key",
          "segment_id",
          "combination",
          "replicate",
          "n_samples",
          "n_mothers",
          "n_fathers",
          "true_sd_ma",
          "true_sd_pa",
          "true_sd_env",
          "true_phi",
          setdiff(
            names(result),
            c(
              "key",
              "segment_id",
              "combination",
              "replicate",
              "n_samples",
              "n_mothers",
              "n_fathers",
              "true_sd_ma",
              "true_sd_pa",
              "true_sd_env",
              "true_phi"
            )
          )
        )
      )
      
      chromosome_results <- rbind(
        chromosome_results,
        result,
        fill = TRUE
      )
      
      completed <- c(completed, key)
      
      # Checkpoint after every successful fit
      fwrite(
        chromosome_results,
        results_file
      )
      
      cat(
        "Completed:", key,
        "| Total successful fits:", nrow(chromosome_results),
        "\n"
      )
      
    } # combination-replicate loop
    
    # Save after each segment as well
    saveRDS(chromosome_models, models_file)
    
    fwrite(chromosome_results, results_file)
    
    if (nrow(chromosome_errors) > 0) {
      fwrite(chromosome_errors, errors_file)
    }
    
    rm(sim)
    gc(verbose = FALSE)
    
  } # segment loop
  
  cat("\nFinished chromosome:", chr, "\n")
  cat("Successful fits:", nrow(chromosome_results), "\n")
  cat("Errors:", nrow(chromosome_errors), "\n")
  
} # chromosome loop

cat("\nAll assigned chromosomes finished.\n")
