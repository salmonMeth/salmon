
library(data.table)
library(glmmTMB)

#path to simulated methylation
sim_file <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/with_4_new_samples/simulated_with_4_extra_samples_methyl_1000per_chr.rds"
#path to simulated metadata

meta_file <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/with_4_new_samples/simulated_with_4_extra_samples_metadata.csv"

out_dir <-  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/with_4_new_samples/model_results"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


sim_methylation <- readRDS(sim_file)
meta <- fread(meta_file)

required_fields <- c(
  "loci", "sample_ids", "coverage", "Cs", "parameters"
)

missing_fields <- setdiff(
  required_fields,
  names(sim_methylation)
)

if (length(missing_fields) > 0L) {
  stop(
    "Missing fields in simulation object: ",
    paste(missing_fields, collapse = ", ")
  )
}

loci <- as.data.table(sim_methylation$loci)
coverage <- sim_methylation$coverage
Cs <- sim_methylation$Cs
sample_ids <- as.character(sim_methylation$sample_ids)

if (!identical(colnames(coverage), sample_ids)) {
  stop("sample_ids must match colnames(coverage), in the same order.")
}

if (!identical(colnames(Cs), sample_ids)) {
  stop("Cs columns do not match sample_ids.")
}

#get metadata

if (!"sample" %in% names(meta)) {
  if (all(c("Running.ID", "TimePoint") %in% names(meta))) {
    meta[, sample := paste(Running.ID, TimePoint, sep = "_")]
  } else {
    stop("Metadata must contain sample, or Running.ID and TimePoint.")
  }
}

required_meta <- c(
  "sample", "ID_ma", "ID_pa", "Weight.g"
)

missing_meta <- setdiff(required_meta, names(meta))

if (length(missing_meta) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

meta[, sample := as.character(sample)]

if (anyDuplicated(meta$sample)) {
  stop("Metadata contains duplicated sample IDs.")
}

# Match metadata to the matrix columns, retaining exact order.
meta <- meta[match(sample_ids, meta$sample)]

if (anyNA(meta$sample) ||
    !identical(meta$sample, sample_ids)) {
  missing_samples <- setdiff(sample_ids, fread(meta_file)$sample)
  
  stop(
    "Some simulated samples could not be matched to metadata. ",
    "Check sample IDs. Unmatched examples: ",
    paste(head(missing_samples, 10), collapse = ", ")
  )
}

if (anyNA(meta$Weight.g)) {
  stop("Weight.g contains missing values for fitted samples.")
}

if (any(!is.finite(meta$Weight.g))) {
  stop("Weight.g contains non-finite values.")
}


meta[, ID_ma := as.character(ID_ma)]
meta[, ID_pa := as.character(ID_pa)]

meta[, Weight.sc := scale(Weight.g)]

#fit reml=true

fit_one_model <- function(dat) {
  
  warning_messages <- character()
  
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB(
        cbind(meth.Data, unmeth.Data) ~
          Weight.sc + (1 | ID_ma) + (1 | ID_pa),
        data = dat,
        family = betabinomial(link = "logit"),
        REML = TRUE
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
    error = NA_character_,
    warnings = paste(
      unique(warning_messages),
      collapse = " | "
    )
  )
}

# get results

extract_model <- function(obj, locus_id, dat) {
  
  fit <- obj$fit
  
  row <- list(
    locus_id = locus_id,
    model = "REML",
    n_observations = nrow(dat),
    n_mothers = uniqueN(dat$ID_ma),
    n_fathers = uniqueN(dat$ID_pa),
    error = obj$error,
    warnings = obj$warnings,
    convergence_code = NA_integer_,
    pdHess = NA,
    converged = FALSE,
    beta_weight = NA_real_,
    beta_weight_se = NA_real_,
    var_maternal = NA_real_,
    var_paternal = NA_real_,
    sd_maternal = NA_real_,
    sd_paternal = NA_real_,
    dispersion_phi = NA_real_,
    dispersion_rho = NA_real_,
    logLik = NA_real_,
    AIC = NA_real_
  )
  
  if (is.null(fit)) {
    return(as.data.table(row))
  }
  
  row$convergence_code <- fit$fit$convergence
  row$pdHess <- isTRUE(fit$sdr$pdHess)
  
  row$converged <- (
    identical(fit$fit$convergence, 0L) &&
      isTRUE(fit$sdr$pdHess)
  )
  
  beta <- fixef(fit)$cond
  
  if ("Weight.sc" %in% names(beta)) {
    row$beta_weight <- unname(beta["Weight.sc"])
    
    row$beta_weight_se <- tryCatch(
      sqrt(vcov(fit)$cond["Weight.sc", "Weight.sc"]),
      error = function(e) NA_real_
    )
  }
  
  vc <- VarCorr(fit)$cond
  
  if ("ID_ma" %in% names(vc)) {
    row$var_maternal <- as.numeric(vc$ID_ma[1, 1])
    row$sd_maternal <- sqrt(row$var_maternal)
  }
  
  if ("ID_pa" %in% names(vc)) {
    row$var_paternal <- as.numeric(vc$ID_pa[1, 1])
    row$sd_paternal <- sqrt(row$var_paternal)
  }
  
  row$dispersion_phi <- tryCatch(
    as.numeric(sigma(fit)),
    error = function(e) NA_real_
  )
  
  if (is.finite(row$dispersion_phi) &&
      row$dispersion_phi > 0) {
    row$dispersion_rho <- 1 / (row$dispersion_phi + 1)
  }
  
  row$logLik <- tryCatch(
    as.numeric(logLik(fit)),
    error = function(e) NA_real_
  )
  
  row$AIC <- tryCatch(
    AIC(fit),
    error = function(e) NA_real_
  )
  
  as.data.table(row)
}

#get true params

parameters <- sim_methylation$parameters

parameters <- as.data.table(parameters)

# These columns are expected from your simulation function.
truth_cols <- c(
  "combs",
  "empirical_rate",
  "true_sd_ma",
  "true_sd_pa",
  "true_var_ma",
  "true_var_pa",
  "true_beta_weight"
)

missing_truth <- setdiff(truth_cols, names(parameters))

if (length(missing_truth) > 0L) {
  stop(
    "Missing simulation truth columns: ",
    paste(missing_truth, collapse = ", "),
    ". Inspect names(sim_methylation$parameters)."
  )
}

# fit by chr

chromosomes <- unique(as.character(loci$chr))

for (chr in chromosomes) {
  
  message("Starting chromosome ", chr)
  
  idx_chr <- which(as.character(loci$chr) == chr)
  
  chr_fits <- vector("list", length(idx_chr))
  chr_results <- vector("list", length(idx_chr))
  
  for (j in seq_along(idx_chr)) {
    
    i <- idx_chr[j]
    
    locus <- loci[i]
    
    locus_id <- paste(
      as.character(locus$chr),
      locus$start,
      locus$end,
      sep = ":"
    )
    
    cov_i <- as.numeric(coverage[i, ])
    cs_i <- as.numeric(Cs[i, ])
    
    # Keep only valid, positive-coverage observations.
    valid <- (
      !is.na(cov_i) &
        !is.na(cs_i) &
        is.finite(cov_i) &
        is.finite(cs_i) &
        cov_i > 0 &
        cs_i >= 0 &
        cs_i <= cov_i &
        cov_i == floor(cov_i) &
        cs_i == floor(cs_i)
    )
    
    dat <- copy(meta[valid])
    
    dat[, `:=`(
      meth.Data = cs_i[valid],
      unmeth.Data = cov_i[valid] - cs_i[valid]
    )]
    
    # Need enough valid observations and some variation in the data.
    if (nrow(dat) < 5L) {
      
      fit_result <- list(
        fit = NULL,
        error = "Fewer than 5 valid observations at this locus.",
        warnings = ""
      )
      
    } else {
      
      fit_result <- fit_one_model(dat)
    }
    
    chr_fits[[j]] <- list(
      locus_id = locus_id,
      REML = fit_result$fit,
      REML_error = fit_result$error,
      REML_warnings = fit_result$warnings
    )
    
    res <- extract_model(
      fit_result,
      locus_id,
      dat
    )
    
    # Add coordinates and simulation truth.
    truth <- parameters[i]
    
    res[, `:=`(
      chr = as.character(locus$chr),
      start = locus$start,
      end = locus$end,
      true_category = as.character(truth$combs),
      empirical_rate = as.numeric(truth$empirical_rate),
      true_sd_ma = as.numeric(truth$true_sd_ma),
      true_sd_pa = as.numeric(truth$true_sd_pa),
      true_var_ma = as.numeric(truth$true_var_ma),
      true_var_pa = as.numeric(truth$true_var_pa),
      true_beta_weight = as.numeric(truth$true_beta_weight)
    )]
    
    # Estimation errors.
    res[, `:=`(
      beta_weight_error =
        beta_weight - true_beta_weight,
      var_maternal_error =
        var_maternal - true_var_ma,
      var_paternal_error =
        var_paternal - true_var_pa,
      sd_maternal_error =
        sd_maternal - true_sd_ma,
      sd_paternal_error =
        sd_paternal - true_sd_pa
    )]
    
    res[, `:=`(
      beta_weight_abs_error = abs(beta_weight_error),
      var_maternal_abs_error = abs(var_maternal_error),
      var_paternal_abs_error = abs(var_paternal_error),
      sd_maternal_abs_error = abs(sd_maternal_error),
      sd_paternal_abs_error = abs(sd_paternal_error)
    )]
    
    chr_results[[j]] <- res
    
    if (j %% 100 == 0) {
      message(
        "Chromosome ", chr, ": ",
        j, " / ", length(idx_chr),
        " loci completed"
      )
    }
  }
  
  # Combine results for this chromosome.
  chr_results_dt <- rbindlist(
    chr_results,
    fill = TRUE
  )
  
  chr_safe <- gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    chr
  )
  
  csv_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_REML_results.csv")
  )
  
  rds_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_REML_results.rds")
  )
  
  fwrite(chr_results_dt, csv_file)
  
  saveRDS(
    list(
      chromosome = chr,
      results = chr_results_dt,
      fits = chr_fits
    ),
    rds_file
  )
  
  message(
    "Saved chromosome ", chr, "\n",
    "  CSV: ", csv_file, "\n",
    "  RDS: ", rds_file
  )
  
  rm(chr_results_dt, chr_results, chr_fits)
  gc()
}

message("All chromosomes finished.")


#merge the csvs

csv_files <- list.files(
  path = out_dir,
  pattern = "^chr_.*_REML_results\\.csv$",
  full.names = TRUE
)

message("Found ", length(csv_files), " chromosome CSV files.")

all_results <- rbindlist(
  lapply(csv_files, fread),
  use.names = TRUE,
  fill = TRUE
)
setorder(all_results, chr, start, end)

# Save the combined CSV
combined_file <- file.path(
  out_dir,
  "all_chromosomes_REML_results.csv"
)

fwrite(all_results, combined_file)

