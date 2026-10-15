library(data.table)
library(glmmTMB)


fit_one_model <- function(dat, family_col = "Fam.ID") {
  warning_messages <- character()
  
  model_formula <- as.formula(
    paste0(
      "cbind(meth.Data, unmeth.Data) ~ ",
      "Weight.g + (1 | `", family_col, "`)"
    )
  )
  
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB(
        formula = model_formula,
        data = dat,
        REML = TRUE,
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
    error = NA_character_,
    warnings = paste(
      unique(warning_messages),
      collapse = " | "
    )
  )
}


extract_model <- function(
    obj, locus_id, dat,
    family_col = "Fam.ID"
) {
  fit <- obj$fit
  
  row <- list(
    locus_id = locus_id,
    model = "REML",
    n_observations = nrow(dat),
    n_families = uniqueN(dat[[family_col]]),
    error = obj$error,
    warnings = obj$warnings,
    convergence_code = NA_integer_,
    pdHess = NA,
    converged = FALSE,
    beta_weight = NA_real_,
    beta_weight_se = NA_real_,
    var_family = NA_real_,
    sd_family = NA_real_,
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
  
  if ("Weight.g" %in% names(beta)) {
    row$beta_weight <- unname(beta["Weight.g"])
    
    row$beta_weight_se <- tryCatch(
      sqrt(vcov(fit)$cond["Weight.g", "Weight.g"]),
      error = function(e) NA_real_
    )
  }
  
  vc <- VarCorr(fit)$cond
  
  if (family_col %in% names(vc)) {
    row$var_family <- as.numeric(
      vc[[family_col]][1, 1]
    )
    row$sd_family <- sqrt(row$var_family)
  }
  
  row$dispersion_phi <- tryCatch(
    as.numeric(sigma(fit)),
    error = function(e) NA_real_
  )
  
  if (
    is.finite(row$dispersion_phi) &&
    row$dispersion_phi > 0
  ) {
    row$dispersion_rho <- 1 / (row$dispersion_phi + 1)
  }
  
  row$logLik <- as.numeric(logLik(fit))
  row$AIC <- AIC(fit)
  
  as.data.table(row)
}


sim_family <- readRDS(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/simulated_T1_FAM_methylation_1000per_chr.rds"
)

family_col <- sim_family$family_col  # "Fam.ID"

path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]

meta <- meta[match(sim_family$sample_ids, meta$sample)]

out_dir <- file.path(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci",
  "family_model_results_by_chromosome"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

loci <- as.data.table(sim_family$loci)
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
    
    cov_i <- as.numeric(sim_family$coverage[i, ])
    cs_i <- as.numeric(sim_family$Cs[i, ])
    ts_i <- cov_i - cs_i
    
    dat <- copy(meta)
    
    dat[, `:=`(
      meth.Data = cs_i,
      unmeth.Data = ts_i
    )]
    
    fit_reml <- fit_one_model(
      dat,
      family_col = family_col
    )
    
    chr_fits[[j]] <- list(
      locus_id = locus_id,
      REML = fit_reml$fit,
      REML_error = fit_reml$error,
      REML_warnings = fit_reml$warnings
    )
    
    res <- extract_model(
      fit_reml,
      locus_id = locus_id,
      dat = dat,
      family_col = family_col
    )
    
    # get true params
    truth <- sim_family$parameters[i]
    
    res[, `:=`(
      chr = as.character(locus$chr),
      start = locus$start,
      end = locus$end,
      true_category = as.character(truth$combs),
      empirical_methyl_rate = truth$empirical_methyl_rate,
      true_sd_family = truth$true_sd_family,
      true_var_family = truth$true_var_family,
      true_beta_weight = truth$true_beta_weight
    )]
    
    # Estimation errors
    res[, `:=`(
      beta_weight_error =
        beta_weight - true_beta_weight,
      var_family_error =
        var_family - true_var_family,
      sd_family_error =
        sd_family - true_sd_family
    )]
    
    res[, `:=`(
      beta_weight_abs_error = abs(beta_weight_error),
      var_family_abs_error = abs(var_family_error),
      sd_family_abs_error = abs(sd_family_error)
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
  
  # Combine this chromosome's results
  chr_results_dt <- rbindlist(
    chr_results,
    fill = TRUE
  )
  
  chr_safe <- gsub(
    "[^A-Za-z0-9_.-]", "_", chr
  )
  
  csv_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_family_REML_results.csv")
  )
  
  rds_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_family_REML_results.rds")
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
