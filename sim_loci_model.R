fit_one_model <- function(dat, reml) {
  warning_messages <- character()
  
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB(
        REML = reml,
        cbind(meth.Data, unmeth.Data) ~
          Weight.g + (1 | ID_ma) + (1 | ID_pa),
        data = dat,
        family = betabinomial(link = "logit")
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages, conditionMessage(w)
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
  
  list(
    fit = fit,
    error = NA_character_,
    warnings = paste(unique(warning_messages), collapse = " | ")
  )
}


extract_model <- function(obj, model_name, locus_id, dat) {
  
  fit <- obj$fit
  
  row <- list(
    locus_id = locus_id,
    model = model_name,
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
  
  if ("Weight.g" %in% names(beta)) {
    row$beta_weight <- unname(beta["Weight.g"])
    
    row$beta_weight_se <- tryCatch(
      sqrt(vcov(fit)$cond["Weight.g", "Weight.g"]),
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
  
  row$logLik <- as.numeric(logLik(fit))
  row$AIC <- AIC(fit)
  
  as.data.table(row)
}


library(data.table)
library(glmmTMB)



sim_methylation <- readRDS(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/simulated_T1_parental_methylation_1000per_chr.rds"
)


path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]

meta <- meta[match(sim_methylation$sample_ids, meta$sample)]


out_dir <- file.path(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci",
  "model_results_by_chromosome"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

loci <- as.data.table(sim_methylation$loci)
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
    
    cov_i <- as.numeric(sim_methylation$coverage[i, ])
    cs_i <- as.numeric(sim_methylation$Cs[i, ])
    ts_i <- cov_i - cs_i
    
    dat <- copy(meta)
    dat[, `:=`(
      meth.Data = cs_i,
      unmeth.Data = ts_i
    )]
    
    
    truth <- sim_methylation$parameters[i]
    
    # Store both models, even if one fails.
    fit_reml <- fit_one_model(dat, reml = TRUE)
    fit_ml <- fit_one_model(dat, reml = FALSE)
    
    chr_fits[[j]] <- list(
      locus_id = locus_id,
      REML = fit_reml$fit,
      ML = fit_ml$fit,
      REML_error = fit_reml$error,
      ML_error = fit_ml$error,
      REML_warnings = fit_reml$warnings,
      ML_warnings = fit_ml$warnings
    )
    
    r1 <- extract_model(fit_reml, "REML", locus_id, dat)
    r2 <- extract_model(fit_ml, "ML", locus_id, dat)
    
    res <- rbindlist(list(r1, r2), fill = TRUE)
    
    res[, `:=`(
      chr = as.character(locus$chr),
      start = locus$start,
      end = locus$end,
      true_category = as.character(truth$combs),
      empirical_rate = truth$empirical_rate,
      true_sd_ma = truth$true_sd_ma,
      true_sd_pa = truth$true_sd_pa,
      true_var_ma = truth$true_var_ma,
      true_var_pa = truth$true_var_pa,
      true_beta_weight = truth$true_beta_weight
    )]
    
    # Estimation errors
    res[, `:=`(
      beta_weight_error = beta_weight - true_beta_weight,
      var_maternal_error = var_maternal - true_var_ma,
      var_paternal_error = var_paternal - true_var_pa,
      sd_maternal_error = sd_maternal - true_sd_ma,
      sd_paternal_error = sd_paternal - true_sd_pa
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
        j, " / ", length(idx_chr), " loci completed"
      )
    }
  }
  
  # Combine results for THIS chromosome only
  chr_results_dt <- rbindlist(chr_results, fill = TRUE)
  
  # One CSV and one RDS per chromosome
  chr_safe <- gsub("[^A-Za-z0-9_.-]", "_", chr)
  
  csv_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_model_results.csv")
  )
  
  rds_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_model_results.rds")
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
  
  # Release memory before the next chromosome
  rm(chr_results_dt, chr_results, chr_fits)
  gc()
}

message("All chromosomes finished.")

