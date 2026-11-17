#dharma to the loci
library(data.table)
library(DHARMa)
library(glmmTMB)

model_dir="/scratch/project_2010912/ezel/salmon/data_files/beta_binomial_for_actual_loci_data"
out_dir <- file.path(model_dir, "DHARMa_T1_sampled_loci")
plot_dir <- file.path(out_dir, "plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

model_files <- list.files(
  model_dir,
  pattern = "^chr_NC_.*_REML_results\\.rds$",
  full.names = TRUE
)

stopifnot(length(model_files) > 0)

extract_test <- function(test_fun, sim_res) {
  tryCatch({
    test <- test_fun(sim_res)
    
    statistic <- if (!is.null(test$statistic)) {
      as.numeric(test$statistic[1])
    } else {
      NA_real_
    }
    
    pvalue <- if (!is.null(test$p.value)) {
      as.numeric(test$p.value[1])
    } else {
      NA_real_
    }
    
    alternative <- if (!is.null(test$alternative)) {
      as.character(test$alternative[1])
    } else {
      NA_character_
    }
    
    list(
      statistic = statistic,
      pvalue = pvalue,
      alternative = alternative,
      status = "OK"
    )
  }, error = function(e) {
    list(
      statistic = NA_real_,
      pvalue = NA_real_,
      alternative = NA_character_,
      status = paste0("ERROR: ", conditionMessage(e))
    )
  })
}


results <- list()

for (file in model_files) {
  
  models <- readRDS(file)
  chr_name <- models$chromosome
  
  cat("\nProcessing", chr_name, "\n")
  
  chr_fits <- models$fits
  
  for (i in seq_along(chr_fits)) {
    
    model_entry <- chr_fits[[i]]
    
    locus_id <- model_entry$locus_id
    model <- model_entry$fit
    
    coords <- strsplit(locus_id, ":", fixed = TRUE)[[1]]
    
    chr <- coords[1]
    start <- as.integer(coords[2])
    end <- as.integer(coords[3])
    
    cat(
      "  Locus", i, "of", length(chr_fits),
      ":", locus_id, "\n"
    )
    
    res <- data.table(
      chr = chr,
      start = start,
      end = end,
      locus_id = locus_id,
      model_status = "NOT_RUN",
      uniformity_statistic = NA_real_,
      uniformity_pvalue = NA_real_,
      uniformity_status = "NOT_RUN",
      dispersion_statistic = NA_real_,
      dispersion_pvalue = NA_real_,
      dispersion_status = "NOT_RUN",
      outliers_statistic = NA_real_,
      outliers_pvalue = NA_real_,
      outliers_status = "NOT_RUN",
      plot_status = "NOT_RUN"
    )
    
    if (is.null(model)) {
      
      res$model_status <- paste0(
        "MODEL_FAILED: ",
        model_entry$error
      )
      
      results[[length(results) + 1]] <- res
      next
    }
    
    if (
      is.null(model$sdr) ||
      !isTRUE(model$sdr$pdHess)
    ) {
      
      res$model_status <- "NON_POSITIVE_HESSIAN"
      
      results[[length(results) + 1]] <- res
      next
    }
    
    res$model_status <- "OK"
    dat <- model_entry$data
    
    if (!is.null(dat) && nrow(dat) > 0) {
      coverage <- dat$meth.Data + dat$unmeth.Data
      methylation <- dat$meth.Data / coverage
      
      valid <- is.finite(methylation) & coverage > 0
      
      res[, `:=`(
        n_samples = nrow(dat),
        n_mothers = uniqueN(dat$ID_ma),
        n_fathers = uniqueN(dat$ID_pa),
        
        total_coverage = sum(coverage, na.rm = TRUE),
        mean_coverage = mean(coverage, na.rm = TRUE),
        median_coverage = median(coverage, na.rm = TRUE),
        sd_coverage = sd(coverage, na.rm = TRUE),
        min_coverage = min(coverage, na.rm = TRUE),
        max_coverage = max(coverage, na.rm = TRUE),
        
        mean_sample_methylation = mean(methylation[valid], na.rm = TRUE),
        median_sample_methylation = median(methylation[valid], na.rm = TRUE),
        sd_sample_methylation = sd(methylation[valid], na.rm = TRUE),
        min_sample_methylation = min(methylation[valid], na.rm = TRUE),
        max_sample_methylation = max(methylation[valid], na.rm = TRUE),
        sample_methylation_range =
          max(methylation[valid], na.rm = TRUE) -
          min(methylation[valid], na.rm = TRUE),
        
        weighted_methylation =
          sum(dat$meth.Data, na.rm = TRUE) /
          sum(coverage, na.rm = TRUE),
        
        n_zero_coverage = sum(coverage == 0, na.rm = TRUE),
        n_valid_methylation_samples = sum(valid)
      )]
    }
    
    sim_res <- tryCatch({
      
      set.seed(123 + i)
      
      simulateResiduals(
        fittedModel = model,
        n = 1000
      )
      
    }, error = function(e) {
      NULL
    })
    
    if (is.null(sim_res)) {
      
      res$model_status <- "DHARMa_SIMULATION_FAILED"
      
      results[[length(results) + 1]] <- res
      next
    }
    
    uniformity <- extract_test(
      testUniformity,
      sim_res
    )
    
    res[, `:=`(
      uniformity_statistic = uniformity$statistic,
      uniformity_pvalue = uniformity$pvalue,
      uniformity_status = uniformity$status
    )]
    
    dispersion <- extract_test(testDispersion, sim_res)
    
    res[, `:=`(
      dispersion_statistic = dispersion$statistic,
      dispersion_pvalue = dispersion$pvalue,
      dispersion_alternative = dispersion$alternative,
      dispersion_status = dispersion$status
    )]
    
    outliers <- extract_test(
      testOutliers,
      sim_res
    )
    
    res[, `:=`(
      outliers_statistic = outliers$statistic,
      outliers_pvalue = outliers$pvalue,
      outliers_status = outliers$status
    )]
    
    plot_file <- file.path(
      plot_dir,
      paste0(gsub(":", "_", locus_id), "_DHARMa.pdf")
    )
    
    res$plot_status <- tryCatch({
      
      pdf(
        plot_file,
        width = 10,
        height = 8
      )
      
      plot(sim_res)
      
      dev.off()
      
      "OK"
      
    }, error = function(e) {
      
      if (dev.cur() > 1) {
        dev.off()
      }
      
      paste0("ERROR: ", conditionMessage(e))
    })
    
    results[[length(results) + 1]] <- res
    
    if (i %% 100 == 0) {
      cat(
        "  Completed", i, "of", length(chr_fits),
        "loci\n"
      )
    }
  }
  
  rm(models, chr_fits)
  gc()
}

results <- rbindlist(
  results,
  fill = TRUE
)

setorder(
  results,
  chr, start, end
)
dharma_results=results
alpha <- 0.05

dharma_results[, `:=`(
  uniformity_flag =
    !is.na(uniformity_pvalue) & uniformity_pvalue < alpha,
  
  outlier_flag =
    !is.na(outliers_pvalue) & outliers_pvalue < alpha,
  
  dispersion_flag =
    !is.na(dispersion_pvalue) &
    !is.na(dispersion_statistic) &
    dispersion_pvalue < alpha &
    (dispersion_statistic > 1.2 | dispersion_statistic < 0.8),
  
  overdispersed =
    !is.na(dispersion_pvalue) &
    !is.na(dispersion_statistic) &
    dispersion_pvalue < alpha &
    dispersion_statistic > 1.2,
  
  underdispersed =
    !is.na(dispersion_pvalue) &
    !is.na(dispersion_statistic) &
    dispersion_pvalue < alpha &
    dispersion_statistic < 0.8
)]

dharma_results[, diagnostic_class := fcase(
  model_status != "OK", "model_or_simulation_failed",
  is.na(uniformity_pvalue) &
    is.na(dispersion_pvalue) &
    is.na(outliers_pvalue), "no_valid_tests",
  uniformity_flag | dispersion_flag | outlier_flag, "flagged",
  default = "unflagged"
)]

chromosome_summary <- dharma_results[
  , .(
    n_loci = .N,
    n_models_ok = sum(model_status == "OK"),
    n_models_failed = sum(model_status != "OK"),
    
    n_uniformity_tested = sum(!is.na(uniformity_pvalue)),
    n_uniformity_flagged = sum(uniformity_flag, na.rm = TRUE),
    pct_uniformity_flagged = if (sum(!is.na(uniformity_pvalue)) > 0)
      100 * mean(uniformity_flag[!is.na(uniformity_pvalue)])
    else NA_real_,
    
    n_dispersion_tested = sum(
      !is.na(dispersion_pvalue) & !is.na(dispersion_statistic)
    ),
    n_overdispersed = sum(overdispersed, na.rm = TRUE),
    n_underdispersed = sum(underdispersed, na.rm = TRUE),
    n_dispersion_flagged = sum(dispersion_flag, na.rm = TRUE),
    
    pct_overdispersed = if (sum(
      !is.na(dispersion_pvalue) & !is.na(dispersion_statistic)
    ) > 0)
      100 * mean(overdispersed[
        !is.na(dispersion_pvalue) & !is.na(dispersion_statistic)
      ])
    else NA_real_,
    
    pct_underdispersed = if (sum(
      !is.na(dispersion_pvalue) & !is.na(dispersion_statistic)
    ) > 0)
      100 * mean(underdispersed[
        !is.na(dispersion_pvalue) & !is.na(dispersion_statistic)
      ])
    else NA_real_,
    
    n_outlier_tested = sum(!is.na(outliers_pvalue)),
    n_outlier_flagged = sum(outlier_flag, na.rm = TRUE),
    
    mean_coverage = mean(mean_coverage, na.rm = TRUE),
    median_coverage = median(median_coverage, na.rm = TRUE),
    mean_weighted_methylation = mean(weighted_methylation, na.rm = TRUE),
    mean_methylation_sd = mean(sd_sample_methylation, na.rm = TRUE),
    mean_methylation_range = mean(sample_methylation_range, na.rm = TRUE)
  ),
  by = chr
]

setorder(chromosome_summary, chr)

fwrite(
  chromosome_summary,
  file.path(out_dir, "DHARMa_T1_chromosome_summary.csv")
)

overall_summary <- dharma_results[
  , .(
    n_loci = .N,
    n_models_ok = sum(model_status == "OK"),
    n_models_failed = sum(model_status != "OK"),
    
    n_uniformity_tested = sum(!is.na(uniformity_pvalue)),
    n_uniformity_flagged = sum(uniformity_flag, na.rm = TRUE),
    
    n_dispersion_tested = sum(
      !is.na(dispersion_pvalue) & !is.na(dispersion_statistic)
    ),
    n_overdispersed = sum(overdispersed, na.rm = TRUE),
    n_underdispersed = sum(underdispersed, na.rm = TRUE),
    n_dispersion_flagged = sum(dispersion_flag, na.rm = TRUE),
    
    n_outlier_tested = sum(!is.na(outliers_pvalue)),
    n_outlier_flagged = sum(outlier_flag, na.rm = TRUE),
    
    mean_coverage = mean(mean_coverage, na.rm = TRUE),
    median_coverage = median(median_coverage, na.rm = TRUE),
    sd_coverage = sd(mean_coverage, na.rm = TRUE),
    
    mean_weighted_methylation = mean(weighted_methylation, na.rm = TRUE),
    sd_weighted_methylation = sd(weighted_methylation, na.rm = TRUE),
    
    mean_methylation_sd = mean(sd_sample_methylation, na.rm = TRUE),
    mean_methylation_range = mean(sample_methylation_range, na.rm = TRUE)
  )
]

fwrite(
  overall_summary,
  file.path(out_dir, "DHARMa_T1_overall_summary.csv")
)

fwrite(
  results,
  file.path(out_dir, "DHARMa_T1_sampled_loci.csv")
)

cat("\n-----------------------------\n")
cat("DHARMa analysis completed\n")
cat("-----------------------------\n")

cat("Loci processed:", nrow(results), "\n")

cat(
  "Models successfully tested:",
  sum(results$model_status == "OK"),
  "\n"
)

cat(
  "Plots saved:",
  sum(results$plot_status == "OK"),
  "\n"
)

cat(
  "Uniformity significant:",
  sum(results$uniformity_pvalue < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Dispersion significant:",
  sum(results$dispersion_pvalue < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Outlier tests significant:",
  sum(results$outliers_pvalue < 0.05, na.rm = TRUE),
  "\n"
)

cat("\nOutput:", out_dir, "\n")