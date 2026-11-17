
library(data.table)
library(brms)
library(posterior)

# --------------------------------------------------
# Directories and input files
# --------------------------------------------------

data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/animal_model"

input_files <- c(
  "simulation_results_31-62_BINOMIAL_MODEL.rds",
  "simulation_results_63-100_BINOMIAL_MODEL.rds",
  "simulation_results_900_segments_BINOMIAL_MODEL.rds"
)

output_file <- file.path(
  data_dir,
  "animal_model_parameter_convergence_all.csv"
)

# Start a fresh output file
if (file.exists(output_file)) {
  file.remove(output_file)
}

# --------------------------------------------------
# Helper: safely extract a parameter estimate
# --------------------------------------------------

get_parameter <- function(tab, parameter, column) {
  
  tab <- as.data.frame(tab)
  
  if (is.null(rownames(tab)) || !parameter %in% rownames(tab)) {
    return(NA_real_)
  }
  
  if (!column %in% colnames(tab)) {
    return(NA_real_)
  }
  
  as.numeric(tab[parameter, column])
}

# --------------------------------------------------
# Process input files
# --------------------------------------------------

for (file in input_files) {
  
  input_path <- file.path(data_dir, file)
  
  cat("\n========================================\n")
  cat("Loading:", file, "\n")
  cat("========================================\n")
  
  results <- readRDS(input_path)
  
  cat("Models:", length(results), "\n")
  
  for (i in seq_along(results)) {
    
    cat("Processing", i, "of", length(results), "\n")
    
    segment_id <- results[[i]]$segment_id
    fit <- results[[i]]$fit
    truth <- results[[i]]$true_parameters
    
    # ----------------------------------------------
    # Model summary
    # ----------------------------------------------
    
    sm <- summary(fit)
    
    fixed <- sm$fixed
    animal <- sm$random$animal
    
    # Fixed effects
    intercept_est <- get_parameter(fixed, "Intercept", "Estimate")
    intercept_low <- get_parameter(fixed, "Intercept", "l-95% CI")
    intercept_high <- get_parameter(fixed, "Intercept", "u-95% CI")
    
    weight_est <- get_parameter(fixed, "Weight_sc", "Estimate")
    weight_low <- get_parameter(fixed, "Weight_sc", "l-95% CI")
    weight_high <- get_parameter(fixed, "Weight_sc", "u-95% CI")
    
    # Animal SD
    animal <- as.data.frame(animal)
    
    animal_sd_est <- animal$Estimate[1]
    animal_sd_low <- animal$`l-95% CI`[1]
    animal_sd_high <- animal$`u-95% CI`[1]
    
    # ----------------------------------------------
    # Full convergence diagnostics
    # Includes individual animal effects
    # ----------------------------------------------
    
    draws <- posterior::as_draws_df(fit)
    
    diagnostics <- posterior::summarise_draws(
      draws,
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    )
    
    max_rhat <- max(diagnostics$rhat, na.rm = TRUE)
    min_bulk <- min(diagnostics$ess_bulk, na.rm = TRUE)
    min_tail <- min(diagnostics$ess_tail, na.rm = TRUE)
    
    n_bad_rhat <- sum(diagnostics$rhat > 1.01, na.rm = TRUE)
    n_low_bulk <- sum(diagnostics$ess_bulk < 400, na.rm = TRUE)
    n_low_tail <- sum(diagnostics$ess_tail < 400, na.rm = TRUE)
    
    # ----------------------------------------------
    # Divergent transitions
    # ----------------------------------------------
    
    np <- brms::nuts_params(fit)
    
    divergences <- sum(
      np$Parameter == "divergent__" & np$Value == 1
    )
    
    # ----------------------------------------------
    # True simulation parameters
    # ----------------------------------------------
    
    VA_true <- as.numeric(truth$VA_true[1])
    SD_true <- sqrt(VA_true)
    
    beta0_true <- as.numeric(truth$beta0_true[1])
    beta1_true <- as.numeric(truth$beta1_true[1])
    
    phi_true <- if ("phi_true" %in% names(truth)) {
      as.numeric(truth$phi_true[1])
    } else {
      NA_real_
    }
    
    mean_methylation <- as.numeric(
      truth$mean_methylation[1]
    )
    
    # ----------------------------------------------
    # Model family
    # ----------------------------------------------
    
    fam <- family(fit)
    
    # ----------------------------------------------
    # One row per segment
    # ----------------------------------------------
    
    out <- data.table(
      
      segment_id = segment_id,
      source_file = file,
      
      model_family = fam$family,
      link = fam$link,
      
      # Animal SD
      SD_true = SD_true,
      VA_true = VA_true,
      SD_estimated = animal_sd_est,
      SD_lower = animal_sd_low,
      SD_upper = animal_sd_high,
      
      SD_absolute_error = abs(animal_sd_est - SD_true),
      
      SD_CI_contains_truth =
        SD_true >= animal_sd_low &
        SD_true <= animal_sd_high,
      
      # Intercept
      beta0_true = beta0_true,
      beta0_estimated = intercept_est,
      beta0_lower = intercept_low,
      beta0_upper = intercept_high,
      
      # Weight effect
      beta1_true = beta1_true,
      beta1_estimated = weight_est,
      beta1_lower = weight_low,
      beta1_upper = weight_high,
      
      # Other simulation information
      phi_true = phi_true,
      mean_methylation = mean_methylation,
      
      # Convergence
      n_parameters = nrow(diagnostics),
      max_Rhat = max_rhat,
      min_bulk_ESS = min_bulk,
      min_tail_ESS = min_tail,
      
      n_Rhat_above_1.01 = n_bad_rhat,
      n_bulk_ESS_below_400 = n_low_bulk,
      n_tail_ESS_below_400 = n_low_tail,
      
      divergent_transitions = divergences
    )
    
    # ----------------------------------------------
    # Append results immediately
    # ----------------------------------------------
    
    fwrite(
      out,
      output_file,
      append = file.exists(output_file)
    )
    
    # Release memory
    rm(
      fit, sm, fixed, animal, draws,
      diagnostics, np, out, truth
    )
    
    gc(verbose = FALSE)
  }
  
  rm(results)
  gc(verbose = FALSE)
  
  cat("Finished:", file, "\n")
}

cat("\n========================================\n")
cat("ALL FILES PROCESSED\n")
cat("Output:", output_file, "\n")
cat("========================================\n")
