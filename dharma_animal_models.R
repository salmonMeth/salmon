library(brms)
library(DHARMa)
library(data.table)


model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/animal_model/final_models/animal_model_T1_all"

output_dir <- file.path(model_dir, "DHARMa_results")

n_simulations <- 500

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "DHARMa_objects"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"),
           recursive = TRUE, showWarnings = FALSE)

results_file <- file.path(output_dir, "DHARMa_results.csv")
errors_file  <- file.path(output_dir, "DHARMa_errors.csv")


rds_files <- list.files(
  path = model_dir,
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

rds_files <- sort(rds_files)

# Avoid accidentally processing DHARMa output files
rds_files <- rds_files[
  !grepl("/DHARMa_results/", rds_files, fixed = TRUE)
]

cat("Number of model files found:", length(rds_files), "\n")

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

# Extract a numeric statistic or p-value from a DHARMa test
extract_number <- function(x, field) {
  
  value <- x[[field]]
  
  if (is.null(value) || length(value) == 0) {
    return(NA_real_)
  }
  
  as.numeric(value[1])
}


# Run an individual DHARMa test safely
run_dharma_test <- function(test_function, dh) {
  
  tryCatch({
    
    test_result <- test_function(dh)
    
    list(
      statistic = extract_number(test_result, "statistic"),
      p_value = extract_number(test_result, "p.value"),
      error = NA_character_
    )
    
  }, error = function(e) {
    
    list(
      statistic = NA_real_,
      p_value = NA_real_,
      error = conditionMessage(e)
    )
    
  })
}


# Append one result row to a CSV
append_result <- function(row, filename) {
  
  fwrite(
    row,
    filename,
    append = file.exists(filename),
    col.names = !file.exists(filename)
  )
}


for (i in seq_along(rds_files)) {
  
  model_file <- rds_files[i]
  
  model_name <- tools::file_path_sans_ext(
    basename(model_file)
  )
  
  cat("\n--------------------------------------------\n")
  cat("Processing model", i, "of", length(rds_files), "\n")
  cat(model_name, "\n")
  cat("--------------------------------------------\n")
  
  start_time <- Sys.time()
  
  tryCatch({
    
    # --------------------------------------------------------
    # Load model
    # --------------------------------------------------------
    
    fit <- readRDS(model_file)
    
    cat("Model loaded\n")
    
    # Check response data
    if (is.null(fit$data$meth) ||
        is.null(fit$data$coverage)) {
      stop("Model data does not contain meth and coverage")
    }
    
    observed <- fit$data$meth
    
    # --------------------------------------------------------
    # Generate posterior predictive simulations
    # --------------------------------------------------------
    
    cat("Generating", n_simulations,
        "posterior predictive simulations...\n")
    
    sim_matrix <- posterior_predict(
      fit,
      ndraws = n_simulations
    )
    
    # DHARMa expects observations in rows and simulations
    # in columns for createDHARMa()
    sim_matrix <- t(sim_matrix)
    
    # Check dimensions
    if (nrow(sim_matrix) != length(observed)) {
      stop(
        "Simulation dimensions do not match observed response: ",
        nrow(sim_matrix), " simulated observations versus ",
        length(observed), " observed observations"
      )
    }
    
    cat("Simulation dimensions:",
        nrow(sim_matrix), "observations x",
        ncol(sim_matrix), "simulations\n")
    
    # --------------------------------------------------------
    # Create DHARMa object
    # --------------------------------------------------------
    
    dh <- createDHARMa(
      simulatedResponse = sim_matrix,
      observedResponse = observed,
      fittedPredicted = rowMeans(sim_matrix),
      integerResponse = TRUE
    )
    
    # --------------------------------------------------------
    # Run diagnostic tests
    # --------------------------------------------------------
    
    cat("Running DHARMa tests...\n")
    
    uniformity_result <- run_dharma_test(
      testUniformity, dh
    )
    
    dispersion_result <- run_dharma_test(
      testDispersion, dh
    )
    
    outlier_result <- run_dharma_test(
      testOutliers, dh
    )

    
    object_file <- file.path(
      output_dir,
      "DHARMa_objects",
      paste0(model_name, "_DHARMa.rds")
    )
    
    saveRDS(dh, object_file)

    
    plot_file <- file.path(
      output_dir,
      "plots",
      paste0(model_name, "_DHARMa.png")
    )
    
    png(
      filename = plot_file,
      width = 1800,
      height = 1400,
      res = 180
    )
    
    plot(dh)
    dev.off()

    
    elapsed <- as.numeric(
      difftime(Sys.time(), start_time, units = "mins")
    )
    
    result_row <- data.table(
      model_file = model_file,
      model_name = model_name,
      n_observations = length(observed),
      n_simulations = ncol(sim_matrix),
      
      uniformity_statistic = uniformity_result$statistic,
      uniformity_p = uniformity_result$p_value,
      
      dispersion_statistic = dispersion_result$statistic,
      dispersion_p = dispersion_result$p_value,
      
      outlier_statistic = outlier_result$statistic,
      outlier_p = outlier_result$p_value,
      
      elapsed_minutes = elapsed,
      status = "success"
    )
    
    append_result(result_row, results_file)
    
    # Save individual test errors, if any
    test_errors <- c(
      uniformity = uniformity_result$error,
      dispersion = dispersion_result$error,
      outliers = outlier_result$error
    )
    
    failed_tests <- names(test_errors)[
      !is.na(test_errors)
    ]
    
    if (length(failed_tests) > 0) {
      
      error_row <- data.table(
        model_name = model_name,
        stage = paste(failed_tests, collapse = ";"),
        error = paste(
          test_errors[failed_tests],
          collapse = " | "
        )
      )
      
      append_result(error_row, errors_file)
    }
    
    cat("SUCCESS:", model_name, "\n")
    cat("Uniformity p-value:",
        uniformity_result$p_value, "\n")
    cat("Dispersion p-value:",
        dispersion_result$p_value, "\n")
    cat("Outlier p-value:",
        outlier_result$p_value, "\n")
    cat("Elapsed minutes:", round(elapsed, 2), "\n")

    
    rm(fit, sim_matrix, dh)
    gc()
    
  }, error = function(e) {
    
    error_message <- conditionMessage(e)
    
    cat("ERROR:", model_name, "\n")
    cat(error_message, "\n")
    
    error_row <- data.table(
      model_name = model_name,
      stage = "model_processing",
      error = error_message
    )
    
    append_result(error_row, errors_file)
    
    # Clean up and continue
    if (exists("fit")) rm(fit)
    if (exists("sim_matrix")) rm(sim_matrix)
    if (exists("dh")) rm(dh)
    
    gc()
  })
}


cat("\n============================================\n")
cat("DHARMa analysis finished\n")
cat("Models found:", length(rds_files), "\n")
cat("Results:", results_file, "\n")
cat("Errors:", errors_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("============================================\n")