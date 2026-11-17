#dharma for the fucked cov 
library(DHARMa)
library(data.table)

base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/data_with_lowered_coverage"

model_dir <- file.path(base_dir, "model_results")
plot_dir <- file.path(base_dir, "dharma_results")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

model_files <- list.files(
  model_dir,
  pattern = "^models_NC_[0-9]+\\.1_T1\\.rds$",
  full.names = TRUE
)

if (length(model_files) == 0) {
  stop("No chromosome model files found in: ", model_dir)
}

for (model_file in model_files) {
  
  chr_accession <- sub(
    "^models_(NC_[0-9]+\\.1)_T1\\.rds$",
    "\\1",
    basename(model_file)
  )
  
  cat("\nProcessing", chr_accession, "\n")
  
  models <- readRDS(model_file)
  results <- list()
  
  for (model_id in names(models)) {
    
    cat("Processing", model_id, "\n")
    
    model <- models[[model_id]]
    
    if (is.null(model$sdr$pdHess) || !model$sdr$pdHess) {
      results[[length(results) + 1]] <- data.table(
        model_id = model_id,
        status = "non_positive_definite_hessian"
      )
      next
    }
    
    tryCatch({
      
      d <- model$frame
      
      sim_res <- DHARMa::simulateResiduals(
        fittedModel = model,
        n = 250
      )
      
      test_uniformity <- DHARMa::testUniformity(sim_res)
      test_dispersion <- DHARMa::testDispersion(sim_res)
      test_outliers <- DHARMa::testOutliers(sim_res)
      
      result <- data.table(
        model_id = model_id,
        status = "success",
        uniformity_statistic = as.numeric(test_uniformity$statistic),
        uniformity_p = test_uniformity$p.value,
        dispersion_statistic = as.numeric(test_dispersion$statistic),
        dispersion_p = test_dispersion$p.value,
        outlier_statistic = as.numeric(test_outliers$statistic),
        outlier_p = test_outliers$p.value
      )
      
      pdf(
        file.path(plot_dir, paste0(model_id, "_DHARMa.pdf"))
      )
      plot(sim_res)
      dev.off()
      
      results[[length(results) + 1]] <- result
      
    }, error = function(e) {
      
      if (dev.cur() != 1) dev.off()
      
      results[[length(results) + 1]] <<- data.table(
        model_id = model_id,
        status = "error",
        error = conditionMessage(e)
      )
    })
  }
  
  results_chr <- rbindlist(results, fill = TRUE)
  
  fwrite(
    results_chr,
    file.path(
      plot_dir,
      paste0("DHARMa_diagnostics_", chr_accession, "_T1.csv")
    )
  )
  
  cat(
    "Finished", chr_accession,
    "| Models:", length(models),
    "| Results:", nrow(results_chr),
    "\n"
  )
}

cat("\nAll chromosomes processed.\n")