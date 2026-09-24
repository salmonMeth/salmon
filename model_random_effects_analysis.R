library(data.table)
library(glmmTMB)

model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_beta_binom_ID_fam"

out_dir <- file.path(model_dir, "model_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(
  model_dir,
  pattern = "^model_.*\\.rds$",
  full.names = TRUE
)

varcorr_results <- list()
counter <- 1

for (f in rds_files) {
  
  cat("Reading:", basename(f), "\n")
  
  x <- readRDS(f)
  
  for (region in names(x)) {
    
    for (tp in c("T1", "T2")) {
      
      model <- x[[region]][[paste0("model_", tp)]]
      
      if (is.null(model) || !inherits(model, "glmmTMB")) {
        next
      }
      
      vc <- VarCorr(model)$cond
      
      if (length(vc) == 0) {
        next
      }
      
      for (re_name in names(vc)) {
        
        mat <- vc[[re_name]]
        
        varcorr_results[[counter]] <- data.table(
          model_file = basename(f),
          region = region,
          timepoint = tp,
          
          random_effect_group = re_name,
          
          variance = as.numeric(mat[1, 1]),
          
          SD = attr(mat, "stddev")[1]
        )
        
        counter <- counter + 1
      }
    }
  }
}

varcorr_results <- rbindlist(
  varcorr_results,
  fill = TRUE
)

fwrite(
  varcorr_results,
  file.path(
    out_dir,
    "random_effect_variance.csv"
  )
)