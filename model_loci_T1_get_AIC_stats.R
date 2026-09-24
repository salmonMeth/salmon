#work with the lcoi based models
path_loci_models="/scratch/project_2010912/ezel/mixed_models/meth_MixMod_T1"
files <- list.files(
  path = path_loci_models,
  full.names = TRUE
)

env <- new.env()
for (f in files) {
  load(f, envir = env)
}

##
f <- files[1]
loaded_name <- load(f, envir = env)
loaded_name
x <- env[[loaded_name]]
#structure of x
#x[[1]]  = location
#x[[2]]  = mod1.remlT.list
#x[[3]]  = mod1.remlF.list

#x[[4]]  = location
#x[[5]]  = mod1.remlT.list
#x[[6]]  = mod1.remlF.list
#etc

extract_and_save_AIC <- function(x, original_file, output_dir) {
    # Create output directory if needed
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # Original file name without path/extension
    base_name <- basename(original_file)
    base_name <- sub("\\.rds$", "", base_name)
    base_name <- sub("\\.RData$", "", base_name)
    base_name <- sub("\\.Rdata$", "", base_name)
    
    # Check that the structure is what we expect
    if (length(x) %% 3 != 0) {
      stop(
        "Unexpected structure: length(x) is not divisible by 3."
      )
    }
    
    # Number of genomic regions
    n_regions <- length(x) / 3
    
    # Extract results
    aic_df <- data.frame(
      
      region_index = seq_len(n_regions),
      
      chr = sapply(seq_len(n_regions), function(i) {
        x[[3*i - 2]][1]
      }),
      
      start = sapply(seq_len(n_regions), function(i) {
        as.numeric(x[[3*i - 2]][2])
      }),
      
      end = sapply(seq_len(n_regions), function(i) {
        as.numeric(x[[3*i - 2]][3])
      }),
      
      # REML = TRUE model
      # --------------------------------------------------
      
      AIC_remlT = sapply(seq_len(n_regions), function(i) {
        
        z <- x[[3*i - 1]]
        
        if (is.list(z) && length(z) >= 5) {
          as.numeric(z[[5]])
        } else {
          NA_real_
        }
      }),
      
      convergence_remlT = sapply(seq_len(n_regions), function(i) {
        
        z <- x[[3*i - 1]]
        
        if (is.list(z) && length(z) >= 4) {
          as.numeric(z[[4]])
        } else {
          NA_real_
        }
      }),
      
      message_remlT = sapply(seq_len(n_regions), function(i) {
        
        z <- x[[3*i - 1]]
        
        if (is.list(z) && length(z) >= 3) {
          as.character(z[[3]])
        } else {
          NA_character_
        }
      }),
      

      # REML = FALSE model
     
      
      AIC_remlF = sapply(seq_len(n_regions), function(i) {
        
        z <- x[[3*i]]
        
        if (is.list(z) && length(z) >= 5) {
          as.numeric(z[[5]])
        } else {
          NA_real_
        }
      }),
      
      convergence_remlF = sapply(seq_len(n_regions), function(i) {
        
        z <- x[[3*i]]
        
        if (is.list(z) && length(z) >= 4) {
          as.numeric(z[[4]])
        } else {
          NA_real_
        }
      }),
      
      message_remlF = sapply(seq_len(n_regions), function(i) {
        
        z <- x[[3*i]]
        
        if (is.list(z) && length(z) >= 3) {
          as.character(z[[3]])
        } else {
          NA_character_
        }
      }),
      
      stringsAsFactors = FALSE
    )

    
    aic_df$AIC_remlT_available <-
      !is.na(aic_df$AIC_remlT)
    
    aic_df$AIC_remlF_available <-
      !is.na(aic_df$AIC_remlF)
    
    aic_df$remlT_converged <-
      aic_df$convergence_remlT == 0
    
    aic_df$remlF_converged <-
      aic_df$convergence_remlF == 0

    
    output_rds <- file.path(
      output_dir,
      paste0(base_name, "_AIC.rds")
    )
    
    output_csv <- file.path(
      output_dir,
      paste0(base_name, "_AIC.csv")
    )
    
    
    
    saveRDS(aic_df, output_rds)
    
    write.csv(
      aic_df,
      output_csv,
      row.names = FALSE
    )
    
    message("Saved RDS: ", output_rds)
    message("Saved CSV: ", output_csv)
    
    return(aic_df)
  }

output_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_loci_T1"

env <- new.env() 
all_aic <- list() 
for (f in files) { cat("Processing:", f, "\n") 
  loaded_name <- load(f, envir = env) 
  x <- env[[loaded_name]] 
  all_aic[[loaded_name]] <- extract_and_save_AIC( x = x, original_file = f, output_dir = output_dir ) }