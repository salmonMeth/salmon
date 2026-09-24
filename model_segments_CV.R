#CV for region based models, the basic model 

library ( glmmTMB )
library(data.table)

#we are modeling the methyl&unmethyl CpG count using weight, mom and dad IDs

path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

#this model separates T1 and T2 samples.
base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_cs_cov_3_samples_removed"
out_dir <- ("/scratch/project_2010912/ezel/salmon/data_files/cross_val")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cs_files <- list.files(
  base_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)
for(cs_file in cs_files){
  
  chr_name <- sub(
    "_segment_cs\\.csv$",
    "",
    basename(cs_file)
  )
  
  cov_file <- file.path(
    base_dir,
    paste0(chr_name, "_segment_coverage.csv")
  )
  
  if(!file.exists(cov_file)){
    cat("Missing coverage file for", chr_name, "\n")
    next
  }
  
  cat("Processing", chr_name, "\n")
  
  meth <- fread(cs_file)
  cov <- fread(cov_file)
  
  meth_long <- melt(
    meth,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "coverage"
  )
  
  dat <- merge(
    meth_long,
    cov_long,
    by = c("chr", "start", "end", "sample")
  )
  
  dat <- merge(
    dat,
    meta,
    by = "sample",
    all.x = TRUE
  )
  
  dat[, unmeth_reads := coverage - meth_reads]
  
  regions <- unique(
    dat[, .(chr, start, end)]
  )
  #cv model fitting
  fit_region_sc_cv <- function(i, K = 5, seed = 123) {
    
    d <- dat[
      chr == regions$chr[i] &
        start == regions$start[i] &
        end == regions$end[i]
    ]
    
    d_T1 <- d[TimePoint == "T1"]
    d_T2 <- d[TimePoint == "T2"]
    
    out <- data.table(
      chr = regions$chr[i],
      start = regions$start[i],
      end = regions$end[i],
      
      n_T1 = nrow(d_T1),
      n_T2 = nrow(d_T2),
      
      RMSE_T1 = NA_real_,
      MAE_T1 = NA_real_,
      R2_T1 = NA_real_,
      status_T1 = "NOT_RUN",
      
      RMSE_T2 = NA_real_,
      MAE_T2 = NA_real_,
      R2_T2 = NA_real_,
      status_T2 = "NOT_RUN"
    )
    
    run_cv <- function(d_time) {
      
      if (
        nrow(d_time) < 10 ||
        length(unique(d_time$Weight.g)) <= 1
      ) {
        return(list(
          RMSE = NA_real_,
          MAE = NA_real_,
          R2 = NA_real_,
          status = "INSUFFICIENT_DATA"
        ))
      }
      
      # Need enough parental pairs for K folds
      d_time[, family := interaction(
        ID_ma,
        ID_pa,
        drop = TRUE,
        lex.order = TRUE
      )]
      
      n_families <- uniqueN(d_time$family)
      
      if (n_families < K) {
        return(list(
          RMSE = NA_real_,
          MAE = NA_real_,
          R2 = NA_real_,
          status = "TOO_FEW_FAMILIES"
        ))
      }
      
      #### everyone from a family goes to the same fold
      
      set.seed(seed)
      
      families <- unique(d_time$family)
      
      fold_id <- sample(
        rep(seq_len(K), length.out = length(families))
      )
      
      family_folds <- data.table(
        family = families,
        fold = fold_id
      )
      
      d_time <- merge(
        d_time,
        family_folds,
        by = "family",
        all.x = TRUE
      )
      
      predictions <- vector("list", K)
      
      
      ##cv
      
      for (k in seq_len(K)) {
        
        train <- d_time[fold != k]
        test  <- d_time[fold == k]
        
        # Scale using TRAINING data only
        train_mean <- mean(train$Weight.g, na.rm = TRUE)
        train_sd <- sd(train$Weight.g, na.rm = TRUE)
        
        if (
          is.na(train_sd) ||
          train_sd == 0
        ) {
          next
        }
        
        train[, Weight_sc :=
                (Weight.g - train_mean) / train_sd]
        
        test[, Weight_sc :=
               (Weight.g - train_mean) / train_sd]
        
        
        ##fit the model, get stats
        
        model <- tryCatch({
          
          glmmTMB(
            cbind(
              meth_reads,
              unmeth_reads
            ) ~
              Weight_sc +
              (1 | ID_ma) +
              (1 | ID_pa),
            data = train,
            family = betabinomial(link = "logit")
          )
          
        }, error = function(e) {
          NULL
        })
        
        
        if (is.null(model)) {
          next
        }
        
        pred <- tryCatch({
          
          predict(
            model,
            newdata = test,
            type = "response",
            allow.new.levels = TRUE
          )
          
        }, error = function(e) {
          NULL
        })
        
        
        if (is.null(pred)) {
          next
        }

        
        observed <- test$meth_reads /
          (test$meth_reads + test$unmeth_reads)
        
        predictions[[k]] <- data.table(
          observed = observed,
          predicted = as.numeric(pred)
        )
      }
      
      
      ###
      
      predictions <- rbindlist(
        predictions,
        fill = TRUE
      )
      
      if (nrow(predictions) == 0) {
        return(list(
          RMSE = NA_real_,
          MAE = NA_real_,
          R2 = NA_real_,
          status = "CV_FAILED"
        ))
      }
      
      error <- predictions$predicted -
        predictions$observed
      
      RMSE <- sqrt(
        mean(error^2, na.rm = TRUE)
      )
      
      MAE <- mean(
        abs(error),
        na.rm = TRUE
      )
      
      # R2
      SSE <- sum(
        error^2,
        na.rm = TRUE
      )
      
      SST <- sum(
        (predictions$observed -
           mean(predictions$observed, na.rm = TRUE))^2,
        na.rm = TRUE
      )
      
      R2 <- 1 - SSE / SST
      
      
      return(list(
        RMSE = RMSE,
        MAE = MAE,
        R2 = R2,
        status = "OK"
      ))
    }
    ##T1
    
    tryCatch({
      
      cv_T1 <- run_cv(d_T1)
      
      out$RMSE_T1 <- cv_T1$RMSE
      out$MAE_T1 <- cv_T1$MAE
      out$R2_T1 <- cv_T1$R2
      out$status_T1 <- cv_T1$status
      
    }, error = function(e) {
      
      out$status_T1 <<-
        paste0("ERROR: ", conditionMessage(e))
      
    })
    
    ##T2
    
    tryCatch({
      
      cv_T2 <- run_cv(d_T2)
      
      out$RMSE_T2 <- cv_T2$RMSE
      out$MAE_T2 <- cv_T2$MAE
      out$R2_T2 <- cv_T2$R2
      out$status_T2 <- cv_T2$status
      
    }, error = function(e) {
      
      out$status_T2 <<-
        paste0("ERROR: ", conditionMessage(e))
      
    })
    
    return(out)
  }

  
  cv_results <- rbindlist(
    lapply(
      seq_len(nrow(regions)),
      fit_region_sc_cv
    ),
    fill = TRUE
  )
  
  out_file <- file.path(
    out_dir,
    paste0("model_", chr_name, ".csv")
  )
  
  fwrite(
    results,
    out_file
  )
  
  cat(
    "Saved",
    nrow(results),
    "regions to",
    out_file,
    "\n\n"
  )
}

cat("Done.\n")
