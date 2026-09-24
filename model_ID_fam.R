library ( glmmTMB )
library(data.table)


#this model separates T1 and T2 samples.

path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]


base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_cs_cov_3_samples_removed"
out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_beta_binom_ID_fam"
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
  #this is for the scaled weight variable
  fit_region_sc <- function(i){
    
    d <- dat[
      chr == regions$chr[i] &
        start == regions$start[i] &
        end == regions$end[i]
    ]
    
    d_T1 <- d[TimePoint == "T1"]
    d_T2 <- d[TimePoint == "T2"]
    
    # Region name
    region_name <- paste0(
      regions$chr[i], "_",
      regions$start[i], "_",
      regions$end[i]
    )
    
    # Objects to store
    model_T1 <- NULL
    model_T2 <- NULL
    
    ####################################################
    # T1
    
    status_T1 <- "NOT_RUN"
    
    if(
      nrow(d_T1) >= 10 &&
      length(unique(d_T1$Weight.g)) > 1
    ){
      
      tryCatch({
        
        d_T1$Weight_sc <- as.numeric(
          scale(d_T1$Weight.g)
        )
        
        model_T1 <- glmmTMB(
          cbind(
            meth_reads,
            unmeth_reads
          ) ~
            Weight_sc +
            (1 | ID_family),
          data = d_T1,
          family = betabinomial(link = "logit")
        )
        
        if(
          !is.null(model_T1$sdr) &&
          isTRUE(model_T1$sdr$pdHess)
        ){
          status_T1 <- "OK"
        } else {
          status_T1 <- "NON_POSITIVE_HESSIAN"
        }
        
      }, error = function(e){
        
        status_T1 <<- paste0(
          "ERROR: ",
          conditionMessage(e)
        )
        
      })
      
    } else {
      
      status_T1 <- "INSUFFICIENT_DATA"
      
    }
    
    ####################################################
    # T2
    
    status_T2 <- "NOT_RUN"
    
    if(
      nrow(d_T2) >= 10 &&
      length(unique(d_T2$Weight.g)) > 1
    ){
      
      tryCatch({
        
        d_T2$Weight_sc <- as.numeric(
          scale(d_T2$Weight.g)
        )
        
        model_T2 <- glmmTMB(
          cbind(
            meth_reads,
            unmeth_reads
          ) ~
            Weight_sc +
            (1 | ID_family),
          data = d_T2,
          family = betabinomial(link = "logit")
        )
        
        if(
          !is.null(model_T2$sdr) &&
          isTRUE(model_T2$sdr$pdHess)
        ){
          status_T2 <- "OK"
        } else {
          status_T2 <- "NON_POSITIVE_HESSIAN"
        }
        
      }, error = function(e){
        
        status_T2 <<- paste0(
          "ERROR: ",
          conditionMessage(e)
        )
        
      })
      
    } else {
      
      status_T2 <- "INSUFFICIENT_DATA"
      
    }
    
    
    ####################################################
    # Return nested structure
    ####################################################
    
    return(
      list(
        location = list(
          chr = regions$chr[i],
          start = regions$start[i],
          end = regions$end[i]
        ),
        
        model_T1 = model_T1,
        model_T2 = model_T2,
        
        status_T1 = status_T1,
        status_T2 = status_T2
      )
    )
  }
  
  
  results <- lapply(
    seq_len(nrow(regions)),
    fit_region_sc
  )
  
  # Give each region its genomic name
  names(results) <- paste0(
    regions$chr,
    "_",
    regions$start,
    "_",
    regions$end
  )
  
  out_file <- file.path(
    out_dir,
    paste0("model_", chr_name, ".rds")
  )
  
  saveRDS(
    results,
    out_file
  )
  
  cat(
    "Saved",
    length(results),
    "regions to",
    out_file,
    "\n\n"
  )
}

cat("Done.\n")

#alt overdispersion parameter function from https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


##################################################
##################################################
# assess the model
library(data.table)
library(glmmTMB)
library(DHARMa)


model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated_beta_binom"

out_dir <- file.path(model_dir, "model_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(
  out_dir,
  "model_fit_diagnostics.csv"
)

rds_files <- list.files(
  model_dir,
  pattern = "^model_.*\\.rds$",
  full.names = TRUE
)

cat("Found", length(rds_files), "RDS files\n")

extract_model_stats <- function(model, model_file, region, timepoint) {
  
  out <- data.table(
    model_file = basename(model_file),
    region = region,
    timepoint = timepoint,
    
    n_obs = NA_integer_,
    n_fixed_parameters = NA_integer_,
    
    AIC = NA_real_,
    BIC = NA_real_,
    logLik = NA_real_,
    
    convergence = NA_integer_,
    pdHess = NA,
    
    fixed_effect_estimate = NA_real_,
    fixed_effect_SE = NA_real_,
    fixed_effect_z = NA_real_,
    fixed_effect_pvalue = NA_real_,
    
    status = "failed",
    error = NA_character_
  )
  
  # -------------------------------------------------------
  # Check model
  # -------------------------------------------------------
  
  if (is.null(model)) {
    out$error <- "NULL model"
    return(out)
  }
  
  if (!inherits(model, "glmmTMB")) {
    out$error <- paste(
      "Unexpected class:",
      paste(class(model), collapse = ", ")
    )
    return(out)
  }
  
  # -------------------------------------------------------
  # Extract information
  # -------------------------------------------------------
  
  tryCatch({
    
    out$n_obs <- nobs(model)
    
    out$n_fixed_parameters <- length(
      fixef(model)$cond
    )
    
    out$AIC <- AIC(model)
    
    out$BIC <- BIC(model)
    
    out$logLik <- as.numeric(
      logLik(model)
    )
    
    out$convergence <- model$fit$convergence
    
    out$pdHess <- if (
      !is.null(model$sdr)
    ) {
      isTRUE(model$sdr$pdHess)
    } else {
      NA
    }
    
    # -----------------------------------------------------
    # Fixed effect: Weight_sc
    # -----------------------------------------------------
    
    coefs <- summary(model)$coefficients$cond
    
    if ("Weight_sc" %in% rownames(coefs)) {
      
      out$fixed_effect_estimate <- coefs[
        "Weight_sc",
        "Estimate"
      ]
      
      out$fixed_effect_SE <- coefs[
        "Weight_sc",
        "Std. Error"
      ]
      
      out$fixed_effect_z <- coefs[
        "Weight_sc",
        "z value"
      ]
      
      out$fixed_effect_pvalue <- coefs[
        "Weight_sc",
        "Pr(>|z|)"
      ]
    }
    
    out$status <- "OK"
    
  }, error = function(e) {
    
    out$error <<- conditionMessage(e)
    
  })
  
  return(out)
}



#####loop
all_results <- list()

counter <- 1

for (f in rds_files) {
  
  cat(
    "\nReading:",
    basename(f),
    "\n"
  )
  
  x <- readRDS(f)
  
  for (region in names(x)) {
    
    # T1
    all_results[[counter]] <- extract_model_stats(
      model = x[[region]]$model_T1,
      model_file = f,
      region = region,
      timepoint = "T1"
    )
    
    counter <- counter + 1
    
    # T2
    all_results[[counter]] <- extract_model_stats(
      model = x[[region]]$model_T2,
      model_file = f,
      region = region,
      timepoint = "T2"
    )
    
    counter <- counter + 1
  }
}

model_statistics <- rbindlist(
  all_results,
  fill = TRUE
)

fwrite(
  model_statistics,
  out_file
)

cat(
  "\nSaved",
  nrow(model_statistics),
  "model records to:\n",
  out_file,
  "\n"
)