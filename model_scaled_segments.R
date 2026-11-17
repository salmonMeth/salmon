library ( glmmTMB )
library(data.table)

#we are modeling the methyl&unmethyl CpG count using weight, mom and dad IDs

path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

scaled_dir <- "/scratch/project_2010912/ezel/salmon/data_files/beta_binom_actual_SCALED_DOWN_COV"

cs_files <- list.files(
  scaled_dir,
  pattern = "_segment_cs\\_scaled.csv$",
  full.names = TRUE
)


results <- list()

for (cs_file in list.files(
  scaled_dir,
  pattern = "_segment_cs_scaled\\.csv$",
  full.names = TRUE
)) {
  
  chr_name <- sub(
    "_segment_cs_scaled\\.csv$",
    "",
    basename(cs_file)
  )
  
  cov_file <- file.path(
    scaled_dir,
    paste0(chr_name, "_segment_coverage_scaled.csv")
  )
  
  if (!file.exists(cov_file)) {
    cat("Missing coverage file for", chr_name, "\n")
    next
  }
  
  cat("\nProcessing", chr_name, "\n")
  
  meth <- fread(cs_file)
  cov <- fread(cov_file)
  
  coord_cols <- c("chr", "start", "end")
  
  sample_cols <- setdiff(names(meth), coord_cols)
  
  stopifnot(
    identical(
      sample_cols,
      setdiff(names(cov), coord_cols)
    ),
    identical(
      meth[, ..coord_cols],
      cov[, ..coord_cols]
    )
  )
  
  meth_long <- melt(
    meth,
    id.vars = coord_cols,
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = coord_cols,
    variable.name = "sample",
    value.name = "coverage"
  )
  
  dat <- merge(
    meth_long,
    cov_long,
    by = c(coord_cols, "sample")
  )
  
  dat <- merge(
    dat,
    meta,
    by = "sample",
    all.x = TRUE
  )
  
  dat[, unmeth_reads := coverage - meth_reads]
  
  dat <- dat[
    is.finite(meth_reads) &
      is.finite(coverage) &
      coverage > 0 &
      meth_reads >= 0 &
      meth_reads <= coverage &
      !is.na(Weight.g) &
      !is.na(ID_ma) &
      !is.na(ID_pa)
  ]
  
  regions <- unique(dat[, ..coord_cols])
  
  models_chr <- list()
  chr_results <- vector("list", nrow(regions))
  
  for (i in seq_len(nrow(regions))) {
    
    chr <- regions$chr[i]
    start <- regions$start[i]
    end <- regions$end[i]
    
    chr_i <- regions$chr[i]
    start_i <- regions$start[i]
    end_i <- regions$end[i]
    
    d <- dat[
      chr == chr_i &
        start == start_i &
        end == end_i
    ]
    
    model_id <- paste(chr, start, end, sep = "_")
    
    out <- data.table(
      chr = chr,
      start = start,
      end = end,
      estimate_T1 = NA_real_,
      SE_T1 = NA_real_,
      z_T1 = NA_real_,
      pvalue_T1 = NA_real_,
      status_T1 = "NOT_RUN",
      estimate_T2 = NA_real_,
      SE_T2 = NA_real_,
      z_T2 = NA_real_,
      pvalue_T2 = NA_real_,
      status_T2 = "NOT_RUN"
    )
    
    for (tp in c("T1", "T2")) {
      
      d_tp <- d[TimePoint == tp]
      
      if (
        nrow(d_tp) < 10 ||
        length(unique(d_tp$Weight.g)) < 2
      ) {
        out[[paste0("status_", tp)]] <- "INSUFFICIENT_DATA"
        next
      }
      
      d_tp[, Weight_sc := as.numeric(scale(Weight.g))]
      
      fit <- tryCatch(
        glmmTMB(
          cbind(meth_reads, unmeth_reads) ~
            Weight_sc +
            (1 | ID_ma) +
            (1 | ID_pa),
          data = d_tp,
          family = betabinomial(link = "logit")
        ),
        error = function(e) e
      )
      
      if (inherits(fit, "error")) {
        out[[paste0("status_", tp)]] <-
          paste0("ERROR: ", conditionMessage(fit))
        next
      }
      
      models_chr[[paste0(model_id, "_", tp)]] <- fit
      
      coefs <- summary(fit)$coefficients$cond
      
      out[[paste0("estimate_", tp)]] <-
        unname(coefs["Weight_sc", "Estimate"])
      
      out[[paste0("SE_", tp)]] <-
        unname(coefs["Weight_sc", "Std. Error"])
      
      out[[paste0("z_", tp)]] <-
        unname(coefs["Weight_sc", "z value"])
      
      out[[paste0("pvalue_", tp)]] <-
        unname(coefs["Weight_sc", "Pr(>|z|)"])
      
      out[[paste0("status_", tp)]] <- if (
        !is.null(fit$sdr) && isTRUE(fit$sdr$pdHess)
      ) "OK" else "NON_POSITIVE_HESSIAN"
    }
    
    chr_results[[i]] <- out
  }
  
  chr_results <- rbindlist(chr_results, fill = TRUE)
  
  saveRDS(
    models_chr,
    file.path(scaled_dir, paste0("models_", chr_name, ".rds"))
  )
  
  fwrite(
    chr_results,
    file.path(scaled_dir, paste0("model_results_", chr_name, ".csv"))
  )
  
  results[[chr_name]] <- chr_results
  
  cat(
    "Saved", length(models_chr), "models and",
    nrow(chr_results), "segment results for", chr_name, "\n"
  )
}

results_dt <- rbindlist(results, fill = TRUE)

fwrite(
  results_dt,
  file.path(scaled_dir, "model_results.csv")
)

saveRDS(
  results_dt,
  file.path(scaled_dir,"model_results.rds")
)