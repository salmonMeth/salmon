library(data.table)
library(glmmTMB)

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/with_extra_20_samples_per_family"

model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/with_extra_20_samples_per_family/model_results"

# Segment-level truth: true effects and segment parameters
truth_file <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/with_extra_20_samples_per_family/true_sample_parameters.csv"

# Expected columns: sample, Weight_sc, ID_ma, ID_pa
sample_params_file <- truth_file

dir.create(
  model_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

stopifnot(
  file.exists(truth_file),
  file.exists(sample_params_file),
  dir.exists(out_dir)
)


truth <- fread(truth_file)

required_truth_cols <- c("segment","chr","start",
                         "end","mean_methylation","weight_beta","true_sd_ma",
                         "true_sd_pa",
                         "true_var_ma",
                         "true_var_pa"
)

missing_truth_cols <- setdiff(
  required_truth_cols,
  names(truth)
)

if (length(missing_truth_cols) > 0L) {
  stop(
    "Missing columns in segment truth file: ",
    paste(missing_truth_cols, collapse = ", ")
  )
}

sim_parameters <- unique(
  truth[
    ,
    .(
      chr = as.character(chr),
      start = as.integer(start),
      end = as.integer(end),
      segment,
      mean_methylation,
      baseline,
      weight_beta = weight_beta,
      true_sd_ma,
      true_sd_pa,
      true_var_ma,
      true_var_pa
    )
  ]
)

#get the sample params

sample_params <- fread(sample_params_file)

required_sample_cols <- c(
  "sample",
  "Weight_sc",
  "ID_ma",
  "ID_pa"
)

missing_sample_cols <- setdiff(
  required_sample_cols,
  names(sample_params)
)

if (length(missing_sample_cols) > 0L) {
  stop(
    "Missing columns in sample truth file: ",
    paste(missing_sample_cols, collapse = ", ")
  )
}

sample_params <- unique(
  sample_params[
    ,
    .(
      sample = as.character(sample),
      Weight_sc = as.numeric(Weight_sc),
      ID_ma,
      ID_pa
    )
  ]
)


#get cs/cov

cs_files <- list.files(
  out_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)

if (length(cs_files) == 0L) {
  stop("No CS files found in: ", out_dir)
}

sim_list <- list()

for (cs_file in cs_files) {
  
  chr_name <- sub(
    "_segment_cs\\.csv$",
    "",
    basename(cs_file)
  )
  
  cov_file <- file.path(
    out_dir,
    paste0(chr_name, "_segment_coverage.csv")
  )
  
  if (!file.exists(cov_file)) {
    warning("Coverage file missing for ", chr_name)
    next
  }
  
  cs <- fread(cs_file)
  coverage <- fread(cov_file)
  
  key_cols <- c("chr", "start", "end")
  
  if (!all(key_cols %in% names(cs)) ||
      !all(key_cols %in% names(coverage))) {
    stop("Missing genomic coordinate columns for ", chr_name)
  }
  
  cs[, chr := as.character(chr)]
  coverage[, chr := as.character(chr)]
  
  # 
  t1_samples <- grep(
    "_T1$",
    names(cs),
    value = TRUE
  )
  
  t1_samples <- setdiff(
    t1_samples,
    "176_T1"
  )

  missing_cov_samples <- setdiff(
    t1_samples,
    names(coverage)
  )
  
  cs <- cs[
    ,
    c(key_cols, t1_samples),
    with = FALSE
  ]
  
  coverage <- coverage[
    ,
    c(key_cols, t1_samples),
    with = FALSE
  ]
    cs_long <- melt(
    cs,
    id.vars = key_cols,
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  coverage_long <- melt(
    coverage,
    id.vars = key_cols,
    variable.name = "sample",
    value.name = "coverage"
  )
  
  cs_long[, sample := as.character(sample)]
  coverage_long[, sample := as.character(sample)]
  
  d <- merge(
    cs_long,
    coverage_long,
    by = c(key_cols, "sample"),
    all = FALSE,
    sort = FALSE
  )
  
  # merge with truths
  d <- merge(
    d,
    sim_parameters[
      ,
      .(chr, start, end, segment)
    ],
    by = key_cols,
    all.x = TRUE,
    sort = FALSE
  )
  
  if (anyNA(d$segment)) {
    unmatched <- unique(
      d[
        is.na(segment),
        .(chr, start, end)
      ]
    )
    
    stop(
      nrow(unmatched),
      " genomic intervals in ",
      chr_name,
      " do not match the segment truth file."
    )
  }
  
  d[
    ,
    `:=`(
      meth_reads = as.numeric(meth_reads),
      coverage = as.numeric(coverage)
    )
  ]
  
  
  d[
    ,
    unmeth_reads := coverage - meth_reads
  ]
  
  sim_list[[chr_name]] <- d
}

if (length(sim_list) == 0L) {
  stop("No chromosome data were successfully loaded.")
}

sim <- rbindlist(
  sim_list,
  use.names = TRUE,
  fill = TRUE
)

# add the true params

sim[
  ,
  sample := as.character(sample)
]

missing_samples <- setdiff(
  unique(sim$sample),
  sample_params$sample
)


sim <- merge(
  sim,
  sample_params,
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)



sim <- merge(
  sim,
  sim_parameters[
    ,
    .(
      chr,
      start,
      end,
      mean_methylation,
      baseline,
      weight_beta,
      true_sd_ma,
      true_sd_pa,
      true_var_ma,
      true_var_pa
    )
  ],
  by = c("chr", "start", "end"),
  all.x = TRUE,
  sort = FALSE
)



####################################
####################################
#fit the model

chromosomes <- unique(sim$chr)

all_results <- list()
all_models <- list()

for (cc in chromosomes) {
  
  d_chr <- sim[chr == cc]
  
  segments_chr <- unique(d_chr$segment)
  
  chr_results <- list()
  chr_models <- list()
  
  for (seg in segments_chr) {
    
    d <- d_chr[
      segment == seg
    ]
    
    #check
    d <- d[
      is.finite(meth_reads) &
        is.finite(unmeth_reads) &
        is.finite(Weight_sc) &
        meth_reads >= 0 &
        unmeth_reads >= 0
    ]
    
    if (nrow(d) == 0L) {
      next
    }
    
    # Both random effects must have multiple levels.
    fit <- tryCatch(
      glmmTMB(
        cbind(meth_reads, unmeth_reads) ~
          Weight_sc +
          (1 | ID_ma) +
          (1 | ID_pa),
        data = d,
        family = glmmTMB::betabinomial(
          link = "logit"
        ),
        REML = TRUE
      ),
      error = function(e) {
        message(
          "FAILED: chr=", cc,
          " segment=", seg,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(fit)) {
      next
    }
    
    # Fixed effects
    fixed <- fixef(fit)$cond
    
    intercept_est <- if ("(Intercept)" %in% names(fixed)) {
      unname(fixed["(Intercept)"])
    } else {
      NA_real_
    }
    
    weight_est <- if ("Weight_sc" %in% names(fixed)) {
      unname(fixed["Weight_sc"])
    } else {
      NA_real_
    }
    
    # Estimated random-effect SDs
    vc <- VarCorr(fit)
    
    sd_ma <- NA_real_
    sd_pa <- NA_real_
    
    if ("ID_ma" %in% names(vc$cond)) {
      sd_ma <- unname(
        attr(vc$cond$ID_ma, "stddev")[1]
      )
    }
    
    if ("ID_pa" %in% names(vc$cond)) {
      sd_pa <- unname(
        attr(vc$cond$ID_pa, "stddev")[1]
      )
    }
    
    # One record of segment truth
    tr <- unique(
      d[
        ,
        .(
          chr,
          segment,
          start,
          end,
          mean_methylation,
          weight_beta,
          true_sd_ma,
          true_sd_pa,
          true_var_ma,
          true_var_pa
        )
      ]
    )

    
    convergence <- fit$fit$convergence
    
    pd_hess <- if (!is.null(fit$sdr$pdHess)) {
      fit$sdr$pdHess
    } else {
      NA
    }
    
    status <- if (
      identical(as.integer(convergence), 0L) &&
      isTRUE(pd_hess)
    ) {
      "OK"
    } else {
      "CHECK"
    }
    
    n_obs <- nrow(d)
    n_ma <- uniqueN(d$ID_ma)
    n_pa <- uniqueN(d$ID_pa)
    
    total_reads <- sum(
      d$meth_reads + d$unmeth_reads,
      na.rm = TRUE
    )
    
    chr_results[[as.character(seg)]] <- data.table(
      chr = tr$chr,
      segment = seg,
      start = tr$start,
      end = tr$end,
      mean_methylation = tr$mean_methylation,
      
      true_weight_beta = tr$weight_beta,
      estimated_weight_beta = weight_est,
      
      true_sd_ma = tr$true_sd_ma,
      estimated_sd_ma = sd_ma,
      
      true_sd_pa = tr$true_sd_pa,
      estimated_sd_pa = sd_pa,
      
      true_var_ma = tr$true_var_ma,
      estimated_var_ma = sd_ma^2,
      
      true_var_pa = tr$true_var_pa,
      estimated_var_pa = sd_pa^2,
      
      error_weight = weight_est - tr$weight_beta,
      error_sd_ma = sd_ma - tr$true_sd_ma,
      error_sd_pa = sd_pa - tr$true_sd_pa,
      error_var_ma = sd_ma^2 - tr$true_var_ma,
      error_var_pa = sd_pa^2 - tr$true_var_pa,
      
      intercept_estimate = intercept_est,
      
      n_obs = n_obs,
      n_ma = n_ma,
      n_pa = n_pa,
      total_reads = total_reads,
      
      AIC = AIC(fit),
      BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      
      convergence = convergence,
      pdHess = pd_hess,
      status = status
    )
    
    chr_models[[as.character(seg)]] <- fit
  }
  
  if (length(chr_results) > 0L) {
    
    chr_results_dt <- rbindlist(
      chr_results,
      use.names = TRUE,
      fill = TRUE
    )
    
    setorder(chr_results_dt, start)
    
    all_results[[cc]] <- chr_results_dt
    all_models[[cc]] <- chr_models
    
    saveRDS(
      list(
        chr = cc,
        results = chr_results_dt,
        models = chr_models,
        formula =
          cbind(meth_reads, unmeth_reads) ~
          Weight_sc +
          (1 | ID_ma) +
          (1 | ID_pa),
        family = glmmTMB::betabinomial(
          link = "logit"
        ),
        REML = TRUE
      ),
      file = file.path(
        model_dir,
        paste0("model_results_chr", cc, ".rds")
      )
    )
    
    cat(
      "\nFitted ",
      length(chr_results),
      " / ",
      length(segments_chr),
      " segments on ",
      cc,
      "\n",
      sep = ""
    )
    
  } else {
    cat(
      "\nNO SUCCESSFUL MODELS on ",
      cc,
      "\n",
      sep = ""
    )
  }
}


if (length(all_results) == 0L) {
  
  warning("No models were successfully fitted.")
  
  results <- data.table()
  
} else {
  
  results <- rbindlist(
    all_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  setorder(results, chr, start)
}

cat(
  "\nSuccessful model fits: ",
  nrow(results),
  "\n",
  sep = ""
)

if (nrow(results) > 0L) {
  
  cat(
    "OK models: ",
    sum(results$status == "OK"),
    "\n",
    sep = ""
  )
}

saveRDS(
  results,
  file = file.path(
    model_dir,
    "all_model_results_REML_T.rds"
  )
)

fwrite(
  results,
  file = file.path(
    model_dir,
    "all_model_results_REML_T.csv"
  )
)

cat(
  "\nResults saved to:\n",
  model_dir,
  "\n"
)
