library(data.table)
library(glmmTMB)


out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/cs_cov"

model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/model_results_REML_T"

dir.create(
  model_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

truth_file <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/model_results/all_model_results.csv"

# Metadata
path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"


res <- fread(truth_file)

sim_parameters <- res[
  ,
  .(
    segment,
    chr,
    start,
    end,
    mean_methylation,
    weight_beta = true_weight_beta,
    true_sd_ma,
    true_sd_pa,
    true_var_ma,
    true_var_pa
  )
]

sim_parameters[
  ,
  chr_name := sprintf(
    "NC_%06d.1",
    59441L + as.integer(chr)
  )
]

sim_parameters[
  ,
  baseline := qlogis(mean_methylation)
]



cs_files <- list.files(
  out_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)


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
  
  
  cs <- fread(cs_file)
  coverage <- fread(cov_file)
  
  cs[
    ,
    segment := seq_len(.N)
  ]
  
  coverage[
    ,
    segment := seq_len(.N)
  ]
  
  
  t1_samples <- grep(
    "_T1$",
    names(cs),
    value = TRUE
  )
  
  # Remove 176_T1
  t1_samples <- setdiff(
    t1_samples,
    "176_T1"
  )
  
  cs <- cs[
    ,
    c(
      "chr",
      "start",
      "end",
      "segment",
      t1_samples
    ),
    with = FALSE
  ]
  
  coverage <- coverage[
    ,
    c(
      "chr",
      "start",
      "end",
      "segment",
      t1_samples
    ),
    with = FALSE
  ]

  
  cs_long <- melt(
    cs,
    id.vars = c(
      "chr",
      "start",
      "end",
      "segment"
    ),
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  
  coverage_long <- melt(
    coverage,
    id.vars = c(
      "chr",
      "start",
      "end",
      "segment"
    ),
    variable.name = "sample",
    value.name = "coverage"
  )
  
  
  d <- merge(
    cs_long,
    coverage_long,
    by = c(
      "chr",
      "start",
      "end",
      "segment",
      "sample"
    ),
    all = FALSE
  )

  
  d[
    ,
    unmeth_reads := coverage - meth_reads
  ]
  
  sim_list[[chr_name]] <- d
}

sim <- rbindlist(
  sim_list,
  use.names = TRUE,
  fill = TRUE
)

meta <- fread(path_ids)

meta[
  ,
  sample := paste(
    Running.ID,
    TimePoint,
    sep = "_"
  )
]

# Keep T1 and remove 176_T1
meta_T1 <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1",
  .(
    sample,
    Weight.g,
    ID_ma,
    ID_pa
  )
]


meta_T1[
  ,
  Weight_sc := as.numeric(scale(Weight.g))
]



sim <- merge(
  sim,
  meta_T1[
    ,
    .(
      sample,
      Weight_sc,
      ID_ma,
      ID_pa
    )
  ],
  by = "sample",
  all.x = TRUE
)

sim <- merge(
  sim,
  sim_parameters[
    ,
    .(
      chr_name,
      start,
      end,
      mean_methylation,
      weight_beta,
      true_sd_ma,
      true_sd_pa,
      true_var_ma,
      true_var_pa
    )
  ],
  by.x = c(
    "chr",
    "start",
    "end"
  ),
  by.y = c(
    "chr_name",
    "start",
    "end"
  ),
  all.x = TRUE
)


chromosomes <- unique(sim$chr)

print(chromosomes)


all_results <- list()
all_models <- list()

for (cc in chromosomes) {
  
  d_chr <- sim[
    chr == cc
  ]
  
  segments <- unique(
    d_chr$segment
  )
  
  chr_results <- list()
  chr_models <- list()

  
  for (seg in segments) {
    
    d <- d_chr[
      segment == seg
    ]

    
    fit <- tryCatch(
      
      glmmTMB(
        cbind(
          meth_reads,
          unmeth_reads
        ) ~
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
          "FAILED: chr=",
          cc,
          " segment=",
          seg,
          " : ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    if (is.null(fit)) {
      next
    }
    
    # get fixed effect estimates
    fixed <- fixef(fit)$cond
    
    intercept_est <- unname(
      fixed["(Intercept)"]
    )
    
    weight_est <- unname(
      fixed["Weight_sc"]
    )
    #get random effects estimates
    
    vc <- VarCorr(fit)
    
    sd_ma <- NA_real_
    sd_pa <- NA_real_
    
    if ("ID_ma" %in% names(vc$cond)) {
      sd_ma <- attr(
        vc$cond$ID_ma,
        "stddev"
      )[1]
    }
    
    if ("ID_pa" %in% names(vc$cond)) {
      sd_pa <- attr(
        vc$cond$ID_pa,
        "stddev"
      )[1]
    }
    
    # add true params
    
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

    
    #diagnostics
    
    convergence <- fit$fit$convergence
    
    pd_hess <- if (
      !is.null(fit$sdr$pdHess)
    ) {
      fit$sdr$pdHess
    } else {
      NA
    }
    
    n_obs <- nrow(d)
    
    n_ma <- uniqueN(
      d$ID_ma
    )
    
    n_pa <- uniqueN(
      d$ID_pa
    )
    
    total_reads <- sum(
      d$meth_reads +
        d$unmeth_reads,
      na.rm = TRUE
    )
#stats
    
    status <- if (
      convergence == 0 &&
      isTRUE(pd_hess)
    ) {
      "OK"
    } else {
      "CHECK"
    }
    chr_results[[as.character(seg)]] <- data.table(
      
      chr = tr$chr,
      
      segment = seg,
      
      start = tr$start,
      
      end = tr$end,
      
      mean_methylation =
        tr$mean_methylation,
      
      true_weight_beta =
        tr$weight_beta,
      
      estimated_weight_beta =
        weight_est,
      
      true_sd_ma =
        tr$true_sd_ma,
      
      estimated_sd_ma =
        sd_ma,
      
      true_sd_pa =
        tr$true_sd_pa,
      
      estimated_sd_pa =
        sd_pa,
      
      true_var_ma =
        tr$true_var_ma,
      
      true_var_pa =
        tr$true_var_pa,
      
      estimated_var_ma =
        sd_ma^2,
      
      estimated_var_pa =
        sd_pa^2,
      
      error_weight =
        weight_est -
        tr$weight_beta,
      
      error_sd_ma =
        sd_ma -
        tr$true_sd_ma,
      
      error_sd_pa =
        sd_pa -
        tr$true_sd_pa,
      
      error_var_ma =
        sd_ma^2 -
        tr$true_var_ma,
      
      error_var_pa =
        sd_pa^2 -
        tr$true_var_pa,
      
      intercept_estimate =
        intercept_est,
      
      n_obs =
        n_obs,
      
      n_ma =
        n_ma,
      
      n_pa =
        n_pa,
      
      total_reads =
        total_reads,
      
      AIC =
        AIC(fit),
      
      BIC =
        BIC(fit),
      
      logLik =
        as.numeric(
          logLik(fit)
        ),
      
      convergence =
        convergence,
      
      pdHess =
        pd_hess,
      
      status =
        status
    )
    chr_models[[as.character(seg)]] <- fit
    
  } # end segment loop

  
  if (length(chr_results) > 0L) {
    
    chr_results_dt <- rbindlist(
      chr_results,
      fill = TRUE
    )
    
    setorder(
      chr_results_dt,
      start
    )
    all_results[[cc]] <- chr_results_dt
    
    
    all_models[[cc]] <- chr_models
    
    saveRDS(
      list(
        chr = cc,
        
        results =
          chr_results_dt,
        
        models =
          chr_models,
        
        formula =
          cbind(
            meth_reads,
            unmeth_reads
          ) ~
          Weight_sc +
          (1 | ID_ma) +
          (1 | ID_pa),
        
        family =
          glmmTMB::betabinomial(
            link = "logit"
          ),
        
        REML = TRUE
      ),
      file = file.path(
        model_dir,
        paste0(
          "model_results_chr",
          cc,
          ".rds"
        )
      )
    )
    
    
    cat(
      "\nFitted ",
      length(chr_results),
      " / ",
      length(segments),
      " segments on ",
      cc,
      "\n"
    )
    
  } else {
    
    cat(
      "\nNO SUCCESSFUL MODELS on ",
      cc,
      "\n"
    )
  }
  
} # end chromosome loop


if (length(all_results) == 0L) {
  
  warning(
    "No models were successfully fitted."
  )
  
  results <- data.table()
  
} else {
  
  results <- rbindlist(
    all_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  setorder(
    results,
    chr,
    start
  )
}



cat(
  "Successful model fits:",
  nrow(results),
  "\n"
)

if (nrow(results) > 0L) {
  
  cat(
    "Chromosomes represented:",
    uniqueN(results$chr),
    "\n"
  )
  
  cat(
    "OK models:",
    sum(results$status == "OK"),
    "\n"
  )
  
  cat(
    "CHECK models:",
    sum(results$status == "CHECK"),
    "\n"
  )
  
  cat(
    "Mean maternal SD error:",
    mean(results$error_sd_ma, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Mean paternal SD error:",
    mean(results$error_sd_pa, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Mean weight-beta error:",
    mean(results$error_weight, na.rm = TRUE),
    "\n"
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

cat("\n")
cat(
  "Results saved to:\n",
  model_dir,
  "\n"
)