#models the matrices filtered for methylation variability
# 950 segments

library(data.table)
library(glmmTMB)


cov_matrix= readRDS("/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/filtered_segment_cov_T1_and_T2.rds")
cs_matrix=readRDS("/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/filtered_segment_cs_T1_and_T2.rds")

setDT(cs_matrix)
setDT(cov_matrix)

# Sample columns
coord_cols <- c("chr", "start", "end")
sample_cols <- setdiff(names(cs_matrix), coord_cols)


setcolorder(cov_matrix, c(coord_cols, sample_cols))
setcolorder(cs_matrix, c(coord_cols, sample_cols))

#need metadata
path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]
meta <- as.data.table(meta)
setkey(meta, sample)

all_results <- list()
all_random_effects <- list()

# Loop through genomic segments
fit_region <- function(i) {
  
  chr_i <- cs_matrix$chr[i]
  start_i <- cs_matrix$start[i]
  end_i <- cs_matrix$end[i]
  
  # Extract methylated reads and coverage for this segment
  meth_reads <- as.numeric(cs_matrix[i, ..sample_cols])
  coverage <- as.numeric(cov_matrix[i, ..sample_cols])
  
  # Construct sample-level dataset
  d <- data.table(
    sample = sample_cols,
    meth_reads = meth_reads,
    coverage = coverage
  )
  
  # Add metadata
  d <- merge(d, meta, by = "sample", all.x = TRUE)
  
  # Calculate unmethylated reads
  d[, unmeth_reads := coverage - meth_reads]
  
  # Remove invalid observations
  d <- d[
    !is.na(meth_reads) &
      !is.na(coverage) &
      coverage > 0 &
      meth_reads >= 0 &
      meth_reads <= coverage &
      !is.na(ID_ma) &
      !is.na(ID_pa) &
      !is.na(Weight.g)
  ]
  
  # Fit model
  tryCatch({
    
    model <- glmmTMB(
      cbind(meth_reads, unmeth_reads) ~
        Weight.g +
        (1 | ID_ma) +
        (1 | ID_pa),
      data = d,
      family = betabinomial(link = "logit")
    )
    
    # Fixed effects
    coefs <- summary(model)$coefficients$cond
    
    # Random-effect variance components
    vc <- VarCorr(model)$cond
    
    var_ma <- as.numeric(vc$ID_ma[1, 1])
    var_pa <- as.numeric(vc$ID_pa[1, 1])
    
    sd_ma <- sqrt(var_ma)
    sd_pa <- sqrt(var_pa)
    
    # Model diagnostics
    convergence <- model$fit$convergence
    pdHess <- model$sdr$pdHess
    
    # Fixed effect estimate
    if ("Weight.g" %in% rownames(coefs)) {
      
      estimate <- coefs["Weight.g", "Estimate"]
      SE <- coefs["Weight.g", "Std. Error"]
      z <- coefs["Weight.g", "z value"]
      pvalue <- coefs["Weight.g", "Pr(>|z|)"]
      
    } else {
      
      estimate <- SE <- z <- pvalue <- NA_real_
    }
    
    # Summary results
    result <- data.table(
      chr = chr_i,
      start = start_i,
      end = end_i,
      
      n_observations = nrow(d),
      n_mothers = uniqueN(d$ID_ma),
      n_fathers = uniqueN(d$ID_pa),
      
      beta_weight = estimate,
      beta_weight_se = SE,
      z_value = z,
      pvalue = pvalue,
      
      var_maternal = var_ma,
      sd_maternal = sd_ma,
      
      var_paternal = var_pa,
      sd_paternal = sd_pa,
      
      dispersion_phi = sigma(model),
      
      convergence = convergence,
      pdHess = pdHess,
      converged = convergence == 0 && isTRUE(pdHess),
      error = NA_character_
    )
    
    # Extract individual maternal and paternal random-effect estimates
    re <- ranef(model)$cond
    
    re_ma <- as.data.table(re$ID_ma, keep.rownames = "ID_ma")
    re_pa <- as.data.table(re$ID_pa, keep.rownames = "ID_pa")
    
    setnames(re_ma, "(Intercept)", "random_effect")
    setnames(re_pa, "(Intercept)", "random_effect")
    
    re_ma[, parent := "maternal"]
    re_pa[, parent := "paternal"]
    
    re_ma[, `:=`(
      chr = chr_i,
      start = start_i,
      end = end_i
    )]
    
    re_pa[, `:=`(
      chr = chr_i,
      start = start_i,
      end = end_i
    )]
    
    random_effects <- rbindlist(
      list(re_ma, re_pa),
      fill = TRUE
    )
    
    list(
      result = result,
      random_effects = random_effects
    )
    
  }, error = function(e) {
    
    list(
      result = data.table(
        chr = chr_i,
        start = start_i,
        end = end_i,
        n_observations = nrow(d),
        convergence = NA_integer_,
        pdHess = NA,
        converged = FALSE,
        error = conditionMessage(e)
      ),
      random_effects = data.table()
    )
  })
}


# Fit all segments
for (i in seq_len(nrow(cs_matrix))) {
  
  if (i %% 100 == 0) {
    cat("Processing segment", i, "of", nrow(cs_matrix), "\n")
  }
  
  fit <- fit_region(i)
  
  all_results[[i]] <- fit$result
  all_random_effects[[i]] <- fit$random_effects
}


# Combine results
results <- rbindlist(all_results, fill = TRUE)

random_effects <- rbindlist(
  all_random_effects,
  fill = TRUE
)

fwrite(
  results,
  "/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/model_results/all_segments_model_results.csv"
)

fwrite(
  random_effects,
  "/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/model_results/all_segments_random_effect_estimates.csv"
)

cat("\nFinished!\n")
cat("Total segments:", nrow(results), "\n")
cat("Successfully converged:",
    sum(results$converged, na.rm = TRUE), "\n")
cat("Failed models:",
    sum(!results$converged, na.rm = TRUE), "\n")
cat("Random effect estimates:", nrow(random_effects), "\n")

