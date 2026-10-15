
library(data.table)
library(methylKit)
library(VGAM)



simulate_family_methylation <- function(
    db,
    simulated_coverage,
    metadata,
    timepoint = "T1",
    family_col = "Fam.ID",
    comb_probs = c(
      none = 0.25,
      family_effect = 0.75
    ),
    rho = 0.02,
    beta_weight_mean = 0,
    beta_weight_sd = 0.3,
    family_sd_min = 0.12,
    family_sd_max = 0.30,
    seed = 123
) {
  
  set.seed(seed)

  
  required_meta <- c("sample", family_col, "Weight.g")
  missing_cols <- setdiff(required_meta, names(metadata))

  allowed_combs <- c("none", "family_effect")
  
  comb_probs <- comb_probs / sum(comb_probs)

  
  coverage <- as.matrix(simulated_coverage$coverage)
  loci <- data.table::as.data.table(
    simulated_coverage$loci
  )
  sample_ids <- simulated_coverage$sample_ids
  
 
  coverage <- coverage[, sample_ids, drop = FALSE]
  
  
  meta <- data.table::as.data.table(
    data.table::copy(metadata)
  )
  
  meta <- meta[match(sample_ids, meta$sample)]
  
  weight_sd <- sd(meta$Weight.g, na.rm = TRUE)
  meta[, Weight_sc := as.numeric(scale(Weight.g))]
  weight <- meta$Weight_sc

  
  family_id <- as.character(meta[[family_col]])
  
  missing_family <- is.na(family_id) | trimws(family_id) == ""
  
  # Missing family IDs become unique per sample
  family_id[missing_family] <- paste0(
    "unknown_family_",
    meta$sample[missing_family]
  )
  
  unique_families <- unique(family_id)
  
  n_families <- length(unique_families)
#get mean methy
  
  original <- data.table::as.data.table(
    methylKit::getData(db)
  )
  
  cs_cols_all <- names(original)[db@numCs.index]
  ts_cols_all <- names(original)[db@numTs.index]
  
  original_samples <- db@sample.ids
  sample_idx <- match(sample_ids, original_samples)
  
  if (anyNA(sample_idx)) {
    stop("Some simulated samples are not present in db@sample.ids.")
  }
  
  cs_cols <- cs_cols_all[sample_idx]
  ts_cols <- ts_cols_all[sample_idx]
  
  key <- function(x) {
    paste(
      as.character(x$chr),
      x$start,
      x$end,
      sep = "\r"
    )
  }
  
  original_rows <- match(key(loci), key(original))
  
  original_cs <- as.matrix(
    original[original_rows, ..cs_cols]
  )
  
  original_ts <- as.matrix(
    original[original_rows, ..ts_cols]
  )
  
  storage.mode(original_cs) <- "numeric"
  storage.mode(original_ts) <- "numeric"
  
  n_loci <- nrow(loci)
  n_samples <- length(sample_ids)
  
  empirical_rate <- rep(NA_real_, n_loci)
  
  for (i in seq_len(n_loci)) {
    
    cs <- original_cs[i, ]
    ts <- original_ts[i, ]
    
    valid <- is.finite(cs) &
      is.finite(ts) &
      cs >= 0 &
      ts >= 0
    
    if (any(valid)) {
      
      total_cs <- sum(cs[valid])
      total_reads <- sum(cs[valid] + ts[valid])
      
      if (total_reads > 0) {
        empirical_rate[i] <- total_cs / total_reads
      }
    }
  }

  
  combinations <- sample(
    names(comb_probs),
    size = n_loci,
    replace = TRUE,
    prob = comb_probs
  )
  
  parameters <- data.table::data.table(
    chr = loci$chr,
    start = loci$start,
    end = loci$end,
    empirical_rate = empirical_rate,
    combs = combinations
  )
  
  # Weight coefficient for each locus
  parameters[, true_beta_weight := rnorm(
    .N,
    mean = beta_weight_mean,
    sd = beta_weight_sd
  )]
  
  # One true family SD per locus
  parameters[, true_sd_family := 0]
  
  parameters[
    combs == "family_effect",
    true_sd_family := runif(
      .N,
      min = family_sd_min,
      max = family_sd_max
    )
  ]
  
  parameters[, true_var_family := true_sd_family^2]

  
  simulated_probability <- matrix(
    NA_real_,
    nrow = n_loci,
    ncol = n_samples,
    dimnames = list(NULL, sample_ids)
  )
  
  simulated_Cs <- matrix(
    NA_integer_,
    nrow = n_loci,
    ncol = n_samples,
    dimnames = list(NULL, sample_ids)
  )
  
  family_effects <- vector("list", n_loci)
  
  # simulate probs
  
  for (i in seq_len(n_loci)) {
    
    # Empirical baseline, bounded to avoid infinite logits
    baseline <- empirical_rate[i]
    
    if (!is.finite(baseline)) {
      baseline <- NA_real_
    } else {
      baseline <- min(max(baseline, 1e-6), 1 - 1e-6)
    }
    
    sd_family <- parameters$true_sd_family[i]
    
    # One random effect for each family at this locus
    u_family <- stats::rnorm(
      n_families,
      mean = 0,
      sd = sd_family
    )
    
    names(u_family) <- unique_families
    
    beta_i <- parameters$true_beta_weight[i]

    if (!is.finite(baseline)) {
      
      p <- rep(NA_real_, n_samples)
      
    } else {
      
      # Samples in the same family share the same effect
      eta <- qlogis(baseline) +
        beta_i * weight +
        unname(u_family[family_id])
      
      p <- plogis(eta)
    }
    
    simulated_probability[i, ] <- p
    
    # Coverage handling: preserve missing and zero coverage
    cov_i <- coverage[i, ]
    
    valid <- !is.na(cov_i)
    positive <- valid & cov_i > 0
    zero <- valid & cov_i == 0
    
    simulated_Cs[i, zero] <- 0L
    
    if (any(positive)) {
      
      if (rho == 0) {
        
        counts <- stats::rbinom(
          n = sum(positive),
          size = as.integer(cov_i[positive]),
          prob = p[positive]
        )
        
      } else {
        
        counts <- VGAM::rbetabinom(
          n = sum(positive),
          size = as.integer(cov_i[positive]),
          prob = p[positive],
          rho = rho
        )
      }
      
      simulated_Cs[i, positive] <- as.integer(counts)
    }
    
    # Save family effects for this locus
    family_effects[[i]] <- data.table::data.table(
      family = unique_families,
      effect = unname(u_family),
      locus_row = i,
      chr = loci$chr[i],
      start = loci$start[i],
      end = loci$end[i],
      true_sd_family = sd_family
    )
  }
  

  
  parameters[, `:=`(
    n_samples = n_samples,
    n_families = n_families,
    true_mean_probability = NA_real_
  )]
  
  for (i in seq_len(n_loci)) {
    
    if (any(is.finite(simulated_probability[i, ]))) {
      
      parameters$true_mean_probability[i] <- mean(
        simulated_probability[i, ],
        na.rm = TRUE
      )
    }
  }
  
  # ------------------------------------------------------------
  # 10. Return simulation results
  # ------------------------------------------------------------
  
  list(
    loci = loci,
    sample_ids = sample_ids,
    coverage = coverage,
    Cs = simulated_Cs,
    methylation_probability = simulated_probability,
    parameters = parameters,
    family_ids = family_id,
    family_col = family_col,
    family_effects = data.table::rbindlist(
      family_effects,
      fill = TRUE
    ),
    beta_weight_mean = beta_weight_mean,
    beta_weight_sd = beta_weight_sd,
    family_sd_min = family_sd_min,
    family_sd_max = family_sd_max,
    rho = rho,
    seed = seed
  )
}



path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]

db=readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds")

ch <- readRDS(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/simulated_T1_locus_coverage_1000per_chr.rds"
)

simulated_coverage <- readRDS(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/simulated_T1_locus_coverage_1000per_chr.rds"
)


sim_methylation <- simulate_family_methylation(
  db,
  simulated_coverage,
  meta,
  timepoint = "T1",
  family_col = "Fam.ID",
  comb_probs = c(
    none = 0.25,
    family_effect = 0.75
  ),
  rho = 0.02,
  beta_weight_mean = 0,
  beta_weight_sd = 0.3,
  family_sd_min = 0.05,
  family_sd_max = 0.50,
  seed = 123
)

saveRDS(
  sim_methylation,
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/simulated_T1_FAM_methylation_1000per_chr.rds"
)
