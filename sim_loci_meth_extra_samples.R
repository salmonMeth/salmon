library(data.table)
library(methylKit)
library(VGAM)

simulate_parental_methylation <- function(
    db,
    simulated_coverage,
    metadata,
    timepoint = "T1",
    comb_probs = c(
      none = 0.25,
      maternal_dominant = 0.25,
      paternal_dominant = 0.25,
      both = 0.25
    ),
    rho = 0.02,
    beta_weight_mean = 0,
    beta_weight_sd = 0.3,
    seed = 123
) {
  
  set.seed(seed)
  
  meta <- as.data.table(copy(metadata))
  
  required_meta <- c(
    "sample", "TimePoint", "ID_ma", "ID_pa", "Weight.g"
  )
  
  missing_cols <- setdiff(required_meta, names(meta))
  
  
  comb_probs <- comb_probs[
    c("none", "maternal_dominant",
      "paternal_dominant", "both")
  ]
  
  comb_probs <- comb_probs / sum(comb_probs)
  
  
  coverage <- as.matrix(simulated_coverage$coverage)
  loci <- as.data.table(copy(simulated_coverage$loci))
  sample_ids <- as.character(simulated_coverage$sample_ids)
  
  
  coverage <- coverage[, sample_ids, drop = FALSE]
  
  # Match every simulated sample to the combined metadata
  meta_idx <- match(sample_ids, meta$sample)
  
  if (anyNA(meta_idx)) {
    stop(
      "Missing metadata for coverage sample IDs: ",
      paste(sample_ids[is.na(meta_idx)], collapse = ", ")
    )
  }
  
  meta <- meta[meta_idx]
  
  
  # Scale weights across all samples included in this simulation
  meta[, Weight_sc := as.numeric(scale(Weight.g))]
  
  # If all weights are identical, scale() produces NaN
  meta[!is.finite(Weight_sc), Weight_sc := 0]
  
  weight <- meta$Weight_sc
  
  mother_id <- as.character(meta$ID_ma)
  father_id <- as.character(meta$ID_pa)
  
  missing_mother <- is.na(mother_id) | mother_id == ""
  missing_father <- is.na(father_id) | father_id == ""
  
  mother_id[missing_mother] <- paste0(
    "unknown_mother_", meta$sample[missing_mother]
  )
  
  father_id[missing_father] <- paste0(
    "unknown_father_", meta$sample[missing_father]
  )
  
  #get actual methyl counts
  
  original <- as.data.table(methylKit::getData(db))
  
  original_samples <- db@sample.ids
  cs_cols_all <- names(original)[db@numCs.index]
  ts_cols_all <- names(original)[db@numTs.index]
  
  if (length(cs_cols_all) != length(original_samples) ||
      length(ts_cols_all) != length(original_samples)) {
    stop("Original methylation count columns do not match sample IDs.")
  }
  
  observed_idx <- which(
    grepl(paste0("_", timepoint, "$"), original_samples)
  )
  
  if (length(observed_idx) == 0L) {
    stop("No original database samples found for ", timepoint)
  }
  
  cs_cols <- cs_cols_all[observed_idx]
  ts_cols <- ts_cols_all[observed_idx]
  
  # Match selected loci to the original database by coordinates
  key <- function(x) {
    paste(
      as.character(x$chr),
      x$start,
      x$end,
      sep = "\r"
    )
  }
  
  original_rows <- match(key(loci), key(original))
  
  if (anyNA(original_rows)) {
    stop(
      sum(is.na(original_rows)),
      " simulated loci could not be matched to original database coordinates."
    )
  }
  
  original_cs <- as.matrix(
    original[original_rows, ..cs_cols]
  )
  
  original_ts <- as.matrix(
    original[original_rows, ..ts_cols]
  )
  
  storage.mode(original_cs) <- "numeric"
  storage.mode(original_ts) <- "numeric"
  
  #mean methyl
  
  n_loci <- nrow(loci)
  n_samples <- length(sample_ids)
  
  empirical_rate <- rep(NA_real_, n_loci)
  
  for (i in seq_len(n_loci)) {
    
    cs <- original_cs[i, ]
    ts <- original_ts[i, ]
    
    valid <- is.finite(cs) & is.finite(ts) &
      cs >= 0 & ts >= 0
    
    if (any(valid)) {
      
      total_cs <- sum(cs[valid])
      total_reads <- sum(cs[valid] + ts[valid])
      
      if (total_reads > 0) {
        empirical_rate[i] <- total_cs / total_reads
      }
    }
  }
  #get random effects
  
  combinations <- sample(
    names(comb_probs),
    size = n_loci,
    replace = TRUE,
    prob = comb_probs
  )
  
  parameters <- data.table(
    chr = loci$chr,
    start = loci$start,
    end = loci$end,
    empirical_rate = empirical_rate,
    combs = combinations
  )
  
  parameters[, true_beta_weight := rnorm(
    .N,
    mean = beta_weight_mean,
    sd = beta_weight_sd
  )]
  
  parameters[, `:=`(
    true_sd_ma = 0,
    true_sd_pa = 0
  )]
  
  parameters[combs == "none", `:=`(
    true_sd_ma = 1e-9,
    true_sd_pa = 1e-9
  )]
  
  parameters[combs == "maternal_dominant", `:=`(
    true_sd_ma = runif(.N, 0.12, 0.30),
    true_sd_pa = runif(.N, 0.00, 0.05)
  )]
  
  parameters[combs == "paternal_dominant", `:=`(
    true_sd_ma = runif(.N, 0.00, 0.05),
    true_sd_pa = runif(.N, 0.12, 0.30)
  )]
  
  parameters[combs == "both", `:=`(
    true_sd_ma = runif(.N, 0.10, 0.30),
    true_sd_pa = runif(.N, 0.10, 0.30)
  )]
  
  parameters[, `:=`(
    true_var_ma = true_sd_ma^2,
    true_var_pa = true_sd_pa^2
  )]

  
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
  
  maternal_effects <- vector("list", n_loci)
  paternal_effects <- vector("list", n_loci)
  
  unique_mothers <- unique(mother_id)
  unique_fathers <- unique(father_id)
  
  for (i in seq_len(n_loci)) {
    
    baseline <- empirical_rate[i]
    
    # No empirical baseline: leave probabilities and Cs as NA
    if (!is.finite(baseline)) {
      next
    }
    
    baseline <- pmin(pmax(baseline, 1e-6), 1 - 1e-6)
    
    sd_ma <- parameters$true_sd_ma[i]
    sd_pa <- parameters$true_sd_pa[i]
    
    # Each parent has one shared effect at this locus
    u_mother <- rnorm(
      length(unique_mothers),
      mean = 0,
      sd = sd_ma
    )
    names(u_mother) <- unique_mothers
    
    u_father <- rnorm(
      length(unique_fathers),
      mean = 0,
      sd = sd_pa
    )
    names(u_father) <- unique_fathers
    
    beta_i <- parameters$true_beta_weight[i]
    
    eta <- qlogis(baseline) +
      unname(u_mother[mother_id]) +
      unname(u_father[father_id]) +
      beta_i * weight
    
    p <- plogis(eta)
    
    simulated_probability[i, ] <- p
    
    cov_i <- coverage[i, ]
    
    # Preserve NA coverage as NA Cs
    missing_cov <- is.na(cov_i)
    
    # Preserve zero coverage as zero Cs
    zero_cov <- !missing_cov & cov_i == 0
    simulated_Cs[i, zero_cov] <- 0L
    
    # Simulate counts only for positive, finite coverage
    positive <- is.finite(cov_i) & cov_i > 0 & is.finite(p)
    
    if (any(positive)) {
      
      size_i <- as.integer(round(cov_i[positive]))
      
      if (rho == 0) {
        counts <- rbinom(
          n = sum(positive),
          size = size_i,
          prob = p[positive]
        )
      } else {
        counts <- VGAM::rbetabinom(
          n = sum(positive),
          size = size_i,
          prob = p[positive],
          rho = rho
        )
      }
      
      simulated_Cs[i, positive] <- as.integer(counts)
    }
    
    maternal_effects[[i]] <- data.table(
      parent = unique_mothers,
      effect = unname(u_mother),
      type = "maternal",
      locus_row = i
    )
    
    paternal_effects[[i]] <- data.table(
      parent = unique_fathers,
      effect = unname(u_father),
      type = "paternal",
      locus_row = i
    )
  }

  
  parameters[, `:=`(
    n_samples = n_samples,
    true_mean_probability = NA_real_
  )]
  
  for (i in seq_len(n_loci)) {
    if (any(is.finite(simulated_probability[i, ]))) {
      parameters$true_mean_probability[i] <-
        mean(simulated_probability[i, ], na.rm = TRUE)
    }
  }
  
  # check
  
  stopifnot(
    all(is.na(simulated_Cs[is.na(coverage)])),
    all(simulated_Cs[!is.na(simulated_Cs)] >= 0),
    all(
      simulated_Cs[!is.na(simulated_Cs)] ==
        floor(simulated_Cs[!is.na(simulated_Cs)])
    )
  )
  
  observed_pairs <- !is.na(simulated_Cs)
  
  if (any(
    simulated_Cs[observed_pairs] >
    coverage[observed_pairs]
  )) {
    stop("Some simulated Cs counts exceed coverage.")
  }
  
  
  list(
    loci = loci,
    sample_ids = sample_ids,
    metadata = meta,
    coverage = coverage,
    Cs = simulated_Cs,
    methylation_probability = simulated_probability,
    parameters = parameters,
    beta_weight_mean = beta_weight_mean,
    beta_weight_sd = beta_weight_sd,
    maternal_effects = rbindlist(maternal_effects, fill = TRUE),
    paternal_effects = rbindlist(paternal_effects, fill = TRUE),
    rho = rho,
    seed = seed
  )
}

#simulated metadata file
path_ids <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/with_20_new_samples/simulated_with_20_extra_samples_metadata.csv"

meta <- fread(path_ids)

#actual db file 
db <- readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds")

#simulated coverage file
ch <- readRDS("/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/with_20_new_samples/simulated_with_20_extra_samples_coverage_1000per_chr.rds")

sim_methylation <- simulate_parental_methylation(
  db = db,
  simulated_coverage = ch,
  metadata = meta,
  timepoint = "T1",
  beta_weight_mean = 0,
  beta_weight_sd = 0.3,
  rho = 0.02,
  seed = 123
)

saveRDS(
  sim_methylation,"/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/with_20_new_samples/simulated_with_20_extra_samples_methyl_1000per_chr.rds"
)
