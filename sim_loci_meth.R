
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
  
  required_meta <- c("sample", "ID_ma", "ID_pa")
  
  
  comb_probs <- comb_probs / sum(comb_probs)
  
  allowed_combs <- c(
    "none",
    "maternal_dominant",
    "paternal_dominant",
    "both"
  )

  
  coverage <- as.matrix(simulated_coverage$coverage)
  loci <- as.data.table(simulated_coverage$loci)
  sample_ids <- simulated_coverage$sample_ids

  
  coverage <- coverage[, sample_ids, drop = FALSE]
  
  meta <- as.data.table(copy(metadata))
  meta <- meta[match(sample_ids, meta$sample)]
  
  meta[, Weight_sc := as.numeric(scale(Weight.g))]
  
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
  
  # get empirical mean
  
  original <- as.data.table(methylKit::getData(db))
  
  cs_cols <- names(original)[db@numCs.index]
  ts_cols <- names(original)[db@numTs.index]
  
  original_samples <- db@sample.ids
  sample_idx <- match(sample_ids, original_samples)

  
  cs_cols <- cs_cols[sample_idx]
  ts_cols <- ts_cols[sample_idx]
  
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
  
  empirical_rate <- rep(NA_real_, nrow(loci))
  
  for (i in seq_len(nrow(loci))) {
    
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

  
  n_loci <- nrow(loci)
  n_samples <- length(sample_ids)
  
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
  #adjust the true params if needed
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
  
#simulate effects
  
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
  
  for (i in seq_len(n_loci)) {
    
   
    # safe log
    baseline <- pmin(
      pmax(empirical_rate[i], 1e-6),
      1 - 1e-6
    )
    
    sd_ma <- parameters$true_sd_ma[i]
    sd_pa <- parameters$true_sd_pa[i]
    
    unique_mothers <- unique(mother_id)
    unique_fathers <- unique(father_id)
    
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

    # shared parents get shared effects
    eta <- qlogis(baseline) +
      unname(u_mother[mother_id]) +
      unname(u_father[father_id]) + beta_i * weight 
    
    p <- plogis(eta)
    
    simulated_probability[i, ] <- p
        cov_i <- coverage[i, ]
    #preserve NA coverages
    valid <- !is.na(cov_i)
    positive <- valid & cov_i > 0
    zero <- valid & cov_i == 0
    
    #preserve 0 coverages
    simulated_Cs[i, zero] <- 0L
    
    if (any(positive)) {
      
      if (rho == 0) {
        counts <- rbinom(
          sum(positive),
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
  
  
  list(
    loci = loci,
    sample_ids = sample_ids,
    coverage = coverage,
    Cs = simulated_Cs,
    methylation_probability = simulated_probability,
    parameters = parameters,
    beta_weight_mean = beta_weight_mean,
    beta_weight_sd = beta_weight_sd,
    maternal_effects = rbindlist(
      maternal_effects,
      fill = TRUE
    ),
    paternal_effects = rbindlist(
      paternal_effects,
      fill = TRUE
    ),
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

sim_methylation <- simulate_parental_methylation(
  db = db,
  simulated_coverage = ch,
  metadata = meta,
  timepoint = "T1",     beta_weight_mean = 0,
  beta_weight_sd = 0.3,
  seed = 123
)

saveRDS(
  sim_methylation,
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/simulated_T1_parental_methylation_1000per_chr.rds"
)