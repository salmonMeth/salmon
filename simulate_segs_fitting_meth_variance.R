
library(data.table)
library(glmmTMB)

# ============================================================
# 1. PATHS
# ============================================================

meth_file <- paste0(
  "/scratch/project_2010912/ezel/salmon/data_files/",
  "beta_binom_models_for_actual_data_split_segments/",
  "meth_T1_all_1107_segments.rds"
)

unmeth_file <- paste0(
  "/scratch/project_2010912/ezel/salmon/data_files/",
  "beta_binom_models_for_actual_data_split_segments/",
  "unmeth_T1_all_1107_segments.rds"
)

metadata_file <- paste0(
  "/scratch/project_2010912/ezel/",
  "hannu_methylation.selectedIDs_10Dec24.txt"
)

model_dir <- paste0(
  "/scratch/project_2010912/ezel/salmon/data_files/",
  "beta_binom_models_for_actual_data_split_segments/",
  "check_for_dharma/model_results_T1_T2_separated/models_T1"
)

output_dir <- paste0(
  "/scratch/project_2010912/ezel/salmon/data_files/",
  "simulated_data/sim_segs_complex"
)

checkpoint_dir <- file.path(output_dir, "checkpoints")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. SETTINGS
# ============================================================

min_samples <- 42
n_replicates <- 30
checkpoint_every <- 100

sd_values <- seq(0, 0.2, by = 0.02)

sd_grid <- CJ(
  sd_ma = sd_values,
  sd_pa = sd_values,
  sd_env = sd_values
)

set.seed(456)

cat("Number of SD combinations:", nrow(sd_grid), "\n")
cat("Replicates per combination:", n_replicates, "\n")
cat(
  "Simulations per segment:",
  nrow(sd_grid) * n_replicates, "\n"
)


# ============================================================
# 3. LOAD COUNT MATRICES
# ============================================================

cat("\nLoading count matrices...\n")

meth <- readRDS(meth_file)
unmeth <- readRDS(unmeth_file)

stopifnot(
  is.matrix(meth),
  is.matrix(unmeth),
  identical(dim(meth), dim(unmeth)),
  identical(rownames(meth), rownames(unmeth)),
  identical(colnames(meth), colnames(unmeth))
)

# Ensure numeric matrices
storage.mode(meth) <- "numeric"
storage.mode(unmeth) <- "numeric"

cat("Methylated matrix:", dim(meth), "\n")
cat("Unmethylated matrix:", dim(unmeth), "\n")


# ============================================================
# 4. LOAD AND MATCH SAMPLE METADATA
# ============================================================

cat("\nLoading sample metadata...\n")

meta <- fread(metadata_file)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

# Identify T1 columns
t1_cols <- grep(
  "^[0-9]+_T1$",
  colnames(meth),
  value = TRUE
)

if (length(t1_cols) == 0) {
  stop("No T1 sample columns found in the methylated matrix.")
}

# Match metadata to matrix column order
sample_map <- meta[
  match(t1_cols, sample),
  .(sample, ID_ma, ID_pa)
]

if (anyNA(sample_map$sample)) {
  stop("Some T1 samples could not be matched to metadata.")
}

if (anyNA(sample_map$ID_ma) || anyNA(sample_map$ID_pa)) {
  stop("Some T1 samples have missing maternal or paternal IDs.")
}

stopifnot(identical(sample_map$sample, t1_cols))

cat("Matched T1 samples:", length(t1_cols), "\n")

# Restrict matrices to T1 samples
meth_mat <- meth[, t1_cols, drop = FALSE]
unmeth_mat <- unmeth[, t1_cols, drop = FALSE]

# Calculate coverage
cov_mat <- meth_mat + unmeth_mat

cat("Final count matrix dimensions:", dim(meth_mat), "\n")


# ============================================================
# 5. PREPARE SEGMENT TARGETS FROM ROW NAMES
# ============================================================

cat("\nPreparing segment coordinates...\n")

segment_targets <- data.table(
  segment_id = rownames(meth_mat)
)

# Parse IDs of the form:
# NC_059442.1_105452687_105548398
#
# The chromosome itself contains an underscore, so use
# a regular expression rather than splitting on every "_".

parsed <- tstrsplit(
  segment_targets$segment_id,
  "_",
  fixed = TRUE
)

# Parse using a regex that preserves the chromosome accession
segment_targets[, chr := sub(
  "^((NC_[0-9]+\\.[0-9]+))_([0-9]+)_([0-9]+)$",
  "\\1",
  segment_id
)]

segment_targets[, start := as.integer(sub(
  "^NC_[0-9]+\\.[0-9]+_([0-9]+)_([0-9]+)$",
  "\\1",
  segment_id
))]

segment_targets[, end := as.integer(sub(
  "^NC_[0-9]+\\.[0-9]+_([0-9]+)_([0-9]+)$",
  "\\2",
  segment_id
))]

if (anyNA(segment_targets$chr) ||
    anyNA(segment_targets$start) ||
    anyNA(segment_targets$end)) {
  stop("Some segment IDs could not be parsed correctly.")
}

# Confirm all 1,107 segments are represented
stopifnot(nrow(segment_targets) == 1107)

cat("Total segments:", nrow(segment_targets), "\n")
cat("Chromosomes:", uniqueN(segment_targets$chr), "\n")

# All segments are included
segments_to_run <- copy(segment_targets)


# ============================================================
# 6. SIMULATION FUNCTION
# ============================================================

simulate_segment <- function(
    dat,
    sd_ma,
    sd_pa,
    sd_env,
    phi,
    target_mean
) {
  
  n <- nrow(dat)
  
  # Resample observed coverage distribution
  coverage <- sample(
    dat$coverage,
    size = n,
    replace = TRUE
  )
  
  # Generate maternal and paternal effects
  mothers <- unique(dat$ID_ma)
  fathers <- unique(dat$ID_pa)
  
  u_ma <- rnorm(length(mothers), 0, sd_ma)
  u_pa <- rnorm(length(fathers), 0, sd_pa)
  
  maternal_effect <- u_ma[
    match(dat$ID_ma, mothers)
  ]
  
  paternal_effect <- u_pa[
    match(dat$ID_pa, fathers)
  ]
  
  # Individual environmental effects
  env_effect <- rnorm(n, 0, sd_env)
  
  eta <- maternal_effect +
    paternal_effect +
    env_effect
  
  # Calibrate intercept to target mean latent probability
  alpha <- uniroot(
    function(a) {
      mean(plogis(a + eta)) - target_mean
    },
    interval = c(-20, 20)
  )$root
  
  p <- plogis(alpha + eta)
  
  # Beta-binomial sampling
  beta_p <- rbeta(
    n,
    shape1 = p * phi,
    shape2 = (1 - p) * phi
  )
  
  simulated_meth <- rbinom(
    n,
    size = coverage,
    prob = beta_p
  )
  
  simulated_unmeth <- coverage - simulated_meth
  
  simulated_methylation <- simulated_meth / coverage
  
  list(
    mean = mean(simulated_methylation),
    variance = var(simulated_methylation),
    coverage = coverage,
    meth = simulated_meth,
    unmeth = simulated_unmeth
  )
}


# ============================================================
# 7. CHECK INPUTS BEFORE SIMULATION
# ============================================================

stopifnot(
  all(segments_to_run$segment_id %in% rownames(meth_mat)),
  all(segments_to_run$segment_id %in% rownames(unmeth_mat))
)

cat(
  "\nSegments available in both count matrices:",
  nrow(segments_to_run), "\n"
)

# Remove old checkpoints so every segment starts from
# combination 1 with the current matrices and settings.
old_checkpoints <- list.files(
  checkpoint_dir,
  pattern = "_checkpoint\\.rds$",
  full.names = TRUE
)

if (length(old_checkpoints) > 0) {
  file.remove(old_checkpoints)
  cat("Removed old checkpoints:", length(old_checkpoints), "\n")
}


# ============================================================
# 8. SIMULATION LOOP
# ============================================================

summary_file <- file.path(
  output_dir,
  "simulation_summary.csv"
)

for (chr in unique(segments_to_run$chr)) {
  
  chr_segments <- segments_to_run[segments_to_run$chr == chr]
  
  cat(
    "\nProcessing chromosome:", chr,
    "| segments:", nrow(chr_segments), "\n"
  )
  
  model_file <- file.path(
    model_dir,
    paste0("models_", chr, "_T1.rds")
  )
  
  if (!file.exists(model_file)) {
    cat("No model file for:", chr, "\n")
    next
  }
  
  chr_models <- readRDS(model_file)
  
  missing_models <- setdiff(
    chr_segments$segment_id,
    names(chr_models)
  )
  
  cat(
    "Matching models:",
    nrow(chr_segments) - length(missing_models),
    "/", nrow(chr_segments), "\n"
  )
  
  if (length(missing_models) > 0) {
    cat("Missing models:", length(missing_models), "\n")
  }
  
  
  # ----------------------------------------------------------
  # SEGMENT LOOP
  # ----------------------------------------------------------
  
  for (i in seq_len(nrow(chr_segments))) {
    
    seg <- chr_segments[i]
    segment_id <- seg$segment_id
    
    cat("\nSegment:", segment_id, "\n")
    
    output_file <- file.path(
      output_dir,
      paste0(segment_id, "_simulations.rds")
    )
    
    
    # --------------------------------------------------------
    # 8.1 GET MODEL AND PHI
    # --------------------------------------------------------
    
    fitted_model <- chr_models[[segment_id]]
    
    if (is.null(fitted_model)) {
      cat("No fitted model for:", segment_id, "\n")
      next
    }
    
    phi <- as.numeric(sigma(fitted_model))
    
    if (!is.finite(phi) || phi <= 0) {
      cat("Invalid phi for:", segment_id, "\n")
      next
    }
    
    
    # --------------------------------------------------------
    # 8.2 EXTRACT OBSERVED COUNTS
    # --------------------------------------------------------
    
    meth_counts <- as.numeric(
      meth_mat[segment_id, ]
    )
    
    unmeth_counts <- as.numeric(
      unmeth_mat[segment_id, ]
    )
    
    coverage <- meth_counts + unmeth_counts
    
    valid <- is.finite(meth_counts) &
      is.finite(unmeth_counts) &
      is.finite(coverage) &
      coverage > 0 &
      meth_counts >= 0 &
      unmeth_counts >= 0 &
      meth_counts <= coverage
    
    if (sum(valid) < min_samples) {
      cat(
        "Too few valid samples:",
        segment_id,
        "|", sum(valid), "\n"
      )
      next
    }
    
    dat <- data.table(
      sample = t1_cols[valid],
      ID_ma = sample_map$ID_ma[valid],
      ID_pa = sample_map$ID_pa[valid],
      meth = meth_counts[valid],
      unmeth = unmeth_counts[valid],
      coverage = coverage[valid]
    )
    
    if (anyNA(dat$ID_ma) || anyNA(dat$ID_pa)) {
      cat("Missing parental IDs:", segment_id, "\n")
      next
    }
    
    dat[, observed_methylation := meth / coverage]
    
    
    # --------------------------------------------------------
    # 8.3 OBSERVED TARGET STATISTICS
    # --------------------------------------------------------
    
    target_mean <- mean(dat$observed_methylation)
    target_variance <- var(dat$observed_methylation)
    
    if (
      !is.finite(target_mean) ||
      target_mean <= 0 ||
      target_mean >= 1 ||
      !is.finite(target_variance) ||
      target_variance <= 0
    ) {
      cat(
        "Invalid target statistics:",
        segment_id,
        "| mean:", target_mean,
        "| variance:", target_variance, "\n"
      )
      next
    }
    
    mean_lower <- target_mean * 0.95
    mean_upper <- target_mean * 1.05
    
    var_lower <- target_variance * 0.95
    var_upper <- target_variance * 1.05
    
    cat(
      "Samples:", nrow(dat),
      "| Mean:", round(target_mean, 4),
      "| Variance:", signif(target_variance, 4),
      "| Phi:", round(phi, 3), "\n"
    )
    
    
    # --------------------------------------------------------
    # 8.4 CHECKPOINT
    # --------------------------------------------------------
    
    checkpoint_file <- file.path(
      checkpoint_dir,
      paste0(segment_id, "_checkpoint.rds")
    )
    
    if (file.exists(checkpoint_file)) {
      
      checkpoint <- readRDS(checkpoint_file)
      
      next_combination <- checkpoint$next_combination
      accepted_chunks <- checkpoint$accepted_chunks
      n_accepted <- checkpoint$n_accepted
      
      cat(
        "Resuming", segment_id,
        "from combination", next_combination, "\n"
      )
      
    } else {
      
      next_combination <- 1L
      accepted_chunks <- list()
      n_accepted <- 0L
    }
    
    
    # --------------------------------------------------------
    # 8.5 SIMULATE
    # --------------------------------------------------------
    
    if (next_combination <= nrow(sd_grid)) {
      
      for (g in seq.int(next_combination, nrow(sd_grid))) {
        
        pars <- sd_grid[g, ]
        
        for (r in seq_len(n_replicates)) {
          
          sim <- simulate_segment(
            dat = dat,
            sd_ma = pars$sd_ma,
            sd_pa = pars$sd_pa,
            sd_env = pars$sd_env,
            phi = phi,
            target_mean = target_mean
          )
          
          mean_ok <- sim$mean >= mean_lower &&
            sim$mean <= mean_upper
          
          variance_ok <- sim$variance >= var_lower &&
            sim$variance <= var_upper
          
          if (mean_ok && variance_ok) {
            
            sim_dt <- copy(dat)
            
            sim_dt[, `:=`(
              coverage = sim$coverage,
              meth = sim$meth,
              unmeth = sim$unmeth,
              methylation = sim$meth / sim$coverage,
              sd_ma = pars$sd_ma,
              sd_pa = pars$sd_pa,
              sd_env = pars$sd_env,
              phi = phi,
              combination = g,
              replicate = r,
              simulated_mean = sim$mean,
              simulated_variance = sim$variance
            )]
            
            accepted_chunks[[length(accepted_chunks) + 1L]] <- sim_dt
            
            n_accepted <- n_accepted + 1L
          }
        }
        
        
        # ----------------------------------------------------
        # 8.6 SAVE CHECKPOINT
        # ----------------------------------------------------
        
        if (
          g %% checkpoint_every == 0 ||
          g == nrow(sd_grid)
        ) {
          
          saveRDS(
            list(
              next_combination = g + 1L,
              accepted_chunks = accepted_chunks,
              n_accepted = n_accepted
            ),
            checkpoint_file
          )
          
          cat(
            segment_id,
            "| Combination:", g, "/", nrow(sd_grid),
            "| Accepted:", n_accepted, "\n"
          )
          
          gc(verbose = FALSE)
        }
      }
    }
    
    
    # --------------------------------------------------------
    # 8.7 SAVE COMPLETED SEGMENT
    # --------------------------------------------------------
    
    accepted_dt <- if (length(accepted_chunks) > 0) {
      rbindlist(
        accepted_chunks,
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      data.table()
    }
    
    saveRDS(accepted_dt, output_file)
    
    if (file.exists(checkpoint_file)) {
      file.remove(checkpoint_file)
    }
    
    cat(
      "Completed:", segment_id,
      "| Accepted simulations:", n_accepted, "\n"
    )
    
    
    # --------------------------------------------------------
    # 8.8 UPDATE SUMMARY
    # --------------------------------------------------------
    
    if (file.exists(summary_file)) {
      
      summary_dt <- fread(summary_file)
      
      summary_dt <- summary_dt[
        segment_id != seg$segment_id
      ]
      
    } else {
      
      summary_dt <- data.table()
    }
    
    new_summary <- data.table(
      segment_id = segment_id,
      target_mean = target_mean,
      target_variance = target_variance,
      phi = phi,
      n_samples = nrow(dat),
      n_combinations = nrow(sd_grid),
      n_attempted = nrow(sd_grid) * n_replicates,
      n_accepted = n_accepted,
      acceptance_rate = n_accepted /
        (nrow(sd_grid) * n_replicates)
    )
    
    summary_dt <- rbindlist(
      list(summary_dt, new_summary),
      use.names = TRUE,
      fill = TRUE
    )
    
    fwrite(summary_dt, summary_file)
    
    rm(accepted_chunks, accepted_dt)
    gc(verbose = FALSE)
  }
  
  rm(chr_models)
  gc(verbose = FALSE)
}

cat("\nAll chromosome processing finished.\n")
