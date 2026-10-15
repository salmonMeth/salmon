library(data.table)
library(methylKit)

simulate_empirical_locus_coverage <- function(
    db,
    timepoint = "T1",
    n_bases = NULL,
    output_dir = "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci",
    seed = 123,
    noise_sd = 0.05
) {
  
  if (!is.null(n_bases)) {
    if (length(n_bases) != 1L ||
        !is.finite(n_bases) ||
        n_bases < 1 ||
        n_bases != floor(n_bases)) {
      stop("n_bases must be NULL or a positive integer.")
    }
    n_bases <- as.integer(n_bases)
  }
  
  if (length(timepoint) != 1L ||
      !grepl("^[A-Za-z0-9]+$", timepoint)) {
    stop("timepoint must be a single valid timepoint label, e.g. 'T1'.")
  }
  
  if (length(noise_sd) != 1L ||
      !is.finite(noise_sd) ||
      noise_sd < 0) {
    stop("noise_sd must be a non-negative number.")
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  set.seed(seed)
  
  
  locus_data <- as.data.table(methylKit::getData(db))
  
  required_coords <- c("chr", "start", "end")
  
  if (!all(required_coords %in% names(locus_data))) {
    stop("The database data must contain chr, start and end columns.")
  }
  
  # Identify the coverage columns using methylBaseDB metadata
  sample_ids <- db@sample.ids
  coverage_indices <- db@coverage.index
  
  if (length(sample_ids) != length(coverage_indices)) {
    stop("Sample IDs and coverage indices do not have matching lengths.")
  }
  
  if (any(coverage_indices < 1L |
          coverage_indices > ncol(locus_data))) {
    stop("Coverage indices do not match the retrieved data columns.")
  }
  
  coverage_cols <- names(locus_data)[coverage_indices]
  
  
  keep_samples <- grepl(
    paste0("_", timepoint, "$"),
    sample_ids
  )
  
  if (!any(keep_samples)) {
    stop("No samples found for timepoint ", timepoint, ".")
  }
  
  selected_sample_ids <- sample_ids[keep_samples]
  selected_coverage_cols <- coverage_cols[keep_samples]
  
  # pick loci
  
  locus_data[, .original_row := .I]
  
  if (!is.null(n_bases)) {
    setorderv(locus_data, c("chr", "start", "end"))
    
    selected_rows <- locus_data[
      ,
      {
        n_select <- min(.N, n_bases)
        .SD[sample.int(.N, size = n_select)]
      },
      by = chr,
      .SDcols = names(locus_data)
    ]$ .original_row
    
    locus_data <- locus_data[
      match(selected_rows, .original_row)
    ]
  }
  
  # Restore genomic ordering
  setorder(locus_data, chr, start, end)
  
  locus_coords <- locus_data[
    ,
    c("chr", "start", "end"),
    with = FALSE
  ]
  
  # -----------------------------------------------
  # 5. Extract observed coverage
  # -----------------------------------------------
  
  coverage_observed <- as.matrix(
    locus_data[, ..selected_coverage_cols]
  )
  
  storage.mode(coverage_observed) <- "numeric"
  colnames(coverage_observed) <- selected_sample_ids
  
  n_loci <- nrow(coverage_observed)
  n_samples <- length(selected_sample_ids)
  
  # -----------------------------------------------
  # 6. Simulate empirical coverage
  # -----------------------------------------------
  
  simulated_coverage <- matrix(
    NA_real_,
    nrow = n_loci,
    ncol = n_samples,
    dimnames = list(NULL, selected_sample_ids)
  )
  
  for (i in seq_len(n_loci)) {
    
    observed <- coverage_observed[i, ]
    
    # Resample only from observed, non-missing coverage values.
    # Keep the original NA positions unchanged.
    valid <- is.finite(observed)
    
    if (!any(valid)) {
      next
    }
    
    n_valid <- sum(valid)
    
    sampled <- sample(
      observed[valid],
      size = n_valid,
      replace = TRUE
    )
    
    noise <- rnorm(
      n_valid,
      mean = 0,
      sd = noise_sd
    )
    
    simulated <- sampled
    
    # Add noise only to positive coverage values.
    positive <- sampled > 0
    
    simulated[positive] <- pmax(
      round(sampled[positive] * (1 + noise[positive])),
      1
    )
    
    # Preserve zero coverage and leave missing positions as NA.
    simulated_coverage[i, valid] <- simulated
  }
  # -----------------------------------------------
  # 7. Save result
  # -----------------------------------------------
  
  selection_label <- if (is.null(n_bases)) {
    "all"
  } else {
    paste0(n_bases, "per_chr")
  }
  
  output_file <- file.path(
    output_dir,
    paste0(
      "simulated_", timepoint,
      "_locus_coverage_", selection_label,
      ".rds"
    )
  )
  
  result <- list(
    timepoint = timepoint,
    sample_ids = selected_sample_ids,
    loci = locus_coords,
    coverage = simulated_coverage,
    n_loci = n_loci,
    n_samples = n_samples,
    n_bases_per_chromosome = n_bases,
    noise_sd = noise_sd,
    seed = seed
  )
  
  saveRDS(result, output_file)
  
  message("Timepoint: ", timepoint)
  message("Samples: ", n_samples)
  message("Selected loci: ", n_loci)
  message("Output: ", output_file)
  
  invisible(result)
}


##############simulate Cs values, keep the ones with NA coverage NA


simulate_empirical_Cs <- function(
    db,
    simulated_coverage,
    noise_sd = 0.02,
    seed = 123
) {
  set.seed(seed)
  
  # Simulated coverage object: result returned by your coverage function
  coverage_mat <- as.matrix(simulated_coverage$coverage)
  loci <- as.data.table(simulated_coverage$loci)
  sample_ids <- simulated_coverage$sample_ids
  
  if (is.null(colnames(coverage_mat))) {
    stop("The simulated coverage matrix must have sample IDs as column names.")
  }
  
  # Align columns and check dimensions
  if (!all(sample_ids %in% colnames(coverage_mat))) {
    stop("Sample IDs do not match the simulated coverage matrix.")
  }
  
  coverage_mat <- coverage_mat[, sample_ids, drop = FALSE]
  
  # Retrieve original methylated-read counts from methylBaseDB
  original <- as.data.table(methylKit::getData(db))
  db_samples <- db@sample.ids
  cs_cols <- names(original)[db@numCs.index]
  
  # Find corresponding T1 samples
  sample_idx <- match(sample_ids, db_samples)
  
  if (anyNA(sample_idx)) {
    stop("Some simulated samples are not present in the original database.")
  }
  
  original_cs <- as.matrix(
    original[, cs_cols[sample_idx], with = FALSE]
  )
  storage.mode(original_cs) <- "numeric"
  
  # Match simulated loci to the original database
  original_key <- original[, .(
    chr,
    start,
    end,
    original_row = .I
  )]
  
  loci_key <- copy(loci)
  loci_key[, simulated_row := .I]
  
  matches <- merge(
    loci_key,
    original_key,
    by = c("chr", "start", "end"),
    all.x = TRUE,
    sort = FALSE
  )
  
  setorder(matches, simulated_row)
  
  if (anyNA(matches$original_row)) {
    stop(
      sum(is.na(matches$original_row)),
      " simulated loci could not be matched to db by coordinates."
    )
  }
  
  n_loci <- nrow(coverage_mat)
  n_samples <- ncol(coverage_mat)
  
  if (nrow(loci) != n_loci) {
    stop("The number of loci and coverage rows do not match.")
  }
  
  # Initialize with NA; only fill positions with simulated coverage
  simulated_Cs <- matrix(
    NA_integer_,
    nrow = n_loci,
    ncol = n_samples,
    dimnames = list(NULL, sample_ids)
  )
  
  diagnostics <- data.table(
    chr = loci$chr,
    start = loci$start,
    end = loci$end,
    original_mean_Cs = NA_real_,
    simulated_mean_Cs = NA_real_,
    n_observed = 0L,
    mean_matched = FALSE
  )
  
  for (i in seq_len(n_loci)) {
    
    cov_i <- coverage_mat[i, ]
    cs_i <- original_cs[matches$original_row[i], ]
    
    # NA coverage must imply NA Cs
    valid_cov <- is.finite(cov_i) & cov_i >= 0
    
    # Only original Cs values that are observed and valid
    valid_cs <- is.finite(cs_i) & cs_i >= 0
    
    if (!any(valid_cov) || !any(valid_cs)) {
      next
    }
    
    original_mean <- mean(cs_i[valid_cs])
    
    # If there are fewer original observed values than samples,
    # resample with replacement; otherwise shuffle the observed values.
    values <- cs_i[valid_cs]
    n_valid <- sum(valid_cov)
    
    shuffled <- sample(
      values,
      size = n_valid,
      replace = length(values) < n_valid
    )
    
    # Add slight multiplicative noise
    noisy <- round(
      shuffled * (1 + rnorm(n_valid, mean = 0, sd = noise_sd))
    )
    
    noisy <- pmax(noisy, 0)
    
    # Respect each sample's simulated coverage
    
    # Respect each sample's simulated coverage
    capacity <- floor(cov_i[valid_cov])
    noisy <- pmin(pmax(noisy, 0), capacity)
    
    # Target sum based on original mean Cs
    target_sum <- round(original_mean * n_valid)
    
    # Exact target may be impossible if coverage is too low
    target_sum <- min(target_sum, sum(capacity))
    
    current_sum <- sum(noisy)
    
    if (current_sum < target_sum) {
      while (current_sum < target_sum) {
        can_increase <- which(noisy < capacity)
        if (length(can_increase) == 0L) break
        
        j <- sample(can_increase, 1)
        noisy[j] <- noisy[j] + 1L
        current_sum <- current_sum + 1
      }
    } else if (current_sum > target_sum) {
      while (current_sum > target_sum) {
        can_decrease <- which(noisy > 0)
        if (length(can_decrease) == 0L) break
        
        j <- sample(can_decrease, 1)
        noisy[j] <- noisy[j] - 1L
        current_sum <- current_sum - 1
      }
    }
    
    # Final bounds check before assignment
    noisy <- as.integer(pmin(pmax(noisy, 0), capacity))
    
    simulated_Cs[i, valid_cov] <- noisy
    
    diagnostics$n_observed[i] <- n_valid
    diagnostics$original_mean_Cs[i] <- original_mean
    diagnostics$simulated_mean_Cs[i] <- mean(noisy)
    diagnostics$mean_matched[i] <-
      isTRUE(all.equal(mean(noisy), original_mean, tolerance = 1 / n_valid))
  }
  
  # Final consistency checks
  stopifnot(
    all(is.na(simulated_Cs[is.na(coverage_mat)])),
    all(simulated_Cs[!is.na(simulated_Cs)] ==
          floor(simulated_Cs[!is.na(simulated_Cs)])),
    all(simulated_Cs[!is.na(simulated_Cs)] >= 0)
  )
  
  observed_pairs <- !is.na(simulated_Cs)
  
  if (any(simulated_Cs[observed_pairs] > coverage_mat[observed_pairs])) {
    stop("Some simulated Cs values exceed simulated coverage.")
  }
  
  list(
    loci = loci,
    sample_ids = sample_ids,
    coverage = coverage_mat,
    Cs = simulated_Cs,
    diagnostics = diagnostics,
    noise_sd = noise_sd,
    seed = seed
  )
}