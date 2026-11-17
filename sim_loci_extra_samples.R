#we need to simulate data with new samples for the loci
#we can extract the new samples we created for the segment based analysis
library(data.table)

simulate_new_metadata <- function(
    metadata,
    timepoint = "T1",
    n_new_per_family = 0,
    weight_noise_sd = 0.05,
    seed = 123
) {
  set.seed(seed)
  
  metadata <- as.data.table(copy(metadata))
  
  required <- c("sample", "TimePoint", "ID_ma", "ID_pa", "Weight.g")
  missing_cols <- setdiff(required, names(metadata))
  
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  if (n_new_per_family < 0 ||
      n_new_per_family != floor(n_new_per_family)) {
    stop("n_new_per_family must be a non-negative integer.")
  }
  
  # Existing samples at the requested timepoint
  existing <- metadata[TimePoint == timepoint]
  
  if (nrow(existing) == 0L) {
    stop("No metadata rows found for timepoint ", timepoint)
  }
  
  if (n_new_per_family == 0L) {
    return(existing)
  }
  
  # One row per distinct maternal-paternal pair
  families <- unique(
    existing[
      !is.na(ID_ma) & !is.na(ID_pa),
      .(ID_ma, ID_pa)
    ]
  )
  
  all_weights <- existing[
    is.finite(Weight.g) & Weight.g >= 0,
    Weight.g
  ]
  
  
  new_rows <- vector("list", nrow(families) * n_new_per_family)
  k <- 0L
  
  for (i in seq_len(nrow(families))) {
    
    ma <- families$ID_ma[i]
    pa <- families$ID_pa[i]
    
    for (j in seq_len(n_new_per_family)) {
      k <- k + 1L
      
      sampled_weight <- sample(all_weights, 1L)
      
      # multiplicative noise
      new_weight <- sampled_weight *
        (1 + rnorm(1L, mean = 0, sd = weight_noise_sd))
      
      new_rows[[k]] <- data.table(
        sample = paste0("SIM", sprintf("%03d", k), "_", timepoint),
        TimePoint = timepoint,
        ID_ma = ma,
        ID_pa = pa,
        Weight.g = max(0, new_weight)
      )
    }
  }
  
  new_metadata <- rbindlist(new_rows)
  
  # Preserve all original columns; additional columns for simulated
  # samples are NA unless they are explicitly assigned above.
  extra_cols <- setdiff(names(metadata), names(new_metadata))
  
  for (col in extra_cols) {
    new_metadata[, (col) := metadata[[col]][NA_integer_]]
  }
  
  setcolorder(new_metadata, names(metadata))
  
  combined <- rbindlist(
    list(existing, new_metadata),
    use.names = TRUE,
    fill = TRUE
  )
  
  combined[, Weight_sc := as.numeric(scale(Weight.g))]
  
  combined[]
}


path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"

meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]

meta[, sample := as.character(sample)]

meta_sim <- simulate_new_metadata(
  metadata = meta,
  timepoint = "T1",
  n_new_per_family = 20,
  weight_noise_sd = 0.05,
  seed = 123
)


simulate_empirical_locus_coverage <- function(
    db,
    combined_meta,
    timepoint = "T1",
    n_new_per_family = 0,
    n_bases = NULL,
    output_dir = "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/with_20_new_samples",
    seed = 123,
    noise_sd = 0.05,
    weight_noise_sd = 0.05
) {
  
  library(data.table)
  
  set.seed(seed)

  simulated_sample_ids <- combined_meta$sample
  
  # get coverage
  
  locus_data <- as.data.table(methylKit::getData(db))
  
  sample_ids_db <- db@sample.ids
  coverage_indices <- db@coverage.index

  coverage_cols <- names(locus_data)[coverage_indices]
  
  observed_ids <- sample_ids_db[
    grepl(paste0("_", timepoint, "$"), sample_ids_db)
  ]
  
  observed_cols <- coverage_cols[
    match(observed_ids, sample_ids_db, nomatch = 0L)
  ]

  
  locus_data[, .original_row := .I]
  
  if (!is.null(n_bases)) {
    
    setorderv(locus_data, c("chr", "start", "end"))
    
    selected_rows <- locus_data[
      ,
      {
        n_select <- min(.N, n_bases)
        .SD[sample.int(.N, n_select)]$ .original_row
      },
      by = chr
    ]$V1
    
    locus_data <- locus_data[
      match(selected_rows, .original_row)
    ]
  }
  
  setorder(locus_data, chr, start, end)
  
  loci <- locus_data[, .(chr, start, end)]
  
  coverage_observed <- as.matrix(
    locus_data[, ..observed_cols]
  )
  
  storage.mode(coverage_observed) <- "numeric"
  
  n_loci <- nrow(coverage_observed)
  n_samples <- length(simulated_sample_ids)
  
  simulated_coverage <- matrix(
    NA_real_,
    nrow = n_loci,
    ncol = n_samples,
    dimnames = list(NULL, simulated_sample_ids)
  )
  
  for (i in seq_len(n_loci)) {
    
    observed <- coverage_observed[i, ]
    
    # we retain NAs
    sampled <- sample(
      observed,
      size = n_samples,
      replace = TRUE
    )
    
    # Add noise only to finite, positive coverage values
    positive <- is.finite(sampled) & sampled > 0
    
    sampled[positive] <- pmax(
      round(
        sampled[positive] *
          (1 + rnorm(
            sum(positive),
            mean = 0,
            sd = noise_sd
          ))
      ),
      1
    )
    
    # Zeros remain zero; sampled NAs remain NA
    simulated_coverage[i, ] <- sampled
  }
  
  selection_label <- if (is.null(n_bases)) {
    "all"
  } else {
    paste0(n_bases, "per_chr")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  output_file <- file.path(
    output_dir,
    paste0(
      "simulated_with_", n_new_per_family,
      "_extra_samples_coverage_", selection_label, ".rds"
    )
  )
  
  result <- list(
    timepoint = timepoint,
    sample_ids = simulated_sample_ids,
    loci = loci,
    coverage = simulated_coverage,
    metadata = combined_meta,
    n_loci = n_loci,
    n_samples = n_samples,
    n_new_samples = nrow(combined_meta),
    n_new_per_family = n_new_per_family,
    noise_sd = noise_sd,
    seed = seed
  )
  
  saveRDS(result, output_file)
  
  fwrite(
    combined_meta,
    file.path(
      output_dir,
      paste0("simulated_with_", n_new_per_family, "_extra_samples_metadata.csv")
    )
  )
  
  message("Timepoint: ", timepoint)
  message("Loci: ", n_loci)
  message("Total simulated samples: ", n_samples)
  message("New samples: ", nrow(combined_meta))
  message("Saved: ", output_file)
  
  invisible(result)
}

db=readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds")


coverage_sim <- simulate_empirical_locus_coverage(
  db = db,
  combined_meta = meta_sim,
  timepoint = "T1",
  n_new_per_family = 20,
  n_bases = 1000,
  noise_sd = 0.05,
  weight_noise_sd = 0.05,
  seed = 123
)
