#vary rho values

#we simulate data and also generate new samples if we want to


library(data.table)
library(VGAM)

set.seed(1293)

timepoint <- "T1"

#change accordingly
extra_samples <- 0L

data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_cs_cov_3_samples_removed"


path_ids <-"/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"

bed_dir <-"/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/beds"

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/data_with_varying_rho/"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  out_dir,
  paste0("simulated_", timepoint, ".rds")
)

stopifnot(
  length(extra_samples) == 1L,
  !is.na(extra_samples),
  extra_samples >= 0,
  extra_samples == as.integer(extra_samples)
)

extra_samples <- as.integer(extra_samples)


# get coverage

chromosomes <- sprintf(
  "NC_%06d.1",
  59442:59470
)

coverage_files <- file.path(
  data_dir,
  paste0(chromosomes, "_segment_coverage.csv")
)

missing_files <- coverage_files[!file.exists(coverage_files)]

if (length(missing_files) > 0L) {
  stop(
    "Missing coverage files:\n",
    paste(missing_files, collapse = "\n")
  )
}

coverage <- rbindlist(
  lapply(coverage_files, function(file) {
    
    x <- fread(file)
    
    sample_cols <- grep(
      paste0("_", timepoint, "$"),
      names(x),
      value = TRUE
    )
    
    x[
      ,
      c("chr", "start", "end", sample_cols),
      with = FALSE
    ]
  }),
  use.names = TRUE,
  fill = TRUE
)

coverage[
  ,
  `:=`(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end)
  )
]

sample_cols_original <- grep(
  paste0("_", timepoint, "$"),
  names(coverage),
  value = TRUE
)


segments <- unique(
  coverage[, .(chr, start, end)]
)

setorder(segments, chr, start, end)

segments[, segment := .I]

n_segments <- nrow(segments)

meta <- fread(path_ids)

meta[
  ,
  sample := paste(Running.ID, TimePoint, sep = "_")
]

meta <- meta[
  TimePoint == timepoint &
    sample != "176_T1"
]

if (anyDuplicated(meta$sample)) {
  stop("Original metadata contains duplicate sample names.")
}

# Make sure each original sample has one family assignment
family_check <- meta[
  ,
  .(
    n_maternal = uniqueN(ID_ma),
    n_paternal = uniqueN(ID_pa)
  ),
  by = sample
]


# Create extra offspring for every distinct parental pair
if (extra_samples > 0L) {
  
  family_templates <- meta[
    ,
    .SD[1L],
    by = .(ID_ma, ID_pa)
  ]
  
  extra_meta <- family_templates[
    rep(seq_len(nrow(family_templates)),
        each = extra_samples)
  ]
  
  extra_meta[
    ,
    Running.ID := paste0(
      "EXTRA",
      sprintf("%05d", seq_len(.N))
    )
  ]
  
  extra_meta[
    ,
    Weight.g := sample(
      meta$Weight.g,
      size = .N,
      replace = TRUE
    )
  ]
  
  extra_meta[
    ,
    sample := paste(Running.ID, TimePoint, sep = "_")
  ]
  
  meta <- rbindlist(
    list(meta, extra_meta),
    use.names = TRUE,
    fill = TRUE
  )
  
  message(
    "Added ", nrow(extra_meta),
    " offspring across ",
    nrow(family_templates),
    " distinct parental pairs."
  )
}

if (anyDuplicated(meta$sample)) {
  stop("Duplicate sample names after adding offspring.")
}

# Scale weights using the combined sample set
meta[
  ,
  Weight_sc := as.numeric(scale(as.numeric(Weight.g)))
]


parent_info <- unique(
  meta[, .(sample, ID_ma, ID_pa)]
)

family_counts <- meta[
  ,
  .(n_offspring = .N),
  by = .(ID_ma, ID_pa)
]


#get coverage
coverage_long <- melt(
  coverage,
  id.vars = c("chr", "start", "end"),
  measure.vars = sample_cols_original,
  variable.name = "sample",
  value.name = "coverage"
)

coverage_long[, sample := as.character(sample)]

coverage_long <- merge(
  coverage_long,
  segments,
  by = c("chr", "start", "end"),
  all.x = TRUE,
  sort = FALSE
)

observed_coverage <- coverage_long[
  is.finite(coverage) & coverage > 0,
  .(observed = list(coverage)),
  by = segment
]

sample_names <- as.character(meta$sample)
n_samples <- length(sample_names)

simulated_coverage <- matrix(
  NA_real_,
  nrow = n_segments,
  ncol = n_samples,
  dimnames = list(NULL, sample_names)
)

for (i in seq_len(n_segments)) {
  
  obs <- observed_coverage[
    segment == i,
    observed[[1L]]
  ]
  
  if (length(obs) == 0L) {
    next
  }
  
  sampled <- sample(
    obs,
    size = n_samples,
    replace = TRUE
  )
  
  # Add 5% proportional noise
  noise <- rnorm(
    n_samples,
    mean = 0,
    sd = 0.05
  )
  
  simulated_coverage[i, ] <- pmax(
    round(sampled * (1 + noise)),
    1
  )
}

simulated_coverage <- as.data.table(
  simulated_coverage
)

simulated_segments <- cbind(
  segments,
  simulated_coverage
)

coverage_sim_long <- melt(
  simulated_segments,
  id.vars = c("segment", "chr", "start", "end"),
  measure.vars = sample_names,
  variable.name = "sample",
  value.name = "coverage"
)

coverage_sim_long[, sample := as.character(sample)]

coverage_rds <- file.path(
  out_dir,
  paste0("coverage_", extra_samples,"_extra_samples_per_fam" ,".rds")
)

saveRDS(
  simulated_segments,
  coverage_rds
)


library(data.table)
library(VGAM)
n_segments <- nrow(segments)

# Output files
output_file <- file.path(out_dir, "simulated_data.rds")

# Sample names
sample_names <- meta$sample


combs <- sample(
  c(
    "none",
    "maternal_dominant",
    "paternal_dominant",
    "both"
  ),
  size = n_segments,
  replace = TRUE,
  prob = c(0.25, 0.25, 0.25, 0.25)
)

weight_betas <- rnorm(
  n_segments,
  mean = 0,
  sd = 0.05
)

sim_parameters <- data.table(
  segment = segments$segment,
  combs = combs,
  weight_beta = weight_betas
)

sim_parameters[
  ,
  `:=`(
    true_sd_ma = 0,
    true_sd_pa = 0
  )
]

# No heritability
sim_parameters[
  combs == "none",
  `:=`(
    true_sd_ma = 1e-9,
    true_sd_pa = 1e-9
  )
]

# Maternal dominant
sim_parameters[
  combs == "maternal_dominant",
  `:=`(
    true_sd_ma = runif(.N, 0.12, 0.30),
    true_sd_pa = runif(.N, 0.00, 0.05)
  )
]

# Paternal dominant
sim_parameters[
  combs == "paternal_dominant",
  `:=`(
    true_sd_ma = runif(.N, 0.00, 0.05),
    true_sd_pa = runif(.N, 0.12, 0.30)
  )
]

# Both parents contribute
sim_parameters[
  combs == "both",
  `:=`(
    true_sd_ma = runif(.N, 0.10, 0.30),
    true_sd_pa = runif(.N, 0.10, 0.30)
  )
]

# Calculate variances
sim_parameters[
  ,
  `:=`(
    true_var_ma = true_sd_ma^2,
    true_var_pa = true_sd_pa^2
  )
]

# ============================================================
# 3. GENERATE VARYING RHO VALUES
# ============================================================

#sim_parameters[,rho := ifelse(
    #runif(.N) < 0.80,rbeta(.N, shape1 = 1.5, shape2 = 3) * 0.035,runif(.N, 0.025, 0.04))]                 #ifelse(runif(.N) < 0.8,rbeta(.N, shape1 = 1.5, shape2 = 3) * 0.04, rbeta(.N, shape1 = 1.5, shape2 = 0.8) * 0.50) ]

 #rho := rlnorm(.N, meanlog = log(0.008), sdlog = 1.2)]
 #rho := ifelse(runif(.N) < 0.8,runif(.N, 0, 0.008),runif(.N, 0.008, 0.026))]
 #rho := rbeta(.N, shape1 = 2, shape2 = 8) * 0.026]
 #rho := runif(.N, 0, .026)]


# Alternatively, use specific discrete rho values:
#
#rho_values <- c(0, 0.01,.03,.05,.07,.09)
#rho_values = seq(0,.03,length.out=10)
#rho_values=c(0,.002,0.003,.004,.005,.006,.008,.01,.012,.013,.015,.018,.02,.023,.027)
#seq(0,.01,length.out=10)
#
#alternatively, we force there to be a certain number of large values

sim_parameters[, rho := ifelse(
  runif(.N) < 0, NA_real_, NA_real_
)]

high_rho_idx <- sample(seq_len(nrow(sim_parameters)), 200)

sim_parameters[
  ,rho := ifelse(seq_len(.N) %in% high_rho_idx,runif(.N, 0.2, 0.5),
                 ifelse(runif(.N) < 0.80,rbeta(.N, shape1 = 2, shape2 = 1.4) * 0.2,runif(.N, 0.2, 0.4)))]

# Inspect rho distribution
print(
  sim_parameters[
    ,
    .N,
    by = rho
  ][order(rho)]
)

bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

all_segments <- rbindlist(
  lapply(bed_files, function(file) {
    
    x <- fread(file, header = FALSE)
    
    data.table(
      chr = as.character(x[[1L]]),
      start = as.integer(x[[2L]]),
      end = as.integer(x[[3L]]),
      mean_methylation = as.numeric(x[[5L]])
    )
  }),
  use.names = TRUE
)

# Convert chromosome names
all_segments[
  ,
  chr := sprintf(
    "NC_%06d.1",
    59441L + as.integer(chr)
  )
]

# Average methylation across BED files
all_segments <- all_segments[
  ,
  .(
    mean_methylation = mean(
      mean_methylation,
      na.rm = TRUE
    )
  ),
  by = .(chr, start, end)
]

all_segments[
  !is.finite(mean_methylation),
  mean_methylation := NA_real_
]

# Match segments
segment_info <- merge(
  segments,
  all_segments,
  by = c("chr", "start", "end"),
  all.x = TRUE,
  sort = FALSE
)

# Calculate baseline on logit scale
segment_info[
  ,
  baseline := qlogis(
    pmin(
      pmax(mean_methylation, 1e-6),
      1 - 1e-6
    )
  )
]

# Add segment information to parameters
sim_parameters <- merge(
  sim_parameters,
  segment_info[
    ,
    .(
      segment,
      chr,
      start,
      end,
      mean_methylation,
      baseline
    )
  ],
  by = "segment",
  all.x = TRUE,
  sort = FALSE
)

setorder(sim_parameters, segment)

if (sim_parameters[, anyNA(baseline)]) {
  stop("Missing segment baselines.")
}

maternal_re <- CJ(
  segment = sim_parameters$segment,
  ID_ma = unique(meta$ID_ma),
  unique = TRUE
)

maternal_re <- merge(
  maternal_re,
  sim_parameters[
    ,
    .(segment, true_sd_ma)
  ],
  by = "segment",
  all.x = TRUE
)

maternal_re[
  ,
  u_ma := rnorm(
    .N,
    mean = 0,
    sd = true_sd_ma
  )
]


paternal_re <- CJ(
  segment = sim_parameters$segment,
  ID_pa = unique(meta$ID_pa),
  unique = TRUE
)

paternal_re <- merge(
  paternal_re,
  sim_parameters[
    ,
    .(segment, true_sd_pa)
  ],
  by = "segment",
  all.x = TRUE
)

paternal_re[
  ,
  u_pa := rnorm(
    .N,
    mean = 0,
    sd = true_sd_pa
  )
]


sim <- CJ(
  segment = sim_parameters$segment,
  sample = sample_names,
  unique = TRUE
)

# Add sample metadata
sim <- merge(
  sim,
  meta,
  by = "sample",
  all.x = TRUE
)

# Add segment parameters
sim <- merge(
  sim,
  sim_parameters[
    ,
    .(
      segment,
      chr,
      start,
      end,
      mean_methylation,
      baseline,
      weight_beta,
      true_sd_ma,
      true_sd_pa,
      true_var_ma,
      true_var_pa,
      combs,
      rho
    )
  ],
  by = "segment",
  all.x = TRUE
)

# Add maternal effects
sim <- merge(
  sim,
  maternal_re[
    ,
    .(segment, ID_ma, u_ma)
  ],
  by = c("segment", "ID_ma"),
  all.x = TRUE
)

# Add paternal effects
sim <- merge(
  sim,
  paternal_re[
    ,
    .(segment, ID_pa, u_pa)
  ],
  by = c("segment", "ID_pa"),
  all.x = TRUE
)


sim[
  ,
  eta := baseline +
    weight_beta * Weight_sc +
    u_ma +
    u_pa
]

# Transform to probability scale
sim[
  ,
  methylation_probability := plogis(eta)
]


sim[, sample := as.character(sample)]

sim <- merge(
  sim,
  coverage_sim_long[
    ,
    .(segment, sample, coverage)
  ],
  by = c("segment", "sample"),
  all.x = TRUE
)

sim[, coverage := as.integer(coverage)]


sim[, meth_reads := NA_integer_]

idx <- which(
  !is.na(sim$coverage) &
    sim$coverage > 0 &
    is.finite(sim$methylation_probability) &
    !is.na(sim$methylation_probability) &
    !is.na(sim$rho)
)

# Each observation uses its segment-specific rho
sim[
  idx,
  meth_reads := as.integer(
    VGAM::rbetabinom(
      n = length(idx),
      size = coverage[idx],
      prob = methylation_probability[idx],
      rho = rho[idx]
    )
  )
]

# Unmethylated reads
sim[
  ,
  unmeth_reads := coverage - meth_reads
]


setorder(
  sim,
  chr,
  start,
  end,
  sample
)

saveRDS(
  sim,
  output_file
)

# Save segment-level simulation parameters separately
fwrite(
  sim_parameters,
  file.path(
    out_dir,
    "simulation_segment_parameters.csv"
  )
)

# cs/coverage

for (cc in unique(sim$chr)) {
  
  d <- sim[chr == cc]
  
  # Methylated reads
  cs <- dcast(
    d,
    chr + start + end ~ sample,
    value.var = "meth_reads"
  )
  
  setorder(cs, start, end)
  
  fwrite(
    cs,
    file.path(
      out_dir,
      paste0(cc, "_segment_cs.csv")
    )
  )
  
  # Coverage
  cov_dt <- dcast(
    d,
    chr + start + end ~ sample,
    value.var = "coverage"
  )
  
  setorder(cov_dt, start, end)
  
  fwrite(
    cov_dt,
    file.path(
      out_dir,
      paste0(cc, "_segment_coverage.csv")
    )
  )
}


true_sample_parameters <- sim[
  ,
  .(
    segment,
    sample,
    ID_ma,
    ID_pa,
    Weight.g,
    Weight_sc,
    chr,
    start,
    end,
    mean_methylation,
    baseline,
    weight_beta,
    true_sd_ma,
    true_sd_pa,
    true_var_ma,
    true_var_pa,
    rho,
    combs,
    u_ma,
    u_pa,
    eta,
    methylation_probability
  )
]

fwrite(
  true_sample_parameters,
  file.path(
    out_dir,
    "true_sample_parameters.csv"
  )
)


cat("\nSimulation completed!\n")

cat(
  "Total segments:",
  uniqueN(sim$segment),
  "\n"
)

cat(
  "Total observations:",
  nrow(sim),
  "\n"
)

cat(
  "Rho range:",
  range(sim_parameters$rho),
  "\n"
)

cat("\nSegments by parental configuration:\n")
print(sim_parameters[, .N, by = combs])

cat("\nSegments by rho:\n")
print(
  sim_parameters[
    ,
    .(
      n = .N,
      min_rho = min(rho),
      max_rho = max(rho),
      mean_rho = mean(rho)
    )
  ]
)



#fit the model



library(data.table)
library(glmmTMB)


data_dir <- out_dir

model_dir <- file.path(data_dir, "model_results")

dir.create(
  model_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


sim_parameters <- fread(
  file.path(
    data_dir,
    "simulation_segment_parameters.csv"
  )
)

# Keep one row per segment
sim_parameters <- unique(
  sim_parameters,
  by = "segment"
)

# Check required columns
required_cols <- c(
  "segment",
  "chr",
  "start",
  "end",
  "weight_beta",
  "true_sd_ma",
  "true_sd_pa",
  "true_var_ma",
  "true_var_pa",
  "rho",
  "mean_methylation"
)

missing_cols <- setdiff(
  required_cols,
  names(sim_parameters)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing columns in simulation parameters: ",
    paste(missing_cols, collapse = ", ")
  )
}

path_ids <-"/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"


meta <- fread(path_ids)

meta <- as.data.table(meta)


meta[
  ,
  sample := paste(Running.ID, TimePoint, sep = "_")
]

meta <- meta[
  TimePoint == timepoint &
    sample != "176_T1"
]

meta[, sample := as.character(sample)]

meta[
  ,
  Weight_sc := as.numeric(scale(as.numeric(Weight.g)))
]


cs_files <- list.files(
  data_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)

all_results <- list()
all_models <- list()


for (cs_file in cs_files) {
  
  cc <- sub(
    "_segment_cs\\.csv$",
    "",
    basename(cs_file)
  )
  
  cov_file <- file.path(
    data_dir,
    paste0(cc, "_segment_coverage.csv")
  )
  
  if (!file.exists(cov_file)) {
    warning("Missing coverage file: ", cov_file)
    next
  }
  
  cat("\nProcessing chromosome:", cc, "\n")
  
  
  cs <- fread(cs_file)
  cov <- fread(cov_file)
  
  coord_cols <- c(
    "chr",
    "start",
    "end"
  )
  
  sample_cols <- intersect(
    setdiff(names(cs), coord_cols),
    setdiff(names(cov), coord_cols)
  )
  
  if (length(sample_cols) == 0) {
    warning("No shared sample columns: ", cc)
    next
  }
  
  # Ensure matching coordinates and ordering
  if (
    nrow(cs) != nrow(cov) ||
    !identical(
      cs[, ..coord_cols],
      cov[, ..coord_cols]
    )
  ) {
    warning("CS and coverage coordinates differ: ", cc)
    next
  }
  
  
  cs_long <- melt(
    cs,
    id.vars = coord_cols,
    measure.vars = sample_cols,
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = coord_cols,
    measure.vars = sample_cols,
    variable.name = "sample",
    value.name = "coverage"
  )
  
  d <- merge(
    cs_long,
    cov_long,
    by = c(
      "chr",
      "start",
      "end",
      "sample"
    ),
    all = FALSE
  )
  
  
  d[, sample := as.character(sample)]
  
  d <- merge(
    d,
    meta,
    by = "sample",
    all.x = TRUE
  )
  d <- merge(
    d,
    sim_parameters[
      ,
      .(
        chr,
        start,
        end,
        segment,
        weight_beta,
        true_sd_ma,
        true_sd_pa,
        true_var_ma,
        true_var_pa,
        rho,
        mean_methylation
      )
    ],
    by = c("chr", "start", "end"),
    all.x = TRUE
  )
  
  
  d[
    ,
    unmeth_reads := coverage - meth_reads
  ]
  
  d[
    ,
    meth_reads := as.integer(meth_reads)
  ]
  
  d[
    ,
    unmeth_reads := as.integer(unmeth_reads)
  ]
  
  # Remove invalid observations
  d <- d[
    !is.na(meth_reads) &
      !is.na(unmeth_reads) &
      coverage > 0 &
      !is.na(ID_ma) &
      !is.na(ID_pa) &
      is.finite(Weight_sc)
  ]
  
  
  segments_chr <- unique(
    d[
      ,
      .(segment, chr, start, end)
    ]
  )
  
  setorder(segments_chr, start)
  
  chr_results <- list()
  chr_models <- list()
  
  
  for (i in seq_len(nrow(segments_chr))) {
    
    seg <- segments_chr$segment[i]
    
    cat(
      "Fitting segment:",
      seg,
      "/",
      nrow(segments_chr),
      "\n"
    )
    
    current_chr <- segments_chr$chr[i]
    current_start <- segments_chr$start[i]
    current_end <- segments_chr$end[i]
    
    dd <- d[
      segment == seg &
        chr == current_chr &
        start == current_start &
        end == current_end
    ]
    
    # Both random effects must have multiple levels
    if (
      uniqueN(dd$ID_ma) < 2 ||
      uniqueN(dd$ID_pa) < 2
    ) {
      next
    }
    
    
    
    fit <- tryCatch(
      glmmTMB(
        cbind(meth_reads, unmeth_reads) ~
          Weight_sc +
          (1 | ID_ma) +
          (1 | ID_pa),
        data = dd,
        family = glmmTMB::betabinomial(
          link = "logit"
        ),
        REML = TRUE
      ),
      error = function(e) {
        message(
          "FAILED: chr=", current_chr,
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
    
    # fixed eff
    
    fixed <- fixef(fit)$cond
    
    intercept_est <- if (
      "(Intercept)" %in% names(fixed)
    ) {
      unname(fixed["(Intercept)"])
    } else {
      NA_real_
    }
    
    weight_est <- if (
      "Weight_sc" %in% names(fixed)
    ) {
      unname(fixed["Weight_sc"])
    } else {
      NA_real_
    }
    
    #rand effc
    
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
    
    var_ma <- sd_ma^2
    var_pa <- sd_pa^2
    
    
    tr <- unique(
      dd[
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
          true_var_pa,
          rho
        )
      ]
    )
    
    if (nrow(tr) != 1) {
      warning(
        "Multiple or missing parameter records for segment ",
        seg
      )
      next
    }

    
    # glmmTMB beta-binomial dispersion:
    # phi = precision parameter
    # rho = 1 / (phi + 1)
    
    phi_est <- tryCatch(
      {
        unname(
          predict(
            fit,
            type = "disp"
          )[1]
        )
      },
      error = function(e) {
        NA_real_
      }
    )
    
    rho_est <- if (
      is.finite(phi_est) &&
      phi_est > 0
    ) {
      1 / (phi_est + 1)
    } else {
      NA_real_
    }
    
    
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
    
    n_obs <- nrow(dd)
    
    n_ma <- uniqueN(dd$ID_ma)
    n_pa <- uniqueN(dd$ID_pa)
    
    total_reads <- sum(
      dd$meth_reads + dd$unmeth_reads,
      na.rm = TRUE
    )
    
    
    chr_results[[as.character(seg)]] <- data.table(
      
      chr = tr$chr,
      segment = seg,
      start = tr$start,
      end = tr$end,
      
      mean_methylation = tr$mean_methylation,
      
      # Fixed effects
      true_weight_beta = tr$weight_beta,
      estimated_weight_beta = weight_est,
      
      # Maternal random effects
      true_sd_ma = tr$true_sd_ma,
      estimated_sd_ma = sd_ma,
      
      true_var_ma = tr$true_var_ma,
      estimated_var_ma = var_ma,
      
      # Paternal random effects
      true_sd_pa = tr$true_sd_pa,
      estimated_sd_pa = sd_pa,
      
      true_var_pa = tr$true_var_pa,
      estimated_var_pa = var_pa,
      
      # Dispersion
      true_rho = tr$rho,
      estimated_rho = rho_est,
      estimated_dispersion_phi = phi_est,
      
      # Errors
      error_weight = weight_est - tr$weight_beta,
      
      error_sd_ma = sd_ma - tr$true_sd_ma,
      error_sd_pa = sd_pa - tr$true_sd_pa,
      
      error_var_ma = var_ma - tr$true_var_ma,
      error_var_pa = var_pa - tr$true_var_pa,
      
      error_rho = rho_est - tr$rho,
      
      # Intercept
      intercept_estimate = intercept_est,
      
      # Sample information
      n_obs = n_obs,
      n_ma = n_ma,
      n_pa = n_pa,
      total_reads = total_reads,
      
      # Model fit
      AIC = AIC(fit),
      BIC = BIC(fit),
      logLik = as.numeric(logLik(fit)),
      
      # Diagnostics
      convergence = convergence,
      pdHess = pd_hess,
      status = status
    )
    
    chr_models[[as.character(seg)]] <- fit
    
  }
  
  
  if (length(chr_results) > 0) {
    
    chr_results_dt <- rbindlist(
      chr_results,
      use.names = TRUE,
      fill = TRUE
    )
    
    setorder(chr_results_dt, start)
    
    all_results[[cc]] <- chr_results_dt
    all_models[[cc]] <- chr_models
    
    # Save chromosome results
    fwrite(
      chr_results_dt,
      file.path(
        model_dir,
        paste0("model_results_chr", cc, ".csv")
      )
    )
    
    # Save fitted models
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
      "\nSaved ",
      nrow(chr_results_dt),
      " segments for ",
      cc,
      "\n",
      sep = ""
    )
    
  } else {
    
    cat(
      "\nNo successful models for ",
      cc,
      "\n",
      sep = ""
    )
    
  }
}


if (length(all_results) > 0) {
  
  final_results <- rbindlist(
    all_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  setorder(
    final_results,
    chr,
    start
  )
  
  fwrite(
    final_results,
    file.path(
      model_dir,
      "all_chromosomes_REML_results.csv"
    )
  )
  
  saveRDS(
    all_models,
    file.path(
      model_dir,
      "all_chromosomes_REML_models.rds"
    )
  )
  
  cat("\n====================================\n")
  cat("ALL MODELS COMPLETED\n")
  cat("====================================\n")
  
  cat(
    "Total fitted segments:",
    nrow(final_results),
    "\n"
  )
  
  cat(
    "Successful convergence:",
    sum(final_results$status == "OK", na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Mean true rho:",
    mean(final_results$true_rho, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Mean estimated rho:",
    mean(final_results$estimated_rho, na.rm = TRUE),
    "\n"
  )
  
  cat(
    "Mean estimated phi:",
    mean(final_results$estimated_dispersion_phi, na.rm = TRUE),
    "\n"
  )
  
} else {
  
  stop("No models were successfully fitted.")
  
}


# Read CSV containing the segments to keep
keep <- fread("/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/segments_variable_T1_or_T2.csv")

# Keep only matching coordinates
results_filtered <- merge(
  final_results,
  unique(keep[, .(chr, start, end)]),
  by = c("chr", "start", "end"),
  all = FALSE
)

# Remove potential duplicate rows
setorder(results_filtered, chr, start, end)

# Save
fwrite(
  results_filtered,
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/vary_rho_values_long_tail_beta_mixture_with_large_values/model_results/filtered_model_results.csv"
)

cat("Original rows:", nrow(final_results), "\n")
cat("Filtered rows:", nrow(results_filtered), "\n")
