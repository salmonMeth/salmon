#we simulate data and also generate new samples if we want to


library(data.table)
library(VGAM)

set.seed(123)

timepoint <- "T1"

#change accordingly
extra_samples <- 10L

data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_cs_cov_3_samples_removed"


path_ids <-"/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"

bed_dir <-"/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/beds"

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/with_extra_10_samples_per_family"

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

####################################################
####################################################
#get ma&pa random effects

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
#create the sd values we use for the parental effects distribution
#we later will pick one value per parent per segment from these distributions
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

sim_parameters[
  ,
  `:=`(
    true_var_ma = true_sd_ma^2,
    true_var_pa = true_sd_pa^2
  )
]


# get mean segment methylation

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

all_segments[
  ,
  chr := sprintf(
    "NC_%06d.1",
    59441L + as.integer(chr)
  )
]

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

segment_info <- merge(
  segments,
  all_segments,
  by = c("chr", "start", "end"),
  all.x = TRUE,
  sort = FALSE
)


segment_info[
  ,
  baseline := qlogis(
    pmin(
      pmax(mean_methylation, 1e-6),
      1 - 1e-6
    )
  )
]

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
  stop("Missing segment baselines after merging BED information.")
}

#generate ma/pa random effects

maternal_re <- CJ(
  segment = sim_parameters$segment,
  ID_ma = unique(meta$ID_ma),
  unique = TRUE
)

maternal_re <- merge(
  maternal_re,
  sim_parameters[, .(segment, true_sd_ma)],
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
  sim_parameters[, .(segment, true_sd_pa)],
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

# Add metadata
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
      combs
    )
  ],
  by = "segment",
  all.x = TRUE
)
# add ma/pa random effects

sim <- merge(
  sim,
  maternal_re[, .(segment, ID_ma, u_ma)],
  by = c("segment", "ID_ma"),
  all.x = TRUE
)

sim <- merge(
  sim,
  paternal_re[, .(segment, ID_pa, u_pa)],
  by = c("segment", "ID_pa"),
  all.x = TRUE
)

#get the logit
sim[
  ,
  eta := baseline +
    weight_beta * Weight_sc +
    u_ma +
    u_pa
]

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
#simulate using beta binomial

rho <- 0.02

sim[, meth_reads := NA_integer_]

idx <- which(
  sim$coverage > 0 &
    is.finite(sim$methylation_probability) &
    !is.na(sim$methylation_probability)
)

sim[
  idx,
  meth_reads := as.integer(
    VGAM::rbetabinom(
      n = length(idx),
      size = coverage[idx],
      prob = methylation_probability[idx],
      rho = rho
    )
  )
]

sim[
  ,
  unmeth_reads := coverage - meth_reads
]

sim[, chr_name := as.character(chr)]

setorder(sim, chr_name, start, end, sample)



saveRDS(
  sim,
  output_file
)

#create cs/cov files we use for modeling

for (cc in unique(sim$chr_name)) {
  
  d <- sim[chr_name == cc]

  cs <- dcast(
    d,
    chr_name + start + end ~ sample,
    value.var = "meth_reads"
  )
  
  setnames(cs, "chr_name", "chr")
  setorder(cs, start, end)
  
  fwrite(
    cs,
    file.path(
      out_dir,
      paste0(cc, "_segment_cs.csv")
    )
  )

  cov_dt <- dcast(
    d,
    chr_name + start + end ~ sample,
    value.var = "coverage"
  )
  
  setnames(cov_dt, "chr_name", "chr")
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
    u_ma,
    u_pa,
    eta,
    methylation_probability
  )
]

fwrite(
  true_sample_parameters,
  file.path(out_dir, "true_sample_parameters.csv")
)