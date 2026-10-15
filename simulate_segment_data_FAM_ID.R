#this simulates data with fam id instead of ma and pa ids


library(data.table)
library(VGAM)


cov <- readRDS(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/coverage.rds"
)

cov <- as.data.table(cov)
n_segments <- nrow(cov)

# Check coverage matrix dimensions and values
print(dim(cov))
print(head(cov[, 1:min(5, ncol(cov)), with = FALSE]))



weight_betas <- rnorm(n_segments, 0, 0.05)

combs <- sample(
  c("none", "family"),
  size = n_segments,
  replace = TRUE,
  prob = c(0.25, 0.75)
)

sim_parameters <- data.table(
  segment = seq_len(n_segments),
  combs = combs
)

sim_parameters[, true_fam_sd := 0]

sim_parameters[
  combs == "none",
  true_fam_sd := 1e-9
]

sim_parameters[
  combs == "family",
  true_fam_sd := runif(.N, 0.1, 0.5)
]

sim_parameters[, true_fam_var := true_fam_sd^2]

sim_parameters[, weight_beta := weight_betas]



path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"

meta <- fread(path_ids)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]

meta[, sample := as.character(sample)]

meta[, Weight_sc := as.numeric(scale(Weight.g))]

# Check sample identifiers
stopifnot(!anyDuplicated(meta$sample))


family_ids <- unique(meta$ID_family[!is.na(meta$ID_family)])

fam_re <- CJ(
  segment = seq_len(n_segments),
  ID_family = family_ids
)

fam_re <- merge(
  fam_re,
  sim_parameters[, .(segment, true_fam_sd)],
  by = "segment",
  all.x = TRUE
)

set.seed(123)

fam_re[, u_fam := rnorm(
  .N,
  mean = 0,
  sd = true_fam_sd
)]



bed_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/beds"

bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

stopifnot(length(bed_files) > 0)

all_segments <- rbindlist(
  lapply(bed_files, function(file) {
    
    x <- fread(file, header = FALSE)
    
    data.table(
      chr = x[[1]],
      start = x[[2]],
      end = x[[3]],
      mean_methylation = as.numeric(x[[5]])
    )
  })
)

# IMPORTANT: This assumes the BED segment order matches
# the row order in coverage.rds.
if (nrow(all_segments) != n_segments) {
  stop(
    "The number of BED segments does not match the number ",
    "of rows in coverage.rds. Verify segment correspondence."
  )
}

all_segments[, segment := seq_len(.N)]

# Avoid infinite baseline values if means equal exactly 0 or 1
eps <- 1e-6

all_segments[
  ,
  baseline := qlogis(
    pmin(pmax(mean_methylation, eps), 1 - eps)
  )
]

sim_parameters <- merge(
  sim_parameters,
  all_segments[, .(
    segment, chr, start, end,
    mean_methylation, baseline
  )],
  by = "segment",
  sort = FALSE
)


sim <- CJ(
  segment = seq_len(n_segments),
  sample = meta$sample
)

sim <- merge(
  sim,
  meta,
  by = "sample",
  all.x = TRUE
)

sim <- merge(
  sim,
  sim_parameters[, .(
    segment,
    chr,
    start,
    end,
    mean_methylation,
    baseline,
    weight_beta,
    true_fam_sd,
    true_fam_var,
    combs
  )],
  by = "segment",
  all.x = TRUE
)

# Merge family effects using the SAME name, u_fam
sim <- merge(
  sim,
  fam_re[, .(segment, ID_family, u_fam)],
  by = c("segment", "ID_family"),
  all.x = TRUE
)

# Individuals without a family ID have no family random effect
sim[is.na(u_fam), u_fam := 0]

# get methylation probs

sim[
  ,
  eta := baseline +
    weight_beta * Weight_sc +
    u_fam
]

sim[, methylation_probability := plogis(eta)]


coverage_dt <- copy(as.data.table(cov))

coverage_dt[, segment := seq_len(.N)]

metadata_cols <- c(
  "chr", "start", "end",
  "segment_id", "distribution", "segment"
)

sample_cols <- setdiff(names(coverage_dt), metadata_cols)


coverage_long <- melt(
  coverage_dt,
  id.vars = "segment",
  measure.vars = sample_cols,
  variable.name = "sample",
  value.name = "coverage"
)

coverage_long[, sample := as.character(sample)]

# make sure coverage is integer
coverage_long[
  ,
  coverage := as.numeric(coverage)
]

coverage_long[
  !is.na(coverage) & is.finite(coverage),
  coverage := as.integer(pmax(0, round(coverage)))
]

sim[, sample := as.character(sample)]

sim <- merge(
  sim,
  coverage_long,
  by = c("segment", "sample"),
  all.x = TRUE
)

sim[, cov := as.integer(coverage)]


rho <- 0.02

sim[
  !is.na(cov) & cov == 0L,
  `:=`(
    meth_reads = 0L,
    unmeth_reads = 0L
  )
]

idx <- which(
  !is.na(sim$cov) &
    sim$cov > 0L &
    is.finite(sim$methylation_probability) &
    !is.na(sim$methylation_probability) &
    sim$methylation_probability >= 0 &
    sim$methylation_probability <= 1
)

set.seed(456)

sim[
  idx,
  meth_reads := as.integer(
    VGAM::rbetabinom(
      n = length(idx),
      size = cov,
      prob = methylation_probability,
      rho = rho
    )
  )
]

sim[
  !is.na(meth_reads) & !is.na(cov),
  unmeth_reads := as.integer(cov - meth_reads)
]
