meth=read.csv("/scratch/project_2010912/ezel/salmon/data_files/beta_binom_models_for_actual_data_split_segments/check_for_dharma/meth.csv")
unmeth=read.csv("/scratch/project_2010912/ezel/salmon/data_files/beta_binom_models_for_actual_data_split_segments/check_for_dharma/unmeth.csv")
head(meth)
head(unmeth)
# Identify T1 sample columns
t1_cols <- grep("^X?[0-9]+_T1$", names(meth), value = TRUE)
# Keep segment identifiers and T1 samples
id_cols <- c("chr", "start", "end")
meth <- meth[, c(id_cols, t1_cols), with = FALSE]
meth <- meth[, c(id_cols, t1_cols), drop = FALSE]
unmeth <- unmeth[, c(id_cols, t1_cols), drop = FALSE]
identical(names(meth), names(unmeth))
grep("_T2$", names(meth), value = TRUE)
dim(meth)
dim(unmeth)
# Coverage matrix
cov <- meth
for (s in t1_cols) {
  cov[[s]] <- meth[[s]] + unmeth[[s]]
}
library(data.table)
cov_dt <- as.data.table(cov)
coverage_summary <- cov_dt[, {
  x <- as.numeric(unlist(.SD))
  .(
    n_samples = sum(is.finite(x)),
    mean_coverage = mean(x, na.rm = TRUE),
    median_coverage = median(x, na.rm = TRUE),
    sd_coverage = sd(x, na.rm = TRUE),
    min_coverage = min(x, na.rm = TRUE),
    max_coverage = max(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE)
  )
}, by = .(chr, start, end), .SDcols = t1_cols]
head(coverage_summary)
# Example: resample coverage for one segment
set.seed(123)


meth_mat <- as.matrix(meth[, t1_cols])
unmeth_mat <- as.matrix(unmeth[, t1_cols])

cov_mat <- meth_mat + unmeth_mat

meth_pct_mat <- meth_mat / cov_mat
meth_pct_mat[!is.finite(meth_pct_mat)] <- NA_real_

segment_targets <- data.table(
  chr = meth$chr,
  start = meth$start,
  end = meth$end,
  mean_methylation = rowMeans(meth_pct_mat, na.rm = TRUE),
  var_methylation = apply(meth_pct_mat, 1, var, na.rm = TRUE),
  sd_methylation = apply(meth_pct_mat, 1, sd, na.rm = TRUE)
)

sample_map <- data.table(
  sample = sub("^X", "", t1_cols)
)

sample_map <- merge(
  sample_map,
  meta[, .(
    sample = paste(Running.ID, TimePoint, sep = "_"),
    ID_ma,
    ID_pa
  )],
  by = "sample",
  all.x = TRUE,
  sort = FALSE
)

resample_segment_coverage <- function(coverage, n = length(coverage)) {
  
  coverage <- coverage[
    is.finite(coverage) & coverage > 0
  ]
  
  sample(
    coverage,
    size = n,
    replace = TRUE
  )
}

set.seed(123)

sim_coverage <- lapply(seq_len(nrow(cov_mat)), function(j) {
  
  resample_segment_coverage(
    cov_mat[j, ],
    n = ncol(cov_mat)
  )
  
})

names(sim_coverage) <- paste(
  segment_targets$chr,
  segment_targets$start,
  segment_targets$end,
  sep = "_"
)

library(data.table)

# Count valid coverage values per segment
n_valid_coverage <- rowSums(
  is.finite(cov_mat) & cov_mat > 0
)

# Build diagnostic table
segment_validity <- data.table(
  segment_index = seq_len(nrow(segment_targets)),
  chr = segment_targets$chr,
  start = segment_targets$start,
  end = segment_targets$end,
  mean_methylation = segment_targets$mean_methylation,
  var_methylation = segment_targets$var_methylation,
  n_valid_methylation = rowSums(is.finite(meth_pct_mat)),
  n_valid_coverage = n_valid_coverage
)

# Individual validity flags
segment_validity[, `:=`(
  invalid_mean = !is.finite(mean_methylation),
  mean_zero = is.finite(mean_methylation) & mean_methylation <= 0,
  mean_one_or_more = is.finite(mean_methylation) & mean_methylation >= 1,
  invalid_variance = !is.finite(var_methylation),
  insufficient_methylation = n_valid_methylation < 2,
  incomplete_coverage = n_valid_coverage < ncol(cov_mat),
  no_valid_coverage = n_valid_coverage == 0
)]

# Overall validity: same criteria as the original selection
segment_validity[, valid := (
  is.finite(mean_methylation) &
    is.finite(var_methylation) &
    mean_methylation > 0 &
    mean_methylation < 1 &
    n_valid_coverage == ncol(cov_mat)
)]

# Count reasons (these overlap if a segment has multiple issues)
reason_summary <- data.table(
  reason = c(
    "Non-finite mean methylation",
    "Mean methylation = 0",
    "Mean methylation >= 1",
    "Non-finite variance",
    "Fewer than 2 valid methylation values",
    "Incomplete coverage across samples",
    "No valid coverage"
  ),
  n_segments = c(
    sum(segment_validity$invalid_mean),
    sum(segment_validity$mean_zero),
    sum(segment_validity$mean_one_or_more),
    sum(segment_validity$invalid_variance),
    sum(segment_validity$insufficient_methylation),
    sum(segment_validity$incomplete_coverage),
    sum(segment_validity$no_valid_coverage)
  )
)

min_samples <- 42

segment_validity[, `:=`(
  enough_coverage = n_valid_coverage >= min_samples,
  valid = (
    is.finite(mean_methylation) &
      is.finite(var_methylation) &
      mean_methylation > 0 &
      mean_methylation < 1 &
      n_valid_coverage >= min_samples
  )
)]

model_dir <- paste0(
  "/scratch/project_2010912/ezel/salmon/data_files/",
  "beta_binom_models_for_actual_data_split_segments/",
  "check_for_dharma/model_results_T1_T2_separated/models_T1"
)

test_chr <- as.character(test_segment$chr)

model_file <- file.path(
  model_dir,
  paste0("models_", test_chr, "_T1.rds")
)


simulate_segment <- function(
    dat,
    sd_ma,
    sd_pa,
    sd_env,
    phi,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(dat)
  
  # Shared parental effects
  mothers <- unique(dat$ID_ma)
  fathers <- unique(dat$ID_pa)
  
  u_ma <- rnorm(length(mothers), 0, sd_ma)
  u_pa <- rnorm(length(fathers), 0, sd_pa)
  
  maternal_effect <- u_ma[match(dat$ID_ma, mothers)]
  paternal_effect <- u_pa[match(dat$ID_pa, fathers)]
  
  # Independent environmental effects
  env_effect <- rnorm(n, 0, sd_env)
  
  eta <- maternal_effect + paternal_effect + env_effect
  
  # Calibrate intercept to target mean probability
  alpha <- uniroot(
    function(a) mean(plogis(a + eta)) - target_mean,
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
    size = dat$coverage,
    prob = beta_p
  )
  
  simulated_unmeth <- dat$coverage - simulated_meth
  
  simulated_pct <- simulated_meth / dat$coverage
  
  list(
    meth = simulated_meth,
    unmeth = simulated_unmeth,
    coverage = dat$coverage,
    methylation = simulated_pct,
    mean = mean(simulated_pct),
    variance = var(simulated_pct),
    maternal_effect = maternal_effect,
    paternal_effect = paternal_effect,
    environmental_effect = env_effect
  )
}

tolerance <- 0.05

mean_lower <- target_mean * (1 - tolerance)
mean_upper <- target_mean * (1 + tolerance)

var_lower <- target_var * (1 - tolerance)
var_upper <- target_var * (1 + tolerance)

data.table(
  statistic = c("Mean methylation", "Variance"),
  target = c(target_mean, target_var),
  lower = c(mean_lower, var_lower),
  upper = c(mean_upper, var_upper)
)
check_simulation <- function(sim, target_mean, target_var, tolerance = 0.05) {
  
  mean_ok <- abs(sim$mean - target_mean) <=
    target_mean * tolerance
  
  variance_ok <- abs(sim$variance - target_var) <=
    target_var * tolerance
  
  mean_ok && variance_ok
}
test_accepted <- (
  test_sim$mean >= mean_lower &&
    test_sim$mean <= mean_upper &&
    test_sim$variance >= var_lower &&
    test_sim$variance <= var_upper
)

simulate_segment <- function(
    dat,
    sd_ma,
    sd_pa,
    sd_env,
    phi,
    target_mean
) {
  
  n <- nrow(dat)
  
  # Resample segment-specific coverage
  coverage <- dat$coverage[
    sample.int(n, size = n, replace = TRUE)
  ]
  
  # One effect per unique parent
  mothers <- unique(dat$ID_ma)
  fathers <- unique(dat$ID_pa)
  
  u_ma <- rnorm(length(mothers), 0, sd_ma)
  u_pa <- rnorm(length(fathers), 0, sd_pa)
  
  maternal_effect <- u_ma[match(dat$ID_ma, mothers)]
  paternal_effect <- u_pa[match(dat$ID_pa, fathers)]
  
  # Independent environmental effects
  env_effect <- rnorm(n, 0, sd_env)
  
  eta <- maternal_effect + paternal_effect + env_effect
  
  # Calibrate latent mean
  alpha <- uniroot(
    function(a) mean(plogis(a + eta)) - target_mean,
    interval = c(-20, 20)
  )$root
  
  p <- plogis(alpha + eta)
  
  # Beta-binomial read generation
  beta_p <- rbeta(
    n,
    shape1 = p * phi,
    shape2 = (1 - p) * phi
  )
  
  simulated_meth <- rbinom(n, coverage, beta_p)
  simulated_pct <- simulated_meth / coverage
  
  list(
    mean = mean(simulated_pct),
    variance = var(simulated_pct),
    coverage = coverage,
    meth = simulated_meth,
    unmeth = coverage - simulated_meth
  )
}
sd_grid <- CJ(
  sd_ma = seq(0, 0.2, by = 0.02),
  sd_pa = seq(0, 0.2, by = 0.02),
  sd_env = seq(0, 0.2, by = 0.02)
)


###########now we simulate
library(data.table)

min_samples <- 42

model_dir <- paste0(
  "/scratch/project_2010912/ezel/salmon/data_files/",
  "beta_binom_models_for_actual_data_split_segments/",
  "check_for_dharma/model_results_T1_T2_separated/models_T1"
)

segment_inputs <- vector("list", nrow(segment_targets))

for (i in seq_len(nrow(segment_targets))) {
  
  seg <- segment_targets[i]
  
  segment_id <- paste(
    seg$chr, seg$start, seg$end, sep = "_"
  )
  
  # Locate original chromosome model
  model_file <- file.path(
    model_dir,
    paste0("models_", seg$chr, "_T1.rds")
  )
  
  if (!file.exists(model_file)) next
  
  chr_models <- readRDS(model_file)
  fitted_model <- chr_models[[segment_id]]
  
  if (is.null(fitted_model)) next
  
  phi <- as.numeric(sigma(fitted_model))
  
  # Locate segment in methylation matrices
  idx <- which(
    meth$chr == seg$chr &
      meth$start == seg$start &
      meth$end == seg$end
  )
  
  if (length(idx) != 1) next
  
  coverage <- as.numeric(cov_mat[idx, ])
  methylation <- as.numeric(meth_pct_mat[idx, ])
  
  valid <- is.finite(coverage) &
    coverage > 0 &
    is.finite(methylation)
  
  if (sum(valid) < min_samples) next
  
  samples <- t1_cols[valid]
  
  dat <- sample_map[
    match(sub("^X", "", samples), sample)
  ]
  
  if (anyNA(dat$ID_ma) || anyNA(dat$ID_pa)) next
  
  segment_inputs[[i]] <- data.table(
    segment_id = segment_id,
    chr = seg$chr,
    start = seg$start,
    end = seg$end,
    n_samples = sum(valid),
    target_mean = mean(methylation[valid]),
    target_variance = var(methylation[valid]),
    phi = phi,
    data = list(data.table(
      sample = sub("^X", "", samples),
      ID_ma = dat$ID_ma,
      ID_pa = dat$ID_pa,
      coverage = coverage[valid],
      observed_methylation = methylation[valid]
    ))
  )
}

segment_inputs <- rbindlist(
  Filter(Negate(is.null), segment_inputs)
)

cat("Segments prepared:", nrow(segment_inputs), "\n")
cat("Segments excluded:", nrow(segment_targets) - nrow(segment_inputs), "\n")