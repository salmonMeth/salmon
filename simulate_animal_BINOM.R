library(data.table)
library(brms)
library(posterior)

data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/animal_model"
# Load saved data
dat_T1 <- readRDS(
  file.path(data_dir, "keep_segments_model_data_T1.rds")
)
A_T1 <- readRDS(
  file.path(data_dir, "A_T1.rds")
)
set.seed(123)
# Select 10 segments
  # Select the first 900 segments in their existing order
segment_ids <- unique(dat_T1$segment_id)[63:100]

set.seed(123)

VA_values <- runif(
  length(segment_ids),
  min = 0,
  max = 2
)
beta1_values <- rnorm(
  length(segment_ids),0,.2)

phi_values <- runif(length(segment_ids),1,40)


simulation_parameters <- data.table(
  segment_id = segment_ids,
  VA_true = VA_values,
  phi_true = beta1_values,
  beta1_true = -phi_values
)

simulation_parameters

simulate_and_fit <- function(seg_id, dat, A,
                             VA_true,
                             phi_true = 40,
                             beta1_true = -1.14,
                             seed = 123) {
  
  set.seed(seed)
  
  # 1. Extract original segment data
  d <- as.data.frame(dat[dat$segment_id == seg_id, ])
  
  ids <- sort(unique(as.character(d$animal)))
  
  # Match relationship matrix to observed animals
  A_seg <- A[ids, ids, drop = FALSE]
  
  d$animal <- factor(
    as.character(d$animal),
    levels = ids
  )
  
  # 2. Original coverage-weighted mean methylation
  mean_observed <- sum(d$meth) / sum(d$coverage)
  
  # 3. Generate correlated animal effects
  L <- t(chol(A_seg))
  
  a_true <- sqrt(VA_true) *
    drop(L %*% rnorm(length(ids)))
  
  names(a_true) <- ids
  
  true_effects <- data.table(
    animal = ids,
    a_true = a_true
  )
  
  # 4. Calibrate intercept to preserve expected mean methylation
  beta0_true <- uniroot(
    function(b0) {
      
      mu <- plogis(
        b0 +
          beta1_true * d$Weight_sc +
          a_true[as.character(d$animal)]
      )
      
      sum(d$coverage * mu) / sum(d$coverage) -
        mean_observed
    },
    interval = c(-20, 20)
  )$root
  
  # 5. Generate methylation probabilities
  eta <- beta0_true +
    beta1_true * d$Weight_sc +
    a_true[as.character(d$animal)]
  
  mu <- plogis(eta)

  d$meth <- rbinom(
    n = nrow(d),
    size = d$coverage,
    prob = mu
  )
  
  d$unmeth <- d$coverage - d$meth
  
  # 7. Fit animal model
  fit <- brm(
    meth | trials(coverage) ~
      1 + Weight_sc +
      (1 | gr(animal, cov = A)),
    data = d,
    data2 = list(A = A_seg),
    family = stats::binomial(link = "logit"),
    chains = 4,
    cores = 2,
    iter = 6000,
    warmup = 3000,
    thin=5,
    seed = seed + 1,
    control = list(
      adapt_delta = 0.99,
      max_treedepth = 15
    ),
    refresh = 0
  )
  
  # 8. Return results
  list(
    segment_id = seg_id,
    fit = fit,
    true_parameters = data.table(
      VA_true = VA_true,
      phi_true = phi_true,
      beta0_true = beta0_true,
      beta1_true = beta1_true,
      mean_methylation = mean_observed
    ),
    true_effects = true_effects,
    simulated_data = d
  )
}

results <- vector("list", length(segment_ids))

for (i in seq_along(segment_ids)) {
  
  cat("Fitting segment", i, "of", length(segment_ids), "\n")
  
  results[[i]] <- simulate_and_fit(
    seg_id = segment_ids[i],
    dat = dat_T1,
    A = A_T1,
    VA_true = simulation_parameters$VA_true[i],
    phi_true = 40,
    beta1_true = -1.14,
    seed = 1000 + i
  )
}

saveRDS(
  results,
  file.path(data_dir, "simulation_results_63-100_BINOMIAL_MODEL.rds")
)