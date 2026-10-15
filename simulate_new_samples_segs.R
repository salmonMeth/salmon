#generate new samples
#we add them to the existing fams, and get their coverage values from pre-existing coverage values +noise

library(dplyr)
library(tidyr)
library(data.table)
library(VGAM)

path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]

meta[, sample := as.character(sample)]

meta[, Weight_sc := as.numeric(scale(Weight.g))]


cov <- readRDS(
  "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/coverage.rds"
)

cov <- as.data.table(cov)
n_segments <- nrow(cov)

family_pairs <- meta %>%
  distinct(ID_ma, ID_pa)

nrow(family_pairs)

#we can adjust this if we want to add not 2 but more samples per family
n_new_per_family <- 2  

new_meta <- family_pairs %>%
  uncount(n_new_per_family, .id = "offspring_number") %>%
  mutate(
    sample = paste0(
      "SIM", sprintf("%03d", row_number()), "_T1"
    ),
    TimePoint = "T1"
  )

# checks
cat("Number of families:", nrow(family_pairs), "\n")
cat("New individuals per family:", n_new_per_family, "\n")
cat("Total new individuals:", nrow(new_meta), "\n")


#gotta generate weight values for the new samples
#we just resample+add some noise
weight_sd <- sd(meta$Weight.g, na.rm = TRUE)
noise_sd <- 0.05 * weight_sd

source_idx <- sample(
  which(!is.na(meta$Weight.g)),
  size = nrow(new_meta),
  replace = TRUE
)

new_weights <- meta$Weight.g[source_idx] +
  rnorm(nrow(new_meta), mean = 0, sd = noise_sd)

data.table::set(new_meta, j = "Weight.g", value = new_weights)

data.table::set(
  new_meta,
  j = "Weight_sc",
  value = as.numeric(scale(new_meta$Weight.g))
)


fwrite(
  new_meta,
  file = "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/new_samples_covs/new_meta_simulated_T1.csv"
)
#get cov

coord_cols <- c("chr", "start", "end", "segment_id")

sample_cols <- setdiff(names(cov), coord_cols)

existing_samples <- intersect(meta$sample, sample_cols)

stopifnot(length(existing_samples) > 0)

new_cov <- cov[, ..coord_cols]
existing_samples <- intersect(
  existing_samples,
  names(cov)
)

set.seed(897)
#resample from the existing coverage distribution+noise
for (i in seq_len(nrow(new_meta))) {
  
  new_id <- new_meta$sample[i]
  
  source_sample <- sample(existing_samples, size = 1)
  source_cov <- as.numeric(cov[[source_sample]])
  
  # add noise
  noise <- rnorm(
    length(source_cov),
    mean = 0,
    sd = 0.05
  )
  
  sampled_cov <- round(source_cov * (1 + noise))
  
  # no negatives
  sampled_cov[!is.na(sampled_cov)] <- pmax(
    0,
    sampled_cov[!is.na(sampled_cov)]
  )
  
  new_cov[[new_id]] <- as.integer(sampled_cov)
}


saveRDS(new_cov, "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/new_samples_covs/coverage.rds" )



#########start here

cov=new_cov
n_segments=nrow(cov)
weight_betas=rnorm(n_segments,0,.05)



#what combinations of heritability we will have in the model
combs <- sample(
  c(
    "none",
    "some"
  ),
  size = n_segments,
  replace = TRUE,
  prob = c(0.25, 0.75)
)

sim_parameters <- data.table(
  segment = seq_len(n_segments),
  combs = combs
)

#init
sim_parameters[, `:=`(
  true_sd_ma = 0,
  true_sd_pa = 0
)]

#no heritability
sim_parameters[
  combs == "none",
  `:=`(
    true_sd_ma = 1e-9,
    true_sd_pa = 1e-9
  )
]

#mostly maternal
sim_parameters[
  combs == "some",
  `:=`(
    true_sd_ma = runif(.N, 0.05, 0.30),
    true_sd_pa = runif(.N, 0.05, 0.3)
  )
]

# get var
sim_parameters[, `:=`(
  true_var_ma = true_sd_ma^2,
  true_var_pa = true_sd_pa^2
)]


parent_info <- unique(
  new_meta[
    grepl("_T1$", sample),
    .(
      sample,
      ID_ma,
      ID_pa
    )
  ]
)



#need the mean segment methylation values
bed_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/beds"

bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

all_segments <- rbindlist(
  lapply(bed_files, function(file) {
    
    x <- fread(
      file,
      header = FALSE
    )
    
    data.table(
      chr = x[[1]],
      start = x[[2]],
      end = x[[3]],
      mean_methylation = x[[5]]
    )
  })
)
all_segments[, segment := seq_len(.N)]

sim_parameters <- merge(
  sim_parameters,
  all_segments[, .(
    segment,
    chr,
    start,
    end,
    mean_methylation
  )],
  by = "segment",
  sort = FALSE
)

sim_parameters[, baseline := qlogis(mean_methylation)]


maternal_re <- CJ(
  segment = seq_len(n_segments),
  ID_ma = unique(parent_info$ID_ma)
)

maternal_re <- merge(
  maternal_re,
  sim_parameters[, .(
    segment,
    true_sd_ma
  )],
  by = "segment",
  all.x = TRUE
)


maternal_re[, u_ma := rnorm(
  .N,
  mean = 0,
  sd = true_sd_ma
)]


paternal_re <- CJ(
  segment = seq_len(n_segments),
  ID_pa = unique(parent_info$ID_pa)
)

paternal_re <- merge(
  paternal_re,
  sim_parameters[, .(
    segment,
    true_sd_pa
  )],
  by = "segment",
  all.x = TRUE
)

paternal_re[, u_pa := rnorm(
  .N,
  mean = 0,
  sd = true_sd_pa
)]


####

sim_parameters[
  ,
  weight_beta := weight_betas
]
# Create all segment-sample combinations
sim <- CJ(
  segment = sim_parameters$segment,
  sample = new_meta$sample
)

# Add sample metadata first (including ID_ma and ID_pa)
sim <- merge(
  sim,
  new_meta,
  by = "sample",
  all.x = TRUE
)

# Now add maternal random effects
sim <- merge(
  sim,
  maternal_re[, .(segment, ID_ma, u_ma)],
  by = c("segment", "ID_ma"),
  all.x = TRUE
)

# Add paternal random effects
sim <- merge(
  sim,
  paternal_re[, .(segment, ID_pa, u_pa)],
  by = c("segment", "ID_pa"),
  all.x = TRUE
)
##################

# Add segment-level parameters
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
      true_var_pa
    )
  ],
  by = "segment",
  all.x = TRUE
)

sim[
  ,
  eta :=
    baseline +
    weight_beta * Weight_sc +
    u_ma +
    u_pa
]


sim[
  ,
  methylation_probability := plogis(eta)
]

coverage_dt <- as.data.table(cov)

coverage_dt[
  ,
  segment := seq_len(.N)
]

coverage_long <- melt(
  coverage_dt,
  id.vars = "segment",
  variable.name = "sample",
  value.name = "coverage"
)

coverage_long[
  ,
  sample := as.character(sample)
]

sim[
  ,
  sample := as.character(sample)
]

# Add coverage
sim <- merge(
  sim,
  coverage_long,
  by = c("segment", "sample"),
  all.x = TRUE
)


## we assume the data comes from a beta binomial distribution with parameterrho=0.02
sim[, coverage := as.numeric(coverage)]

rho <- 0.02
sim[, meth_reads := 0L]

idx <- which(
  sim$cov > 0 &
    is.finite(sim$methylation_probability) &
    !is.na(sim$methylation_probability)
)

sim[
  idx,
  meth_reads := VGAM::rbetabinom(
    n = length(idx),
    size = sim$cov[idx],
    prob = sim$methylation_probability[idx],
    rho = rho
  )
]

sim[, unmeth_reads := sim$cov - meth_reads]

saveRDS(sim,"/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/new_samples_covs/sim_params.rds")


#save the data in the format we used to fit the previous models
library(data.table)

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/new_samples_covs/cs_cov"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

sim[
  ,
  chr_name := sprintf(
    "NC_%06d.1",
    59441L + as.integer(chr)
  )
]

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

