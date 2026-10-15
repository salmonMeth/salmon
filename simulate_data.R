library(data.table)

set.seed(123)

# Change this to "T2" later
timepoint <- "T1"

data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_cs_cov_3_samples_removed"

output_file <- paste0(
  "/scratch/project_2010912/ezel/salmon/data_files/segmentation/",
  "simulated_", timepoint, "_coverage_empirical.rds"
)
chromosomes <- sprintf(
  "NC_%06d.1",
  59442:59470
)

coverage_files <- file.path(
  data_dir,
  paste0(chromosomes, "_segment_coverage.csv")
)

length(coverage_files)
coverage <- rbindlist(
  lapply(seq_along(coverage_files), function(i) {
    #remove the T2 samples
    x <- fread(coverage_files[i])
    sample_cols <- grep(
      paste0("_", timepoint, "$"),
      names(x),
      value = TRUE
    )
    
    x <- x[
      ,
      c("chr", "start", "end", sample_cols),
      with = FALSE
    ]
    
    x
    
  }),
  use.names = TRUE
  
)

sample_cols <- grep(
  paste0("_", timepoint, "$"),
  names(coverage),
  value = TRUE
)

length(sample_cols)
coverage_long <- melt(
  coverage,
  id.vars = c("chr", "start", "end"),
  measure.vars = sample_cols,
  variable.name = "sample",
  value.name = "coverage"
)
setorder(
  segments,
  chr,
  start,
  end
)

segments[, segment_id := .I]

nrow(segments)
n_segments <- nrow(segments)
n_samples <- length(sample_cols)

simulated_coverage <- matrix(
  NA_real_,
  nrow = n_segments,
  ncol = n_samples
)

colnames(simulated_coverage) <- sample_cols


for (i in seq_len(n_segments)) {
  
  #get the actual coverage values
  observed <- coverage_long[
    chr == segments$chr[i] &
      start == segments$start[i] &
      end == segments$end[i],
    coverage
  ]
  
  
  observed <- observed[
    is.finite(observed) &
      observed > 0
  ]
  
  if (length(observed) > 0) {
    
    # resample from the observed values
    sampled <- sample(
      observed,
      size = n_samples,
      replace = TRUE
    )
    
    #add noise
    noise <- rnorm(
      n_samples,
      mean = 0,
      sd = 0.05
    )
    #noise is proportional to the simulated coverage
    simulated_coverage[i, ] <-
      sampled * (1 + noise)
    
    simulated_coverage[i, ] <-
      pmax(
        simulated_coverage[i, ],
        1
      )
  }
}

simulated_coverage <- as.data.table(
  simulated_coverage
)

simulated_segments <- cbind(
  segments,
  simulated_coverage
)

############################
############################
#simulate methylation
#meth ~weight*beta_w+maternal+paternal
# for every segment we need one maternal and one paternal distribution
#lets assume 25% have 0 heritability, the remanining will come from various distributions, some just maternal, some just paternal, some both

cov=readRDS("/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/coverage.rds")
n_segments=nrow(cov)
weight_betas=rnorm(n_segments,0,.05)



#what combinations of heritability we will have in the model
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
  combs == "maternal_dominant",
  `:=`(
    true_sd_ma = runif(.N, 0.12, 0.30),
    true_sd_pa = runif(.N, 0.00, 0.05)
  )
]
#mostly paternal
sim_parameters[
  combs == "paternal_dominant",
  `:=`(
    true_sd_ma = runif(.N, 0.00, 0.05),
    true_sd_pa = runif(.N, 0.12, 0.30)
  )
]

#both
sim_parameters[
  combs == "both",
  `:=`(
    true_sd_ma = runif(.N, 0.10, 0.30),
    true_sd_pa = runif(.N, 0.10, 0.30)
  )
]

# get var
sim_parameters[, `:=`(
  true_var_ma = true_sd_ma^2,
  true_var_pa = true_sd_pa^2
)]

#need metadata 
path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]



parent_info <- unique(
  meta[
    grepl("_T1$", sample),
    .(
      sample,
      ID_ma,
      ID_pa
    )
  ]
)


parent_info

#now we need to generate one parent value per parent per segment using the sd values we generated

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

####now we generate one methylation value per segmet per sample
meta[
  ,
  Weight_sc := as.numeric(scale(Weight.g))
]

sim_parameters[
  ,
  weight_beta := weight_betas
]


sim <- CJ(
  segment = sim_parameters$segment,
  sample = meta$sample
)

# add the metadata
sim <- merge(
  sim,
  meta,
  by = "sample",
  all.x = TRUE
)

#add segment info
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

###gotta make sure ppl with the same parent get the same parental effects

sim <- merge(
  sim,
  maternal_re[
    ,
    .(
      segment,
      ID_ma,
      u_ma
    )
  ],
  by = c("segment", "ID_ma"),
  all.x = TRUE
)

sim <- merge(
  sim,
  paternal_re[
    ,
    .(
      segment,
      ID_pa,
      u_pa
    )
  ],
  by = c("segment", "ID_pa"),
  all.x = TRUE
)

# now we construct the data
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

###add the coverage


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

#save the data in the format we used to fit the previous models
library(data.table)

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_data/cs_cov"

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
