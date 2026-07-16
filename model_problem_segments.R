library(data.table)

model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated"

seg_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated/problem_segments"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

model_files <- list.files(
  model_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

all_bad_cov <- list()
all_bad_cs <- list()

for(model_file in model_files){
  
  chr_name <- sub(
    "^model_(NC_[0-9\\.]+)\\.csv$",
    "\\1",
    basename(model_file)
  )
  
  res <- fread(model_file)
  
  bad <- res[
    status_T1 != "OK" |
      status_T2 != "OK"
  ]
  
  if(nrow(bad) == 0)
    next
  
  cov_file <- file.path(
    seg_dir,
    paste0(chr_name, "_segment_coverage.csv")
  )
  
  cs_file <- file.path(
    seg_dir,
    paste0(chr_name, "_segment_cs.csv")
  )
  
  cov_dt <- fread(cov_file)
  cs_dt <- fread(cs_file)
  
  bad_cov <- merge(
    bad,
    cov_dt,
    by = c("chr","start","end"),
    all.x = TRUE
  )
  
  bad_cs <- merge(
    bad,
    cs_dt,
    by = c("chr","start","end"),
    all.x = TRUE
  )
  
  all_bad_cov[[chr_name]] <- bad_cov
  all_bad_cs[[chr_name]] <- bad_cs
}

all_bad_cov <- rbindlist(all_bad_cov, fill = TRUE)
all_bad_cs <- rbindlist(all_bad_cs, fill = TRUE)

fwrite(
  all_bad_cov,
  file.path(
    out_dir,
    "all_problem_segments_coverage.csv"
  )
)

fwrite(
  all_bad_cs,
  file.path(
    out_dir,
    "all_problem_segments_cs.csv"
  )
)
############stats regarding the coverage of the problematic segments
library(data.table)

bad_cov <- fread(
  "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated/problem_segments/all_problem_segments_coverage.csv"
)

bad_cs <- fread(
  "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated/problem_segments/all_problem_segments_cs.csv"
)

meta_cols <- c(
  "chr","start","end",
  "estimate_T1","SE_T1","z_T1","pvalue_T1","status_T1",
  "estimate_T2","SE_T2","z_T2","pvalue_T2","status_T2"
)

cov_cols <- setdiff(colnames(bad_cov), meta_cols)

segment_cov_summary <- bad_cov[
  ,
  {
    x <- as.numeric(.SD)
    
    .(
      mean_cov = mean(x, na.rm = TRUE),
      median_cov = median(x, na.rm = TRUE),
      sd_cov = sd(x, na.rm = TRUE),
      min_cov = min(x, na.rm = TRUE),
      q25_cov = quantile(x, 0.25, na.rm = TRUE),
      q75_cov = quantile(x, 0.75, na.rm = TRUE),
      max_cov = max(x, na.rm = TRUE),
      n_zero_cov = sum(x == 0, na.rm = TRUE),
      prop_zero_cov = mean(x == 0, na.rm = TRUE)
    )
  },
  by = .(
    chr,
    start,
    end,
    status_T1,
    status_T2
  ),
  .SDcols = cov_cols
]

cs_cols <- setdiff(colnames(bad_cs), meta_cols)

segment_cs_summary <- bad_cs[
  ,
  {
    x <- as.numeric(.SD)
    
    .(
      mean_cs = mean(x, na.rm = TRUE),
      median_cs = median(x, na.rm = TRUE),
      sd_cs = sd(x, na.rm = TRUE),
      min_cs = min(x, na.rm = TRUE),
      q25_cs = quantile(x, 0.25, na.rm = TRUE),
      q75_cs = quantile(x, 0.75, na.rm = TRUE),
      max_cs = max(x, na.rm = TRUE),
      n_zero_cs = sum(x == 0, na.rm = TRUE),
      prop_zero_cs = mean(x == 0, na.rm = TRUE)
    )
  },
  by = .(
    chr,
    start,
    end,
    status_T1,
    status_T2
  ),
  .SDcols = cs_cols
]

segment_problem_summary <- merge(
  segment_cov_summary,
  segment_cs_summary,
  by = c(
    "chr",
    "start",
    "end",
    "status_T1",
    "status_T2"
  )
)

fwrite(
  segment_problem_summary,
  "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated/problem_segments/problem_segment_summary.csv"
)
######################################################
######################################################
#save some of the error throwing model summaries

library(data.table)
library(glmmTMB)

model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated"
seg_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"
out_dir <- paste0(
  model_dir,
  "/problem_segments"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


model_files <- list.files(
  model_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

all_results <- rbindlist(
  lapply(model_files, fread),
  fill = TRUE
)

bad <- all_results[
  status_T1 != "OK" |
    status_T2 != "OK"
]

both_bad <- bad[
  status_T1 != "OK" &
    status_T2 != "OK"
][1]

t1_bad <- bad[
  status_T1 != "OK" &
    status_T2 == "OK"
][1]

t2_bad <- bad[
  status_T1 == "OK" &
    status_T2 != "OK"
][1]

both_ok <- all_results[
  status_T1 == "OK" &
    status_T2 == "OK"
][1]

path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[  ,sample := paste(Running.ID, TimePoint,sep = "_")]

#####################################################
#helpers

load_chr_dat <- function(chr_name){
  
  meth <- fread(
    file.path(
      seg_dir,
      paste0(chr_name, "_segment_cs.csv")
    )
  )
  
  cov <- fread(
    file.path(
      seg_dir,
      paste0(chr_name, "_segment_coverage.csv")
    )
  )
  
  meth_long <- melt(
    meth,
    id.vars = c("chr","start","end"),
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = c("chr","start","end"),
    variable.name = "sample",
    value.name = "coverage"
  )
  
  dat <- merge(
    meth_long,
    cov_long,
    by = c(
      "chr",
      "start",
      "end",
      "sample"
    )
  )
  
  dat <- merge(
    dat,
    meta,
    by = "sample"
  )
  
  dat[
    ,
    unmeth_reads := coverage - meth_reads
  ]
  
  dat
}

inspect_segment <- function(seg, model_fun) {
  
  if (
    is.null(seg) ||
    nrow(seg) == 0
  ) {
    cat("No segment found\n")
    return(NULL)
  }
  
  cat("\n=====================================\n")
  
  print(seg)
  
  dat <- load_chr_dat(seg$chr)
  
  d <- dat[
    chr == seg$chr &
      start == seg$start &
      end == seg$end
  ]
  
  cat(
    "\nSamples:",
    nrow(d),
    "\n"
  )
  
  cat("\nCoverage summary:\n")
  print(summary(d$coverage))
  
  cat(
    "\nZero coverage:",
    sum(d$coverage == 0),
    "\n"
  )
  
  ###################################################
  # T1
  
  d_t1 <- d[
    grepl("_T1$", sample)
  ]
  cat("\n========== T1 ==========\n")
  
  try(
    print(
      summary(
        model_fun(d_t1)
      )
    ),
    silent = FALSE
  )
  
  ###################################################
  # T2
  
  d_t2 <- d[
    grepl("_T2$", sample)
  ]
  cat("\n========== T2 ==========\n")
  try(
    print(
      summary(
        model_fun(d_t2)
      )
    ),
    silent = FALSE
  )
  invisible(NULL)
}
#########################
##########################
#define a few models
bb_model <- function(data) {
  glmmTMB(
    cbind(meth_reads, unmeth_reads) ~
      Weight.g +
      (1 | ID_ma) +
      (1 | ID_pa),
    data = data,
    family = betabinomial(link = "logit")
  )
}

bin_model <- function(data) {
  glmmTMB(
    cbind(meth_reads, unmeth_reads) ~
      Weight.g +
      (1 | ID_ma) +
      (1 | ID_pa),
    data = data,
    family = binomial(link = "logit")
  )
}
bin_model_weight_scaled <- function(data) {
  data$Weight.g.scaled <- scale(data$Weight.g)
  
  glmmTMB(
    cbind(meth_reads, unmeth_reads) ~
      Weight.g.scaled +
      (1 | ID_ma) +
      (1 | ID_pa),
    data = data,
    family = binomial(link = "logit")
  )
}
#####################################################
####

#save

writeLines(
  capture.output(
    inspect_segment(both_bad)
  ),
  file.path(
    out_dir,
    "both_bad_summary.txt"
  )
)

writeLines(
  capture.output(
    inspect_segment(t1_bad)
  ),
  file.path(
    out_dir,
    "t1_bad_summary.txt"
  )
)

writeLines(
  capture.output(
    inspect_segment(t2_bad)
  ),
  file.path(
    out_dir,
    "t2_bad_summary.txt"
  )
)

writeLines(
  capture.output(
    inspect_segment(both_ok)
  ),
  file.path(
    out_dir,
    "both_ok_summary.txt"
  )
)

cat("Done\n")

############################################
############################################
#compare the problematic parts of the two models
library(data.table)

old <- fread(
  "/scratch/project_2010912/ezel/salmon/data_files/model_results/problem_segments/problem_segment_summary.csv"
)

new <- fread(
  "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated/problem_segments/problem_segment_summary.csv"
)
old[
  ,
  seg_id := paste(
    chr,
    start,
    end,
    sep = "_"
  )
]

new[
  ,
  seg_id := paste(
    chr,
    start,
    end,
    sep = "_"
  )
]
old_set <- unique(old$seg_id)
new_set <- unique(new$seg_id)

shared <- intersect(old_set, new_set)

old_only <- setdiff(old_set, new_set)

new_only <- setdiff(new_set, old_set)
comparison_summary <- data.table(
  old_problem_segments = length(old_set),
  new_problem_segments = length(new_set),
  shared_problem_segments = length(shared),
  old_only_problem_segments = length(old_only),
  new_only_problem_segments = length(new_only)
)

comparison_summary

cat(
  "Fraction of old problematic segments still problematic:",
  round(
    length(shared) /
      length(old_set),
    3
  ),
  "\n"
)

cat(
  "Fraction of new problematic segments already problematic before:",
  round(
    length(shared) /
      length(new_set),
    3
  ),
  "\n"
)