#we gotta model randomly picked loci

library(data.table)
library(methylKit)
library(data.table)
library(DHARMa)
library(glmmTMB)

path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
out_dir <- 
  "/scratch/project_2010912/ezel/salmon/data_files/model_loci_for_dharma_2"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- fread(path_ids)

meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

meta <- meta[
  grepl("_T1$", sample) &
    sample != "176_T1"
]

meta[, Weight.sc := as.numeric(scale(Weight.g))]
meta[, sample := as.character(sample)]

db <- readRDS(
  "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds"
)

locus_data <- as.data.table(methylKit::getData(db))

sample_ids_db <- as.character(db@sample.ids)

t1_idx <- which(grepl("_T1$", sample_ids_db))
sample_ids_t1 <- sample_ids_db[t1_idx]

coverage_cols <- names(locus_data)[db@coverage.index][t1_idx]
cs_cols <- names(locus_data)[db@numCs.index][t1_idx]
ts_cols <- names(locus_data)[db@numTs.index][t1_idx]

meta <- meta[match(sample_ids_t1, sample)]

stopifnot(
  !anyNA(meta$sample),
  identical(meta$sample, sample_ids_t1),
  length(coverage_cols) == nrow(meta),
  length(cs_cols) == nrow(meta),
  length(ts_cols) == nrow(meta)
)

locus_data[, chr_num := as.integer(
  sub("NC_059", "", chr)
) - 441]

meth_mat <- as.matrix(
  locus_data[, ..cs_cols] /
    locus_data[, ..coverage_cols]
)

meth_range <- apply(meth_mat, 1, function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  max(x) - min(x)
})

locus_data[, meth_range := meth_range]

qualifying_loci <- locus_data[
  chr_num %in% 1:29 & meth_range > 0.1
]

counts <- qualifying_loci[
  ,
  .(n_qualifying = .N),
  by = chr_num
]

missing_chr <- setdiff(1:29, counts$chr_num)

if (length(missing_chr) > 0) {
  counts <- rbind(
    counts,
    data.table(
      chr_num = missing_chr,
      n_qualifying = 0L
    )
  )
}

counts <- counts[order(chr_num)]

print(counts)

if (any(counts$n_qualifying < 5000)) {
  stop("Some chromosomes have fewer than 1,000 qualifying loci.")
}

set.seed(123)

diagnostics <- fread("/scratch/project_2010912/ezel/salmon/data_files/beta_binomial_for_actual_loci_data/DHARMa_T1_sampled_loci/DHARMa_T1_sampled_loci.csv")

existing_loci <- unique(diagnostics[, .(chr, start, end)])

qualifying_loci <- qualifying_loci[
  !existing_loci,
  on = .(chr, start, end)
]

sampled_loci <- qualifying_loci[
  ,
  .SD[sample(.N, min(.N, 3000))],
  by = chr_num
]

stopifnot(
  !anyDuplicated(sampled_loci[, .(chr, start, end)])
)

loci <- sampled_loci
chromosomes <- unique(as.character(loci$chr))

fit_one_model <- function(dat) {
  warning_messages <- character()
  
  fit <- tryCatch(
    withCallingHandlers(
      glmmTMB(
        REML = TRUE,
        cbind(meth.Data, unmeth.Data) ~
          Weight.sc + (1 | ID_ma) + (1 | ID_pa),
        data = dat,
        family = betabinomial(link = "logit")
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      error = conditionMessage(fit),
      warnings = paste(unique(warning_messages), collapse = " | ")
    ))
  }
  
  list(
    fit = fit,
    error = NA_character_,
    warnings = paste(unique(warning_messages), collapse = " | ")
  )
}


extract_model <- function(obj, model_name, locus_id, dat) {
  
  fit <- obj$fit
  
  row <- list(
    locus_id = locus_id,
    model = model_name,
    n_observations = nrow(dat),
    n_mothers = uniqueN(dat$ID_ma),
    n_fathers = uniqueN(dat$ID_pa),
    error = obj$error,
    warnings = obj$warnings,
    convergence_code = NA_integer_,
    pdHess = NA,
    converged = FALSE,
    beta_weight = NA_real_,
    beta_weight_se = NA_real_,
    var_maternal = NA_real_,
    var_paternal = NA_real_,
    sd_maternal = NA_real_,
    sd_paternal = NA_real_,
    dispersion_phi = NA_real_,
    dispersion_rho = NA_real_,
    logLik = NA_real_,
    AIC = NA_real_
  )
  
  if (is.null(fit)) {
    return(as.data.table(row))
  }
  
  row$convergence_code <- fit$fit$convergence
  row$pdHess <- isTRUE(fit$sdr$pdHess)
  
  row$converged <- (
    identical(fit$fit$convergence, 0L) &&
      isTRUE(fit$sdr$pdHess)
  )
  
  beta <- fixef(fit)$cond
  
  if ("Weight.sc" %in% names(beta)) {
    row$beta_weight <- unname(beta["Weight.sc"])
    
    row$beta_weight_se <- tryCatch(
      sqrt(vcov(fit)$cond["Weight.sc", "Weight.sc"]),
      error = function(e) NA_real_
    )
  }
  
  vc <- VarCorr(fit)$cond
  
  if ("ID_ma" %in% names(vc)) {
    row$var_maternal <- as.numeric(vc$ID_ma[1, 1])
    row$sd_maternal <- sqrt(row$var_maternal)
  }
  
  if ("ID_pa" %in% names(vc)) {
    row$var_paternal <- as.numeric(vc$ID_pa[1, 1])
    row$sd_paternal <- sqrt(row$var_paternal)
  }
  
  row$dispersion_phi <- tryCatch(
    as.numeric(sigma(fit)),
    error = function(e) NA_real_
  )
  
  if (is.finite(row$dispersion_phi) &&
      row$dispersion_phi > 0) {
    row$dispersion_rho <- 1 / (row$dispersion_phi + 1)
  }
  
  row$logLik <- tryCatch(
    as.numeric(logLik(fit)),
    error = function(e) NA_real_
  )
  
  row$AIC <- tryCatch(
    AIC(fit),
    error = function(e) NA_real_
  )
  
  as.data.table(row)
}


for (chr in chromosomes) {
  
  message("Starting chromosome ", chr)
  
  idx_chr <- which(as.character(loci$chr) == chr)
  
  chr_fits <- vector("list", length(idx_chr))
  chr_results <- vector("list", length(idx_chr))
  
  for (j in seq_along(idx_chr)) {
    
    i <- idx_chr[j]
    
    locus <- loci[i]
    
    locus_id <- paste(
      as.character(locus$chr),
      locus$start,
      locus$end,
      sep = ":"
    )
    cov_i <- as.numeric(locus[1, ..coverage_cols])
    cs_i <- as.numeric(locus[1, ..cs_cols])
    ts_i <- as.numeric(locus[1, ..ts_cols])
    
    dat <- copy(meta)
    
    dat[, `:=`(
      meth.Data = cs_i,
      unmeth.Data = ts_i
    )]
    
    fit_reml <- fit_one_model(dat)
    
    chr_fits[[j]] <- list(
      locus_id = locus_id,
      fit = fit_reml$fit,
      data = copy(dat),
      error = fit_reml$error,
      warnings = fit_reml$warnings
    )
    
    res <- extract_model(
      fit_reml,
      "REML",
      locus_id,
      dat
    )
    
    res[, `:=`(
      chr = as.character(locus$chr),
      start = locus$start,
      end = locus$end,
      meth_range = locus$meth_range
    )]
    
    chr_results[[j]] <- res
    
    if (j %% 100 == 0) {
      message(
        "Chromosome ", chr, ": ",
        j, " / ", length(idx_chr), " loci completed"
      )
    }
  }
  
  chr_results_dt <- rbindlist(chr_results, fill = TRUE)
  
  chr_safe <- gsub("[^A-Za-z0-9_.-]", "_", chr)
  
  csv_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_REML_results.csv")
  )
  
  rds_file <- file.path(
    out_dir,
    paste0("chr_", chr_safe, "_REML_results.rds")
  )
  
  fwrite(chr_results_dt, csv_file)
  
  saveRDS(
    list(
      chromosome = chr,
      results = chr_results_dt,
      fits = chr_fits
    ),
    rds_file
  )
  
  message(
    "Saved chromosome ", chr, "\n",
    "  CSV: ", csv_file, "\n",
    "  RDS: ", rds_file
  )
  
  rm(chr_results_dt, chr_results, chr_fits)
  gc()
}