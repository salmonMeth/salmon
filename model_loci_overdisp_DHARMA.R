library(methylKit)
library(data.table)
library(glmmTMB)
library(GenomicRanges)
library(DHARMa)

db_file <- "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds"
path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]


out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_loci_results/chr3"
dharma_dir <- file.path(out_dir, "DHARMa")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(dharma_dir, recursive=TRUE, showWarnings=FALSE)

# Load data
db <- readRDS(db_file)

dat <- methylKit::getData(db)
chr <- "NC_059444.1"

db_1000 <- dat[dat$chr == chr, ]
# Select NC_059442.1 and first 1000 CpGs

#db_1000 <- db_chr[1:10]

#db_1000 <- db_chr[1:nrow(db_chr)]

sample_ids <- db@sample.ids

# Convert methylKit object to long format
dat <- rbindlist(lapply(seq_along(sample_ids), function(i)
  data.table(
    Running.ID=sample_ids[i],
    chr=db_1000$chr,
    start=db_1000$start,
    end=db_1000$end,
    meth_reads=db_1000[[paste0("numCs", i)]],
    unmeth_reads=db_1000[[paste0("numTs", i)]]
  )
))

# Add metadata
meta_sub <- meta[, .(sample, TimePoint, Weight.g, ID_ma, ID_pa)]
dat <- merge(dat, meta_sub, by.x="Running.ID", by.y="sample", all.x=TRUE)

# Check metadata
stopifnot(!anyNA(dat$TimePoint), !anyNA(dat$Weight.g),
          !anyNA(dat$ID_ma), !anyNA(dat$ID_pa))

# Standardize weight across all samples
dat[, Weight_sc := as.numeric(scale(Weight.g))]

# GLMM fitting function
fit_model <- function(d) {
  if (nrow(d) < 3)
    return(list(model=NULL, sigma=NA_real_, coef=NULL,
                convergence=NA_integer_, pdHess=FALSE,
                AIC=NA_real_, error="Too few observations"))
  
  tryCatch({
    m <- glmmTMB(
      cbind(meth_reads, unmeth_reads) ~ Weight_sc + (1|ID_ma) + (1|ID_pa),
      data=d, family=betabinomial(link="logit")
    )
    co <- summary(m)$coefficients$cond
    list(
      model=m,
      sigma=sigma(m),
      coef=co["Weight_sc",],
      convergence=m$fit$convergence,
      pdHess=isTRUE(m$sdr$pdHess),
      AIC=AIC(m),
      error=NA_character_
    )
  }, error=function(e)
    list(model=NULL, sigma=NA_real_, coef=NULL,
         convergence=NA_integer_, pdHess=FALSE,
         AIC=NA_real_, error=conditionMessage(e))
  )
}

# DHARMa function
run_dharma <- function(m) {
  if (is.null(m))
    return(list(uniformity=NULL, dispersion=NULL,
                outliers=NULL, error="No model"))
  
  tryCatch({
    set.seed(123)
    sim <- simulateResiduals(m, n=1000, plot=FALSE)
    list(
      uniformity=testUniformity(sim),
      dispersion=testDispersion(sim),
      outliers=testOutliers(sim),
      error=NA_character_
    )
  }, error=function(e)
    list(uniformity=NULL, dispersion=NULL,
         outliers=NULL, error=conditionMessage(e))
  )
}

get_p <- function(x) if (is.null(x)) NA_real_ else x$p.value
get_stat <- function(x) {if (is.null(x)) return(NA_real_)
  unname(x$statistic)
}
loci <- unique(dat[, .(chr, start, end)])

glmm_results <- vector("list", nrow(loci))
dharma_results <- vector("list", nrow(loci))

# Fit models
for (i in seq_len(nrow(loci))) {
  
  loc <- loci[i]
  d <- dat[chr==loc$chr & start==loc$start & end==loc$end]
  
  fit_T1 <- fit_model(d[TimePoint=="T1"])
  fit_T2 <- fit_model(d[TimePoint=="T2"])
  
  # GLMM results only
  glmm_results[[i]] <- data.table(
    chr=loc$chr, start=loc$start, end=loc$end,
    sigma_T1=fit_T1$sigma,
    estimate_T1=if(!is.null(fit_T1$coef)) fit_T1$coef["Estimate"] else NA_real_,
    SE_T1=if(!is.null(fit_T1$coef)) fit_T1$coef["Std. Error"] else NA_real_,
    z_T1=if(!is.null(fit_T1$coef)) fit_T1$coef["z value"] else NA_real_,
    p_T1=if(!is.null(fit_T1$coef)) fit_T1$coef["Pr(>|z|)"] else NA_real_,
    convergence_T1=fit_T1$convergence,
    pdHess_T1=fit_T1$pdHess,
    AIC_T1=fit_T1$AIC,
    model_error_T1=fit_T1$error,
    sigma_T2=fit_T2$sigma,
    estimate_T2=if(!is.null(fit_T2$coef)) fit_T2$coef["Estimate"] else NA_real_,
    SE_T2=if(!is.null(fit_T2$coef)) fit_T2$coef["Std. Error"] else NA_real_,
    z_T2=if(!is.null(fit_T2$coef)) fit_T2$coef["z value"] else NA_real_,
    p_T2=if(!is.null(fit_T2$coef)) fit_T2$coef["Pr(>|z|)"] else NA_real_,
    convergence_T2=fit_T2$convergence,
    pdHess_T2=fit_T2$pdHess,
    AIC_T2=fit_T2$AIC,
    model_error_T2=fit_T2$error
  )
  
  # DHARMa
  dh_T1 <- run_dharma(fit_T1$model)
  dh_T2 <- run_dharma(fit_T2$model)
  
  dharma_results[[i]] <- data.table(
    chr   = loc$chr,
    start = loc$start,
    end   = loc$end,
    
    # T1
    DHARMa_uniformity_statistic_T1 =get_stat(dh_T1$uniformity),
    
    DHARMa_uniformity_pvalue_T1 =get_p(dh_T1$uniformity),
    
    DHARMa_dispersion_statistic_T1 =
      get_stat(dh_T1$dispersion),
    
    DHARMa_dispersion_pvalue_T1 =
      get_p(dh_T1$dispersion),
    
    DHARMa_outlier_statistic_T1 =
      get_stat(dh_T1$outliers),
    
    DHARMa_outlier_pvalue_T1 =
      get_p(dh_T1$outliers),
    
    DHARMa_error_T1 =
      dh_T1$error,
    # T2
    DHARMa_uniformity_statistic_T2 =
      get_stat(dh_T2$uniformity),
    
    DHARMa_uniformity_pvalue_T2 =
      get_p(dh_T2$uniformity),
    
    DHARMa_dispersion_statistic_T2 =
      get_stat(dh_T2$dispersion),
    
    DHARMa_dispersion_pvalue_T2 =
      get_p(dh_T2$dispersion),
    
    DHARMa_outlier_statistic_T2 =
      get_stat(dh_T2$outliers),
    
    DHARMa_outlier_pvalue_T2 =
      get_p(dh_T2$outliers),
    
    DHARMa_error_T2 =
      dh_T2$error
  )
  
  message(i, "/", nrow(loci), " ", loc$chr, ":", loc$start)
rm(fit_T1, fit_T2, dh_T1, dh_T2, d)

if (i %% 100 == 0)
    gc()
}



# Combine results
glmm_results_dt <- rbindlist(glmm_results, fill=TRUE)
dharma_results_dt <- rbindlist(dharma_results, fill=TRUE)

# Save GLMM outputs
saveRDS(glmm_results_dt, file.path(out_dir, "glmm_chr3_results.rds"))
fwrite(glmm_results_dt, file.path(out_dir, "glmm_chr3_results.tsv"), sep="\t")

# Save data
saveRDS(dat, file.path(out_dir, "dat_chr3.rds"))

# Save DHARMa outputs separately
saveRDS(dharma_results_dt, file.path(dharma_dir, "DHARMa_chr3_results.rds"))
fwrite(dharma_results_dt, file.path(dharma_dir, "DHARMa_chr3_results.tsv"), sep="\t")

message("Finished.")
message("GLMM results: ", out_dir)
message("DHARMa results: ", dharma_dir)