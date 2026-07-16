library(glmmTMB)
library(data.table)
library(lmtest)

#################
path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]
###
base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"
out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results/binom_vs_betabinom"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

###function to fit the models
fit_region_compare <- function(i){
  
  d <- dat[
    chr == regions$chr[i] &
      start == regions$start[i] &
      end == regions$end[i]
  ]
  
  out <- data.table(
    chr = regions$chr[i],
    start = regions$start[i],
    end = regions$end[i],
    
    AIC_binom = NA_real_,
    AIC_bb = NA_real_,
    
    BIC_binom = NA_real_,
    BIC_bb = NA_real_,
    
    logLik_binom = NA_real_,
    logLik_bb = NA_real_,
    
    LRT_stat = NA_real_,
    LRT_df = NA_real_,
    LRT_p = NA_real_,
    
    phi = NA_real_,
    status = "NOT_RUN"
  )
  
  if(
    nrow(d) < 10 ||
    length(unique(d$Weight.g)) < 2
  ){
    out$status <- "INSUFFICIENT_DATA"
    return(out)
  }
  
  tryCatch({
    
    d$Weight_sc <- as.numeric(scale(d$Weight.g))
    
    m_binom <- glmmTMB(
      cbind(meth_reads, unmeth_reads) ~
        Weight_sc +
        (1 | ID_ma) +
        (1 | ID_pa),
      data = d,
      family = binomial()
    )
    
    m_bb <- glmmTMB(
      cbind(meth_reads, unmeth_reads) ~
        Weight_sc +
        (1 | ID_ma) +
        (1 | ID_pa),
      data = d,
      family = betabinomial(link = "logit")
    )
    
    out$AIC_binom <- AIC(m_binom)
    out$AIC_bb <- AIC(m_bb)
    
    out$BIC_binom <- BIC(m_binom)
    out$BIC_bb <- BIC(m_bb)
    
    out$logLik_binom <- as.numeric(logLik(m_binom))
    out$logLik_bb <- as.numeric(logLik(m_bb))
    
    lrt <- anova(m_binom, m_bb)
    
    out$LRT_stat <- lrt$Chisq[2]
    out$LRT_df <- lrt$Df[2]
    out$LRT_p <- lrt$`Pr(>Chisq)`[2]
    
    out$phi <- sigma(m_bb)
    
    if(
      !is.null(m_bb$sdr) &&
      isTRUE(m_bb$sdr$pdHess)
    ){
      out$status <- "OK"
    } else {
      out$status <- "NON_POSITIVE_HESSIAN"
    }
    
    out
    
  }, error = function(e){
    
    out$status <- paste0(
      "ERROR: ",
      conditionMessage(e)
    )
    
    out
  })
}

cs_files <- list.files(
  base_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)



files <- list.files(
  out_dir,
  pattern = "_model_comparison\\.csv$",
  full.names = TRUE
)
for(cs_file in cs_files){
  
  chr_name <- sub(
    "_segment_cs\\.csv$",
    "",
    basename(cs_file)
  )
  
  cov_file <- file.path(
    base_dir,
    paste0(chr_name, "_segment_coverage.csv")
  )
  
  if(!file.exists(cov_file)){
    cat("Missing coverage file for", chr_name, "\n")
    next
  }
  
  cat("Processing", chr_name, "\n")
  
  meth <- fread(cs_file)
  cov <- fread(cov_file)
  
  meth_long <- melt(
    meth,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "coverage"
  )
  
  dat <- merge(
    meth_long,
    cov_long,
    by = c("chr", "start", "end", "sample")
  )
  
  dat <- merge(
    dat,
    meta,
    by = "sample",
    all.x = TRUE
  )
  
  dat[, unmeth_reads := coverage - meth_reads]
  
  regions <- unique(
    dat[, .(chr, start, end)]
  )
  
  results <- rbindlist(
    lapply(seq_len(nrow(regions)), fit_region_compare),
    fill = TRUE
  )
  
  out_file <- file.path(
    out_dir,
    paste0("model_", chr_name, ".csv")
  )
  
  fwrite(
    results,
    out_file
  )
  
  cat(
    "Saved",
    nrow(results),
    "regions to",
    out_file,
    "\n\n"
  )
}

cat("Done.\n")



######################
files <- list.files(
  out_dir,
  pattern = "^model_.*\\.csv$",
  full.names = TRUE
)
all_results <- rbindlist(
  lapply(files, fread),
  fill = TRUE
)
#################################
#################################
#summary stats across all chromosomes
summary_stats <- data.table(
  
  total_segments = nrow(all_results),
  
  mean_AIC_binom_minus_bb=
    mean(
      all_results$AIC_binom -
        all_results$AIC_bb,
      na.rm = TRUE
    ),
  
  median_AIC_binom_minus_bb  =
    median(
      all_results$AIC_binom -
        all_results$AIC_bb,
      na.rm = TRUE
    ),
  
  mean_BIC_binom_minus_bb  =
    mean(
      all_results$BIC_binom -
        all_results$BIC_bb,
      na.rm = TRUE
    ),
  
  median_BIC_binom_minus_bb =
    median(
      all_results$BIC_binom -
        all_results$BIC_bb,
      na.rm = TRUE
    ),
  fraction_segments_BB_better_AIC =
    mean(
      all_results$AIC_bb <
        all_results$AIC_binom,
      na.rm = TRUE
    ),
  
  fraction_segments_BB_better_BIC =
    mean(
      all_results$BIC_bb <
        all_results$BIC_binom,
      na.rm = TRUE
    ),
  
  fraction_LRT_significant =
    mean(
      all_results$LRT_p < 0.05,
      na.rm = TRUE
    )
)

print(summary_stats)