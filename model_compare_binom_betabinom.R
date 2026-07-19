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
  
  d <- dat[chr == regions$chr[i] &
             start == regions$start[i] &
             end == regions$end[i]]
  
  out <- data.table(
    chr = regions$chr[i], start = regions$start[i], end = regions$end[i],
    
    AIC_binom_T1 = NA_real_, AIC_bb_T1 = NA_real_,
    BIC_binom_T1 = NA_real_, BIC_bb_T1 = NA_real_,
    logLik_binom_T1 = NA_real_, logLik_bb_T1 = NA_real_,
    LRT_stat_T1 = NA_real_, LRT_df_T1 = NA_real_, LRT_p_T1 = NA_real_,
    phi_T1 = NA_real_, status_T1 = "NOT_RUN",
    
    AIC_binom_T2 = NA_real_, AIC_bb_T2 = NA_real_,
    BIC_binom_T2 = NA_real_, BIC_bb_T2 = NA_real_,
    logLik_binom_T2 = NA_real_, logLik_bb_T2 = NA_real_,
    LRT_stat_T2 = NA_real_, LRT_df_T2 = NA_real_, LRT_p_T2 = NA_real_,
    phi_T2 = NA_real_, status_T2 = "NOT_RUN"
  )
  
  fit_one_group <- function(dd, suffix){
    
    if(nrow(dd) < 10 || length(unique(dd$Weight.g)) < 2){
      out[[paste0("status_", suffix)]] <<- "INSUFFICIENT_DATA"
      return()
    }
    
    tryCatch({
      
      dd[, Weight_sc := as.numeric(scale(Weight.g))]
      
      m_binom <- glmmTMB(
        cbind(meth_reads, unmeth_reads) ~
          Weight_sc + (1|ID_ma) + (1|ID_pa),
        data = dd,
        family = binomial()
      )
      
      m_bb <- glmmTMB(
        cbind(meth_reads, unmeth_reads) ~
          Weight_sc + (1|ID_ma) + (1|ID_pa),
        data = dd,
        family = betabinomial(link = "logit")
      )
      
      lrt <- anova(m_binom, m_bb)
      
      out[[paste0("AIC_binom_", suffix)]] <<- AIC(m_binom)
      out[[paste0("AIC_bb_", suffix)]] <<- AIC(m_bb)
      
      out[[paste0("BIC_binom_", suffix)]] <<- BIC(m_binom)
      out[[paste0("BIC_bb_", suffix)]] <<- BIC(m_bb)
      
      out[[paste0("logLik_binom_", suffix)]] <<- as.numeric(logLik(m_binom))
      out[[paste0("logLik_bb_", suffix)]] <<- as.numeric(logLik(m_bb))
      
      out[[paste0("LRT_stat_", suffix)]] <<- lrt$Chisq[2]
      out[[paste0("LRT_df_", suffix)]] <<- lrt$Df[2]
      out[[paste0("LRT_p_", suffix)]] <<- lrt$`Pr(>Chisq)`[2]
      
      out[[paste0("phi_", suffix)]] <<- sigma(m_bb)
      
      out[[paste0("status_", suffix)]] <<-
        if(!is.null(m_bb$sdr) && isTRUE(m_bb$sdr$pdHess))
          "OK"
      else
        "NON_POSITIVE_HESSIAN"
      
    }, error = function(e){
      
      out[[paste0("status_", suffix)]] <<-
        paste0("ERROR: ", conditionMessage(e))
    })
  }
  
  fit_one_group(d[TimePoint == "T1"], "T1")
  fit_one_group(d[TimePoint == "T2"], "T2")
  
  out
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
  
  mean_AIC_advantage_BB_T1 =
    mean(all_results$AIC_binom_T1 - all_results$AIC_bb_T1, na.rm=TRUE),
  
  median_AIC_advantage_BB_T1 =
    median(all_results$AIC_binom_T1 - all_results$AIC_bb_T1, na.rm=TRUE),
  
  mean_BIC_advantage_BB_T1 =
    mean(all_results$BIC_binom_T1 - all_results$BIC_bb_T1, na.rm=TRUE),
  
  median_BIC_advantage_BB_T1 =
    median(all_results$BIC_binom_T1 - all_results$BIC_bb_T1, na.rm=TRUE),
  
  fraction_segments_BB_better_AIC_T1 =
    mean(all_results$AIC_bb_T1 < all_results$AIC_binom_T1, na.rm=TRUE),
  
  fraction_segments_BB_better_BIC_T1 =
    mean(all_results$BIC_bb_T1 < all_results$BIC_binom_T1, na.rm=TRUE),
  
  fraction_segments_LRT_significant_T1 =
    mean(all_results$LRT_p_T1 < 0.05, na.rm=TRUE),
  
  
  
  mean_AIC_advantage_BB_T2 =
    mean(all_results$AIC_binom_T2 - all_results$AIC_bb_T2, na.rm=TRUE),
  
  median_AIC_advantage_BB_T2 =
    median(all_results$AIC_binom_T2 - all_results$AIC_bb_T2, na.rm=TRUE),
  
  mean_BIC_advantage_BB_T2 =
    mean(all_results$BIC_binom_T2 - all_results$BIC_bb_T2, na.rm=TRUE),
  
  median_BIC_advantage_BB_T2 =
    median(all_results$BIC_binom_T2 - all_results$BIC_bb_T2, na.rm=TRUE),
  
  fraction_segments_BB_better_AIC_T2 =
    mean(all_results$AIC_bb_T2 < all_results$AIC_binom_T2, na.rm=TRUE),
  
  fraction_segments_BB_better_BIC_T2 =
    mean(all_results$BIC_bb_T2 < all_results$BIC_binom_T2, na.rm=TRUE),
  
  fraction_segments_LRT_significant_T2 =
    mean(all_results$LRT_p_T2 < 0.05, na.rm=TRUE)
  
)

print(summary_stats)

summary_counts <- data.table(
  
  segments_with_valid_T1 =
    sum(!is.na(all_results$AIC_binom_T1) &
          !is.na(all_results$AIC_bb_T1)),
  
  segments_BB_better_AIC_T1 =
    sum(all_results$AIC_bb_T1 <
          all_results$AIC_binom_T1,
        na.rm=TRUE),
  
  segments_BB_better_BIC_T1 =
    sum(all_results$BIC_bb_T1 <
          all_results$BIC_binom_T1,
        na.rm=TRUE),
  
  segments_with_valid_T2 =
    sum(!is.na(all_results$AIC_binom_T2) &
          !is.na(all_results$AIC_bb_T2)),
  
  segments_BB_better_AIC_T2 =
    sum(all_results$AIC_bb_T2 <
          all_results$AIC_binom_T2,
        na.rm=TRUE),
  
  segments_BB_better_BIC_T2 =
    sum(all_results$BIC_bb_T2 <
          all_results$BIC_binom_T2,
        na.rm=TRUE)
  
)

print(summary_counts)

fwrite(summary_stats,"/scratch/project_2010912/ezel/salmon/data_files/model_results/binom_vs_betabinom/summary_stats_AIC_BIC.csv")