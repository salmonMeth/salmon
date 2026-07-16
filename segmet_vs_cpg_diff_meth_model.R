
library(limma)
library(ggplot2)
library(qvalue)
#library(dendextend)
library ( glmmTMB )
#library(snow)
############
#modeling the segmented methylation
path_segment_matrix="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/segment_matrix.rds"
path_cpg_matrix="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/cpg_meth_matrix.rds"


########################
########################
########################

limma_methyl_plot <- function(path_segment_matrix) {
  segment_matrix <- readRDS(path_segment_matrix)
  sample_names <- colnames(segment_matrix)
  group <- factor(
    ifelse(grepl("_T1$", sample_names), "T1", "T2")
  )
#model
  design <- model.matrix(~group)
  fit <- lmFit(segment_matrix, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = 2, number = Inf)
  sig <- res[res$adj.P.Val < 0.05 & abs(res$logFC) > 0.02, ]
  
  # p&q-values
  pvals <- fit$p.value[,2]
  cat("\nSummary of raw p-values:\n")
  print(summary(pvals))
  cat("\nSummary of q-values:\n")
  print(summary(qvalue(pvals)))
  cat("\nNumber of significant segments (FDR < 0.05 & |logFC| > 0.02): ",
      nrow(sig), "\n")
  
  # Volcano plot
  res$Significant <- res$adj.P.Val < 0.05
  volcano <- ggplot(res,
                    aes(logFC,
                        -log10(adj.P.Val),
                        color = Significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("lightpink", "#4B0082")) +
    theme_bw() +
    labs(
      title = "Differential methylation",
      x = "logFC",
      y = "-log10(FDR)"
    )
  print(volcano)
  
  # Histogram of p-values
  hist_plot <- ggplot(res, aes(P.Value)) +
    geom_histogram(bins = 40) +
    theme_bw() +
    labs(
      title = "Distribution of p-values",
      x = "P-value",
      y = "Count"
    )
  
  print(hist_plot)
  
  return(list(
    fit = fit,
    results = res,
    pvalues = pvals
  ))
}
#######
#fit the models
res_cpg= limma_methyl_plot(path_cpg_matrix)
res_segment= limma_methyl_plot(path_segment_matrix)

#########################################
#########################################
#compare the results

compare_limma_results <- function(res1, res2,
                                  name1 = "CpG",
                                  name2 = "Segment",
                                  thresholds = c(0.01, 0.05, 0.1)) {
  
  library(qvalue)
  
  summarize_result <- function(x) {
    
    p <- x$results$P.Value
    bh <- x$results$adj.P.Val
    bonf <- p.adjust(p, method = "bonferroni")
    qv <- qvalue(p)$qvalues
    
    data.frame(
      Threshold = thresholds,
      `P-value` = sapply(thresholds, function(t) mean(p < t)),
      `BH (FDR)` = sapply(thresholds, function(t) mean(bh < t)),
      Bonferroni = sapply(thresholds, function(t) mean(bonf < t)),
      `Q-value` = sapply(thresholds, function(t) mean(qv < t)),
      check.names = FALSE
    )
  }
  
  tab1 <- summarize_result(res1)
  tab2 <- summarize_result(res2)
  
  cat("\n========================================\n")
  cat(name1, "\n")
  cat("Number of tests:", nrow(res1$results), "\n\n")
  print(tab1, row.names = FALSE)
  
  cat("\n========================================\n")
  cat(name2, "\n")
  cat("Number of tests:", nrow(res2$results), "\n\n")
  print(tab2, row.names = FALSE)
  invisible(list(
    first = tab1,
    second = tab2
  ))
}

compare_limma_results(
  res_cpg,
  res_segment,
  name1 = "CpG",
  name2 = "Segments"
)


############model
path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
path_chr_1_cs="/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed/NC_059442.1_segment_cs.csv"
path_chr_1_cov="/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed/NC_059442.1_segment_coverage.csv"
meth <- fread(path_chr_1_cs)
cov=fread(path_chr_1_cov)
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

#meth matrix
meth_long <- melt(
  meth,
  id.vars = c("chr", "start", "end"),
  variable.name = "sample",
  value.name = "meth_reads"
)

#get cov matrix
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
  by = "sample"
)

dat$unmeth_reads <- dat$coverage - dat$meth_reads
#try fitting the model to 1 region first
regions <- unique(dat[, .(chr, start, end)])
fit_region <- function(i) {
  
  d <- dat[
    chr == regions$chr[i] &
      start == regions$start[i] &
      end == regions$end[i]
  ]
  
  tryCatch({
    
    model <- glmmTMB(
      cbind(meth_reads, unmeth_reads) ~
        Weight.g +
        (1 | ID_ma) +
        (1 | ID_pa),
      data = d,
      family = betabinomial(link = "logit")
    )
    
    coefs <- summary(model)$coefficients$cond
    
    data.table(
      chr = regions$chr[i],
      start = regions$start[i],
      end = regions$end[i],
      estimate = unname(coefs["Weight.g", "Estimate"]),
      SE = unname(coefs["Weight.g", "Std. Error"]),
      z = unname(coefs["Weight.g", "z value"]),
      pvalue = unname(coefs["Weight.g", "Pr(>|z|)"]),)
  }, error = function(e) {
    NULL
  })
}
################################
################################
#loop over all the regions in chr1

results <- rbindlist(
  lapply(seq_len(nrow(regions)), fit_region),
  fill = TRUE
)

path_results="/scratch/project_2010912/ezel/salmon/data_files/segment_methylation_cov_cs/model_NC_059442.1.csv"
# save
fwrite(
  results,
  path_results
)

