#power analysis for the simulated model
library(data.table)


#accuracy analysis at different thresholds
results <-read.csv("/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/model_results_REML_T/all_model_results_REML_T.csv")
setDT(results)

accuracy_thresholds <- c(0.02, 0.05, 0.08)

accuracy_table <- rbindlist(lapply(
  c("ma", "pa"),
  function(parent) {
    true_col <- paste0("true_sd_", parent)
    est_col  <- paste0("estimated_sd_", parent)
    
    rbindlist(lapply(accuracy_thresholds, function(tol) {
      valid <- !is.na(results[[true_col]]) &
        !is.na(results[[est_col]])
      
      error <- abs(results[[est_col]][valid] -
                     results[[true_col]][valid])
      
      data.table(
        parent = if (parent == "ma") "Mother" else "Father",
        tolerance = tol,
        n_segments = sum(valid),
        n_accurate = sum(error <= tol),
        accuracy_percent = 100 * mean(error <= tol)
      )
    }))
  }
))

print(accuracy_table)



####
#

detection_threshold <- 0.05

detection_table <- rbindlist(lapply(
  c("ma", "pa"),
  function(parent) {
    true_col <- paste0("true_sd_", parent)
    est_col  <- paste0("estimated_sd_", parent)
    
    valid <- !is.na(results[[true_col]]) &
      !is.na(results[[est_col]])
    
    true_present <- results[[true_col]][valid] >
      detection_threshold
    
    estimated_present <- results[[est_col]][valid] >
      detection_threshold
    
    TP <- sum(true_present & estimated_present)
    FP <- sum(!true_present & estimated_present)
    TN <- sum(!true_present & !estimated_present)
    FN <- sum(true_present & !estimated_present)
    
    data.table(
      parent = if (parent == "ma") "Mother" else "Father",
      threshold = detection_threshold,
      TP = TP,
      FP = FP,
      TN = TN,
      FN = FN,
      power = if (TP + FN > 0) TP / (TP + FN) else NA_real_,
      true_negative_detection = if (TN + FP > 0) TN / (TN + FP) else NA_real_,
      false_positive_rate = if (FP + TN > 0) FP / (FP + TN) else NA_real_,
      false_negative_rate = if (FN + TP > 0) FN / (FN + TP) else NA_real_
    )
  }
))

print(detection_table)



#roc analysis for different thresholds 
library(pROC)
library(ggplot2)
library(dplyr)

# model results file
data <- read.csv("/scratch/project_2010912/ezel/salmon/data_files/simulated_data/no_extras/model_results/all_model_results_REML_T.csv")

thresholds <- c(0.001,.01, .05,.08,.1,.12,.14)

calculate_rocs <- function(data, true_var, estimated_var, thresholds) {
  
  results <- list()
  
  for (t in thresholds) {
    
    true_class <- ifelse(data[[true_var]] > t, 1, 0)
    
    roc_obj <- roc(
      response = true_class,
      predictor = data[[estimated_var]],
      quiet = TRUE
    )
    
    results[[as.character(t)]] <- data.frame(
      threshold = t,
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      AUC = as.numeric(auc(roc_obj))
    )
  }
  
  bind_rows(results)
}

# for ma
roc_results_ma <- calculate_rocs(
  data, "true_sd_ma", "estimated_sd_ma", thresholds
)


p=ggplot(roc_results_ma,
       aes(x = FPR, y = TPR,
           color = factor(threshold))) +
  
  geom_line(linewidth = 1) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey"
  ) +
  
  coord_equal() +
  
  labs(
    title = "ROC curves, Ma SD- All segments",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "True sd threshold"
  ) +
  
  theme_minimal()

#ggsave("my_plot.png", plot = p, width = 8, height = 6, dpi = 300)


#repeat for the variable segments

both_T1_T2_segs=read.csv("/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/segments_variable_both_T1_T2.csv")

library(dplyr)

data_filtered <- data %>%
  semi_join(
    both_T1_T2_segs,
    by = c("chr", "start", "end")
  )

roc_filtered<- calculate_rocs(
  data, "true_sd_ma", "estimated_sd_ma", thresholds
)

#ROC curves, Ma SD- All segments

p_filtered=ggplot(roc_filtered,
         aes(x = FPR, y = TPR,
             color = factor(threshold))) +
  
  geom_line(linewidth = 1) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    color = "grey"
  ) +
  
  coord_equal() +
  
  labs(
    title = "ROC curves, Ma SD- Variable segments",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "True sd threshold"
  ) +
  
  theme_minimal()


#compare with the loci

loci_results=read.csv("/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/model_results_by_chromosome/stats/all_chromosomes_combined.csv")