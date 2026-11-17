
#calculate the true&false positives etc

res=read.csv("/scratch/project_2010912/ezel/salmon/data_files/simulated_data/second_try/model_results_REML_T/all_model_results_REML_T.csv")
results=res
# Define whether the true and estimated maternal effects are present
results$true_effect <- results$true_sd_ma > 0.05
results$estimated_effect <- results$estimated_sd_ma > 0.05
# Create the four data frames
true_positive <- results[
  results$true_effect & results$estimated_effect, ]
true_negative <- results[
  !results$true_effect & !results$estimated_effect, ]
false_positive <- results[
  !results$true_effect & results$estimated_effect, ]
false_negative <- results[
  results$true_effect & !results$estimated_effect, ]
dim(true_positive)
TP <- nrow(true_positive)
TN <- nrow(true_negative)
FP <- nrow(false_positive)
FN <- nrow(false_negative)
power <- TP / (TP + FN)
specificity <- TN / (TN + FP)
false_positive_rate <- FP / (FP + TN)
false_negative_rate <- FN / (FN + TP)
data.frame(
  TP = TP,
  TN = TN,
  FP = FP,
  FN = FN,
  power = power,
  specificity = specificity,
  false_positive_rate = false_positive_rate,
  false_negative_rate = false_negative_rate
)