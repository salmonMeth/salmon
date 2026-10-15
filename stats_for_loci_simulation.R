#extract some stats etc


library(data.table)
library(ggplot2)


results_dir <- "/scratch/project_2010912/ezel/salmon/data_files/simulated_loci/model_results_by_chromosome"

out_dir <- file.path(results_dir, "stats")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  results_dir,
  pattern = "^chr_.*_model_results\\.csv$",
  full.names = TRUE
)

message("Found ", length(files), " chromosome CSV files.")


results_list <- lapply(files, function(f) {
  x <- fread(f)
  x[, source_file := basename(f)]
  x
})

results <- rbindlist(results_list, fill = TRUE)

# Ensure convergence is interpreted as logical
results[, converged := as.logical(converged)]

# Save combined table
fwrite(results, file.path(out_dir, "all_chromosomes_combined.csv"))

print(results[, .N, by = model])

# --------------------------------------------------
# 3. Overall convergence statistics
# --------------------------------------------------

# Each locus should normally have one REML row and one ML row.
convergence_overall <- results[, .(
  n_attempted = .N,
  n_converged = sum(converged %in% TRUE),
  convergence_percent = 100 * mean(converged %in% TRUE),
  n_fit_errors = sum(!is.na(error) & nzchar(error)),
  n_optimizer_nonzero = sum(
    !is.na(convergence_code) & convergence_code != 0
  ),
  n_pdHess_false = sum(pdHess %in% FALSE)
), by = model]

print(convergence_overall)

fwrite(
  convergence_overall,
  file.path(out_dir, "convergence_overall.csv")
)

# Convergence by chromosome and model
convergence_by_chr <- results[, .(
  n_attempted = .N,
  n_converged = sum(converged %in% TRUE),
  convergence_percent = 100 * mean(converged %in% TRUE),
  n_fit_errors = sum(!is.na(error) & nzchar(error))
), by = .(chr, model)]

fwrite(
  convergence_by_chr,
  file.path(out_dir, "convergence_by_chromosome.csv")
)

results_with_means <- merge(
  results,
  locus_means[, .(locus_id, mean_methylation)],
  by = "locus_id",
  all.x = TRUE
)
results_with_means[, convergence_status := fifelse(
  converged %in% TRUE, "Converged", "Non-converged"
)]


reml_ma <- results_with_means[
  model == "REML" &
    is.finite(true_sd_ma) &
    is.finite(var_maternal) &
    var_maternal >= 0
]

# Convert estimated variance to estimated SD
reml_ma[, estimated_sd_ma := sqrt(var_maternal)]

p_ma <- ggplot(
  reml_ma,
  aes(x = true_sd_ma, y = estimated_sd_ma)
) +
  geom_point(alpha = 0.25, size = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "red",
    linewidth = 0.8
  ) +
  labs(
    title = "Estimated vs. true maternal SD (REML)",
    x = "True maternal SD",
    y = "Estimated maternal SD"
  ) +
  theme_bw(base_size = 14) + coord_cartesian(ylim = c(0, 2))


print(p_ma)
ggsave(
  file.path(out_dir, "REML_estimated_vs_true_maternal_SD.png"),
  p_ma,
  width = 8,
  height = 7,
  dpi = 300
)

reml_pa <- results_with_means[
  model == "REML" &
    is.finite(true_sd_pa) &
    is.finite(var_paternal) &
    var_paternal >= 0
]

# Convert estimated variance to estimated SD
reml_pa[, estimated_sd_pa := sqrt(var_paternal)]

p_pa <- ggplot(
  reml_pa,
  aes(x = true_sd_pa, y = estimated_sd_pa)
) +
  geom_point(alpha = 0.25, size = 1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "red",
    linewidth = 0.8
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Estimated vs. true paternal SD (REML)",
    x = "True paternal SD",
    y = "Estimated paternal SD"
  ) +
  theme_bw(base_size = 14)

print(p_pa)

ggsave(
  file.path(out_dir, "REML_estimated_vs_true_maternal_SD.png"),
  p_ma,
  width = 8,
  height = 7,
  dpi = 300
)

library(data.table)

# Use converged REML models only
reml_converged <- results_with_means[
  model == "REML" & converged %in% TRUE
]

# Calculate estimated SDs from estimated variances
reml_converged[, estimated_sd_ma := fifelse(
  is.finite(var_maternal) & var_maternal >= 0,
  sqrt(var_maternal), NA_real_
)]

reml_converged[, estimated_sd_pa := fifelse(
  is.finite(var_paternal) & var_paternal >= 0,
  sqrt(var_paternal), NA_real_
)]

# Bin each parent's true SD
sd_bin <- function(x) {
  fcase(
    !is.finite(x) | x < 0, NA_character_,
    x <= 0.1, "0-0.1",
    x <= 0.2, ">0.1-0.2",
    x > 0.2, ">0.2"
  )
}

reml_converged[, true_ma_bin := sd_bin(true_sd_ma)]
reml_converged[, true_pa_bin := sd_bin(true_sd_pa)]

# Convert to long format so maternal and paternal results
# can be summarized using the same function
ma_long <- reml_converged[
  !is.na(true_ma_bin) & is.finite(estimated_sd_ma),
  .(
    parent = "Maternal",
    true_sd_bin = true_ma_bin,
    true_sd = true_sd_ma,
    estimated_sd = estimated_sd_ma
  )
]

pa_long <- reml_converged[
  !is.na(true_pa_bin) & is.finite(estimated_sd_pa),
  .(
    parent = "Paternal",
    true_sd_bin = true_pa_bin,
    true_sd = true_sd_pa,
    estimated_sd = estimated_sd_pa
  )
]

sd_long <- rbindlist(list(ma_long, pa_long))

# Summary statistics by parent and true-SD bin
summary_stats <- sd_long[
  ,
  .(
    n = .N,
    mean_true_sd = mean(true_sd),
    mean_estimated_sd = mean(estimated_sd),
    sd_estimated_sd = sd(estimated_sd),
    median_estimated_sd = median(estimated_sd),
    q25_estimated_sd = quantile(estimated_sd, 0.25, names = FALSE),
    q75_estimated_sd = quantile(estimated_sd, 0.75, names = FALSE),
    min_estimated_sd = min(estimated_sd),
    max_estimated_sd = max(estimated_sd),
    mean_bias = mean(estimated_sd - true_sd),
    mean_absolute_error = mean(abs(estimated_sd - true_sd)),
    rmse = sqrt(mean((estimated_sd - true_sd)^2))
  ),
  by = .(parent, true_sd_bin)
]

# Order the bins and round numeric columns
summary_stats[, true_sd_bin := factor(
  true_sd_bin,
  levels = c("0-0.1", ">0.1-0.2", ">0.2")
)]

setorder(summary_stats, parent, true_sd_bin)

numeric_cols <- names(summary_stats)[
  vapply(summary_stats, is.numeric, logical(1))
]
summary_stats[, (numeric_cols) := lapply(.SD, round, 4),
              .SDcols = numeric_cols]

print(summary_stats)

