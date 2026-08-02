#we want to see if the diff.meth levels for the x-y pairs of the same homeolog blocks are 
#more correlated than any random x-y pair

#we will use the files we created containing stats of the blocks where both x and y sides have CpGs mapping to them
library(dplyr)
library(ggplot2)

# Load data
blocks <- read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_stats_for_paired_blocks.csv")
# Keep blocks with CpGs on both sides
blocks_filtered <- blocks %>%
  filter(
    n_CpGs_x > 0,
    n_CpGs_y > 0
  )

# Optional: require more reliable estimates
# blocks_filtered <- blocks_filtered %>%
#   filter(
#     n_CpGs_x >= 5,
#     n_CpGs_y >= 5
#   )

# Remove missing methylation values
blocks_filtered <- blocks_filtered %>%
  filter(
    !is.na(mean_meth_diff_x),
    !is.na(mean_meth_diff_y)
  )

#actual corr

observed_cor <- cor(
  blocks_filtered$mean_meth_diff_x,
  blocks_filtered$mean_meth_diff_y,
  method = "pearson"
)

print(observed_cor)

#randomize
set.seed(123)
n_perm <- 10000

random_cor <- replicate(
  n_perm,
  cor(
    blocks_filtered$mean_meth_diff_x,
    sample(blocks_filtered$mean_meth_diff_y),
    method = "pearson"
  )
)

# empirical p-value
p_value <- mean(
  abs(random_cor) >= abs(observed_cor)
)

print(p_value)
#we want to get the null dsitribution(randomly assign xs to ys and check the correlation of those)

plot_df <- data.frame(
  random_cor = random_cor
)

p=ggplot(plot_df, aes(random_cor)) +
  geom_histogram(
    bins = 60,
    fill = "#D8B4FE",
    colour = "white"
  ) +
  geom_vline(
    xintercept = observed_cor,
    colour = "#DB2777",
    linewidth = 1.2
  ) +
  theme_classic() +
  labs(
    title = paste0(
      "Homeolog methylation correlation (MethylKit)\n",
      "Observed r = ",
      round(observed_cor, 3),
      ", p = ",
      signif(p_value, 3)
    ),
    x = "Randomized correlation",
    y = "Count"
  )

ggsave("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/dss/correlation_null_test.png",p)