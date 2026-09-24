library(dplyr)
library(magrittr)
library(purrr)
library(tidyr)
library(ggplot2)
library(stringr)

#calculates the average methylation of homeolog blocks where both sides have cpgs
#and also includes soem stats for the diff methylatoin results(not related to homeologs)
mk_homeo=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_diff_meth_paired.csv")
dss_homeo=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/dss/dss_diff_meth_paired.csv")

#we are looking at the methylation levels of the blocks where both
#x and y sides of the block have cpgs mapping to them
mk_meth_stats <- mk_homeo %>%
  filter(!is.na(homeolog_block_name)) %>%
  group_by(homeolog_block_name, homeolog_side) %>%
  summarise(
    n_CpGs = n(),
    mean_meth_diff = mean(meth.diff, na.rm = TRUE),
    median_meth_diff = median(meth.diff, na.rm = TRUE),
    sd_meth_diff = sd(meth.diff, na.rm = TRUE),
    min_meth_diff = min(meth.diff, na.rm = TRUE),
    max_meth_diff = max(meth.diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = homeolog_side,
    values_from = c(
      n_CpGs,
      mean_meth_diff,
      median_meth_diff,
      sd_meth_diff,
      min_meth_diff,
      max_meth_diff
    ),
    names_sep = "_"
  ) %>%
  filter(
    n_CpGs_x > 0,
    n_CpGs_y > 0
  )

write.csv(mk_meth_stats,"/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/stats_for_paired_blocks.csv")

mk_meth_stats <- dss_homeo %>%
  filter(!is.na(homeolog_block_name)) %>%
  group_by(homeolog_block_name, homeolog_side) %>%
  summarise(
    n_CpGs = n(),
    mean_meth_diff = mean(meth.diff, na.rm = TRUE),
    median_meth_diff = median(meth.diff, na.rm = TRUE),
    sd_meth_diff = sd(meth.diff, na.rm = TRUE),
    min_meth_diff = min(meth.diff, na.rm = TRUE),
    max_meth_diff = max(meth.diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = homeolog_side,
    values_from = c(
      n_CpGs,
      mean_meth_diff,
      median_meth_diff,
      sd_meth_diff,
      min_meth_diff,
      max_meth_diff
    ),
    names_sep = "_"
  ) %>%
  filter(
    n_CpGs_x > 0,
    n_CpGs_y > 0
  )

write.csv(mk_meth_stats,"/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/dss/stats_for_paired_blocks.csv")


######################################
######################################
######################################
##also some additional stats regarding the more general diff-meth cpgs

common_sites <- inner_join(
  mk_sites,
  dss_sites,
  by = c("chr_ncbi", "start")
)

mk_only <- anti_join(
  mk_sites,
  dss_sites,
  by = c("chr_ncbi", "start")
)

dss_only <- anti_join(
  dss_sites,
  mk_sites,
  by = c("chr_ncbi", "start")
)
library(GenomicRanges)

common_gr <- GRanges(
  seqnames = common_sites$chr_ncbi,
  ranges = IRanges(
    start = common_sites$start,
    end = common_sites$start
  )
)

mk_only_gr <- GRanges(
  seqnames = mk_only$chr_ncbi,
  ranges = IRanges(
    start = mk_only$start,
    end = mk_only$start
  )
)

dss_only_gr <- GRanges(
  seqnames = dss_only$chr_ncbi,
  ranges = IRanges(
    start = dss_only$start,
    end = dss_only$start
  )
)
mk_distance <- distanceToNearest(
  mk_only_gr,
  common_gr
)

mk_distances <- data.frame(
  mk_site = mk_only$start[queryHits(mk_distance)],
  chr = as.character(seqnames(mk_only_gr)[queryHits(mk_distance)]),
  distance = mcols(mk_distance)$distance
)
dss_distance <- distanceToNearest(
  dss_only_gr,
  common_gr
)

dss_distances <- data.frame(
  dss_site = dss_only$start[queryHits(dss_distance)],
  chr = as.character(seqnames(dss_only_gr)[queryHits(dss_distance)]),
  distance = mcols(dss_distance)$distance
)
summary(dss_distances$distance)
summary(mk_distances$distance)
########################################
########################################
#get the mean methylation values for the two sides of all the homeolog blocks
mk_homeo=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_diff_meth_paired.csv")
dss_homeo=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/dss/dss_diff_meth_paired.csv")
library(dplyr)
library(tidyr)

mk_block_means <- mk_homeo %>%
  filter(!is.na(homeolog_block_name)) %>%
  group_by(homeolog_block_name, homeolog_side) %>%
  summarise(
    mean_meth_diff = mean(meth.diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = homeolog_side,
    values_from = mean_meth_diff,
    names_prefix = "mean_meth_diff_"
  ) %>%
  arrange(homeolog_block_name)

write.csv(
  mk_block_means,
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_block_mean_methylation.csv",
  row.names = FALSE
)

### t-test
blocks <- read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/dss/dss_block_mean_methylation.csv")

meth_groups <- blocks %>%
  select(homeolog_block_name, mean_meth_diff_x, mean_meth_diff_y) %>%
  mutate(group = ifelse(!is.na(mean_meth_diff_x) & !is.na(mean_meth_diff_y), "Both sides", "One side only")) %>%
  pivot_longer(c(mean_meth_diff_x, mean_meth_diff_y), values_to = "meth_diff", values_drop_na = TRUE) %>%
  select(group, meth_diff)

meth_groups %>% 
  group_by(group) %>% 
  summarise(n=n(), mean=mean(meth_diff), median=median(meth_diff), sd=sd(meth_diff))

t.test(meth_diff ~ group, data=meth_groups)

#################################
#################################
#################################
#unrelated to homeologs
###lets check the spread of cpgs identified by the methods
library(dplyr)
library(GenomicRanges)
library(ggplot2)

# Read data
mk <- read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/methylkit_diff_meth_bases.csv")
dss <- read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/dss_diff_meth_bases.csv")
mk_short=dss[,c("seqnames","start","end","chr_ncbi")]
mk_short

mk_nn <- mk_short %>%
  arrange(chr_ncbi, start) %>%
  group_by(chr_ncbi) %>%
  mutate(
    prev_distance = start - lag(start),
    next_distance = lead(start) - start,
    nearest_distance = case_when(
      is.na(prev_distance) ~ next_distance,
      is.na(next_distance) ~ prev_distance,
      TRUE ~ pmin(prev_distance, next_distance)
    )
  ) %>%
  ungroup() %>%
  select(
    chr_ncbi,
    start,
    nearest_distance,
    prev_distance,
    next_distance
  )

#save the results
sink("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/dss_distance_summary.txt")
mk_nn %>%
  group_by(chr_ncbi) %>%
  group_walk(~{
    
    cat("\n====================================\n")
    cat("Chromosome", .y$chr_ncbi, "\n")
    cat("====================================\n")
    
    cat("Number of significant CpGs:",
        nrow(.x), "\n")
    
    cat("Mean nearest-neighbour distance:",
        round(mean(.x$nearest_distance, na.rm = TRUE), 2),
        "bp\n")
    
    cat("Median nearest-neighbour distance:",
        round(median(.x$nearest_distance, na.rm = TRUE), 2),
        "bp\n")
    
    cat("Standard deviation:",
        round(sd(.x$nearest_distance, na.rm = TRUE), 2),
        "bp\n")
    
    q <- quantile(
      .x$nearest_distance,
      probs = c(0.25, 0.75, 0.90, 0.95),
      na.rm = TRUE
    )
    
    cat("25th percentile:", round(q[1], 2), "bp\n")
    cat("75th percentile:", round(q[2], 2), "bp\n")
    cat("90th percentile:", round(q[3], 2), "bp\n")
    cat("95th percentile:", round(q[4], 2), "bp\n")
  })

sink()
