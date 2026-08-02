library(dplyr)
library(magrittr)
library(purrr)
library(tidyr)
library(ggplot2)
library(stringr)


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