#plots for the segmentation

library(data.table)
library(ggplot2)
library(tidyverse)

###############################################
#PLOT 1
#to compare the number of segment groups
median_consensus <- fread("segmentation_output/summaries/chromosome_summary.csv")
weighted_consensus <- fread("segmentation_output_cov_weighted_consensus/summaries/chromosome_summary.csv")

consensus_dt <- merge(
  median_consensus[, .(chr = chromosome, K_median = n_segment_groups)],
  weighted_consensus[, .(chr = chromosome, K_weighted = n_segment_groups)],
  by = "chr",
  all = TRUE
)

plot_dt <- copy(segment_counts)

dir.create("segmentation_output/T1_vs_T2_plots", recursive = TRUE, showWarnings = FALSE)

# Sort chromosomes
chroms <- sort(unique(plot_dt$chr))

# Split into 4 chunks
chunks <- split(chroms, cut(seq_along(chroms), 4, labels = FALSE))

for (i in seq_along(chunks)) {
  
  chr_chunk <- chunks[[i]]
  
  dt_plot <- plot_dt[chr %in% chr_chunk]
  
  # force discrete treatment
  dt_plot[, group := factor(group)]
  
  y_breaks <- seq(min(dt_plot$K, na.rm = TRUE),
                  max(dt_plot$K, na.rm = TRUE),
                  by = 1)
  
  p <- ggplot(
    dt_plot,
    aes(x = group, y = K, color = group)
  ) +
    geom_point(
      position = position_jitter(width = 0.15, height = 0),
      size = 2,
      alpha = 0.8
    ) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.35,
      fatten = 0
    ) +
    facet_wrap(~chr, ncol = 3, scales = "fixed") +
    scale_y_continuous(breaks = y_breaks) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      plot.title = element_text(size = 11)
    ) +
    labs(
      title = paste0("Total segment group number(K) per chromosome (block ", i, "/4)"),
      x = NULL,
      y = "K (# segment groups)"
    )
  
  # subset consensus
  sub_cons <- consensus_dt[chr %in% chr_chunk]
  
  same_dt <- sub_cons[K_median == K_weighted]
  diff_dt <- sub_cons[K_median != K_weighted]
  
  # identical consensus (single line)
  if (nrow(same_dt) > 0) {
    p <- p +
      geom_hline(
        data = same_dt,
        aes(yintercept = K_median),
        color = "purple",
        linewidth = 0.7,
        inherit.aes = FALSE
      )
  }
  
  # differing consensus (two lines per chromosome)
  if (nrow(diff_dt) > 0) {
    p <- p +
      geom_hline(
        data = diff_dt,
        aes(yintercept = K_median),
        color = "blue",
        linewidth = 0.7,
        inherit.aes = FALSE
      ) +
      geom_hline(
        data = diff_dt,
        aes(yintercept = K_weighted),
        color = "darkorange",
        linewidth = 0.7,
        inherit.aes = FALSE
      )
  }
  
  ggsave(
    filename = paste0("segmentation_output/T1_vs_T2_plots/total_segment_group_num_per_chromosome_", i, ".png"),
    plot = p,
    width = 14,
    height = 8,
    dpi = 300
  )
}

##########################################################################
#############################
#PLOT 2
#plot number of total segments per sample per chromosome

count_bed_segments <- function(file) {
  nrow(fread(file))
}

###
# get num fo segments
base_dir <- "segmentation_output/per_sample"
samples <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)

all_dt <- rbindlist(lapply(samples, function(sample) {
  
  sample_dir <- file.path(base_dir, sample, "beds")
  bed_files <- list.files(sample_dir, pattern = "_segments\\.bed$", full.names = TRUE)
  
  rbindlist(lapply(bed_files, function(f) {
    
    chr <- gsub("_segments\\.bed", "", basename(f))
    
    data.table(
      sample = sample,
      chr = chr,
      burden = count_bed_segments(f)
    )
  }))
}))

# t1 vs t2

all_dt[, group := ifelse(grepl("T2", sample), "T2", "T1")]


chromosomes <- sort(unique(all_dt$chr))

median_counts <- data.table(
  chr = chromosomes,
  burden_median = sapply(chromosomes, function(chr) {
    count_bed_segments(paste0("segmentation_output/beds/", chr, "_segments.bed"))
  })
)

weighted_counts <- data.table(
  chr = chromosomes,
  burden_weighted = sapply(chromosomes, function(chr) {
    count_bed_segments(
      paste0(
        "segmentation_output_cov_weighted_consensus/beds/",
        chr,
        "_segments.bed"
      )
    )
  })
)

consensus_dt <- merge(median_counts, weighted_counts, by = "chr", all = TRUE)

dir.create(
  "segmentation_output/T1_vs_T2_plots",
  recursive = TRUE,
  showWarnings = FALSE
)


chroms <- sort(unique(all_dt$chr))
chunks <- split(chroms, cut(seq_along(chroms), 4, labels = FALSE))

#loop
for (i in seq_along(chunks)) {
  
  chr_chunk <- chunks[[i]]
  dt_plot <- all_dt[chr %in% chr_chunk]
  dt_plot[, group := factor(group)]
  
  p <- ggplot(
    dt_plot,
    aes(x = group, y = burden, color = group)
  ) +
    geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.35) +
    facet_wrap(~chr, ncol = 3, scales = "fixed") +
    

  scale_y_continuous(
    breaks = pretty(dt_plot$burden, n = 4)
  ) +
    
  theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      plot.title = element_text(size = 11),
      
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_line(color = "grey92")
    ) +
    
    labs(
      title = paste0(
        "Segment number per chromosome (block ",
        i,
        "/4)"
      ),
      x = NULL,
      y = "Number of segments (BED rows)"
    )
  
  sub_cons <- consensus_dt[chr %in% chr_chunk]
  
  p <- p +
    geom_hline(
      data = sub_cons,
      aes(yintercept = burden_median),
      color = "blue",
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    geom_hline(
      data = sub_cons,
      aes(yintercept = burden_weighted),
      color = "darkorange",
      linewidth = 0.7,
      inherit.aes = FALSE
    )
  

  ggsave(
    filename = paste0(
      "segmentation_output/T1_vs_T2_plots/total_segment_num_per_chromosome_",
      i,
      ".png"
    ),
    plot = p,
    width = 14,
    height = 8,
    dpi = 300
  )
}

##############################################################
# PLOT 3
# we want to plot the average methylation per sample group per chromsome

dir_path <- "segmentation_output/sample_summaries"

files <- list.files(
  dir_path,
  pattern = "_segment_group_summary\\.csv$",
  full.names = TRUE
)

df <- purrr::map_dfr(files, function(f) {
  
  sample_name <- str_extract(basename(f), "^[0-9]+_T[12]")
  
  readr::read_csv(f, show_col_types = FALSE) %>%
    mutate(
      sample = sample_name,
      condition = ifelse(str_detect(sample_name, "T1"), "T1", "T2")
    )
})


df <- df %>%
  group_by(sample, chromosome) %>%
  mutate(K = max(seg.group, na.rm = TRUE)) %>%
  ungroup()

# Ensure K is integer
df <- df %>%
  mutate(K = as.integer(K))


out_dir <- "segmentation_output/T1_vs_T2_plots/mean_methylation_per_chromosome"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chromosomes <- unique(df$chromosome)

for (chr in chromosomes) {
  
  df_chr <- df %>%
    filter(chromosome == chr) %>%
    
    # IMPORTANT: order K properly within chromosome
    mutate(K = factor(K, levels = sort(unique(K))))
  
  p <- ggplot(df_chr, aes(
    x = factor(seg.group),
    y = mean_meth,
    color = condition
  )) +
    
    geom_point(
      position = position_jitter(width = 0.15),
      size = 2,
      alpha = 0.8
    ) +
    
    stat_summary(
      aes(group = interaction(seg.group, condition)),
      fun = mean,
      geom = "line",
      linewidth = 0.8
    ) +
    
    facet_wrap(~ K, nrow = 1) +
    
    labs(
      title = paste0("Mean methylation per segment group - ", chr),
      x = "Segment group",
      y = "Mean methylation",
      color = "Condition (T1 / T2)"
    ) +
    
    theme_bw() +
    
    theme(
      strip.text = element_text(face = "bold")
    )

  ggsave(
    filename = file.path(out_dir, paste0(chr, ".png")),
    plot = p,
    width = 11,
    height = 6,
    dpi = 300
  )
}

#for some reason, 744_T2(NC_059465.1), 838_T2(NC_059469.1) and 826_T2(1st chr)
#have some issues with their segmentation