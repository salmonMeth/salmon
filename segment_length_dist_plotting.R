base <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation"


dirs <- c(
  original = file.path(base, "segmentation_output_G1_5/beds"),
  minseg20 = file.path(base, "segmentation_output_G1_5_minseg_20/beds"),
  maxint40 = file.path(base, "segmentation_output_G1_5_maxInt_40/beds")
)

library(data.table)
library(ggplot2)

seg_lengths <- rbindlist(lapply(names(dirs), function(method) {
  files <- list.files(
    dirs[[method]],
    pattern = "^NC_[0-9]+\\.1_segments\\.bed$",
    full.names = TRUE
  )
  
  rbindlist(lapply(files, function(f) {
    bed <- fread(f, header = FALSE)
    
    data.table(
      chromosome = as.integer(
        sub(
          "^NC_",
          "",
          sub("\\.1_segments\\.bed$", "", basename(f))
        )
      ) - 59441,
      start = bed[[2]],
      end = bed[[3]],
      length = bed[[3]] - bed[[2]] + 1,
      method = method
    )
  }))
}))

#plot segment length distribution
p=ggplot(seg_lengths, aes(x = length, colour = method)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Segment length (bp)",
    y = "Density",
    colour = "Segmentation"
  )

ggsave(
  "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_length_distribution.png",
  p,
  width = 8,
  height = 6,
  dpi = 300
)
#####################################
####################################
#comparing the legtnhs of segments where total number of segment groups is 2

library(data.table)
library(ggplot2)

seg_data <- rbindlist(lapply(names(dirs), function(method) {
  files <- list.files(
    dirs[[method]],
    pattern = "^NC_[0-9]+\\.1_segments\\.bed$",
    full.names = TRUE
  )
  
  rbindlist(lapply(files, function(f) {
    bed <- fread(f, header = FALSE)
    
    data.table(
      chromosome = as.integer(
        sub(
          "^NC_",
          "",
          sub("\\.1_segments\\.bed$", "", basename(f))
        )
      ) - 59441,
      start = bed[[2]],
      end = bed[[3]],
      length = bed[[3]] - bed[[2]],
      seg_group = bed[[4]],
      mean_meth = bed[[5]],
      method = method
    )
  }))
}))

two_group_chr <- seg_data[
  ,
  .(n_groups = uniqueN(seg_group)),
  by = .(method, chromosome)
][n_groups == 2]

two_group_data <- seg_data[
  two_group_chr,
  on = .(method, chromosome),
  nomatch = 0
]
p <- ggplot(
  two_group_data,
  aes(
    x = methylation_group,
    y = length,
    fill = methylation_group
  )
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  facet_wrap(~method) +
  theme_bw() +
  labs(
    x = "Methylation group",
    y = "Segment length (bp)",
    fill = "Methylation"
  )

ggsave(
  "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_length_high_low_comparison.png",
  p,
  width = 8,
  height = 6,
  dpi = 300
)


#############
############
#for 3 levels

three_group_chr <- seg_data[
  ,
  .(n_groups = uniqueN(seg_group)),
  by = .(method, chromosome)
][n_groups == 3]

three_group_data <- seg_data[
  three_group_chr,
  on = .(method, chromosome),
  nomatch = 0
]

group_levels <- three_group_data[
  ,
  .(group_mean = mean(mean_meth)),
  by = .(method, chromosome, seg_group)
]

group_levels[
  ,
  methylation_group := c("Low", "Medium", "High")[
    order(order(group_mean))
  ],
  by = .(method, chromosome)
]

three_group_data <- group_levels[
  three_group_data,
  on = .(method, chromosome, seg_group)
]

three_group_data[, methylation_group := factor(
  methylation_group,
  levels = c("High", "Medium", "Low")
)]
p <- ggplot(
  three_group_data,
  aes(
    x = methylation_group,
    y = length,
    fill = methylation_group
  )
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.7) +
  scale_y_log10() +
  scale_fill_manual(
    values = c(
      "High" = "#D73027",
      "Medium" = "#FC8D59",
      "Low" = "#4575B4"
    )
  ) +
  facet_wrap(~method) +
  theme_bw() +
  labs(
    x = "Methylation group",
    y = "Segment length (bp)",
    fill = "Methylation"
  )

p



###########################
###########################
#methylation stats
library(data.table)
two_group_chr <- seg_data[
  ,
  .(n_groups = uniqueN(seg_group)),
  by = .(method, chromosome)
][n_groups == 2]

two_group_data <- seg_data[
  two_group_chr,
  on = .(method, chromosome),
  nomatch = 0
]
group_levels_2 <- two_group_data[
  ,
  .(group_mean = mean(mean_meth)),
  by = .(method, chromosome, seg_group)
]
group_levels_2[
  ,
  methylation_group := c("Low", "High")[
    order(order(group_mean))
  ],
  by = .(method, chromosome)
]
two_group_data <- group_levels_2[
  two_group_data,
  on = .(method, chromosome, seg_group)
]
two_group_data[, methylation_group := factor(
  methylation_group,
  levels = c("High", "Low")
)]

p_meth <- ggplot(
  two_group_data,
  aes(
    x = methylation_group,
    y = mean_meth,
    fill = methylation_group
  )
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "High" = "#D73027",
      "Low" = "#4575B4"
    )
  ) +
  facet_wrap(~method) +
  theme_bw() +
  labs(
    x = "Methylation group",
    y = "Mean methylation",
    fill = "Methylation"
  )

p_meth


#################
################
#3 groups

three_group_chr <- seg_data[
  ,
  .(n_groups = uniqueN(seg_group)),
  by = .(method, chromosome)
][n_groups == 3]

three_group_data <- seg_data[
  three_group_chr,
  on = .(method, chromosome),
  nomatch = 0
]

group_levels <- three_group_data[
  ,
  .(group_mean = mean(mean_meth)),
  by = .(method, chromosome, seg_group)
]

group_levels[
  ,
  methylation_group := c("Low", "Medium", "High")[
    order(order(group_mean))
  ],
  by = .(method, chromosome)
]

three_group_data <- group_levels[
  three_group_data,
  on = .(method, chromosome, seg_group)
]

three_group_data[, methylation_group := factor(
  methylation_group,
  levels = c("High", "Medium", "Low")
)]

p_meth_3 <- ggplot(
  three_group_data,
  aes(
    x = methylation_group,
    y = mean_meth,
    fill = methylation_group
  )
) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.7) +
  scale_fill_manual(
    values = c(
      "High" = "#D73027",
      "Medium" = "#FC8D59",
      "Low" = "#4575B4"
    )
  ) +
  facet_wrap(~method) +
  theme_bw() +
  labs(
    x = "Methylation group",
    y = "Mean methylation",
    fill = "Methylation"
  )

p_meth_3
ggsave(
  "/scratch/project_2010912/ezel/salmon/data_files/segmentation/methylation_level_violin_3_levels.png",
  p_meth_3,
  width = 8,
  height = 6,
  dpi = 300
)



#methylation stats
two_group_data[
  seg_group == 2,
  .(
    n_chromosomes = uniqueN(chromosome),
    mean_group_mean = mean(group_mean),
    median_group_mean = median(group_mean),
    sd_group_mean = sd(group_mean),
    min_group_mean = min(group_mean),
    max_group_mean = max(group_mean)
  ),
  by = method
]


#####for 3 methylation levls
three_group_chr <- seg_data[
  ,
  .(n_groups = uniqueN(seg_group)),
  by = .(method, chromosome)
][n_groups == 3]

three_group_data <- seg_data[
  three_group_chr,
  on = .(method, chromosome),
  nomatch = 0
]

group_levels <- three_group_data[
  ,
  .(group_mean = mean(mean_meth)),
  by = .(method, chromosome, seg_group)
]

group_levels[
  ,
  methylation_group := c("Low", "Medium", "High")[
    order(order(group_mean))
  ],
  by = .(method, chromosome)
]

three_group_data <- group_levels[
  three_group_data,
  on = .(method, chromosome, seg_group)
]

three_group_data[, methylation_group := factor(
  methylation_group,
  levels = c("High", "Medium", "Low")
)]

three_group_data[
  ,
  .(
    group_mean = unique(group_mean)
  ),
  by = .(method, chromosome, seg_group, methylation_group)
]

#segment length &mean methylatio plots

library(ggplot2)

p_corr <- ggplot(
  seg_data,
  aes(
    x = length,
    y = mean_meth
  )
) +
  geom_point(
    alpha = 0.25,
    size = 1
  ) +
  scale_x_log10() +
  facet_wrap(~method) +
  theme_bw() +
  labs(
    x = "Segment length (bp)",
    y = "Mean methylation"
  )

p_corr
