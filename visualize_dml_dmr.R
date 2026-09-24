library(karyoploteR)
library(GenomicRanges)
library(Rsamtools)
library(dplyr)

mk =read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_mk.csv")
path_fasta <- "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa"
fai <- scanFaIndex(path_fasta)
seqinfo <- seqinfo(fai)

mk_plot <- mk %>%
  mutate(
    chr = as.character(
      as.integer(
        sub("\\.1$", "",
            sub("^NC_059", "", chr)
        )
      ) - 441
    )
  ) %>%
  filter(chr %in% as.character(1:29))



######################
######################
#loop thru all the chrs, save the plots
library(ggplot2)
library(dplyr)

# --------------------------------------------------
# Chromosome lengths
# --------------------------------------------------

chr_lengths <- data.frame(
  chr = as.character(seqnames(fai)),
  length = width(fai)
) %>%
  filter(chr %in% as.character(1:29)) %>%
  mutate(chr_num = as.integer(chr))


mk_plot <- mk %>%
  mutate(
    chr_num = as.integer(
      sub(
        "\\.1$",
        "",
        sub("^NC_059", "", chr)
      )
    ) - 441
  ) %>%
  filter(chr_num %in% 1:29)

# each plot contains one grp

chromosome_groups <- list(
  1:7,
  8:14,
  15:21,
  22:29
)

# loop

for (i in seq_along(chromosome_groups)) {
  
  chromosomes <- chromosome_groups[[i]]
    chr_plot <- chr_lengths %>%
    filter(chr_num %in% chromosomes)
  
  mk_plot_i <- mk_plot %>%
    filter(chr_num %in% chromosomes)
    max_length <- max(chr_plot$length)
  
  # Create plot
  p <- ggplot() +
    
    geom_segment(
      data = chr_plot,
      aes(
        x = 0,
        xend = length,
        y = chr_num,
        yend = chr_num
      ),
      linewidth = 3,
      colour = "grey70",
      lineend = "round"
    ) +
    
    # Methylation positions
    #we might need to change start to pos etc depending on the csv file we are using

    geom_segment(
      data = mk_plot_i,
      aes(
        #x = start,
        #xend = start,
        x = pos,
        xend = pos,
        y = chr_num - 0.15,
        yend = chr_num + 0.15
      ),
      colour = "#D73027",
      linewidth = 0.7
    ) +
    
    # X axis
    scale_x_continuous(
      name = "Genomic position (Mb)",
      limits = c(0, max_length * 1.05),
      labels = function(x) x / 1e6,
      expand = c(0, 0)
    ) +
    
    # Y axis
    scale_y_continuous(
      name = "Chromosome",
      breaks = chromosomes,
      labels = chromosomes,
      limits = c(
        min(chromosomes) - 0.5,
        max(chromosomes) + 0.5
      )
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
  
  # Save
  ggsave(
    filename = paste0("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/visualizations/mk",
                      "salmon_chromosomes_",
                      min(chromosomes),
                      "_",
                      max(chromosomes),
                      ".pdf"
    ),
    plot = p,
    width = 10,
    height = 6
  )
  
  print(p)
}

####################################
####################################
#for the segments
segments=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_segs.csv")
segments_plot <- segments %>%
  mutate(
    chr_num = as.integer(
      sub(
        "\\.1$",
        "",
        sub("^NC_059", "", chr)
      )
    ) - 441
  )

chromosome <- 4
cluster_gap <- 5000
buffer_multiplier <- 3
min_buffer <- 2500

# Chromosome data
chr_length <- chr_lengths %>%
  filter(chr_num == chromosome) %>%
  pull(length)

segments_chr <- segments_plot %>%
  filter(chr_num == chromosome) %>%
  arrange(start, end)

# Find clusters
segments_chr <- segments_chr %>%
  mutate(
    gap = start - lag(end),
    new_cluster = if_else(
      is.na(gap) | gap > cluster_gap,
      1L,
      0L
    ),
    cluster = cumsum(new_cluster)
  )

# Cluster boundaries
clusters <- segments_chr %>%
  group_by(cluster) %>%
  summarise(
    cluster_start = min(start),
    cluster_end = max(end),
    .groups = "drop"
  ) %>%
  mutate(
    cluster_span = cluster_end - cluster_start,
    buffer = pmax(
      cluster_span * buffer_multiplier,
      min_buffer
    ),
    zoom_start = pmax(
      0,
      cluster_start - buffer
    ),
    zoom_end = pmin(
      chr_length,
      cluster_end + buffer
    )
  )



# Chromosome overview
p_overview <- ggplot() +
  geom_segment(
    aes(
      x = 0,
      xend = chr_length,
      y = 1,
      yend = 1
    ),
    linewidth = 5,
    colour = "grey75",
    lineend = "round"
  ) +
  geom_segment(
    data = segments_chr,
    aes(
      x = start,
      xend = end,
      y = 1,
      yend = 1
    ),
    colour = "#D73027",
    linewidth = 5,
    lineend = "butt"
  ) +
  geom_rect(
    data = clusters,
    aes(
      xmin = zoom_start,
      xmax = zoom_end,
      ymin = 0.72,
      ymax = 1.28
    ),
    fill = "#D73027",
    alpha = 0.12,
    colour = "#D73027"
  ) +
  geom_text(
    data = clusters,
    aes(
      x = (zoom_start + zoom_end) / 2,
      y = 1.4,
      label = paste0("Region ", cluster)
    ),
    size = 3
  ) +
  scale_x_continuous(
    name = "Genomic position (Mb)",
    labels = function(x) x / 1e6,
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    breaks = 1,
    labels = paste("Chromosome", chromosome),
    limits = c(0.5, 1.7)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 10),
    
  )

# Zoom plots
zoom_plots <- vector("list", nrow(clusters))

for (i in seq_len(nrow(clusters))) {
  
  region <- clusters[i, ]
  
  seg_i <- segments_chr %>%
    filter(
      end >= region$zoom_start,
      start <= region$zoom_end
    )
  
  zoom_plots[[i]] <- ggplot() +
    geom_segment(
      aes(
        x = region$zoom_start,
        xend = region$zoom_end,
        y = 1,
        yend = 1
      ),
      linewidth = 3,
      colour = "grey75"
    ) +
    geom_rect(
      data = seg_i,
      aes(
        xmin = start,
        xmax = end,
        ymin = 0.65,
        ymax = 1.35
      ),
      fill = "#D73027",
      colour = "#A50026",
      linewidth = 0.4
    ) +
    scale_x_continuous(
      name = "Genomic position (kb)",
      limits = c(
        region$zoom_start,
        region$zoom_end
      ),
      labels = function(x) round(x / 1000, 1),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = 1,
      labels = paste0("Region ", region$cluster),
      limits = c(0.3, 1.7)
    ) +
    theme_classic(base_size = 13) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      axis.title.x = element_text(size = 9),
      axis.text.x = element_text(size = 7)
      
    )
}

# Combine
p_final <- p_overview

for (i in seq_along(zoom_plots)) {
  p_final <- p_final / zoom_plots[[i]]
}

p_final <- wrap_plots(
  c(list(p_overview), zoom_plots),
  ncol = 1,
  heights = c(1.2, rep(1, length(zoom_plots)))
)

p_final
