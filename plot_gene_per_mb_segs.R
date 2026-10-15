#plot exon or gene per mb of segments
#using /scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/segment_per_mb_gene_exon_info.csv

path="/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/segment_per_mb_gene_exon_info.csv"

library(ggplot2)

# Read data
df <- read.csv(path)

# Make seg_group a factor for categorical coloring
df$seg_group <- factor(df$seg_group)

#remove last 2 chrs
df_filt <- df[!df$chr %in% c(28, 29), ]


#genes
p_genes <- ggplot(
  df,
  aes(
    x = seg_mean,
    y = genes_per_mb,
    color = seg_group
  )
) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    x = "Segment mean methylation",
    y = "Genes per Mb",
    color = "Segment group"
  ) +
  theme_classic()

p_genes

#exons
# -------------------------
# Exons per Mb
# -------------------------

p_exons <- ggplot(
  df,
  aes(
    x = seg_mean,
    y = exons_per_mb,
    color = seg_group
  )
) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    x = "Segment mean methylation",
    y = "Exons per Mb",
    color = "Segment group"
  ) +
  theme_classic()

p_exons