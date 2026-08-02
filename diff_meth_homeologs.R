#we want to check the homeologs of diff methylated bases
library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

mk_homeo=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk_with_homeos.csv")

mk_short <- mk_homeo[, c("seqnames","start","end","chr_ncbi","meth.diff", "gene","in_homeolog_block","homeolog_block_name",
  "partner_chr","partner_start","partner_end")]
blocks <- unique(na.omit(mk_short$homeolog_block_name))

paired_df <- map_dfr(blocks, function(block) {
  
  primary <- mk_short %>%
    filter(homeolog_block_name == block)
  
  partner <- mk_short %>%
    filter(
      chr_ncbi == primary$partner_chr[1],
      start >= primary$partner_start[1],
      start <= primary$partner_end[1]
    )
  
  primary$region <- "primary"
  
  if (nrow(partner) > 0) {
    partner$region <- "partner"
    bind_rows(primary, partner)
  } else {
    primary
  }
})

paired_df <- paired_df %>%
  arrange(homeolog_block_name, region, start)

  
write.csv(
    paired_df,
    "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_diff_meth_paired.csv",
    row.names = FALSE)
  
  
  ###stats 
library(dplyr)
  
library(stringr)

block_summary <- paired_df %>%
  group_by(homeolog_block_name) %>%
  summarise(
    primary_hits = sum(region == "primary"),
    partner_hits = sum(region == "partner"),
    
    primary_chr = first(chr_ncbi[region == "primary"]),
    primary_start = if (any(region == "primary")) min(start[region == "primary"]) else NA_integer_,
    primary_end = if (any(region == "primary")) max(end[region == "primary"]) else NA_integer_,
    
    partner_chr = first(partner_chr[region == "primary"]),
    partner_start = first(partner_start[region == "primary"]),
    partner_end = first(partner_end[region == "primary"]),
    
    .groups = "drop")

# Blocks with hits on both sides
matched_blocks <- block_summary %>%
  filter(primary_hits > 0, partner_hits > 0)

# Blocks with hits on only one side
unmatched_blocks <- block_summary %>%
  filter(primary_hits == 0 | partner_hits == 0)


matched_blocks <- matched_blocks %>%
  mutate(block_order = as.numeric(str_extract(homeolog_block_name, "^[0-9]+"))) %>%
  arrange(block_order, homeolog_block_name) %>%
  select(-block_order)

unmatched_blocks <- unmatched_blocks %>%
  mutate(block_order = as.numeric(str_extract(homeolog_block_name, "^[0-9]+"))) %>%
  arrange(block_order, homeolog_block_name) %>%
  select(-block_order)

# Save
write.csv(
  matched_blocks,
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_homeolog_blocks_matched.csv",
  row.names = FALSE
)

write.csv(
  unmatched_blocks,
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_homeolog_blocks_unmatched.csv",
  row.names = FALSE
)

  
  
###distribution of cpg&blocks
dss_short=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/dss/dss_homeolog_blocks_matched.csv")

mk_short=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk/mk_homeolog_blocks_matched.csv")
library(tidyr)

block_compare <- full_join(
  mk_short %>%
    select(homeolog_block_name,
           mk_primary = primary_hits,
           mk_partner = partner_hits),
  dss_short %>%
    select(homeolog_block_name,
           dss_primary = primary_hits,
           dss_partner = partner_hits),
  by = "homeolog_block_name"
) %>%
  mutate(
    across(
      c(mk_primary, mk_partner, dss_primary, dss_partner),
      ~replace_na(.x, 0)
    )
  )


##plot


common_blocks <- block_compare %>%
  filter(
    homeolog_block_name %in% mk_short$homeolog_block_name,
    homeolog_block_name %in% dss_short$homeolog_block_name
  )

plot_df <- common_blocks %>%
  select(
    homeolog_block_name,
    mk_primary,
    mk_partner,
    dss_primary,
    dss_partner
  ) %>%
  pivot_longer(
    cols = c(mk_primary, mk_partner, dss_primary, dss_partner),
    names_to = c("method", "region"),
    names_pattern = "(mk|dss)_(primary|partner)",
    values_to = "n_cpgs"
  ) %>%
  mutate(
    method = recode(
      method,
      mk = "MethylKit",
      dss = "DSS"
    ),
    region = factor(
      region,
      levels = c("primary", "partner")
    )
  )


plot_df <- plot_df %>%
  mutate(
    homeolog_block_name = reorder(
      homeolog_block_name,
      n_cpgs,
      FUN = sum
    )
  )

p=ggplot(
  plot_df,
  aes(
    x = homeolog_block_name,
    y = n_cpgs,
    color = region,
    shape = method
  )
) +
  geom_point(
    size = 2.2,
    alpha = 0.71,
    position = position_jitter(
      width = 0.04,
      height = 0
    )
  ) +
  scale_color_manual(
    values = c(
      primary = "#2563EB",   # vibrant blue
      partner = "#DB2777"    # vibrant pink
    ),
    labels = c(
      primary = "Primary block",
      partner = "Partner block"
    )
  ) +
  scale_shape_manual(
    values = c(
      MethylKit = 16,  # filled circle
      DSS = 4          # cross
    )
  ) +
  scale_y_log10() +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    ),
    panel.grid.major = element_line(
      colour = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor = element_line(
      colour = "grey92",
      linewidth = 0.2
    ),
    legend.position = "right"
  ) +
  labs(
    x = "Homeolog block",
    y = "Number of significant CpGs",
    color = "Block type",
    shape = "Method"
  )

ggsave(
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/MK_vs_DSS_homeolog_block_CpG_counts.png",
  ,plot=p,width = 12,
  height = 6,
  dpi = 300
)

###############some stats