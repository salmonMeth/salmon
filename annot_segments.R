#we want to see if we get anything meaningful if we annotate the segments

library(data.table)
library(GenomicRanges)
library(IRanges)
library(GenomicFeatures)
library(ggplot2)

#
path_annot= "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
bed_dir <- "segmentation_output/beds"
##
files <- list.files(bed_dir, full.names = TRUE)

seg_list <- lapply(files, function(f) {
  fread(f, header = FALSE)
})

names(seg_list) <- gsub(".*/|\\.bed", "", files)

seg_list <- lapply(seg_list, function(df) {
  
  df <- df[, 1:5]
  
  setnames(df,
           c("chr", "start", "end", "seg_group", "mean_meth"))
  
  df$chr <- as.character(df$chr)
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  df$seg_group <- as.integer(df$seg_group)
  df$mean_meth <- as.numeric(df$mean_meth)
  
  df
})


seg_gr_list <- lapply(seg_list, function(df) {
  
  GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end),
    seg_group = df$seg_group,
    mean_meth = df$mean_meth
  )
})

#check

length(seg_gr_list)                  
length(seg_gr_list[[1]])            
head(seg_gr_list[[1]])
seqnames(seg_gr_list[[1]])
width(seg_gr_list[[1]])[1:5]

###
#so we try it on just one chromosome at first

seg_df <- seg_list[[1]]
seg_gr <- GRanges(
  seqnames = seg_df$chr,
  ranges = IRanges(
    start = seg_df$start,
    end = seg_df$end
  ),
  seg_group = seg_df$seg_group,
  mean_meth = seg_df$mean_meth
)


txdb <- txdbmaker::makeTxDbFromGFF(
  path_annot)
genes_gr <- genes(txdb)
tx_gr <- transcripts(txdb)
promoters_gr <- promoters(txdb, upstream = 1000, downstream = 200)

promoters_chr1 <- promoters_gr[
  as.character(seqnames(promoters_gr)) == "1"
]

length(promoters_chr1)

hits <- findOverlaps(
  seg_gr,
  promoters_chr1
)
#create the df
df <- as.data.frame(seg_gr)

df$is_promoter <- FALSE
df$gene_id <- NA_character_

df$is_promoter[queryHits(hits)] <- TRUE
gene_list <- split(
  mcols(promoters_chr1)$tx_name[subjectHits(hits)],
  queryHits(hits)
)

gene_strings <- sapply(
  gene_list,
  function(x) paste(unique(x), collapse = ";")
)

df$tx_name <- NA_character_

df$tx_name[as.integer(names(gene_strings))] <- gene_strings

######################################################
#loop over all chromosomes


outdir <- "segment_annotations"
dir.create(outdir, showWarnings = FALSE)

for (f in files) {
  
  cat("Processing:", basename(f), "\n")
  
  seg <- fread(f)
  
  colnames(seg)[1:5] <- c(
    "chr",
    "start",
    "end",
    "seg_group",
    "mean_meth"
  )
  
  seg_gr <- GRanges(
    seqnames = as.character(seg$chr),
    ranges = IRanges(
      start = seg$start,
      end = seg$end
    ),
    seg_group = seg$seg_group,
    mean_meth = seg$mean_meth
  )
  
  
  gene_hits <- findOverlaps(
    seg_gr,
    genes_gr,
    ignore.strand = TRUE
  )
  
  n_genes <- tabulate(
    queryHits(gene_hits),
    nbins = length(seg_gr)
  )
  
  gene_bp <- integer(length(seg_gr))
  
  if(length(gene_hits) > 0){
    
    ov_gene <- pintersect(
      seg_gr[queryHits(gene_hits)],
      genes_gr[subjectHits(gene_hits)]
    )
    
    tmp <- tapply(
      width(ov_gene),
      queryHits(gene_hits),
      sum
    )
    
    gene_bp[as.integer(names(tmp))] <- tmp
  }
  
  #promoters
  
  promoter_hits <- findOverlaps(
    seg_gr,
    promoters_gr,
    ignore.strand = TRUE
  )
  
  n_promoters <- tabulate(
    queryHits(promoter_hits),
    nbins = length(seg_gr)
  )
  
  promoter_bp <- integer(length(seg_gr))
  
  if(length(promoter_hits) > 0){
    
    ov_prom <- pintersect(
      seg_gr[queryHits(promoter_hits)],
      promoters_gr[subjectHits(promoter_hits)]
    )
    
    tmp <- tapply(
      width(ov_prom),
      queryHits(promoter_hits),
      sum
    )
    
    promoter_bp[as.integer(names(tmp))] <- tmp
  }
  #genes
  
  gene_ids <- rep("", length(seg_gr))
  
  if(length(gene_hits) > 0){
    
    gene_list <- split(
      genes_gr$gene_id[subjectHits(gene_hits)],
      queryHits(gene_hits)
    )
    
    gene_ids[as.integer(names(gene_list))] <-
      sapply(
        gene_list,
        function(x)
          paste(unique(x), collapse = ";")
      )
  }
    # gene names
  
  gene_names <- rep("", length(seg_gr))
  
  if("Name" %in% colnames(mcols(genes_gr))){
    
    gene_name_list <- split(
      genes_gr$Name[subjectHits(gene_hits)],
      queryHits(gene_hits)
    )
    
    gene_names[as.integer(names(gene_name_list))] <-
      sapply(
        gene_name_list,
        function(x)
          paste(unique(na.omit(x)), collapse = ";")
      )
  }
  
  ##get the info into a df
  
  df <- data.frame(
    seqnames = as.character(seqnames(seg_gr)),
    start = start(seg_gr),
    end = end(seg_gr),
    width = width(seg_gr),
    seg_group = seg$seg_group,
    mean_meth = seg$mean_meth,
    
    n_genes = n_genes,
    gene_bp = gene_bp,
    prop_gene_bp = gene_bp / width(seg_gr),
    genes_per_mb = n_genes / (width(seg_gr) / 1e6),
    
    n_promoters = n_promoters,
    promoter_bp = promoter_bp,
    prop_promoter_bp = promoter_bp / width(seg_gr),
    
    gene_ids = gene_ids,
    gene_names = gene_names
  )
  
  outfile <- file.path(
    outdir,
    paste0(
      sub("\\.bed$", "", basename(f)),
      "_annotated.csv"
    )
  )
  
  write.csv(
    df,
    file = outfile,
    row.names = FALSE,
    quote = TRUE,
    na = "NA"
  )
  
  cat("Saved:", outfile, "\n")
}


#we need to make use of this info
#lets get the toal number of segment groups every chromosome has from the bed files

bed_dir <- "segmentation_output/beds"
files <- list.files(bed_dir, full.names = TRUE)
files <- files[grepl("\\.bed$", files)]

get_chr_name <- function(f) {
  gsub(".*/|_segments\\.bed$", "", f)
}

chr_names <- get_chr_name(files)

K_values <- numeric(length(files))

for (i in seq_along(files)) {
  dt <- fread(files[i], header = FALSE)
  K_values[i] <- max(dt[[4]], na.rm = TRUE)
}

chr_K <- data.table(
  chr = chr_names,
  K = K_values
)

setorder(chr_K, chr)

chr_K
#trying it for k=2 first

chr_k2 <- chr_K[K == 2, chr]
#
in_dir  <- "segment_annotations"
out_dir <- "group_level_stats_k2"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# chromosomes with exactly 2 segment groups
files_k2 <- file.path(in_dir, paste0(chr_k2, "_segments_annotated.csv"))

seg_list <- lapply(files_k2, function(f) {
  dt <- fread(f)
  
  dt[, chr := gsub("_segments_annotated.csv", "", basename(f))]
  
  # safety: ensure numeric
  dt[, seg_group := as.integer(seg_group)]
  
  dt
})

seg_k2 <- rbindlist(seg_list, fill = TRUE)

seg_k2[, seg_group := as.integer(seg_group)]group_stats_k2 <- seg_k2[
  ,
  .(
    n_segments = .N,
    
    mean_meth = mean(mean_meth, na.rm = TRUE),
    sd_meth   = sd(mean_meth, na.rm = TRUE),
    
    mean_genes = mean(n_genes, na.rm = TRUE),
    sd_genes   = sd(n_genes, na.rm = TRUE),
    
    mean_gene_bp = mean(gene_bp, na.rm = TRUE),
    sd_gene_bp   = sd(gene_bp, na.rm = TRUE),
    
    mean_gene_density = mean(genes_per_mb, na.rm = TRUE),
    sd_gene_density   = sd(genes_per_mb, na.rm = TRUE),
    
    mean_promoters = mean(n_promoters, na.rm = TRUE),
    sd_promoters   = sd(n_promoters, na.rm = TRUE),
    
    mean_promoter_bp = mean(promoter_bp, na.rm = TRUE),
    sd_promoter_bp   = sd(promoter_bp, na.rm = TRUE)
  ),
  by = seg_group
]
group_stats_k2[, K := 2]
setcolorder(group_stats_k2, c("K", "seg_group"))

#plots
seg_k2[, seg_group := as.factor(seg_group)]
ggplot(seg_k2, aes(x = seg_group, y = genes_per_mb)) +
  geom_boxplot(fill = "darkgreen") +
  theme_minimal() +
  labs(
    title = "Gene Density by Segment Group (K=2)",
    x = "Segment group",
    y = "Genes per Mb"
  )
ggplot(seg_k2, aes(x = seg_group, y = mean_meth)) +
  geom_boxplot(fill = "steelblue") +
  theme_minimal() +
  labs(
    title = "Mean Methylation by Segment Group (K=2)",
    x = "Segment group",
    y = "Mean methylation"
  )


ggplot(seg_k2, aes(x = seg_group, y = width)) +
  geom_boxplot(fill = "darkorange") +
  scale_y_log10() +  # important: lengths are highly skewed
  theme_minimal() +
  labs(
    title = "Segment Length by Segment Group (K=2)",
    x = "Segment group",
    y = "Segment width (log10 scale)"
  )

#save 
out_dir <- "seg_k2_outputs"
dir.create(out_dir, showWarnings = FALSE)
ggsave(file.path(out_dir, "k2_methylation_boxplot.png"), p_meth, width = 6, height = 4)
ggsave(file.path(out_dir, "k2_gene_density_boxplot.png"), p_gene, width = 6, height = 4)
ggsave(file.path(out_dir, "k2_segment_length_boxplot.png"), p_len, width = 6, height = 4)

fwrite(group_stats_k2,
       file.path(out_dir, "group_stats_k2.csv"))



########
#and the following is the loop


in_dir  <- "segment_annotations"
base_out <- "group_level_stats_by_K"

dir.create(base_out, showWarnings = FALSE, recursive = TRUE)

K_values_unique <- sort(unique(chr_K$K))

for (k in K_values_unique) {
  
  message("Processing K = ", k)
  
  chr_k <- chr_K[K == k, chr]
  
  files_k <- file.path(
    in_dir,
    paste0(chr_k, "_segments_annotated.csv")
  )

  seg_list <- lapply(files_k, function(f) {
    dt <- fread(f)
    
    dt[, chr := gsub("_segments_annotated.csv", "", basename(f))]
    dt[, seg_group := as.integer(seg_group)]
    
    dt
  })
  
  seg_k <- rbindlist(seg_list, fill = TRUE)
  seg_k[, seg_group := as.factor(seg_group)]

  group_stats_k <- seg_k[
    ,
    .(
      n_segments = .N,
      
      mean_meth = mean(mean_meth, na.rm = TRUE),
      sd_meth   = sd(mean_meth, na.rm = TRUE),
      
      mean_genes = mean(n_genes, na.rm = TRUE),
      sd_genes   = sd(n_genes, na.rm = TRUE),
      
      mean_gene_bp = mean(gene_bp, na.rm = TRUE),
      sd_gene_bp   = sd(gene_bp, na.rm = TRUE),
      
      mean_gene_density = mean(genes_per_mb, na.rm = TRUE),
      sd_gene_density   = sd(genes_per_mb, na.rm = TRUE),
      
      mean_promoters = mean(n_promoters, na.rm = TRUE),
      sd_promoters   = sd(n_promoters, na.rm = TRUE),
      
      mean_promoter_bp = mean(promoter_bp, na.rm = TRUE),
      sd_promoter_bp   = sd(promoter_bp, na.rm = TRUE)
    ),
    by = seg_group
  ]
  
  group_stats_k[, K := k]
  setcolorder(group_stats_k, c("K", "seg_group"))

  out_dir <- file.path(base_out, paste0("K_", k))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  p_meth <- ggplot(seg_k, aes(x = seg_group, y = mean_meth)) +
    geom_boxplot(fill = "steelblue") +
    theme_minimal() +
    labs(
      title = paste("Mean Methylation by Segment Group (K=", k, ")", sep=""),
      x = "Segment group",
      y = "Mean methylation"
    )
  
  p_gene <- ggplot(seg_k, aes(x = seg_group, y = genes_per_mb)) +
    geom_boxplot(fill = "darkgreen") +
    theme_minimal() +
    labs(
      title = paste("Gene Density by Segment Group (K=", k, ")", sep=""),
      x = "Segment group",
      y = "Genes per Mb"
    )
  
  p_len <- ggplot(seg_k, aes(x = seg_group, y = width)) +
    geom_boxplot(fill = "darkorange") +
    scale_y_log10() +
    theme_minimal() +
    labs(
      title = paste("Segment Length by Segment Group (K=", k, ")", sep=""),
      x = "Segment group",
      y = "Width (log10)"
    )

  ggsave(file.path(out_dir, paste0("K", k, "_methylation_boxplot.png")), p_meth, width = 6, height = 4)
  ggsave(file.path(out_dir, paste0("K", k, "_gene_density_boxplot.png")), p_gene, width = 6, height = 4)
  ggsave(file.path(out_dir, paste0("K", k, "_segment_length_boxplot.png")), p_len, width = 6, height = 4)

  fwrite(
    group_stats_k,
    file.path(out_dir, paste0("group_stats_K", k, ".csv"))
  )
  
}

fwrite(chr_K, "group_level_stats_by_K/chromosome_K_values.csv")
