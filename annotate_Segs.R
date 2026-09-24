library(data.table)


bed_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/beds"

out_file <- file.path(
  bed_dir,
  "segment_functional_annotation.csv"
)

threshold <- 0.10
###get the gff
gff_file <- "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
gff <- fread(
  gff_file,
  sep = "\t",
  header = FALSE,
  comment.char = "#"
)

setnames(
  gff,
  c(
    "chr",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
  )
)

# Keep chromosomes 1-29 only
gff <- gff[
  chr %in% as.character(1:29)
]

gff[, chr := as.character(chr)]
gff[, start := as.numeric(start)]
gff[, end := as.numeric(end)]


genes <- gff[
  type == "gene",
  .(
    chr,
    start = start - 1,
    end,
    score = 0,
    strand
  )
]

cat("Genes:", nrow(genes), "\n")

exons <- gff[
  type == "exon",
  .(
    chr,
    start = start - 1,
    end,
    score = 0,
    strand
  )
]

cat("Exons:", nrow(exons), "\n")

cds <- gff[
  type == "CDS",
  .(
    chr,
    start = start - 1,
    end,
    score = 0,
    strand
  )
]

cat("CDS:", nrow(cds), "\n")

gene_gff <- gff[type == "gene"]

gene_gff[, tss := ifelse(strand == "+", start, end)]

# Create ±1 kb promoter around TSS
promoters <- gene_gff[
  ,
  .(
    chr = chr,
    start = pmax(1, tss - 1000) - 1,
    end = tss + 1000,
    score = 0,
    strand = strand
  )
]

functional_regions <- list(
  promoter = promoters,
  gene     = genes,
  exon     = exons,
  CDS      = cds
)

functional_regions <- lapply(
  functional_regions,
  function(x) {
    
    x <- as.data.table(x)
    
    x[, chr := as.character(chr)]
    x[, start := as.numeric(start)]
    x[, end := as.numeric(end)]
    
    x[, width := end - start]
    
    x
  }
)



count_boundary_associated <- function(
    seg_chr,
    seg_start,
    seg_end,
    annotation,
    threshold = 0.10
) {
  
  x <- annotation[
    chr == as.character(seg_chr)
  ]
  
  if (nrow(x) == 0) {
    return(0L)
  }
  
  start_tol <- x$width * threshold
  end_tol   <- x$width * threshold
  
  associated <- (
    seg_start >= x$start - start_tol &
      seg_start <= x$start + start_tol &
      seg_end >= x$end - end_tol &
      seg_end <= x$end + end_tol
  )
  
  sum(associated, na.rm = TRUE)
}

count_contained <- function(
    seg_chr,
    seg_start,
    seg_end,
    annotation
) {
  
  x <- annotation[
    chr == as.character(seg_chr)
  ]
  
  if (nrow(x) == 0) {
    return(0L)
  }
  
  contained <- (
    x$start >= seg_start &
      x$end <= seg_end
  )
  
  sum(contained, na.rm = TRUE)
}

bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)


results <- vector(
  "list",
  length(bed_files)
)

for (i in seq_along(bed_files)) {
  
  f <- bed_files[i]
  
  cat(
    "Processing",
    basename(f),
    "\n"
  )
  
  lines <- readLines(f)
  
  if (
    length(lines) > 0 &&
    grepl("^track", lines[1])
  ) {
    
    bed <- fread(
      text = lines[-1],
      header = FALSE
    )
    
  } else {
    
    bed <- fread(
      text = lines,
      header = FALSE
    )
  }
  
  
  # --------------------------------------------------------
  # CHECK BED
  # --------------------------------------------------------
  
  if (ncol(bed) < 9) {
    
    stop(
      "Unexpected number of columns in ",
      basename(f),
      ": ",
      ncol(bed)
    )
  }
  
  
  setnames(
    bed,
    1:9,
    c(
      "chr",
      "start",
      "end",
      "seg_group",
      "seg_mean",
      "num_mark",
      "start_row",
      "end_row",
      "rgb"
    )
  )
  
  
  bed[, chr := as.character(chr)]
  bed[, start := as.numeric(start)]
  bed[, end := as.numeric(end)]
  
  results[[i]] <- bed[
    ,
    {
      
      promoter_boundary_n <- count_boundary_associated(
        chr,
        start,
        end,
        functional_regions$promoter,
        threshold
      )
      
      gene_boundary_n <- count_boundary_associated(
        chr,
        start,
        end,
        functional_regions$gene,
        threshold
      )
      
      exon_boundary_n <- count_boundary_associated(
        chr,
        start,
        end,
        functional_regions$exon,
        threshold
      )
      
      CDS_boundary_n <- count_boundary_associated(
        chr,
        start,
        end,
        functional_regions$CDS,
        threshold
      )
      
      
      # -----------------------------------------------
      # Completely contained
      # -----------------------------------------------
      
      promoter_contained_n <- count_contained(
        chr,
        start,
        end,
        functional_regions$promoter
      )
      
      gene_contained_n <- count_contained(
        chr,
        start,
        end,
        functional_regions$gene
      )
      
      exon_contained_n <- count_contained(
        chr,
        start,
        end,
        functional_regions$exon
      )
      
      CDS_contained_n <- count_contained(
        chr,
        start,
        end,
        functional_regions$CDS
      )
      
      
      list(
        
        # Boundary associations
        promoter_boundary_n = promoter_boundary_n,
        gene_boundary_n     = gene_boundary_n,
        exon_boundary_n     = exon_boundary_n,
        CDS_boundary_n      = CDS_boundary_n,
        
        # Regions contained inside segment
        promoter_contained_n = promoter_contained_n,
        gene_contained_n     = gene_contained_n,
        exon_contained_n     = exon_contained_n,
        CDS_contained_n      = CDS_contained_n
      )
    },
    by = .(
      chr,
      start,
      end,
      seg_group,
      seg_mean,
      num_mark,
      start_row,
      end_row,
      rgb
    )
  ]
  
  
  results[[i]][
    ,
    bed_file := basename(f)
  ]
}

annotation_results <- rbindlist(
  results,
  use.names = TRUE,
  fill = TRUE
)

annotation_results[
  ,
  promoter_boundary := promoter_boundary_n > 0
]

annotation_results[
  ,
  gene_boundary := gene_boundary_n > 0
]

annotation_results[
  ,
  exon_boundary := exon_boundary_n > 0
]

annotation_results[
  ,
  CDS_boundary := CDS_boundary_n > 0
]


annotation_results[
  ,
  promoter_contained := promoter_contained_n > 0
]

annotation_results[
  ,
  gene_contained := gene_contained_n > 0
]

annotation_results[
  ,
  exon_contained := exon_contained_n > 0
]

annotation_results[
  ,
  CDS_contained := CDS_contained_n > 0
]

annotation_results[
  ,
  boundary_threshold := threshold
]


fwrite(
  annotation_results,
  out_file
)

cat(
  "\nSaved:",
  out_file,
  "\n"
)

cat(
  "Total segments:",
  nrow(annotation_results),
  "\n"
)

cols=c("chr","start","end",(colnames(annotation_results))[10:27])

short <- annotation_results[, ..cols]

#plot the number of genes in the segments

gene_counts <- table(short$gene_contained_n)

# -------------------------
# Plot 1: 0–10 genes
# -------------------------
counts_0_10 <- gene_counts[
  as.numeric(names(gene_counts)) >= 0 &
    as.numeric(names(gene_counts)) <= 2
]

barplot(
  counts_0_10,
  xlab = "Number of genes in segment",
  ylab = "Number of segments",
  main = "Segments containing 0–2 genes"
)


# -------------------------
# Plot 2: >10 genes
# -------------------------
counts_over_10 <- gene_counts[
  as.numeric(names(gene_counts)) > 2
]

barplot(
  counts_over_10,
  xlab = "Number of genes in segment",
  ylab = "Number of segments",
  main = "Segments containing >2 genes"
)

##

gene_rich_segments <- short[
  gene_contained_n > 3
]

gene_rich_segments[, segment_length := end - start + 1]

gene_rich_segments[, genes_per_Mb :=
                     gene_contained_n / (segment_length / 1e6)
]

ggplot(
  plot_data,
  aes(
    x = segment_length_mb,
    fill = gene_group
  )
) +
  geom_histogram(
    bins = 40,
    position = "stack"
  ) +
  scale_x_log10() +
  labs(
    x = "Segment length (Mb, log scale)",
    y = "Number of segments",
    fill = "Genes contained"
  ) +
  theme_classic()