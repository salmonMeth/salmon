#here we will see what happens if we create two different consensus profiles for 
#T1 and T2 separately to understand if it actually makes a difference we should be
#worried about

#libs
library(GenomicRanges)
library(dplyr)
library(GenomeInfoDb)
library(methylKit)
#
db <- readRDS("db.rds")
df <- getData(db)

sample_ids <- db@sample.ids

colnames(df) <- c(
  "chr", "start", "end", "strand",
  as.vector(sapply(sample_ids, function(id) {
    c(
      paste0("coverage_", id),
      paste0("numCs_", id),
      paste0("numTs_", id)
    )
  }))
)

M <- as.matrix(df[, grep("^numCs_", colnames(df))])
Cov <- as.matrix(df[, grep("^coverage_", colnames(df))])

colnames(M) <- sample_ids
colnames(Cov) <- sample_ids

T1_samples <- sample_ids[grepl("_T1$", sample_ids)]
T2_samples <- sample_ids[grepl("_T2$", sample_ids)]

#
consensus_T1 <- rowSums(M[, T1_samples, drop = FALSE], na.rm = TRUE) /
  rowSums(Cov[, T1_samples, drop = FALSE], na.rm = TRUE)

consensus_T2 <- rowSums(M[, T2_samples, drop = FALSE], na.rm = TRUE) /
  rowSums(Cov[, T2_samples, drop = FALSE], na.rm = TRUE)

consensus_T1[rowSums(Cov[, T1_samples, drop = FALSE], na.rm = TRUE) == 0] <- NA
consensus_T2[rowSums(Cov[, T2_samples, drop = FALSE], na.rm = TRUE) == 0] <- NA

keep_T1 <- !is.na(consensus_T1)

gr_T1 <- GRanges(
  seqnames = df$chr[keep_T1],
  ranges = IRanges(
    start = df$start[keep_T1],
    end = df$end[keep_T1]
  )
)
mcols(gr_T1)$meth <- consensus_T1[keep_T1]

keep_T2 <- !is.na(consensus_T2)

gr_T2 <- GRanges(
  seqnames = df$chr[keep_T2],
  ranges = IRanges(
    start = df$start[keep_T2],
    end = df$end[keep_T2]
  )
)
mcols(gr_T2)$meth <- consensus_T2[keep_T2]

############################################
#####LOOP
gr= gr_T1
chroms <- seqlevels(gr)

for (i in seq_along(chroms)) {
  chr <- chroms[i]
  message("Processing: ", chr)
  gr_chr <- gr[seqnames(gr) == chr]
  gr_chr <- dropSeqlevels(
    gr_chr,
    setdiff(seqlevels(gr_chr), chr),
    pruning.mode = "coarse"
  )
  out_dir <- "segmentation_output_T1_consensus/diagnostics"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  #we want to capture the diagnostic plots as well
  
  png(paste0("segmentation_output_T1_consensus/diagnostics/", chr, "_diagnostic.png"),
      width = 1200, height = 800)
  
  seg_chr <- methSeg(
    gr_chr,
    diagnostic.plot = TRUE,
    minSeg = 10,
    maxInt = 100
  )
  
  dev.off()
  
  #bed file creation parts
  out_dir <- "segmentation_output_T1_consensus/beds"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  bed_file <- file.path(out_dir, paste0(chr, "_segments.bed"))
  
  methSeg2bed(seg_chr, bed_file)
  
  df <- read.table(
    bed_file,
    sep = "\t",
    header = FALSE,
    skip = 1,
    colClasses = c(
      "character",
      "integer",
      "integer",
      "integer",
      "numeric",
      "character",
      "integer",
      "integer",
      "character"
    )
  )
  
  map <- setNames(seq_along(unique(df$V1)), unique(df$V1))
  df$V1 <- map[df$V1]
  
  write.table(
    df,
    file = bed_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}
#
#these bed files dont have column names to keeo it compatible with the
#IGV browser but we can add them later if we want

add_col_names= function(bed_file) {
  df <- read.table(bed_file, sep = "\t", header = FALSE)
  colnames(df) <- c(
    "chr",
    "start",
    "end",
    "segment_group", #tells us which segment group this segment belongs to
    "segment_mean", #tells us the mean methylartion value of the defined segment
    "strand",
    "thickStart",
    "thickEnd",
    "rgb"
  )
  return(df)
}
###########################################################
######getting the stats from the bed files

bed_dir <- "segmentation_output_T1_consensus/beds"

bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

chr_summary <- list()
group_summary <- list()

for (f in bed_files) {
  
  chr_name <- sub("_segments\\.bed$", "", basename(f))
  
  df <- read.table(
    f,
    sep = "\t",
    header = FALSE
  )
  
  colnames(df) <- c(
    "chr",
    "start",
    "end",
    "seg.group",
    "seg.mean",
    "strand",
    "thickStart",
    "thickEnd",
    "rgb"
  )
  
  df$width <- df$end - df$start
  
  chr_summary[[chr_name]] <- data.frame(
    chromosome = chr_name,
    n_segment_groups = length(unique(df$seg.group)),
    n_segments = nrow(df),
    total_bp = sum(df$width),
    mean_segment_length = mean(df$width),
    median_segment_length = median(df$width),
    sd_segment_length = sd(df$width),
    mean_methylation = mean(df$seg.mean),
    median_methylation = median(df$seg.mean)
  )
  
  group_summary[[chr_name]] <- df %>%
    group_by(seg.group) %>%
    summarise(
      chromosome = chr_name,
      
      n_segments = n(),
      total_bp = sum(width),
      
      mean_length = mean(width),
      median_length = median(width),
      sd_length = sd(width),
      
      min_length = min(width),
      q25_length = quantile(width, 0.25),
      q75_length = quantile(width, 0.75),
      max_length = max(width),
      
      mean_meth = mean(seg.mean),
      median_meth = median(seg.mean),
      sd_meth = sd(seg.mean),
      
      min_meth = min(seg.mean),
      q25_meth = quantile(seg.mean, 0.25),
      q75_meth = quantile(seg.mean, 0.75),
      max_meth = max(seg.mean),
      cv_meth = sd(seg.mean) / mean(seg.mean),
      
      largest_segment = max(width),
      smallest_segment = min(width),
      
      .groups = "drop"
    ) %>%
    mutate(
      prop_genome = total_bp / sum(total_bp),
    )
}

chr_summary <- bind_rows(chr_summary)
group_summary <- bind_rows(group_summary)

chr_summary <- chr_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

group_summary <- group_summary %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))


dir.create("segmentation_output_T1_consensus/summaries",
           recursive = TRUE,
           showWarnings = FALSE)

write.csv(
  chr_summary,
  "segmentation_output_T1_consensus/summaries/chromosome_summary.csv",
  row.names = FALSE
)

write.csv(
  group_summary,
  "segmentation_output_T1_consensus/summaries/segment_group_summary.csv",
  row.names = FALSE
)