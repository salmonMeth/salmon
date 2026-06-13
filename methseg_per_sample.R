#we can try segmenting each sample separately to see if the results are drastically different

library(GenomicRanges)
library(GenomeInfoDb)
library(methylKit)
library(dplyr)
library(tidyr)

#all this from segment_bssq.R
db=readRDS("db.rds")
df <- getData(db)
sample_ids <- db@sample.ids
n_samples <- length(sample_ids)
#we change the colnames of the df to include the sample.id as well
colnames(df) <- c(
  "chr", "start", "end", "strand",
  as.vector(sapply(sample_ids, function(id) {
    c(paste0("coverage_", id),
      paste0("numCs_", id),
      paste0("numTs_", id))
  }))
)
##########################

path_group_m="segmentation_output/segmnet_group_numbers_across_samples.csv"
sample_ids2 <- sample_ids

M <- as.matrix(df[, grep("^numCs_", colnames(df))])
Cov <- as.matrix(df[, grep("^coverage_", colnames(df))])
colnames(M) <- sample_ids2
colnames(Cov) <- sample_ids2

#loop over the samples
#try with a small list first
sample_ids <- sample_ids2[1:2]

for (sample_name in sample_ids) {
    meth <- M[, sample_name] / Cov[, sample_name]
  
  # remove NAs
  keep <- !is.na(meth)
  
  meth <- meth[keep]
  df_sample <- df[keep, ]
  
  # build GRanges
  gr <- GRanges(
    seqnames = df_sample$chr,
    ranges = IRanges(
      start = df_sample$start,
      end = df_sample$end
    )
  )
  
  mcols(gr)$meth <- meth
  
  # output directories
  sample_dir <- file.path(
    "segmentation_output",
    sample_name
  )
  
  diagnostic_dir <- file.path(
    sample_dir,
    "diagnostics"
  )
  
  bed_dir <- file.path(
    sample_dir,
    "beds"
  )
  
  dir.create(
    diagnostic_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    bed_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # chromosome loop
  chroms <- unique(as.character(seqnames(gr)))
  
  for (chr in chroms) {
    
    message("  Processing chromosome: ", chr)
    
    gr_chr <- gr[seqnames(gr) == chr]
    
    gr_chr <- dropSeqlevels(
      gr_chr,
      setdiff(seqlevels(gr_chr), chr),
      pruning.mode = "coarse"
    )
    
    # skip tiny chromosomes
    if (length(gr_chr) < 5000) {
      message("    Skipping (<5000 CpGs)")
      next
    }
    
    # diagnostic plot
    png(
      file.path(
        diagnostic_dir,
        paste0(chr, "_diagnostic.png")
      ),
      width = 1200,
      height = 800
    )
    
    seg_chr <- methSeg(
      gr_chr,
      diagnostic.plot = TRUE,
      minSeg = 10,
      maxInt = 100
    )
    
    dev.off()
    
    # BED export
    bed_file <- file.path(
      bed_dir,
      paste0(chr, "_segments.bed")
    )
    
    methSeg2bed(
      seg_chr,
      bed_file
    )
    
    # convert chromosome names to integers for IGV if desired
    bed_df <- read.table(
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
    
    chr_map <- setNames(
      seq_along(unique(bed_df$V1)),
      unique(bed_df$V1)
    )
    
    bed_df$V1 <- chr_map[bed_df$V1]
    
    write.table(
      bed_df,
      file = bed_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
}


# we want to get some stats from the segments we created for every sample
base_dir <- "segmentation_output"
sample_dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

group_info <- list()

##############################################
#extract info
for (sdir in sample_dirs) {
  
  sample_name <- basename(sdir)
  
  bed_files <- list.files(
    file.path(sdir, "beds"),
    pattern = "_segments\\.bed$",
    full.names = TRUE
  )
  
  for (f in bed_files) {
    
    chr <- sub("_segments\\.bed$", "", basename(f))
    
    bed <- read.table(
      f,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    
    colnames(bed) <- c(
      "chr",
      "start",
      "end",
      "segment_group",
      "segment_mean",
      "strand",
      "thickStart",
      "thickEnd",
      "rgb"
    )
    
    group_info[[length(group_info) + 1]] <- data.frame(
      sample = sample_name,
      chromosome = chr,
      n_groups = length(unique(bed$segment_group))
    )
  }
}

#this has sample id, chromosome id and number of groups
group_info <- bind_rows(group_info)
#

group_matrix <- group_info %>%
  select(chromosome, sample, n_groups) %>%
  pivot_wider(
    names_from = sample,
    values_from = n_groups
  )
group_matrix
readr::write_csv(
  group_matrix,
  path_group_m
)


###########extract some stats we can just use the saved file 
group_matrix <- read.csv(
  path_group_m,
  check.names = FALSE
)

sample_cols <- colnames(group_matrix)[-1]


#we get the consensus number for segmetn_group numbers for each chromosome
get_mode <- function(x) {
  x <- na.omit(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

group_matrix$mode_groups <- apply(
  group_matrix[, sample_cols],
  1,
  get_mode
)
#we want to know if any sample consistently deviates

deviation_summary <- data.frame(
  sample = sample_cols,
  n_deviations = 0,
  net_difference = 0
)

for(i in seq_along(sample_cols)) {
  
  s <- sample_cols[i]
  
  vals <- group_matrix[[s]]
  
  deviation_summary$n_deviations[i] <-
    sum(vals != group_matrix$mode_groups,
        na.rm = TRUE)
  deviation_summary$net_difference[i] <-
    sum(
      group_matrix[[s]] -
        group_matrix$mode_groups,
      na.rm = TRUE
    )

}

deviation_summary <-
  deviation_summary[
    order(
      deviation_summary$n_deviations,
      decreasing = TRUE
    ),
  ]
#
path_deviation= "segmentation_output/deviation_from_consensus_per_sample.csv"

readr::write_csv(
  deviation_summary,
  path_deviation
)