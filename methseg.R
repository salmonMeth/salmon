#using methSeg to segment the CpGs
#we use the consensus levels(so the methylation levels averaged over all samples)
#to determine the methylation ratios and then create the regions

#libs
library(GenomicRanges)


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

M <- as.matrix(df[, grep("^numCs_", colnames(df))])
Cov <- as.matrix(df[, grep("^coverage_", colnames(df))])
colnames(M) <- sample_ids
colnames(Cov) <- sample_ids

meth_mat <- M / Cov

sum(is.na(meth_mat))/(dim(meth_mat)[1]*dim(meth_mat)[2])
#6% NA

#consensus <- rowMeans(meth_mat, na.rm = TRUE)
consensus_median <- apply(meth_mat, 1, median, na.rm = TRUE)

keep <- !is.na(consensus_median)

consensus_median <- consensus_median[keep]

df_cons <- df[keep, ]
gr <- GRanges(
  seqnames = df_cons$chr,
  ranges = IRanges(df_cons$start, df_cons$end)
)

mcols(gr)$meth <- consensus_median

#try it with 1 chromosome
chr <- unique(as.character(seqnames(gr)))[1]
gr_chr <- dropSeqlevels(
  gr_chr,
  setdiff(seqlevels(gr_chr), "NC_059442.1"),
  pruning.mode = "coarse"
)
##
seg_chr <- methSeg(
  gr_chr,
  diagnostic.plot = TRUE,
  minSeg = 10,
  maxInt = 100
)

#we need to fix the name formatting to be able to use this in IGV
methSeg2bed(seg_chr, "chr1_segments.bed")
df <- read.table(
  "chr1_segments.bed",
  sep = "\t",
  header = FALSE,
  skip = 1,
  colClasses = c(
    "character",  # chr (IMPORTANT)
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
head(df$V1)
unique(df$V1)
map <- setNames(seq_along(unique(df$V1)), unique(df$V1))
df$V1 <- map[df$V1]
write.table(
  df,
  file = "chr1_segments_fixed.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
###now we want to do this over all chromosomes as a loop
#####LOOP
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
  out_dir <- "segmentation_output/diagnostics"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  #we want to capture the diagnostic plots as well
  
  png(paste0("segmentation_output/diagnostics/", chr, "_diagnostic.png"),
      width = 1200, height = 800)
  
  seg_chr <- methSeg(
    gr_chr,
    diagnostic.plot = TRUE,
    minSeg = 10,
    maxInt = 100
  )
  
  dev.off()
  
  #bed file creation parts
  out_dir <- "segmentation_output/beds"
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