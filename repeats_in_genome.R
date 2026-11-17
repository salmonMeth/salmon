#work with the repeats file

library(data.table)

repeat_file <- "/scratch/project_2010912/ezel/salmon/data_files/GCF_905237065.1_Ssal_v3.1.sorted.fa.out"

repeats <- fread(
  repeat_file,
  skip = 3,
  header = FALSE,
  fill = TRUE
)

setnames(
  repeats,
  c(
    "score", "divergence", "deletion", "insertion",
    "chr", "start", "end", "query_left",
    "strand", "repeat_name", "repeat_class",
    "repeat_start", "repeat_end", "repeat_left", "ID"
  )
)

head(repeats)


repeat_file_2 <- "/scratch/project_2010912/ezel/salmon/data_files/V2_GCF_905237065.1_Ssal_v3.1.sorted.fa.out"

repeats2 <- data.table::fread(
  repeat_file_2,
  skip = 3,
  header = FALSE,
  fill = TRUE
)


# RepeatModeler annotation columns
setnames(
  repeats2,
  c(
    "score", "divergence", "deletion", "insertion",
    "chr", "start", "end", "query_left",
    "strand", "repeat_name", "repeat_class",
    "repeat_start", "repeat_end", "repeat_left", "ID"
  )
)

# Classes to exclude
exclude_classes <- c(
  "Unknown",
  "tRNA",
  "Satellite",
  "Simple_repeat",
  "rRNA",
  "Low_complexity"
)

# Retain all other repeat annotations
te_repeats <- repeats2[
  !repeat_class %in% exclude_classes
]

# Check dimensions and remaining classes
dim(te_repeats)
sort(table(te_repeats$repeat_class), decreasing = TRUE)

# Save the filtered annotation
saveRDS(
  te_repeats,
  file = "salmon_Ssal_v3.1_TE_annotations.rds"
)


#compare with beds

bed_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_3_samples_removed/beds"

te_file <- "/scratch/project_2010912/ezel/salmon/data_files/salmon_Ssal_v3.1_TE_annotations.rds"

output_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/te_overlap"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load TE annotations
te <- as.data.table(readRDS(te_file))

# Read all segment BED files
bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

stopifnot(length(bed_files) > 0)

segments <- rbindlist(
  lapply(bed_files, function(f) {
    x <- fread(f, header = FALSE)
    
    # Standard BED: chromosome, start, end, optional extra columns
    setnames(x, 1:3, c("chr", "bed_start", "bed_end"))
    
    x[, segment_file := basename(f)]
    x
  }),
  fill = TRUE
)

# Create a unique segment ID
segments[, segment_id := .I]

# bed&repeat masker compatibility
segments[, `:=`(
  start = as.integer(bed_start) + 1L,
  end   = as.integer(bed_end)
)]

#define tolerance
segments[, segment_length := end - start + 1L]
segments[, tolerance := 0.10 * segment_length]

te <- te[
  !is.na(chr) & !is.na(start) & !is.na(end) &
    start <= end
]

segments <- segments[
  !is.na(chr) & !is.na(start) & !is.na(end) &
    start <= end
]


# Map RefSeq accessions to chromosome numbers
chr_map <- data.table(
  accession = paste0("NC_",  sprintf("%06d", 59442:59470), ".1"),
  chr_num = as.character(1:29)
)

# Standardize chromosome names in both datasets
segments[, chr := as.character(chr)]
te[, chr := as.character(chr)]

# Convert TE accessions to chromosome numbers
te <- merge(
  te,
  chr_map,
  by.x = "chr",
  by.y = "accession",
  all = FALSE
)

te[, chr := chr_num]
te[, chr_num := NULL]

# Check the mapping
head(unique(te$chr))
head(unique(segments$chr))

# Recalculate shared chromosomes
common_chr <- intersect(
  unique(segments$chr),
  unique(te$chr)
)

print(common_chr)
length(common_chr)


segments_chr <- segments[chr %in% common_chr]
te_chr <- te[chr %in% common_chr]

library(GenomicRanges)
library(IRanges)

seg_gr <- GRanges(
  seqnames = segments_chr$chr,
  ranges = IRanges(
    start = pmax(1L, segments_chr$start -
                   ceiling(segments_chr$tolerance)),
    end = segments_chr$end +
      ceiling(segments_chr$tolerance)
  ),
  segment_id = segments_chr$segment_id
)

te_gr <- GRanges(
  seqnames = te_chr$chr,
  ranges = IRanges(start = te_chr$start, end = te_chr$end),
  te_row = seq_len(nrow(te_chr))
)

hits <- findOverlaps(seg_gr, te_gr, ignore.strand = TRUE)

# Recover matched rows
candidate_segments <- segments_chr[
  match(mcols(seg_gr)$segment_id[queryHits(hits)], segment_id)
]

candidate_te <- te_chr[
  mcols(te_gr)$te_row[subjectHits(hits)]
]

# Apply the stricter border-to-border criterion
matched <- data.table(
  segment_id = candidate_segments$segment_id,
  segment_chr = candidate_segments$chr,
  segment_start = candidate_segments$start,
  segment_end = candidate_segments$end,
  segment_length = candidate_segments$segment_length,
  segment_file = candidate_segments$segment_file,
  te_start = candidate_te$start,
  te_end = candidate_te$end,
  te_name = candidate_te$repeat_name,
  te_class = candidate_te$repeat_class
)

matched[, tolerance := 0.10 * segment_length]

matched <- matched[
  abs(te_start - segment_start) <= tolerance &
    abs(te_end - segment_end) <= tolerance
]


#compare diff meth loci with TEs

library(GenomicRanges)
library(IRanges)
library(data.table)
mk=read.csv("/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_mk.csv")
# Convert mk and TE annotations to data.tables
mk_dt <- as.data.table(mk)
te_dt <- as.data.table(te)

# Ensure coordinates are numeric
mk_dt[, `:=`(
  start = as.integer(start),
  end   = as.integer(end)
)]

te_dt[, `:=`(
  start = as.integer(start),
  end   = as.integer(end)
)]

# Keep valid coordinates and matching chromosomes
mk_dt <- mk_dt[
  !is.na(chr) & !is.na(start) & !is.na(end) &
    start >= 1 & end >= start
]

te_dt <- te_dt[
  !is.na(chr) & !is.na(start) & !is.na(end) &
    start >= 1 & end >= start
]



te_gr <- GRanges(
  seqnames = te_dt$chr,
  ranges = IRanges(start = te_dt$start, end = te_dt$end)
)

chr_map <- data.table(
  accession = paste0("NC_", sprintf("%06d", 59442:59470), ".1"),
  chr = as.character(1:29)
)

mk_dt[, chr_num := chr_map$chr[
  match(as.character(mk_dt$chr), chr_map$accession)
]]


mk_gr <- GRanges(
  seqnames = mk_dt$chr_num,
  ranges = IRanges(start = mk_dt$start, end = mk_dt$end)
)

hits <- findOverlaps(mk_gr, te_gr, ignore.strand = TRUE)

n_overlap <- length(unique(queryHits(hits)))

library(data.table)
library(GenomicRanges)
library(IRanges)

# Find overlaps (using the chromosome-matched GRanges objects)
hits <- findOverlaps(mk_gr, te_gr, ignore.strand = TRUE)

# Indices of overlapping methylation loci and TEs
mk_idx <- queryHits(hits)
te_idx <- subjectHits(hits)

# Start with the original methylation information
overlapping_loci <- copy(mk_dt[mk_idx])

# Add TE information
overlapping_loci[, `:=`(
  te_chr         = as.character(seqnames(te_gr))[te_idx],
  te_start       = start(te_gr)[te_idx],
  te_end         = end(te_gr)[te_idx],
  te_strand      = as.character(strand(te_gr))[te_idx],
  te_name        = te_dt$repeat_name[te_idx],
  te_class       = te_dt$repeat_class[te_idx],
  te_id          = te_dt$ID[te_idx]
)]

