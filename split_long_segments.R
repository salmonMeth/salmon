#split the long segments into 2
library(data.table)

bed_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segmentation_output_G1_5/beds"

bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

# Read and combine all BED files
bed_list <- lapply(bed_files, function(f) {
  
  dt <- fread(f, header = FALSE)
  
  # Keep chromosome and filename
  dt[, chr := sub("_segments\\.bed$", "", basename(f))]
  
  dt
})

segments <- rbindlist(bed_list, fill = TRUE)

keep_file="/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/segments_variable_T1_or_T2.csv"

keep=fread(keep_file)

if (names(keep)[1] == "" || names(keep)[1] == "V1") {
  keep[, 1 := NULL]
}

setnames(keep, tolower(names(keep)))
# Rename BED coordinate columns
segments[, `:=`(
  start = V2,
  end   = V3
)]

# Clean chromosome strings
segments[, chr := trimws(as.character(chr))]
keep[, chr := trimws(as.character(chr))]
# Keep only segments whose chr, start and end occur in keep
segments_filtered <- merge(
  segments,
  unique(keep[, .(chr, start, end)]),
  by = c("chr", "start", "end"),
  all = FALSE
)

# Keep only genomic coordinates
segments_split <- copy(
  segments_filtered[, .(chr, start, end)]
)

# Ensure integer coordinates
segments_split[, `:=`(
  start = as.integer(start),
  end = as.integer(end)
)]

# Calculate segment lengths
segments_split[, length := end - start]

# Identify the 150 largest segments
top150_indices <- order(-segments_split$length)[1:min(150, .N)]

# Extract segments to split
to_split <- segments_split[top150_indices]

# Split into two approximately equal fragments
first_fragments <- to_split[, .(
  chr,
  start,
  end = as.integer(floor((start + end) / 2))
)]

second_fragments <- to_split[, .(
  chr,
  start = as.integer(floor((start + end) / 2) + 1),
  end
)]

# Retain all other segments unchanged
unchanged <- segments_split[-top150_indices, .(chr, start, end)]

# Combine everything
segments_final <- rbindlist(list(
  unchanged,
  first_fragments,
  second_fragments
))

# Sort and check
setorder(segments_final, chr, start, end)

stopifnot(
  is.integer(segments_final$start),
  is.integer(segments_final$end),
  all(segments_final$end >= segments_final$start)
)

cat("Original segments:", nrow(segments_filtered), "\n")
cat("Final segments:", nrow(segments_final), "\n")
cat("Segments split:", nrow(to_split), "\n")

db = readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db.rds")
# Ensure segments_final is a data.table
setDT(segments_final)

# Ensure integer coordinates
segments_final[, `:=`(
  start = as.integer(start),
  end = as.integer(end)
)]

# Get chromosomes
chr_names <- unique(segments_final$chr)
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

#to remove the corrupted samples
# Samples to remove
remove_samples <- c("746_T2", "176_T1", "861_T2")

# Keep the remaining sample IDs
sample_ids_keep <- setdiff(sample_ids, remove_samples)

# Columns to retain
cols_to_keep <- c(
  "chr", "start", "end", "strand",
  as.vector(sapply(sample_ids_keep, function(id) {
    c(
      paste0("coverage_", id),
      paste0("numCs_", id),
      paste0("numTs_", id)
    )
  }))
)

# Subset the data frame
df <- df[, cols_to_keep]
sample_ids <- sample_ids_keep
n_samples <- length(sample_ids)

####
setDT(df)
setDT(segments_final)

# Ensure chromosome names are compatible
df[, chr := as.character(chr)]
segments_final[, chr := as.character(chr)]

# Check that chromosome names match
setdiff(unique(df$chr), unique(segments_final$chr))
# Segments: convert BED coordinates to 1-based, closed
segments_final[, segment_id := .I]

segments_gr <- GRanges(
  seqnames = segments_final$chr,
  ranges = IRanges(
    start = segments_final$start + 1L,
    end = segments_final$end
  ),
  segment_id = segments_final$segment_id
)
# Methylation loci: already 1-based positions
loci_gr <- GRanges(
  seqnames = df$chr,
  ranges = IRanges(
    start = df$start,
    end = df$end
  )
)

hits <- findOverlaps(
  query = loci_gr,
  subject = segments_gr,
  ignore.strand = TRUE
)

length(hits)

# Get the locus and segment indices for every overlap
query_idx <- queryHits(hits)
subject_idx <- subjectHits(hits)

overlap_dt <- data.table(
  segment_id = mcols(segments_gr)$segment_id[subjectHits(hits)],
  locus_index = queryHits(hits)
)

for (col in cs_cols) {
  overlap_dt[, (col) := df[[col]][locus_index]]
}

head(overlap_dt)

segment_matrix <- merge(
  segments_final,
  segment_counts,
  by = "segment_id",
  all.x = TRUE
)

# Replace missing counts with zero
for (col in cs_cols) {
  set(
    segment_matrix,
    which(is.na(segment_matrix[[col]])),
    col,
    0
  )
}

setorder(segment_matrix, chr, start)
