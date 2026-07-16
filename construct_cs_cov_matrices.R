#we need to have the actual numbers(natural numbers) of the methylated and unmethylated
#counts per segment to fit the models
#so we extract them (for the reduced sample set) from the db object

library(data.table)
library(GenomicRanges)
library(dplyr)
library(GenomeInfoDb)
library(methylKit)

path_db="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db.rds"
#all this from segment_bssq.R
db=readRDS(path_db)
#######
path_m="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/cpg_meth_total_matrix_3_samples_removed.rds"
M <- readRDS(path_m)
path_cov=    "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/cpg_coverage_total_matrix_3_samples_removed.rds"
Cov <- readRDS(path_cov)
#####################
####################
####################
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

cpg_gr <- GRanges(
  seqnames = df$chr,
  ranges = IRanges(
    start = df$start,
    end = df$end
  )
)

colnames(M) <- sample_ids
colnames(Cov) <- sample_ids

#########################################################
#########################################################

# recreate chromosome lookup

unique_chr <- unique(df$chr)

chrom_lookup <- setNames(
  unique_chr,
  as.character(seq_along(unique_chr))
)

#########################################################

bed_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_3_samples_removed/beds"

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bed_files <- list.files(
  bed_dir,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

for(segment_file in bed_files){
  
  file_chr <- sub(
    "_segments\\.bed$",
    "",
    basename(segment_file)
  )
  
  cat("Processing", file_chr, "\n")
  
  seg <- fread(segment_file, header = FALSE)
  
  chr_name <- sub(
    "_segments\\.bed$",
    "",
    basename(segment_file)
  )
  
  
  
  regions <- GRanges(
    seqnames = rep(chr_name, nrow(seg)),
    ranges = IRanges(
      start = seg$V2,
      end = seg$V3
    )
  )
  
  keep_chr <- as.character(seqnames(cpg_gr)) == chr_name
  
  cpg_chr <- cpg_gr[keep_chr]
  
  M_chr <- M[keep_chr, , drop = FALSE]
  
  Cov_chr <- Cov[keep_chr, , drop = FALSE]
  
  cat(
    "CpGs:", nrow(M_chr),
    "Segments:", length(regions),
    "\n"
  )
  
  hits <- findOverlaps(cpg_chr, regions)
  
  cat(
    "Overlaps:", length(hits),
    "\n"
  )
  
  segment_cs <- matrix(
    0,
    nrow = length(regions),
    ncol = ncol(M_chr)
  )
  
  segment_cov <- matrix(
    0,
    nrow = length(regions),
    ncol = ncol(Cov_chr)
  )
  
  for(i in seq_len(length(regions))){
    
    idx <- queryHits(hits)[subjectHits(hits) == i]
    
    if(length(idx) == 0)
      next
    
    segment_cs[i, ] <- colSums(
      M_chr[idx, , drop = FALSE],
      na.rm = TRUE
    )
    
    segment_cov[i, ] <- colSums(
      Cov_chr[idx, , drop = FALSE],
      na.rm = TRUE
    )
  }
  
  colnames(segment_cs) <- colnames(M)
  
  colnames(segment_cov) <- colnames(Cov)
  
  base_dt <- data.table(
    chr = as.character(seqnames(regions)),
    start = start(regions),
    end = end(regions)
  )
  
  cs_dt <- cbind(
    copy(base_dt),
    as.data.table(segment_cs)
  )
  
  cov_dt <- cbind(
    copy(base_dt),
    as.data.table(segment_cov)
  )
  
  fwrite(
    cs_dt,
    file.path(
      out_dir,
      paste0(file_chr, "_segment_cs.csv")
    )
  )
  
  fwrite(
    cov_dt,
    file.path(
      out_dir,
      paste0(file_chr, "_segment_coverage.csv")
    )
  )
  
  cat(
    "Saved",
    nrow(cs_dt),
    "segments\n\n"
  )
}

cat("Done.\n")