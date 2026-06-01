library(methylKit)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(Rsamtools)
library(GenomicRanges)

#
path_genome = "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa"
#
fai <- scanFaIndex(path_genome)
seqlens <- seqlengths(fai)
#
db <- readRDS("db.rds")
df <- getData(db)
sample_ids <- db@sample.ids
##
chr_map <- setNames(
  paste0("NC_", sprintf("%06d", 59442:59470), ".1"),
  as.character(1:29)
)
#get the chromosome lengths
chr_seqlens <- seqlens[1:29]
names(chr_seqlens) <- chr_map[names(chr_seqlens)]
###
#function for tiling
build_tile_matrix <- function(df,
                              sample_ids,
                              chr_seqlens,
                              tile_size = 5000,
                              verbose = TRUE) {

  # 1. Build tiles
  tiles <- tileGenome(
    seqlengths = chr_seqlens,
    tilewidth = tile_size,
    cut.last.tile.in.chrom = TRUE
  )
  
  # 2. CpG GRanges
  cpg_gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(df$start, df$end)
  )
  # 3. Map CpGs to tiles (done once)
  hits <- findOverlaps(cpg_gr, tiles)
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  # 4. Output matrix
  n_tiles <- length(tiles)
  n_samples <- length(sample_ids)
  
  tile_mat <- matrix(
    NA_real_,
    nrow = n_tiles,
    ncol = n_samples,
    dimnames = list(NULL, sample_ids)
  )
  
  # 5. Loop over samples
  for (i in seq_along(sample_ids)) {
    
    numCs <- df[[paste0("numCs", i)]]
    numTs <- df[[paste0("numTs", i)]]
    
    Cs_tile <- tapply(
      numCs[qh],
      sh,
      sum,
      na.rm = TRUE
    )
    
    Ts_tile <- tapply(
      numTs[qh],
      sh,
      sum,
      na.rm = TRUE
    )
    
    meth_tile <- Cs_tile / (Cs_tile + Ts_tile)
    
    tile_mat[
      as.integer(names(meth_tile)),
      i
    ] <- meth_tile
    
    if (verbose && i %% 10 == 0) {
      cat("Finished sample", i, "of", n_samples, "\n")
    }
  }
  
  # 6. Tile metadata
  tile_info <- data.frame(
    chr = as.character(seqnames(tiles)),
    start = start(tiles),
    end = end(tiles)
  )
  return(list(
    matrix = tile_mat,
    tiles = tile_info,
    hits = hits,
    tile_granges = tiles
  ))
}
#arguments of the function:
#df, sample_ids,chr_seqlens, tile_size = 5000)
#get the matrices
tile_5K= (build_tile_matrix(df,sample_ids,chr_seqlens,5e3))$matrix
tile_5K_matrix=tile_5K$matrix
tile_5K_coord=tile_5K$tiles
#10K
tile_10K= build_tile_matrix(df,sample_ids,chr_seqlens,1e4)
tile_10K_matrix=tile_10K$matrix
tile_10K_coord=tile_10K$tiles
#20K
tile_20K= build_tile_matrix(df,sample_ids,chr_seqlens,2e4)
tile_20K_matrix=tile_20K$matrix
tile_20K_coord=tile_20K$tiles

#check the sparsity of the matrices
sum(is.na(tile_5K_matrix))/(dim(tile_5K_matrix)[1]*dim(tile_5K_matrix)[2])
sum(is.na(tile_10K_matrix))/(dim(tile_10K_matrix)[1]*dim(tile_10K_matrix)[2])
sum(is.na(tile_20K_matrix))/(dim(tile_20K_matrix)[1]*dim(tile_20K_matrix)[2])
#
#save the results
saveRDS(
  list(
    mat_5K = tile_5K_matrix,
    coord_5K=tile_5K_coord,
    mat_10K = tile_10K_matrix,
    coord_10K=tile_10K_coord,
    mat_20K=tile_20K_matrix,
    coord_20K=tile_20K_coord,
    sample_ids = sample_ids
  ),
  file = "methylation_tiles_all.rds",
  compress = "xz"
)
T1_mean <- rowMeans(tile_mat_5K[, T1_idx], na.rm=TRUE)
T2_mean <- rowMeans(tile_mat_5K[, T2_idx], na.rm=TRUE)
#######################

#rest might be redundant
#TBD
db <- readRDS("db.rds")
df <- getData(db)
sample_ids <- db@sample.ids

cs_cols <- grep("^numCs", colnames(df))
ts_cols <- grep("^numTs", colnames(df))
cov_cols <- grep("^coverage", colnames(df))

gr_list <- lapply(seq_along(sample_ids), function(i) {
  numCs <- df[[cs_cols[i]]]
  numTs <- df[[ts_cols[i]]]
  cov   <- df[[cov_cols[i]]]
  denom <- numCs + numTs
  ok <- !is.na(denom) & denom > 0
  meth <- rep(NA_real_, length(denom))
  meth[ok] <- numCs[ok] / denom[ok]
  GRanges(
    seqnames = df$chr,
    ranges   = IRanges(df$start, df$end),
    strand   = df$strand,
    methylation = meth,
    coverage = cov
  )
})
#we can make it 500 as well
tile_size <- 1000