library(methylKit)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(Rsamtools)
#
path_genome = "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa"
#
fai <- scanFaIndex(path_genome)
seqlens <- seqlengths(fai)
#
chr_map <- setNames(
  paste0("NC_", sprintf("%06d", 59442:59470), ".1"),
  as.character(1:29)
)
chr_seqlens <- seqlens[1:29]
names(chr_seqlens) <- chr_map[names(chr_seqlens)]
###

library(GenomicRanges)

tiles <- tileGenome(
  seqlengths = chr_seqlens,
  tilewidth = 1000, #we might change this
  cut.last.tile.in.chrom = TRUE
)

tiles_5K <- tileGenome(
  seqlengths = chr_seqlens,
  tilewidth = 5000, #we might change this
  cut.last.tile.in.chrom = TRUE
)

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