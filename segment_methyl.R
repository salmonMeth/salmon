library(data.table)
library(methylKit)
library(GenomicRanges)
# Directories
segment_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_cov_weighted_consensus/beds"
##
db = readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db.rds")
##############
##############
###############################################################first we get the coverage, methylation etc matrices
#from methSeg.R
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
################
################
#for 1 chr

segment_file <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_cov_weighted_consensus/beds/NC_059443.1_segments.bed"

segments <- fread(segment_file, header = FALSE)

chr_name <- sub(
  "_segments\\.bed$",
  "",
  basename(segment_file)
)

segments <- data.table(
  chr = chr_name,
  start = segments$V2,
  end = segments$V3
)


regions <- GRanges(
  seqnames = segments$chr,
  ranges = IRanges(
    start = segments$start,
    end = segments$end
  )
)

#we want to subset M and Cov etc
chr_name <- unique(as.character(seqnames(regions)))
keep_chr <- as.character(seqnames(gr)) == chr_name

gr_chr <- gr[keep_chr]

M_chr <- M[keep_chr, ]
Cov_chr <- Cov[keep_chr, ]
hits <- findOverlaps(gr_chr, regions)
#create the segment *methylation matrix, rows are segments, cols are samples
segment_meth <- matrix(NA_real_,nrow = length(regions),ncol = ncol(M_chr))

for(i in seq_len(length(regions))) {
  
  idx <- queryHits(hits)[subjectHits(hits) == i]
  
  if(length(idx) == 0)
    next
  
  numCs_seg <- colSums(M_chr[idx, , drop = FALSE], na.rm = TRUE)
  cov_seg <- colSums(Cov_chr[idx, , drop = FALSE],na.rm = TRUE)
  segment_meth[i, ] <- numCs_seg / cov_seg
}

colnames(segment_meth) <- sample_ids
###################################
###################################
#now we loop over all the chromosomes
#there is some debug printing commented out
library(GenomeInfoDb)

bed_dir <-  "/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_cov_weighted_consensus/beds"


bed_files <- list.files(bed_dir, pattern = "_segments\\.bed$", full.names = TRUE)

out_dir = "/scratch/project_2010912/ezel/salmon/data_files/segment_methylation"


dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for(segment_file in bed_files){
  
  chr_name <- sub("_segments\\.bed$", "", basename(segment_file))
  
  cat("Processing", chr_name, "\n")
  
  segments <- fread(segment_file, header = FALSE)
  
  segments <- data.table(chr = chr_name,start = segments$V2,end = segments$V3)
  
  regions <- GRanges(seqnames = segments$chr,ranges = IRanges(start = segments$start,
                                      end = segments$end))
  
  this_chr <- sub("_segments\\.bed$", "", basename(segment_file))
  keep_chr <- as.character(seqnames(gr)) == this_chr
  gr_chr <- gr[keep_chr]
  gr_chr <- keepSeqlevels(gr_chr, chr_name, pruning.mode = "coarse")   
  
  #cat("BED chromosome:", this_chr, "\n")
  #cat("GR chromosome :", unique(as.character(seqnames(gr_chr))), "\n")
  
  M_chr <- M[keep_chr, ]
  Cov_chr <- Cov[keep_chr, ]
  
  hits <- findOverlaps(gr_chr, regions)
  
  segment_meth <- matrix(NA_real_, nrow = length(regions), ncol = ncol(M_chr))
  
  for(i in seq_len(length(regions))){
    
    idx <- queryHits(hits)[subjectHits(hits) == i]
    
    if(length(idx) == 0)
      next
    
    numCs_seg <- colSums(M_chr[idx, , drop = FALSE], na.rm = TRUE)
    cov_seg <- colSums(Cov_chr[idx, , drop = FALSE], na.rm = TRUE)
    
    segment_meth[i, ] <- numCs_seg / cov_seg
  }
  
  colnames(segment_meth) <- sample_ids
  
  segment_dt <- data.table(
    chr = as.character(seqnames(regions)),
    start = start(regions),
    end = end(regions)
  )
  
  segment_dt <- cbind(segment_dt, as.data.table(segment_meth))
  
  fwrite(segment_dt,
         file.path(out_dir,
                   paste0(chr_name, "_segment_methylation.csv")))
  
  cat("Saved", nrow(segment_dt), "segments\n")
}

cat("Done.\n")
