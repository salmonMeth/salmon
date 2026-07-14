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

path_segment_matrix="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/cpg_meth_matrix.rds"
saveRDS(segment_matrix,path_segment_matrix)

consensus_meth <- rowSums(M, na.rm = TRUE) / rowSums(Cov, na.rm = TRUE)


keep <- !is.na(consensus_meth)

consensus_meth <- consensus_meth[keep]

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

##############################################
##############################################
#PCA to visualize
segment_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_methylation"
files <- list.files(
  segment_dir,
  pattern = "_segment_methylation\\.csv$",
  full.names = TRUE
)
all_segments <- lapply(files, fread)
# create unique segment IDs
for(i in seq_along(all_segments)){
  dt <- all_segments[[i]]
  dt$segment_id <- paste(
    dt$chr,
    dt$start,
    dt$end,
    sep = "_"
  )
  all_segments[[i]] <- dt
}
all_segments <- rbindlist(all_segments, fill = TRUE)
sample_cols <- setdiff(names(all_segments),c("chr","start","end","segment_id"))
segment_matrix <- as.matrix(all_segments[, ..sample_cols])
rownames(segment_matrix) <- all_segments$segment_id
na_frac <- rowMeans(is.na(segment_matrix))

for(i in seq_len(nrow(segment_matrix))){
  idx <- is.na(segment_matrix[i, ])
  if(any(idx)){
    segment_matrix[i, idx] <- mean(segment_matrix[i, ], na.rm = TRUE)
  }
}
########################
########################
#save the methylation matrix
path_segment_matrix="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/segment_matrix.rds"
saveRDS(segment_matrix,path_segment_matrix)
#########################
#pca
#########################

pca <- prcomp(
  t(segment_matrix),
  center = TRUE,
  scale. = TRUE
)
sample_names <- colnames(segment_matrix)
group <- ifelse(
  grepl("_T1$", sample_names),
  "T1",
  "T2"
)
pca_df <- data.frame(
  Sample = sample_names,
  Group = group,
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)


p <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCA of segment methylation",
    x = paste0(
      "PC1 (",
      round(100 * summary(pca)$importance[2,1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(100 * summary(pca)$importance[2,2], 1),
      "%)"
    )
  )

pca_segment="/scratch/project_2010912/ezel/salmon/data_files/segment_methylation/pca_segment_meth.png"
ggsave(pca_segment,
       plot = p,
       width = 8,
       height = 6,
       dpi = 300
)
