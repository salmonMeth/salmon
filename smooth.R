library(stats)

#matrices of the tiled genomes
matrix_file="methylation_tiles_all.rds"


#we need a function that will take in our tile matrices and a
#filtering function and return the smoothed matrix
#NOTE, this also changes the number of bins, so it reapplies tiling
#SO WE WILL USE ANOTHER FUNCTION
smooth_tiles <- function(
    rds_file,
    tile_level = "5K",
    filter_fun = NULL,
    span = 0.05,
    return_coords = TRUE,
    verbose = TRUE
) {
  
  obj <- readRDS(rds_file)
  
  mat_name <- paste0("mat_", tile_level)
  coord_name <- paste0("coord_", tile_level)
  
  mat <- obj[[mat_name]]
  coords <- obj[[coord_name]]
  sample_ids <- obj$sample_ids
  
  if (verbose) {
    cat("Loaded:", mat_name, "\n")
    cat("Matrix dim:", dim(mat), "\n")
  }
  
  # filtering
  if (!is.null(filter_fun)) {
    
    keep <- filter_fun(mat, coords)
    
    mat <- mat[keep, , drop = FALSE]
    coords <- coords[keep, , drop = FALSE]
    
    if (verbose) {
      cat("After filtering:", sum(keep), "tiles retained\n")
    }
  }
  
  # 3. Prepare output matrix
  smoothed_mat <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat))
  
  x <- coords$start  # genomic coordinate proxy
  
  # 4. LOESS per sample
  for (i in seq_len(ncol(mat))) {
    
    y <- mat[, i]
    
    ok <- !is.na(y)
    
    if (sum(ok) < 10) {
      warning(paste("Skipping sample", sample_ids[i], "- too few points"))
      next
    }
    
    fit <- loess(
      y[ok] ~ x[ok],
      span = span,
      degree = 1,
      family = "gaussian"
    )
    
    smoothed_mat[, i] <- predict(fit, newdata = x)
    
    if (verbose && i %% 10 == 0) {
      cat("Smoothed sample", i, "of", ncol(mat), "\n")
    }
  }
  
  # 5. Return result
  if (return_coords) {
    return(list(
      smoothed_matrix = smoothed_mat,
      coords = coords,
      sample_ids = sample_ids
    ))
  } else {
    return(smoothed_mat)
  }
}

####
#different loess function
smooth_tiles_2 <- function(
    rds_file,
    tile_level,
    filter_fun = NULL,
    span = 0.05,
    return_coords = TRUE,
    verbose = TRUE
) {
  
  # 1. Load object
  obj <- readRDS(rds_file)
  
  mat <- obj[[paste0("mat_", tile_level)]]
  coords <- obj[[paste0("coord_", tile_level)]]
  sample_ids <- obj$sample_ids
  
  if (verbose) {
    cat("Loaded tile level:", tile_level, "\n")
    cat("Matrix dim:", dim(mat), "\n")
  }
  
  # 2. Optional filtering (must preserve row alignment if applied)
  if (!is.null(filter_fun)) {
    
    keep <- filter_fun(mat, coords)
    
    mat <- mat[keep, , drop = FALSE]
    coords <- coords[keep, , drop = FALSE]
    
    if (verbose) {
      cat("After filtering:", sum(keep), "rows retained\n")
    }
  }
  
  # 3. Prepare output (STRICT same shape)
  smoothed <- matrix(
    NA_real_,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )
  
  x <- coords$start
  
  # 4. LOESS per sample (ROW-PRESERVING)
  for (i in seq_len(ncol(mat))) {
    
    y <- mat[, i]
    ok <- !is.na(y)
    
    if (sum(ok) < 10) {
      warning(paste("Skipping sample:", colnames(mat)[i]))
      next
    }
    
    fit <- loess(
      y[ok] ~ x[ok],
      span = span,
      degree = 1,
      family = "gaussian"
    )
    
    # CRITICAL: predict on full grid (preserves row count)
    smoothed[, i] <- predict(fit, newdata = x)
    
    if (verbose && i %% 10 == 0) {
      cat("Smoothed sample", i, "of", ncol(mat), "\n")
    }
  }
  
  # 5. Return
  if (return_coords) {
    return(list(
      smooth_mat = smoothed,
      coords = coords,
      sample_ids = sample_ids
    ))
  } else {
    return(smoothed)
  }
}
###

filter_fun <- function(mat, coords) {rowSums(!is.na(mat)) >= 5}
res <- smooth_tiles_2(matrix_file, tile_level = "20K", NULL,span = 0.05)
smooth_mat <- res$smooth_mat
coords <- res$coords
#save the matrices
saveRDS(
  list(
    smooth_mat = smooth_mat,
    coords=coords
  ),
  file = "smoothed_methylation_tiles.rds",
  compress = "xz"
)
#plot to compare the non-smoothed-tiled and the smoothed-tiled meth profiles
tiled_obj <- readRDS("methylation_tiles_all.rds")
tiled_mat <- tiled_obj$mat_20K
tiled_coords <- tiled_obj$coord_20K
#
smooth_obj <- readRDS("smoothed_methylation_tiles.rds")
smooth_mat <- smooth_obj$smooth_mat
smooth_coords <- smooth_obj$coords
#plot ten random samples
sample_subset <- sample(colnames(tiled_mat), 10)

#df for plotting
df_tiled <- data.frame(
  pos = tiled_coords$start,
  tiled_mat_subset
)
df_smoothed <- data.frame(
  pos = smooth_coords$start,
  smooth_mat_subset
)

###subset
tiled_mat_subset <- tiled_mat[, sample_subset, drop = FALSE]
smooth_mat_subset <- smooth_mat[, sample_subset, drop = FALSE]

coords <- tiled_coords
#try it with 1 sample and filter

i <- colnames(tiled_mat)[1]
df <- data.frame(
  pos = coords$start,
  raw = tiled_mat[, i],
  smooth = smooth_mat[, i]
)
set.seed(1)
idx <- seq(1, nrow(df), by = 10)
df_plot <- df[idx, ]
library(ggplot2)
ggplot(df_plot, aes(x = pos)) +
  geom_line(aes(y = raw), alpha = 0.3) +
  geom_line(aes(y = smooth), color = "red", linewidth = 1) +
  theme_bw() +
  labs(
    title = paste("Raw vs Smoothed methylation:", i),
    x = "Genomic position",
    y = "Methylation"
  )
