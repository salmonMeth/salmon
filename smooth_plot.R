matrix_file="data_files/filtered_methylation_tiles.rds"

obj=readRDS(matrix_file)

f_10_m=obj$filtered_10_m
f_10_c=obj$filtered_10_c
f_20_m=obj$filtered_20_m
f_20_c=obj$filtered_20_c


###
smooth_vector <- function(x, window = 5) {
  n <- length(x)
  out <- rep(NA_real_, n)
  half <- floor(window / 2)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    vals <- x[lo:hi]
    if (all(is.na(vals))) {
      out[i] <- NA_real_
    } else {
      out[i] <- mean(vals, na.rm = TRUE)
    }
  }
  out
}
##
smooth_tile_matrix <- function(mat,
                               coord,
                               window = 5) {
  mat_smooth <- mat
  chromosomes <- unique(coord$chr)
  for (chr in chromosomes) {
    idx <- which(coord$chr == chr)
    # ensure genomic order
    idx <- idx[order(coord$start[idx])]
    for (j in seq_len(ncol(mat))) {
      mat_smooth[idx, j] <-
        smooth_vector(mat[idx, j], window = window)
    }
  }
  mat_smooth
}
####

f_10_m_smooth <- smooth_tile_matrix(
  mat = f_10_m,
  coord = f_10_c,
  window = 10
)

f_20_m_smooth  <- smooth_tile_matrix(
  mat = f_20_m,
  coord = f_20_c,
  window = 10
)

#loess smoothing

smooth_vector_loess <- function(x, pos, span = 0.1) {
  
  ok <- !is.na(x) & !is.na(pos)
    if (sum(ok) < 5) return(x)
  
  fit <- loess(x[ok] ~ pos[ok], span = span, degree = 1)
  
  pred <- rep(NA_real_, length(x))
  pred[ok] <- predict(fit, newdata = pos[ok])
  
  pred
}

smooth_tile_matrix_loess <- function(mat, coord, span = 0.1) {
  
  mat_smooth <- mat
  chromosomes <- unique(coord$chr)
  
  for (chr in chromosomes) {
    idx <- which(coord$chr == chr)
    idx <- idx[order(coord$start[idx])]
    
    pos <- coord$start[idx]
    
    for (j in seq_len(ncol(mat))) {
      
      mat_smooth[idx, j] <-
        smooth_vector_loess(mat[idx, j], pos, span = span)
    }
  }
  
  mat_smooth
}

f_10_m_loess<- smooth_tile_matrix_loess(
  mat = f_10_m,
  coord = f_10_c,
  span = 0.1
)

f_20_m_loess <- smooth_tile_matrix_loess(
  mat = f_20_m,
  coord = f_20_c,
  span = 0.1
)

#plot to see the difference in smoothed vs non smooted values for the 
#i th sample and a given range of values
plot_smoot_vs_reg_1=function(smooth,reg,i, min, max,coords) {
  reg_vals=reg[min:max,i]
  smooth_vals=smooth[min:max,i]
  y_min=min(reg_vals)
  y_max=max(reg_vals)
  plot(loess_1_10[1:100,i],ylim = c(y_min,y_max))
  points(smooth_vals,col="lightblue")
}

plot_2_smooths <- function(smooth1, smooth2, i, min_idx, max_idx, coords) {
  s1 <- smooth1[min_idx:max_idx, i]
  s2 <- smooth2[min_idx:max_idx, i]
  chr <- unique(coords$chr[min_idx:max_idx])
  bp  <- coords$start[min_idx:max_idx]
  ylim <- range(c(s1, s2), na.rm = TRUE)
  plot(s1, type = "l",
       col = "#5E60CE", lwd = 2,
       ylim = ylim,
       main = paste0(chr, ": ", min(bp), "-", max(bp), " bp"),
       xlab = "Bin index", ylab = "Smoothed signal")
  lines(s2, col = "#4CC9F0", lwd = 2)
  legend("bottomright",
         legend = c("Smoothing 1", "Smoothing 2"),
         col = c("#5E60CE", "#4CC9F0"),
         lty = 1,
         lwd = 2,
         bty = "n")
}
plot_2_smooths (f_10_m_loess,f_10_m_smooth,2,500,700,f_10_c)