# we want to check the overlaps between the segments found using the various consensus profiles.
#Specifically here we are going to look at the coverage weighted mean obtained from T1, T2 or both of them together.

library(data.table)

compare_segment_intervals <- function(
    bed_t1,
    bed_t2,
    bed_weighted,
    tol = 100000,
    out_dir = "segment_overlap_results",
    chr_name = NULL
){
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  read_segments <- function(file){
    dt <- fread(file, header = FALSE)
    dt[, .(start = V2, end = V3, meth = V5)]
  }
  
  segment_distance <- function(s1, e1, s2, e2){
    max(abs(s1 - s2), abs(e1 - e2))
  }
  
  A <- read_segments(bed_t1)
  B <- read_segments(bed_t2)
  C <- read_segments(bed_weighted)
  
  used_B <- rep(FALSE, nrow(B))
  used_C <- rep(FALSE, nrow(C))
  
  matches <- list()
  idx <- 1
  
  for(i in seq_len(nrow(A))){
    
    a_start <- A$start[i]
    a_end <- A$end[i]
    
    dB <- if(nrow(B) > 0){
      x <- mapply(segment_distance, a_start, a_end, B$start, B$end)
      x[used_B] <- Inf
      x
    } else Inf
    
    dC <- if(nrow(C) > 0){
      x <- mapply(segment_distance, a_start, a_end, C$start, C$end)
      x[used_C] <- Inf
      x
    } else Inf
    
    hitB <- any(is.finite(dB))
    hitC <- any(is.finite(dC))
    
    jB <- if(hitB) which.min(dB) else NA
    jC <- if(hitC) which.min(dC) else NA
    
    if(hitB && hitC){
      
      starts <- c(a_start, B$start[jB], C$start[jC])
      ends <- c(a_end, B$end[jB], C$end[jC])
      
      start_spread <- max(starts) - min(starts)
      end_spread <- max(ends) - min(ends)
      
      if(start_spread <= tol && end_spread <= tol){
        
        used_B[jB] <- TRUE
        used_C[jC] <- TRUE
        
        matches[[idx]] <- data.table(
          category = "ABC",
          T1_start = a_start,
          T2_start = B$start[jB],
          W_start = C$start[jC],
          T1_end = a_end,
          T2_end = B$end[jB],
          W_end = C$end[jC],
          start_spread = start_spread,
          end_spread = end_spread,
          max_spread = max(start_spread, end_spread)
        )
        
        idx <- idx + 1
        next
      }
    }
    
    if(hitB){
      
      starts <- c(a_start, B$start[jB])
      ends <- c(a_end, B$end[jB])
      
      start_spread <- max(starts) - min(starts)
      end_spread <- max(ends) - min(ends)
      
      if(start_spread <= tol && end_spread <= tol){
        
        used_B[jB] <- TRUE
        
        matches[[idx]] <- data.table(
          category = "AB",
          T1_start = a_start,
          T2_start = B$start[jB],
          W_start = NA,
          T1_end = a_end,
          T2_end = B$end[jB],
          W_end = NA,
          start_spread = start_spread,
          end_spread = end_spread,
          max_spread = max(start_spread, end_spread)
        )
        
        idx <- idx + 1
        next
      }
    }
    
    if(hitC){
      
      starts <- c(a_start, C$start[jC])
      ends <- c(a_end, C$end[jC])
      
      start_spread <- max(starts) - min(starts)
      end_spread <- max(ends) - min(ends)
      
      if(start_spread <= tol && end_spread <= tol){
        
        used_C[jC] <- TRUE
        
        matches[[idx]] <- data.table(
          category = "AC",
          T1_start = a_start,
          T2_start = NA,
          W_start = C$start[jC],
          T1_end = a_end,
          T2_end = NA,
          W_end = C$end[jC],
          start_spread = start_spread,
          end_spread = end_spread,
          max_spread = max(start_spread, end_spread)
        )
        
        idx <- idx + 1
        next
      }
    }
  }
  
  for(i in which(!used_B)){
    
    dC <- mapply(segment_distance, B$start[i], B$end[i], C$start, C$end)
    dC[used_C] <- Inf
    
    if(any(is.finite(dC))){
      
      jC <- which.min(dC)
      
      starts <- c(B$start[i], C$start[jC])
      ends <- c(B$end[i], C$end[jC])
      
      start_spread <- max(starts) - min(starts)
      end_spread <- max(ends) - min(ends)
      
      if(start_spread <= tol && end_spread <= tol){
        
        used_C[jC] <- TRUE
        
        matches[[idx]] <- data.table(
          category = "BC",
          T1_start = NA,
          T2_start = B$start[i],
          W_start = C$start[jC],
          T1_end = NA,
          T2_end = B$end[i],
          W_end = C$end[jC],
          start_spread = start_spread,
          end_spread = end_spread,
          max_spread = max(start_spread, end_spread)
        )
        
        idx <- idx + 1
      }
    }
  }
  
  shared_matches <- rbindlist(matches, fill = TRUE)
  #save 
  
  n_T1 <- nrow(A)
  n_T2 <- nrow(B)
  n_W  <- nrow(C)
  
  abc_n <- sum(shared_matches$category == "ABC", na.rm = TRUE)
  
  global_stats <- data.table(
    chromosome = chr_name,
    T1_segments = n_T1,
    T2_segments = n_T2,
    Weighted_segments = n_W,
    ABC_matches = abc_n
  )

  
  fwrite(
    shared_matches,
    file.path(out_dir, paste0(chr_name, "_matches_tol_", tol, ".csv"))
  )
  
  list(
    matches = shared_matches,
    global_stats = global_stats
  )
}

bed_dir_T1 <- "segmentation_output_T1_consensus/beds"
bed_dir_T2 <- "segmentation_output_T2_consensus/beds"
bed_dir_W  <- "segmentation_output_cov_weighted_consensus/beds"

files_T1 <- list.files(bed_dir_T1, pattern = "_segments\\.bed$", full.names = TRUE)

all_global_stats <- list()

for (f1 in files_T1) {
  
  chr <- sub("_segments\\.bed$", "", basename(f1))
  
  f2 <- file.path(bed_dir_T2, paste0(chr, "_segments.bed"))
  f3 <- file.path(bed_dir_W,  paste0(chr, "_segments.bed"))
  
  if (!file.exists(f2) || !file.exists(f3)) {
    cat("Skipping", chr, "- missing files\n")
    next
  }
  
  res <- compare_segment_intervals(
    bed_t1 = f1,
    bed_t2 = f2,
    bed_weighted = f3,
    tol = 100000,
    out_dir = "segment_overlap_results",
    chr_name = chr
  )
  
  all_global_stats[[chr]] <- res$global_stats
}


all_stats_dt <- rbindlist(all_global_stats)

fwrite(
  all_stats_dt,
  "segment_overlap_results/ALL_CHROMOSOMES_global_stats_tol_100000.csv"
)

