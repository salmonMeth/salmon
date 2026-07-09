#########################to get the matches
# Compare overlaps between segments found using different consensus profiles.
# the segments found using all the samples is the reference
#they are all cov weighted mean

library(data.table)
library(ggplot2)

#paths
dir_T1 <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_T1_consensus/beds"
dir_T2 <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_T2_consensus/beds"
dir_weighted <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_cov_weighted_consensus/beds"
#
out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation_3way_overlap"

#rel_tol is the tolerance we use to declare the segments a match 
# we say it is a match as long as the boundaries of the T1 and T2 segments 
#have a score of <.2, where score considers both the length of the segments and
#the boundaries(we also consider the length of the reference segment when we define how far the other boundaries are allowed
#to be from the reference boundary)

compare_segment_intervals <- function(
    bed_t1,
    bed_t2,
    bed_weighted,
    rel_tol = 0.2,
    out_dir = "segment_overlap_results",
    chr_name = NULL
){
  
  
  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  
  read_segments <- function(file){
    
    dt <- fread(
      file,
      header = FALSE
    )
    
    dt[, .(
      start = as.numeric(V2),
      end = as.numeric(V3)
    )]
  }
  
  # Score 
  segment_score <- function(
    ref_start,
    ref_end,
    cand_start,
    cand_end
  ){
    
    ref_length <- ref_end - ref_start + 1
    cand_length <- cand_end - cand_start + 1
    
    
    if(ref_length <= 0){
      ref_length <- 1
    }
    
    
    max(
      abs(ref_start - cand_start) / ref_length,
      abs(ref_end - cand_end) / ref_length,
      abs(ref_length - cand_length) / ref_length
    )
    
  }
  
  
  
  
  
  # Score agreement among the three segments
  triplet_score <- function(
    weighted_start,
    weighted_end,
    t1_start,
    t1_end,
    t2_start,
    t2_end
  ){
    
    
    weighted_length <- weighted_end - weighted_start + 1
    t1_length <- t1_end - t1_start + 1
    t2_length <- t2_end - t2_start + 1
    
    
    if(weighted_length <= 0){
      weighted_length <- 1
    }
    
    start_spread <- diff(
      range(
        c(
          weighted_start,
          t1_start,
          t2_start
        )
      )
    )
    
    
    end_spread <- diff(
      range(
        c(
          weighted_end,
          t1_end,
          t2_end
        )
      )
    )
    
    
    length_spread <- max(
      abs(weighted_length - t1_length),
      abs(weighted_length - t2_length)
    )
    
    
    score <- max(
      start_spread / weighted_length,
      end_spread / weighted_length,
      length_spread / weighted_length
    )
    
    
    
    list(
      score = score,
      start_spread = start_spread,
      end_spread = end_spread,
      length_spread = length_spread
    )
    
  }
  
  
  ############################################################
  # Get the segments
  
  
  seg_t1 <- read_segments(bed_t1)
  seg_t2 <- read_segments(bed_t2)
  seg_weighted <- read_segments(bed_weighted)
  
  
  
  used_t1 <- rep(FALSE, nrow(seg_t1))
  used_t2 <- rep(FALSE, nrow(seg_t2))
  
  
  
  matches <- list()
  match_id <- 1

  #loop over everything in the reference(cov weighted from all)
  for(i in seq_len(nrow(seg_weighted))){
    
    
    weighted_start <- seg_weighted$start[i]
    weighted_end <- seg_weighted$end[i]

    scores_t1 <- mapply(
      segment_score,
      weighted_start,
      weighted_end,
      seg_t1$start,
      seg_t1$end
    )
    
    scores_t1[used_t1] <- Inf
    
    scores_t2 <- mapply(
      segment_score,
      weighted_start,
      weighted_end,
      seg_t2$start,
      seg_t2$end
    )
    
    scores_t2[used_t2] <- Inf
    
    candidate_t1 <- which(
      scores_t1 <= rel_tol
    )
    
    
    candidate_t2 <- which(
      scores_t2 <= rel_tol
    )
    if(
      length(candidate_t1) == 0 ||
      length(candidate_t2) == 0
    ){
      next
    }

    best_score <- Inf
    best <- NULL

    
    for(j in candidate_t1){
      
      for(k in candidate_t2){
        
        
        triplet <- triplet_score(
          
          weighted_start,
          weighted_end,
          
          seg_t1$start[j],
          seg_t1$end[j],
          
          seg_t2$start[k],
          seg_t2$end[k]
          
        )
        
        
        if(triplet$score < best_score){
          
          best_score <- triplet$score
          
          
          best <- list(
            
            t1_index = j,
            t2_index = k,
            
            score = triplet$score,
            
            start_spread =
              triplet$start_spread,
            
            end_spread =
              triplet$end_spread,
            
            length_spread =
              triplet$length_spread
            
          )
          
        }
        
      }
      
    }

    # Save only ABC matches
    
    if(
      !is.null(best) &&
      best$score <= rel_tol
    ){
      
      
      used_t1[best$t1_index] <- TRUE
      used_t2[best$t2_index] <- TRUE
      
      
      
      
      matches[[match_id]] <- data.table(
        
        category = "ABC",
        
        
        Weighted_start = weighted_start,
        Weighted_end = weighted_end,
        
        
        T1_start =
          seg_t1$start[best$t1_index],
        
        T1_end =
          seg_t1$end[best$t1_index],
        
        
        T2_start =
          seg_t2$start[best$t2_index],
        
        T2_end =
          seg_t2$end[best$t2_index],
        
        
        score = best$score,
        
        
        start_spread =
          best$start_spread,
        
        end_spread =
          best$end_spread,
        
        length_spread =
          best$length_spread
        
      )
      
      
      
      match_id <- match_id + 1
      
    }
    
  }
  
  
  
  
  ##save
  
  
  if(length(matches) > 0){
    
    shared_matches <- rbindlist(matches)
    
  } else {
    
    shared_matches <- data.table()
    
  }
  
  
  
  
  global_stats <- data.table(
    
    chromosome = chr_name,
    
    T1_segments = nrow(seg_t1),
    
    T2_segments = nrow(seg_t2),
    
    Weighted_segments = nrow(seg_weighted),
    
    ABC_matches = nrow(shared_matches)
    
  )
  
  
  
  fwrite(
    
    shared_matches,
    
    file.path(
      out_dir,
      paste0(
        chr_name,
        "_matches_relTol_",
        rel_tol,
        ".csv"
      )
    )
    
  )
  
  
  
  return(
    list(
      matches = shared_matches,
      global_stats = global_stats
    )
  )
  
}

relative_tolerance <- 0.2

##

files_T1 <- list.files(
  
  dir_T1,
  
  pattern = "_segments\\.bed$",
  
  full.names = TRUE
  
)

cat(
  "Found T1 files:",
  length(files_T1),
  "\n"
)

all_global_stats <- list()


for(f1 in files_T1){
  
  
  chr <- sub(
    "_segments\\.bed$",    "",    basename(f1))
  
  
  
  file_T2 <- file.path(
    dir_T2,
    paste0(chr,"_segments.bed")
    
  )
  
  
  
  file_weighted <- file.path(
    
    dir_weighted,    paste0(chr,"_segments.bed")
  )
  
  
  
  if(
    !file.exists(file_T2) ||
    !file.exists(file_weighted)
  ){
    
    cat(
      "Skipping",
      chr,
      "- missing file\n"
    )
    
    next
    
  }
  
  
  
  cat(
    "Processing",
    chr,
    "\n"
  )
  
  
  result <- compare_segment_intervals(
    
    bed_t1 = f1,
    
    bed_t2 = file_T2,
    
    bed_weighted = file_weighted,
    
    rel_tol = relative_tolerance,
    
    out_dir = out_dir,
    
    chr_name = chr
    
  )
  
  
  
  all_global_stats[[chr]] <- result$global_stats
  
}

############################################################
# Save summary stats

if(length(all_global_stats) > 0){
  
  
  all_stats_dt <- rbindlist(
    all_global_stats
  )
  
  
  
  fwrite(all_stats_dt,file.path(out_dir,paste0(
        "ALL_CHROMOSOMES_global_stats_relTol_",relative_tolerance,".csv")))
} else {
  stop(
    "No chromosomes were processed."
  )
  
}


#####################plots for the scores

match_files <- list.files(
  out_dir,
  pattern = "_matches_relTol_0.2\\.csv$",
  full.names = TRUE)

all_matches <- rbindlist(
  lapply(match_files, function(f){
    dt <- fread(f)
    # extract chromosome from filename
    dt[, chromosome := sub(
      "_matches_relTol_0.2\\.csv$",
      "",
      basename(f)
    )]
    
    dt
    
  }),
  fill = TRUE
)


all_matches[, weighted_length :=
              Weighted_end - Weighted_start + 1]

# Separate non-zero matches
nonzero_scores <- all_matches[
  score > 0
]

# Count exact matches per chromosome
zero_counts <- all_matches[
  score == 0,
  .(
    zero_matches = .N
  ),
  by = chromosome
]


# Add chromosomes that might have no zero matches
all_chromosomes <- unique(all_matches$chromosome)

zero_counts <- merge(
  data.table(chromosome = all_chromosomes),
  zero_counts,
  by = "chromosome",
  all.x = TRUE
)

zero_counts[
  is.na(zero_matches),
  zero_matches := 0
]


# Plot

ggplot() +
  
  # non-zero scores
  geom_jitter(
    data = nonzero_scores,
    aes(
      x = chromosome,
      y = score
    ),
    width = 0.2,
    alpha = 0.5,
    size = 1.5,
    color = "steelblue"
  ) +
  
  
  # zero-score bubbles
  geom_point(
    data = zero_counts,
    aes(
      x = chromosome,
      y = 0,
      size = zero_matches
    ),
    shape = 21,
    fill = "orange",
    color = "black",
    alpha = 0.7
  ) +
  
  
  # tolerance line
  geom_hline(
    yintercept = 0.2,
    color = "red",
    linetype = "dashed"
  ) +
  
  
  scale_size_continuous(
    name = "Number of exact matches\n(score = 0)"
  ) +
  
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1
    )
  ) +
  
  labs(
    title = "3way segment similarity by chromosome",
    subtitle = "Blue = non-identical matches, orange bubble = exact boundary matches",
    x = "Chromosome",
    y = "Relative similarity score"
  )

#Save the plot
ggsave(
  file.path(out_dir, "overlap_scores.jpg"),
  width = 8,
  height = 6, dpi=300
)
