#we get the T1, T2 and cumulative beds to find the sections where they disagree
#disagreement is defined as the lack of a 3way match, and we merge the overlapping or adjacent parts

library(data.table)

project_dir <- "/scratch/project_2010912/ezel/rSalmon/salmon-main/data_files"

bed_dir_T1 <- file.path(project_dir,
                        "segmentation_output_T1_consensus/beds")

bed_dir_T2 <- file.path(project_dir,
                        "segmentation_output_T2_consensus/beds")

bed_dir_W <- file.path(project_dir,
                       "segmentation_output_cov_weighted_consensus/beds")

out_dir <- file.path(project_dir,
                     "segment_overlap_results",
                     "disagreement_regions")

dir.create(out_dir,
           recursive = TRUE,
           showWarnings = FALSE)

read_segments <- function(file){
  
  dt <- fread(file, header = FALSE)
  
  dt[, .(
    start = V2,
    end = V3
  )]
}

segment_distance <- function(s1,e1,s2,e2){
  
  max(
    abs(s1-s2),
    abs(e1-e2)
  )
}

merge_intervals <- function(dt){
  
  if(nrow(dt) == 0)
    return(dt)
  
  setorder(dt, start, end)
  
  merged <- list()
  
  cur_start <- dt$start[1]
  cur_end <- dt$end[1]
  
  idx <- 1
  
  for(i in 2:nrow(dt)){
    
    if(dt$start[i] <= cur_end + 1){
      
      cur_end <- max(cur_end,
                     dt$end[i])
      
    } else {
      
      merged[[idx]] <- data.table(
        start = cur_start,
        end = cur_end
      )
      
      idx <- idx + 1
      
      cur_start <- dt$start[i]
      cur_end <- dt$end[i]
    }
  }
  
  merged[[idx]] <- data.table(
    start = cur_start,
    end = cur_end
  )
  
  rbindlist(merged)
}

find_disagreement_regions <- function(
    bed_t1,
    bed_t2,
    bed_w,
    chr_name,
    tol = 100000
){
  
  A <- read_segments(bed_t1)
  B <- read_segments(bed_t2)
  C <- read_segments(bed_w)
  
  used_B <- rep(FALSE,nrow(B))
  used_C <- rep(FALSE,nrow(C))
  
  disagreements <- list()
  k <- 1
  
  for(i in seq_len(nrow(A))){
    
    dB <- mapply(
      segment_distance,
      A$start[i],
      A$end[i],
      B$start,
      B$end
    )
    
    dC <- mapply(
      segment_distance,
      A$start[i],
      A$end[i],
      C$start,
      C$end
    )
    
    dB[used_B] <- Inf
    dC[used_C] <- Inf
    
    if(!any(is.finite(dB)))
      next
    
    if(!any(is.finite(dC)))
      next
    
    jB <- which.min(dB)
    jC <- which.min(dC)
    
    starts <- c(
      A$start[i],
      B$start[jB],
      C$start[jC]
    )
    
    ends <- c(
      A$end[i],
      B$end[jB],
      C$end[jC]
    )
    
    start_spread <- max(starts) - min(starts)
    end_spread <- max(ends) - min(ends)
    
    if(start_spread > tol ||
       end_spread > tol)
      next
    
    used_B[jB] <- TRUE
    used_C[jC] <- TRUE
    
    union_start <- min(starts)
    union_end <- max(ends)
    
    shared_start <- max(starts)
    shared_end <- min(ends)
    
    if(shared_start > union_start){
      
      disagreements[[k]] <- data.table(
        chr = chr_name,
        start = union_start,
        end = shared_start
      )
      
      k <- k + 1
    }
    
    if(union_end > shared_end){
      
      disagreements[[k]] <- data.table(
        chr = chr_name,
        start = shared_end,
        end = union_end
      )
      
      k <- k + 1
    }
  }
  
  if(length(disagreements) == 0){
    
    final_dt <- data.table(
      chr = character(),
      start = integer(),
      end = integer()
    )
    
  } else {
    
    final_dt <- rbindlist(disagreements)
    
    merged <- merge_intervals(
      final_dt[, .(start,end)]
    )
    
    final_dt <- data.table(
      chr = chr_name,
      start = merged$start,
      end = merged$end,
      width = merged$end - merged$start
    )
  }
  
  fwrite(
    final_dt,
    file.path(
      out_dir,
      paste0(chr_name,
             "_disagreement_regions.bed")
    ),
    sep = "\t",
    col.names = FALSE
  )
  
  final_dt
}
###execute
files_T1 <- list.files(
  bed_dir_T1,
  pattern = "_segments\\.bed$",
  full.names = TRUE
)

all_stats <- list()

for(f1 in files_T1){
  
  chr <- sub(
    "_segments\\.bed$",
    "",
    basename(f1)
  )
  
  f2 <- file.path(
    bed_dir_T2,
    paste0(chr,"_segments.bed")
  )
  
  f3 <- file.path(
    bed_dir_W,
    paste0(chr,"_segments.bed")
  )
  
  if(!file.exists(f2) ||
     !file.exists(f3))
    next
  
  dt <- find_disagreement_regions(
    bed_t1 = f1,
    bed_t2 = f2,
    bed_w = f3,
    chr_name = chr,
    tol = 100000
  )
  
  all_stats[[chr]] <- data.table(
    chromosome = chr,
    n_disagreement_regions = nrow(dt),
    total_bp = sum(dt$width)
  )
}

summary_dt <- rbindlist(all_stats)

fwrite(
  summary_dt,
  file.path(
    out_dir,
    "disagreement_region_summary.csv"
  )
)
##########################################
# now we want to annotate these regions
#########
######## we need to create a grrange obj from our beds of the disagreements
bed_dir <- "/scratch/project_2010912/ezel/rSalmon/salmon-main/data_files/segment_overlap_results/disagreement_regions"
path_annot_genes ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"

#
bed_files <- list.files(
  bed_dir,
  pattern = "_disagreement_regions\\.bed$",
  full.names = TRUE
)

all_regions <- rbindlist(lapply(bed_files, function(f){
  
  if (file.info(f)$size == 0) {
    message("Skipping empty file: ", basename(f))
    return(NULL)
  }
  
  dt <- fread(f, header = FALSE)
  
  setnames(dt, c("chr","start","end","width"))
  
  dt
  
}), fill = TRUE)
#####################
####we get the annotations


gene <- txdbmaker::makeTxDbFromGFF(path_annot_genes)

gene_gr <- genes(gene)
exon_gr <- exons(gene)
intron_gr <- unlist(intronsByTranscript(gene))

gene_parts <- list(
  promoters = promoters(gene, upstream = 2000, downstream = 200),
  exons = exon_gr,
  introns = intron_gr
)

proms <- gene_parts$promoters
exs   <- gene_parts$exons
ints  <- gene_parts$introns


p_hits <- findOverlaps(name_dict, proms)
e_hits <- findOverlaps(name_dict, exs)
i_hits <- findOverlaps(name_dict, ints)



final_table <- as.data.frame(name_dict)

final_table$chr <- as.character(seqnames(name_dict))
final_table$start <- start(name_dict)
final_table$end <- end(name_dict)
final_table$width <- width(name_dict)



final_table$is_promoter <- 0
final_table$is_exon <- 0
final_table$is_intron <- 0

final_table$is_promoter[queryHits(p_hits)] <- 1
final_table$is_exon[queryHits(e_hits)] <- 1
final_table$is_intron[queryHits(i_hits)] <- 1

final_table$is_intergenic <-
  ifelse(
    final_table$is_promoter == 0 &
      final_table$is_exon == 0 &
      final_table$is_intron == 0,
    1,
    0
  )

#tss

near_idx <- nearest(name_dict, proms)

final_table$nearest_transcript <- proms$tx_name[near_idx]
final_table$dist_to_tss <- distance(name_dict, proms[near_idx])

#annot

tss <- promoters(proms, upstream = 0, downstream = 1)

near <- distanceToNearest(name_dict, tss)

idx <- subjectHits(near)

site_pos <- start(name_dict)
tss_pos <- start(tss[idx])

strand_dir <- ifelse(as.character(strand(tss[idx])) == "+", 1, -1)

signed_dist <- (site_pos - tss_pos) * strand_dir

final_table$signed_dist_to_tss <- signed_dist
final_table$is_upstream <- signed_dist < 0



write.csv(
  final_table,
  file = "/scratch/project_2010912/ezel/rSalmon/salmon-main/data_files/segment_overlap_results/disagreement_regions/disagree_annot.csv",
  row.names = FALSE
)