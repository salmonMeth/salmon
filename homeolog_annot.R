library(readxl)

#paths
path_output_annot="/scratch/project_2010912/ezel/diff_meth_Full_Annot.csv"

#we get this from the supplement of the paper
#https://pmc.ncbi.nlm.nih.gov/articles/PMC8127823/
homeol_path <- "/scratch/project_2010912/ezel/synteny_AtlanticSalmon_Ssal_v3.1.tsv"
#now we get the homeologs and match them

homeologs <- read.delim(homeol_path)
homeologs <- as.data.frame(homeologs)

dat= read.csv(path_output_annot)

dat$chr_num <- as.integer(dat$chr_ncbi)
##initialize the new cols for our partner chromosome data
dat$in_homeolog_block <- FALSE
dat$homeolog_block_name <- NA

dat$partner_chr <- NA
dat$partner_start <- NA
dat$partner_end <- NA

#we loop thru everything in the homeologs file

for (i in seq_len(nrow(homeologs))) {
  hx <- homeologs$chr_x[i]
  hy <- homeologs$chr_y[i]
  sx <- homeologs$start_x[i]
  ex <- homeologs$end_x[i]
  sy <- homeologs$start_y[i]
  ey <- homeologs$end_y[i]
  bname <- homeologs$homeolog_block[i]
  # --- chr X side ---
  idx_x <- which(
    dat$chr_num == hx &
      dat$start >= sx &
      dat$start <= ex
  )
  if (length(idx_x) > 0) {
    dat$in_homeolog_block[idx_x] <- TRUE
    dat$homeolog_block_name[idx_x] <- bname
    dat$partner_chr[idx_x] <- hy
    dat$partner_start[idx_x] <- sy
    dat$partner_end[idx_x] <- ey
  }
  # --- chr Y side ---
  idx_y <- which(
    dat$chr_num == hy &
      dat$start >= sy &
      dat$start <= ey
  )
  if (length(idx_y) > 0) {
    dat$in_homeolog_block[idx_y] <- TRUE
    dat$homeolog_block_name[idx_y] <- bname
    dat$partner_chr[idx_y] <- hx
    dat$partner_start[idx_y] <- sx
    dat$partner_end[idx_y] <- ex
  }
}

path_output_homeo="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/homeologs/mk_with_homeos.csv"
write.csv(dat, path_output_homeo, row.names = FALSE)