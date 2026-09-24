library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# --------------------------------------------------
# Files
# --------------------------------------------------

bed_file <- "/path/to/your/segments.bed"
fai_file <- "/scratch/project_2010912/ezel/Salmon/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.fai"


# --------------------------------------------------
# Read chromosome lengths from FAI
# --------------------------------------------------

fai <- read.table(
  fai_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(fai)[1:2] <- c("chr", "length")

chr_lengths <- fai %>%
  select(chr, length) %>%
  filter(grepl("^NC_", chr)) %>%
  mutate(
    chr_num = match(
      chr,
      paste0("NC_", 59441 + 1:29, ".1")
    )
  ) %>%
  filter(!is.na(chr_num)) %>%
  select(chr_num, chr, chr_length = length)

chr_lengths
