#associate dmrs to genes if they are within 5000 bp of the genes coding region

#
path_dmr="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_segs.csv"

path_gff="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"

library(rtracklayer)
library(GenomicRanges)
library(dplyr)


dmrs <- read.csv(path_dmr)

gff <- import(path_gff)

# Genes
genes <- gff[gff$type == "gene"]

# CDS
cds <- gff[gff$type == "CDS"]

# Chromosome mapping
chr_map <- setNames(
  paste0("NC_", sprintf("%06d", 59442:59470), ".1"),
  as.character(1:29)
)

# Keep only CDS records on chromosomes 1-29
cds_chr <- cds[
  as.character(seqnames(cds)) %in% names(chr_map)
]

# Get the chromosome levels actually wanted
wanted_levels <- names(chr_map)

# Explicitly restrict seqlevels
seqlevels(cds_chr) <- wanted_levels

seqlevels(cds_chr) <- unname(chr_map)
dmr_gr <- GRanges(
  seqnames = dmrs$chr,
  ranges = IRanges(
    start = dmrs$start,
    end = dmrs$end
  )
)

chr_map <- setNames(
  paste0("NC_", sprintf("%06d", 59442:59470), ".1"),
  as.character(1:29)
)

seqlevels(cds_chr) <- unname(chr_map)

cds_5kb <- cds_chr

start(cds_5kb) <- pmax(1, start(cds_5kb) - 5000)

hits <- findOverlaps(dmr_gr, cds_5kb)
end(cds_5kb) <- end(cds_5kb) + 5000


hits <- findOverlaps(dmr_gr, cds_5kb)

# CDS -> transcript
cds_map <- data.frame(
  cds_index = seq_along(cds_chr),
  transcript_id = as.character(cds_chr$Parent),
  stringsAsFactors = FALSE
)

# Remove the "transcript:" prefix
cds_map$transcript_id <- sub(
  "^transcript:",
  "",
  cds_map$transcript_id
)

mrna <- gff[gff$type == "mRNA"]
mrna <- gff[gff$type == "mRNA"]

tx_map <- data.frame(
  transcript_id = sub("^transcript:", "", mrna$ID),
  gene_id = sub("^gene:", "", as.character(mrna$Parent)),
  stringsAsFactors = FALSE
)

head(tx_map)

cds_map <- data.frame(
  cds_index = seq_along(cds_chr),
  transcript_id = sub(
    "^transcript:",
    "",
    as.character(cds_chr$Parent)
  ),
  stringsAsFactors = FALSE
)

cds_map <- cds_map %>%
  left_join(tx_map, by = "transcript_id")

cds_map <- cds_map %>%
  filter(!is.na(gene_id))

dmr_gene <- data.frame(
  dmr_index = queryHits(hits),
  cds_index = subjectHits(hits)
)
dmr_gene <- dmr_gene %>%
  left_join(cds_map, by = "cds_index")
dmr_gene <- dmr_gene %>%
  mutate(
    chr = as.character(seqnames(dmr_gr))[dmr_index],
    dmr_start = start(dmr_gr)[dmr_index],
    dmr_end = end(dmr_gr)[dmr_index]
  )

dmr_gene <- dmr_gene %>%
  distinct(
    dmr_index,
    gene_id,
    .keep_all = TRUE
  )

dmr_gene$cds_start <- start(cds_chr)[dmr_gene$cds_index]
dmr_gene$cds_end   <- end(cds_chr)[dmr_gene$cds_index]

dmr_gene$distance_to_CDS <- pmax(
  dmr_gene$cds_start - dmr_gene$dmr_end,
  dmr_gene$dmr_start - dmr_gene$cds_end,
  0
)
dmr_gene_final <- dmr_gene %>%
  select(
    dmr_index,
    chr,
    dmr_start,
    dmr_end,
    gene_id,
    transcript_id,
    cds_start,
    cds_end,
    distance_to_CDS
  )