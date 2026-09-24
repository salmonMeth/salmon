find_matches <- function(ref, query, tolerance = 500) {
  
  sapply(seq_len(nrow(ref)), function(i) {
    any(
      query$chr == ref$chr[i] &
        abs(query$start - ref$start[i]) <= tolerance &
        abs(query$end - ref$end[i]) <= tolerance
    )
  })
}
match_40 <- find_matches(
  minseg10_maxint100,
  maxint40
)

match_20 <- find_matches(
  minseg10_maxint100,
  minseg20
)

#####################################################
#####################################################

#check how fragmented the ends of the chromosomes are

library(data.table)

base <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation"
path_fai="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.fai"



fai <- fread(path_fai,header = FALSE)

# FASTA index columns:
# V1 = chromosome
# V2 = chromosome length
chr_lengths <- data.frame(
  chromosome = as.character(fai$V1),
  chromosome_length = fai$V2
)

get_terminal_segments <- function(bed_file, chr_lengths, frac = 0.10) {
  filename <- basename(bed_file)
  accession <- sub("_segments\\.bed$", "", filename)
  
  if (!grepl("^NC_[0-9]+\\.1$", accession))
    return(NULL)
  
  accession_number <- as.integer(
    sub("^NC_", "", sub("\\.1$", "", accession))
  )
  
  chromosome <- accession_number - 59441
  
  if (chromosome < 1 || chromosome > 29)
    return(NULL)
  
  chr_length <- chr_lengths$chromosome_length[
    chr_lengths$chromosome == as.character(chromosome)
  ][1]
  
  if (is.na(chr_length))
    stop("Could not find chromosome ", chromosome, " in chr_lengths")
  
  bed <- fread(bed_file, header = FALSE)
  
  terminal_start <- floor(chr_length * (1 - frac))
  
  terminal <- bed[
    bed[[2]] >= terminal_start &
      bed[[3]] <= chr_length
  ]
  
  data.frame(
    chromosome = chromosome,
    accession = accession,
    chromosome_length = chr_length,
    terminal_start = terminal_start,
    terminal_end = chr_length,
    n_terminal_segments = nrow(terminal),
    n_total_segments = nrow(bed),
    fraction_of_segments_terminal = nrow(terminal) / nrow(bed)
  )
}


dirs <- c(
  original = file.path(base, "segmentation_output_G1_5/beds"),
  minseg20 = file.path(base, "segmentation_output_G1_5_minseg_20/beds"),
  maxint40 = file.path(base, "segmentation_output_G1_5_maxInt_40/beds")
)

results_list <- list()

for (set_name in names(dirs)) {
  files <- list.files(
    dirs[[set_name]],
    pattern = "^NC_[0-9]+\\.1_segments\\.bed$",
    full.names = TRUE
  )
  
  for (f in files) {
    x <- get_terminal_segments(f, chr_lengths)
    x$segmentation <- set_name
    x$file <- basename(f)
    results_list[[length(results_list) + 1]] <- x
  }
}

results <- rbindlist(results_list)
setorder(results, chromosome, segmentation)

write.csv(results,"/scratch/project_2010912/ezel/salmon/data_files/segmentation/chr_end_segments.csv", row.names = FALSE)


p1 <- ggplot(
  subset(results, chromosome <= 15),
  aes(x = chromosome, y = fraction_of_segments_terminal, color = segmentation)
) +
  geom_hline(
    yintercept = 0.10,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    size = 3,
    alpha = 0.7
  )+
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Chromosome",
    y = "Fraction of segments in terminal 10%",
    color = "Segmentation"
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_line(color = "grey90")
  )

p2 <- ggplot(
  subset(results, chromosome >= 16),
  aes(x = chromosome, y = fraction_of_segments_terminal, color = segmentation)
) +
  geom_hline(
    yintercept = 0.10,
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    size = 3,
    alpha = 0.7
  )+
  scale_x_continuous(breaks = 16:29) +
  scale_y_continuous(labels = percent) +
  labs(
    x = "Chromosome",
    y = "Fraction of segments in terminal 10%",
    color = "Segmentation"
  ) +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_line(color = "grey90")
  )

ggsave(
  file.path(base, "segments_in_last10pc.png"),
  p1,
  width = 10,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(base, "segments_in_last10pc_2.png"),
  p2,
  width = 10,
  height = 6,
  dpi = 300
)

###############################################
###############################################
#check the ones where at least 20% of all the segments in the whole chromosome 
#are located in the last 10% of the chromosome
library(dplyr)

res=read.csv("/scratch/project_2010912/ezel/salmon/data_files/segmentation/chr_end_segments.csv", row.names = FALSE)

ge.2=res[res$fraction_of_segments_terminal >.2,]
ge.2_filt <- ge.2 %>%
  distinct(accession, terminal_start, terminal_end, .keep_all = TRUE)

write.csv(ge.2_filt,"/scratch/project_2010912/ezel/salmon/data_files/segmentation/chr_end_segments_filtered.csv", row.names = FALSE)
#we should annotate this, maybe look at the cpg islands too
library(data.table)
library(GenomicRanges)
library(IRanges)
library(GenomicFeatures)
library(ggplot2)
library(rtracklayer)
#

ge.2_filt=read.csv("/scratch/project_2010912/ezel/salmon/data_files/segmentation/chr_end_segments_filtered.csv", row.names = FALSE)

path_annot= "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
txdb <- txdbmaker::makeTxDbFromGFF(path_annot)
exons_gr <- exons(txdb)
introns_gr <- unlist(intronsByTranscript(txdb, use.names = TRUE))
promoters_gr <- promoters(txdb, upstream = 1000, downstream = 200)
genes_gr <- genes(txdb)

cpg_gr <- rtracklayer::import("/scratch/project_2010912/ezel/islands.bed")

cpg_chr <- paste0("NC_", sprintf("%06d", 59441 + 1:29), ".1")

cpg_gr <- cpg_gr[
  as.character(seqnames(cpg_gr)) %in% cpg_chr
]

cpg_gr <- GRanges(
  seqnames = as.character(
    as.integer(
      sub(
        "^NC_",
        "",
        sub("\\.1$", "", as.character(seqnames(cpg_gr)))
      )
    ) - 59441
  ),
  ranges = ranges(cpg_gr),
  name = mcols(cpg_gr)$name
)


#we reduce so that if we have continous exons for example
#150-250 and 200-300, we dont count the length as 100+100 but 300-150, so 150
exons_unique <- reduce(exons_gr, ignore.strand = TRUE)
introns_unique <- reduce(introns_gr, ignore.strand = TRUE)
cpg_unique <- reduce(cpg_gr, ignore.strand = TRUE)

annotate_terminal <- function(chr, terminal_start, terminal_end) {
  
  region <- GRanges(
    seqnames = as.character(chr),
    ranges = IRanges(
      start = terminal_start,
      end = terminal_end
    )
  )
  
  region_bp <- width(region)
  
  gene_hits <- findOverlaps(
    region, genes_gr,
    ignore.strand = TRUE
  )
  
  promoter_hits <- findOverlaps(
    region, promoters_gr,
    ignore.strand = TRUE
  )
  
  exon_hits <- findOverlaps(
    region, exons_gr,
    ignore.strand = TRUE
  )
  
  intron_hits <- findOverlaps(
    region, introns_gr,
    ignore.strand = TRUE
  )
  
  cpg_hits <- findOverlaps(
    region, cpg_gr,
    ignore.strand = TRUE
  )
  
  gene_cov <- findOverlaps(
    region, genes_gr,
    ignore.strand = TRUE
  )
  
  gene_bp <- if (length(gene_cov) > 0) {
    sum(width(reduce(
      pintersect(
        region[queryHits(gene_cov)],
        genes_gr[subjectHits(gene_cov)]
      ),
      ignore.strand = TRUE
    )))
  } else 0
  
  promoter_bp <- if (length(promoter_hits) > 0) {
    sum(width(reduce(
      pintersect(
        region[queryHits(promoter_hits)],
        promoters_gr[subjectHits(promoter_hits)]
      ),
      ignore.strand = TRUE
    )))
  } else 0
  
  exon_region_hits <- findOverlaps(
    region, exons_unique,
    ignore.strand = TRUE
  )
  
  exon_bp <- if (length(exon_region_hits) > 0) {
    sum(width(reduce(
      pintersect(
        region[queryHits(exon_region_hits)],
        exons_unique[subjectHits(exon_region_hits)]
      ),
      ignore.strand = TRUE
    )))
  } else 0
  
  intron_region_hits <- findOverlaps(
    region, introns_unique,
    ignore.strand = TRUE
  )
  
  intron_bp <- if (length(intron_region_hits) > 0) {
    sum(width(reduce(
      pintersect(
        region[queryHits(intron_region_hits)],
        introns_unique[subjectHits(intron_region_hits)]
      ),
      ignore.strand = TRUE
    )))
  } else 0
  
  cpg_region_hits <- findOverlaps(
    region, cpg_unique,
    ignore.strand = TRUE
  )
  
  cpg_bp <- if (length(cpg_region_hits) > 0) {
    sum(width(reduce(
      pintersect(
        region[queryHits(cpg_region_hits)],
        cpg_unique[subjectHits(cpg_region_hits)]
      ),
      ignore.strand = TRUE
    )))
  } else 0
  
  data.frame(
    chromosome = chr,
    terminal_start = terminal_start,
    terminal_end = terminal_end,
    region_bp = region_bp,
    
    n_genes = length(unique(subjectHits(gene_hits))),
    gene_bp = gene_bp,
    prop_gene_bp = gene_bp / region_bp,
    
    n_promoters = length(unique(subjectHits(promoter_hits))),
    promoter_bp = promoter_bp,
    prop_promoter_bp = promoter_bp / region_bp,
    
    n_exons = length(unique(subjectHits(exon_hits))),
    exon_bp = exon_bp,
    prop_exon_bp = exon_bp / region_bp,
    
    n_introns = length(unique(subjectHits(intron_hits))),
    intron_bp = intron_bp,
    prop_intron_bp = intron_bp / region_bp,
    
    n_cpg_islands = length(unique(subjectHits(cpg_hits))),
    cpg_bp = cpg_bp,
    prop_cpg_bp = cpg_bp / region_bp
  )
}

terminal_results <- rbindlist(
  lapply(seq_len(nrow(ge.2_filt)), function(i) {
    annotate_terminal(
      chr = ge.2_filt$chromosome[i],
      terminal_start = ge.2_filt$terminal_start[i],
      terminal_end = ge.2_filt$terminal_end[i]
    )
  }),
  fill = TRUE
)
write.csv(terminal_results,"/scratch/project_2010912/ezel/salmon/data_files/segmentation/chr_end_annotations.csv", row.names = FALSE)
#we can annotate the remaining chrs to compare

chr_lengths <- data.frame(
  chromosome = as.character(fai$V1),
  chromosome_length = fai$V2
)
chr_lengths$chromosome <- as.character(chr_lengths$chromosome)

remaining_chr <- chr_lengths[
  chr_lengths$chromosome %in% as.character(1:29) &
    !chr_lengths$chromosome %in% as.character(ge.2_filt$chromosome),
]

remaining_chr$chromosome <- as.integer(remaining_chr$chromosome)

remaining_chr$terminal_start <- floor(
  remaining_chr$chromosome_length * 0.90
) + 1

remaining_chr$terminal_end <- remaining_chr$chromosome_length

remaining_results <- rbindlist(
  lapply(seq_len(nrow(remaining_chr)), function(i) {
    annotate_terminal(
      chr = remaining_chr$chromosome[i],
      terminal_start = remaining_chr$terminal_start[i],
      terminal_end = remaining_chr$terminal_end[i]
    )
  }),
  fill = TRUE
)
write.csv(remaining_results,"/scratch/project_2010912/ezel/salmon/data_files/segmentation/chr_end_annotations_NOT_fragmented_ends.csv", row.names = FALSE)
