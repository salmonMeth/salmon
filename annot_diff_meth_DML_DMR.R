library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(dplyr)

path_annot <- "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"

txdb <- txdbmaker::makeTxDbFromGFF(path_annot)


########
path_mk_bases= "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_mk.csv"
path_mk_segs ="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/	diff_meth_segs.csv.csv"
path_dss_bases= "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_segs.csv"

mk_bases= read.csv(path_mk_bases)
mk_segs= read.csv(path_mk_segs)
dss_bases=read.csv(path_dss_bases)

################################
################################
################################
# ANNOTATE THE SEGMENTS
# WE STILL NEED THE DATAFRAME WITH THE BASES TO GET THE CHR NUMBER MAPPING(NCIB TO INTEGER FORMAT) RIGHT
## gotta map the NCIB chr names to integers
## Build mapping from the RefSeq chromosome names
## for mk_segs only
annotate_dm_reg= function(mk_bases,mk_segs, path_out){
ref_chr <- sort(unique(mk_bases$chr))
ens_chr <- as.character(seq_along(ref_chr))
chr_map <- setNames(ens_chr, ref_chr)
mk_segs$chr_ens <- chr_map[mk_segs$chr]
#get the annotations
genes_gr <- genes(txdb)
tx_gr <- transcripts(txdb)
promoters_gr <- promoters(
  txdb,
  upstream = 2000,
  downstream = 500
)

tss_gr <- promoters(
  txdb,
  upstream = 0,
  downstream = 1
)

exons_gr <- exons(txdb)

introns_gr <- unlist(
  intronsByTranscript(txdb),
  use.names = FALSE
)

#create the granges objects
#for the mk_segs
dm_seg_gr <- GRanges(
  seqnames = mk_segs$chr,
  ranges = IRanges(
    start = mk_segs$start,
    end   = mk_segs$end
  )
)
mcols(dm_seg_gr)$chr_ens <- chr_map[as.character(seqnames(dm_seg_gr))]
mcols(dm_seg_gr)$meth.diff <- mk_segs$meth.diff
mcols(dm_seg_gr)$qvalue    <- mk_segs$qvalue

if("mean_meth" %in% names(mk_segs))
  mcols(dm_seg_gr)$mean_meth <- mk_segs$mean_meth

if("direction" %in% names(mk_segs))
  mcols(dm_seg_gr)$direction <- mk_segs$direction

##create a new granges obj
dm_seg_gr_ens <- GRanges(
  seqnames = Rle(as.character(dm_seg_gr$chr_ens)),
  ranges = ranges(dm_seg_gr),
  strand = strand(dm_seg_gr)
)

mcols(dm_seg_gr_ens) <- mcols(dm_seg_gr)
##get the hits
gene_hits <- findOverlaps(dm_seg_gr_ens, genes_gr, ignore.strand = TRUE)

promoter_hits <- findOverlaps(dm_seg_gr_ens, promoters_gr, ignore.strand = TRUE)

tss_hits <- findOverlaps(dm_seg_gr_ens, tss_gr, ignore.strand = TRUE)

exon_hits <- findOverlaps(dm_seg_gr_ens, exons_gr, ignore.strand = TRUE)

intron_hits <- findOverlaps(dm_seg_gr_ens, introns_gr, ignore.strand = TRUE)
###fill the results dataframe
dm_annot <- as.data.frame(dm_seg_gr_ens)
dm_annot$gene <- FALSE
dm_annot$promoter <- FALSE
dm_annot$tss <- FALSE
dm_annot$exon <- FALSE
dm_annot$intron <- FALSE

dm_annot$gene[unique(queryHits(gene_hits))] <- TRUE
dm_annot$promoter[unique(queryHits(promoter_hits))] <- TRUE
dm_annot$tss[unique(queryHits(tss_hits))] <- TRUE
dm_annot$exon[unique(queryHits(exon_hits))] <- TRUE
dm_annot$intron[unique(queryHits(intron_hits))] <- TRUE

gene_ids <- rep(NA_character_, length(dm_seg_gr_ens))

gene_list <- split(
  genes_gr$gene_id[subjectHits(gene_hits)],
  queryHits(gene_hits)
)

gene_ids[as.integer(names(gene_list))] <-
  sapply(gene_list, function(x)
    paste(unique(x), collapse = ";")
  )

dm_annot$gene_id <- gene_ids
gene_names <- rep(NA_character_, length(dm_seg_gr_ens))

if ("gene_name" %in% colnames(mcols(genes_gr))) {
  
  gene_name_list <- split(
    genes_gr$gene_name[subjectHits(gene_hits)],
    queryHits(gene_hits)
  )
  
  gene_names[as.integer(names(gene_name_list))] <-
    sapply(gene_name_list, function(x)
      paste(unique(x), collapse = ";")
    )
  
} else if ("Name" %in% colnames(mcols(genes_gr))) {
  
  gene_name_list <- split(
    genes_gr$Name[subjectHits(gene_hits)],
    queryHits(gene_hits)
  )
  
  gene_names[as.integer(names(gene_name_list))] <-
    sapply(gene_name_list, function(x)
      paste(unique(x), collapse = ";")
    )
}

dm_annot$gene_name <- gene_names

write.csv(dm_annot, path_out)
}

########################
########################
########################
#NOW WE ANNOTATE THE BASES, BOTH DSS AND METHYLKIT

library(GenomicRanges)
library(GenomicFeatures)

#function to get the annotations for diff meth bases
#flexible to accomodate both methlkit and dss style diff methy base results
annotate_dm_bases<- function(mk_bases, txdb, path_out){
  tx_chr <- seqlevels(txdb)
  tx_chr_num <- tx_chr[grepl("^[0-9]+$", tx_chr)]
  tx_chr_num <- as.character(sort(as.numeric(tx_chr_num)))
  # RefSeq chromosomes from methylation data
  ref_chr <- unique(mk_bases$chr)
  ref_chr_main <- ref_chr[grepl("^NC_", ref_chr)]
  ref_chr_main <- sort(ref_chr_main)

  # RefSeq -> NCBI mapping
  chr_map <- setNames(
    tx_chr_num,
    ref_chr_main
  )
  ############################
  # Create DML GRanges
  ############################
  if("start" %in% colnames(mk_bases)){
    
    dm_start <- mk_bases$start
    dm_end   <- mk_bases$end
    
  } else if("pos" %in% colnames(mk_bases)){
    
    # DSS output: single CpG position
    dm_start <- mk_bases$pos
    dm_end   <- mk_bases$pos
    
  } else {
    
    stop("Cannot find start/end or pos columns")
    
  }
  
  
  # Strand handling
  if("strand" %in% colnames(mk_bases)){
    dm_strand <- mk_bases$strand
  } else {
    dm_strand <- "*"
  }
  
  
  dm_base_gr <- GRanges(
    seqnames = mk_bases$chr,
    ranges = IRanges(
      start = dm_start,
      end   = dm_end
    ),
    strand = dm_strand
  )#################
  
  
  # Keep original RefSeq chromosome
  # Add NCBI chromosome as metadata
  mcols(dm_base_gr)$chr_ncbi <- chr_map[
    as.character(seqnames(dm_base_gr))
  ]
  
  ############################
  # Add methylation statistics
  ############################
  
  if("pvalue" %in% colnames(mk_bases)){
    mcols(dm_base_gr)$pvalue <- mk_bases$pvalue
  }
  
  if("pval" %in% colnames(mk_bases)){
    mcols(dm_base_gr)$pvalue <- mk_bases$pval
  }
  
  
  if("qvalue" %in% colnames(mk_bases)){
    mcols(dm_base_gr)$qvalue <- mk_bases$qvalue
  }
  
  if("fdr" %in% colnames(mk_bases)){
    mcols(dm_base_gr)$qvalue <- mk_bases$fdr
  }
  
  
  if("meth.diff" %in% colnames(mk_bases)){
    mcols(dm_base_gr)$meth.diff <- mk_bases$meth.diff
  }
  
  if("diff" %in% colnames(mk_bases)){
    mcols(dm_base_gr)$meth.diff <- mk_bases$diff
  }
  ###############
 
  # Check mapping
  if(any(is.na(dm_base_gr$chr_ncbi))){
    warning(
      "Some chromosomes could not be mapped to NCBI chromosome names"
    )
  }
  ############################
  # Create GRanges for overlap
  ############################
  
  dm_base_gr_ncbi <- GRanges(
    seqnames = as.character(mcols(dm_base_gr)$chr_ncbi),
    ranges = ranges(dm_base_gr),
    strand = strand(dm_base_gr)
  )
  
  mcols(dm_base_gr_ncbi) <- mcols(dm_base_gr)

  ############################
  # Get annotations
  ############################
  
  genes_gr <- genes(txdb)
  
  promoters_gr <- promoters(
    txdb,
    upstream = 2000,
    downstream = 500
  )
  
  tss_gr <- promoters(
    txdb,
    upstream = 0,
    downstream = 1
  )
  
  exons_gr <- exons(txdb)
  
  introns_gr <- unlist(
    intronsByTranscript(txdb),
    use.names = FALSE
  )

  ############################
  # Find overlaps
  ############################
  
  gene_hits <- findOverlaps(
    dm_base_gr_ncbi,
    genes_gr,
    ignore.strand = TRUE
  )
  
  promoter_hits <- findOverlaps(
    dm_base_gr_ncbi,
    promoters_gr,
    ignore.strand = TRUE
  )
  
  tss_hits <- findOverlaps(
    dm_base_gr_ncbi,
    tss_gr,
    ignore.strand = TRUE
  )
  
  exon_hits <- findOverlaps(
    dm_base_gr_ncbi,
    exons_gr,
    ignore.strand = TRUE
  )
  
  intron_hits <- findOverlaps(
    dm_base_gr_ncbi,
    introns_gr,
    ignore.strand = TRUE
  )
  ############################
  # Create annotation table
  ############################
  
  dm_annot <- as.data.frame(dm_base_gr)
  
  dm_annot$gene <- FALSE
  dm_annot$promoter <- FALSE
  dm_annot$tss <- FALSE
  dm_annot$exon <- FALSE
  dm_annot$intron <- FALSE
  
  
  dm_annot$gene[
    unique(queryHits(gene_hits))
  ] <- TRUE
  
  dm_annot$promoter[
    unique(queryHits(promoter_hits))
  ] <- TRUE
  
  dm_annot$tss[
    unique(queryHits(tss_hits))
  ] <- TRUE
  
  dm_annot$exon[
    unique(queryHits(exon_hits))
  ] <- TRUE
  
  dm_annot$intron[
    unique(queryHits(intron_hits))
  ] <- TRUE

  ############################
  # Add gene IDs
  ############################
  
  gene_ids <- rep(
    NA_character_,
    length(dm_base_gr)
  )
  gene_list <- split(
    genes_gr$gene_id[subjectHits(gene_hits)],
    queryHits(gene_hits)
  )
  
  
  gene_ids[
    as.integer(names(gene_list))
  ] <- sapply(
    gene_list,
    function(x)
      paste(unique(x), collapse=";")
  )
  
  
  dm_annot$gene_id <- gene_ids
  
  
  
  ############################
  # Add gene names
  ############################
  
  gene_names <- rep(
    NA_character_,
    length(dm_base_gr)
  )
  
  
  if("gene_name" %in% colnames(mcols(genes_gr))){
    
    gene_name_list <- split(
      genes_gr$gene_name[subjectHits(gene_hits)],
      queryHits(gene_hits)
    )
    
    gene_names[
      as.integer(names(gene_name_list))
    ] <- sapply(
      gene_name_list,
      function(x)
        paste(unique(x), collapse=";")
    )
    
    
  } else if("Name" %in% colnames(mcols(genes_gr))){
    
    gene_name_list <- split(
      genes_gr$Name[subjectHits(gene_hits)],
      queryHits(gene_hits)
    )
    
    gene_names[
      as.integer(names(gene_name_list))
    ] <- sapply(
      gene_name_list,
      function(x)
        paste(unique(x), collapse=";")
    )
  }
  
  
  dm_annot$gene_name <- gene_names
  
  
  
  ############################
  # Save output
  ############################
  
  write.csv(
    dm_annot,
    path_out,
    row.names = FALSE
  )
  
  
  return(dm_annot)
}

annotate_dm_bases(mk_bases, txdb, "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/methylkit_diff_meth_bases2.csv")

annotate_dm_bases(dss_bases, txdb, "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/dss_diff_meth_bases.csv")

#ANNOTATION DONE
####################################
####################################
#lets check the results, using the csv files

path_dml_mk="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/methylkit_diff_meth_bases.csv"
path_dml_dss="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/dss_diff_meth_bases.csv"


dml_mk <- read.csv(path_dml_mk,stringsAsFactors = FALSE)

dml_dss <- read.csv(path_dml_dss,stringsAsFactors = FALSE)
feature_summary <- function(x){
  
  features <- c(
    "gene",
    "promoter",
    "tss",
    "exon",
    "intron"
  )
  
  data.frame(
    feature = features,
    n = sapply(
      features,
      function(f) sum(x[[f]])
    ),
    percent = sapply(
      features,
      function(f) mean(x[[f]]) * 100
    ),
    total = nrow(x)
  )
}


summary_mk <- feature_summary(dml_mk)
summary_dss <- feature_summary(dml_dss)


summary_mk$method <- "Methylkit"
summary_dss$method <- "DSS"

feature_comparison <- merge(
  summary_mk,
  summary_dss,
  by = "feature",
  suffixes = c("_Methylkit", "_DSS")
)

##########
##########
#lets compare the breakdown of the genomic features identified by the different methods
#feature_comparison
#good plot
library(tidyr)
library(ggplot2)


plot_data <- feature_comparison %>%
  dplyr::select(
    feature,
    percent_Methylkit, n_Methylkit,
    percent_DSS, n_DSS
  ) %>%
  tidyr::pivot_longer(
    cols = -feature,
    names_to = c(".value", "method"),
    names_sep = "_"
  )


p <- ggplot(
  plot_data,
  aes(
    x = feature,
    y = percent,
    fill = method
  )
) +
  geom_col(position = position_dodge(0.9)) +
  geom_text(
    aes(label = n),
    position = position_dodge2(width = 0.9),
    vjust = 1.5,
    size = 2.8,
    colour = "black",  fontface = "bold"
  )+
  ylab("% of DMLs") +
  xlab("") +
  ggtitle("Genomic feature distribution of differentially methylated bases") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    )
  )


ggsave(
  filename = "/scratch/project_2010912/ezel/salmon/data_files/use_in_thesis/DML_feature_comparison.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)
####
####
common_genes <- intersect(
  genes_mk,
  genes_dss
)
dss_only <- setdiff(
  genes_dss,
  genes_mk
)
methylKit_only <- setdiff(
  genes_mk,
  genes_dss
)

write.csv(
  data.frame(gene_id=common_genes),
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/common_genes.csv",
  row.names=FALSE
)
write.csv(
  data.frame(gene_id=methylKit_only),
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/method1_specific_genes.csv",
  row.names=FALSE
)
write.csv(
  data.frame(gene_id=dss_only),
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/DSS_specific_genes.csv",
  row.names=FALSE
)

########################
########################
#for the segments by methylKit

# Read annotated DMR segments
dmr_seg <- read.csv(
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/methylkit_diff_meth_segs_annotated.csv",
  stringsAsFactors = FALSE
)


# Summarize features
feature_summary <- function(x){
  
  features <- c(
    "gene",
    "promoter",
    "tss",
    "exon",
    "intron"
  )
  
  data.frame(
    feature = features,
    n = sapply(
      features,
      function(f) sum(x[[f]])
    ),
    percent = sapply(
      features,
      function(f) mean(x[[f]]) * 100
    )
  )
}


dmr_summary <- feature_summary(dmr_seg)

p <- ggplot(
  dmr_summary,
  aes(
    x = feature,
    y = percent
  )
) +
  geom_col(fill = "#DA70D6") +
  geom_text(
    aes(label = n),
    vjust = 1.3,       # puts text inside the bars
    colour = "white",
    size = 3,
    fontface = "bold"
  ) +
  ylab("% of DMRs") +
  xlab("") +
  ggtitle("Genomic feature distribution of differentially methylated regions") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    )
  )


ggsave(
  "/scratch/project_2010912/ezel/salmon/data_files/use_in_thesis/DMR_segment_feature_distribution.png",
  plot = p,
  width = 8,
  height = 5,
  dpi = 300
)
