#to get some stats fro the updated homeolog file
#and to get the info about promoters, tss distance etc
library(Rsamtools)
library(data.table)

path="synteny_AtlanticSalmon_Ssal_v3.1.tsv"
path_genome = "/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa"

#
fai <- scanFaIndex(path_genome)
genome_size <- sum(seqlengths(fai), na.rm = TRUE)

df <- read.table(
  path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)
setDT(df)

# block lengths
df[, len_x := end_x - begin_x]
df[, len_y := end_y - begin_y]

# per-chromosome totals
cov_x <- df[, .(bp = sum(len_x, na.rm = TRUE)), by = chromosome_x]
cov_y <- df[, .(bp = sum(len_y, na.rm = TRUE)), by = chromosome_y]

total_x_bp <- sum(cov_x$bp)
total_y_bp <- sum(cov_y$bp)

prop_x <- total_x_bp / genome_size
prop_y <- total_y_bp / genome_size

prop_x
prop_y
#################



gr_x <- GRanges(
  seqnames = df$chromosome_x,
  ranges = IRanges(df$begin_x, df$end_x)
)

gr_y <- GRanges(
  seqnames = df$chromosome_y,
  ranges = IRanges(df$begin_y, df$end_y)
)

# reduce merges overlaps
gr_x_red <- reduce(gr_x)
gr_y_red <- reduce(gr_y)

mapped_x <- sum(width(gr_x_red))
mapped_y <- sum(width(gr_y_red))

prop_x_unique <- mapped_x / genome_size
prop_y_unique <- mapped_y / genome_size


#get the tss and promoters
library(rtracklayer)


path_annot_genes ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
gff <- import(path_annot_genes)

genes <- gff[gff$type == "gene"]

tss <- resize(genes, width = 1, fix = "start")  # strand-aware TSS
###


library(GenomicRanges)
library(GenomeInfoDb)
path_annot_genes ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
final_table=read.csv("diff_meth_updated_homeologs_8thJune.csv")
name_dict=final_table
site_gr <- GRanges(
  seqnames = name_dict$chr_num,
  ranges = IRanges(start = name_dict$start,
                   end = name_dict$end)
)
gene=txdbmaker::makeTxDbFromGFF(path_annot_genes)
gene_gr <- genes(gene)
exon_gr <- exons(gene)
intron_gr <- unlist(intronsByTranscript(gene))

gene_parts <- list(
  promoters = promoters(gene, upstream=2000, downstream=200),
  exons     = exon_gr,
  introns   = intron_gr
)
prom_gr <- gene_parts$promoters

tss <- promoters(prom_gr, upstream = 0, downstream = 1)
promoter_region <- promoters(prom_gr, upstream = 2000, downstream = 200)
olap_prom <- findOverlaps(site_gr, promoter_region)
final_table$in_promoter <- FALSE
final_table$in_promoter[unique(queryHits(olap_prom))] <- TRUE
near_idx <- nearest(site_gr, prom_gr)
final_table$nearest_transcript <- prom_gr$tx_name[near_idx]
final_table$dist_to_tss <- distance(site_gr, prom_gr[near_idx])
near_tss <- distanceToNearest(site_gr, tss)
idx <- subjectHits(near_tss)
site_pos <- start(site_gr)
tss_pos <- start(tss[idx])
strand_dir <- ifelse(as.character(strand(tss[idx])) == "+", 1, -1)
signed_dist <- (site_pos - tss_pos) * strand_dir
final_table$signed_dist_to_tss <- signed_dist
final_table$is_upstream <- signed_dist < 0
hist(signed_dist, breaks = 100)


write.csv(final_table,"diff_meth_with_promoter_10thJune.csv")