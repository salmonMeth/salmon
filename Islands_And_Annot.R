#####libs
BiocManager::install("rtracklayer")
BiocManager::install("Biostrings")
library(Biostrings)
library(rtracklayer)
library(GenomicRanges)
library(readxl)

#paths
#path = "C:/Users/ezele/Downloads/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz"
#path ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz"
#path2 ="C:/Users/ezele/Downloads/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa"
#path2 ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa"

#pathDiff ="C:/Users/ezele/Downloads/diff_meth_results_T1_vs_T2_Full.csv"
pathDiff="/scratch/project_2010912/ezel/diff_meth_results_T1_vs_T2_Full.csv"

pathIslands="/scratch/project_2010912/ezel/islands.bed"
gff_url <- "https://salmobase.org/datafiles/genomes/AtlanticSalmon/Ssal_v3.1/annotations/Ensembl/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
homeol_path <- "C:/Users/ezele/Downloads/homeologs.xlsx"



#top_level_url <- "https://ftp.ensembl.org/pub/release-106/fasta/salmo_salar/dna/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.gz"
# Download the top_level
#download.file(url = top_level_url, destfile = path, mode = "wb")
#R.utils::gunzip(path,remove = FALSE)

#R.utils::gunzip(filename = path, destname = path2, remove = FALSE, overwrite = TRUE)
#this genome is not needed I am pretty sure
#genome <- readDNAStringSet(path2)

download.file(gff_url, destfile = "salmon.gff.gz", mode = "wb")

gff <- import("salmon.gff.gz")
dat= read.csv(pathDiff)
regions <- GRanges(
  seqnames = dat$chr,
  ranges = IRanges(
    start = dat$start,
    end = dat$end))
regions

cpg=import(pathIslands)
str(cpg)
hits <- findOverlaps(regions, cpg)
##add the cpg island stats to the datatable
dat$in_cpg <- FALSE
dat$cpg_id <- NA
q <- queryHits(hits)
s <- subjectHits(hits)
dat$in_cpg[q] <- TRUE
dat$cpg_id[q] <- mcols(cpg)$name[s]

####now we have added a column to check if our DMNs are in CpG islands
#now we will change the chromosome numbering system to match with the annotations
mcols(regions)$chr_num <- seqnames(regions)
dat$chr_num <- dat$chr


####
#save csv
pathOutput="C:/Users/ezele/Downloads/diff_meth_with_cpg_islands.csv"
write.csv(dat, pathOutput, row.names = FALSE)
# we have no diff meth nucleotides in the 22nd chromosome

##we want to annotate the other chromosome segments first
#we get rid of the contigs etc in the annot file
new_seqnames <- as.integer(sub("\\.1$", "",sub("^NC_0*", "", as.character(seqnames(regions))) )) - 59441
mcols(regions)$chr_num <- new_seqnames
dat$chr_num <- new_seqnames

gff_chr <- gff[as.character(seqnames(gff)) %in% as.character(1:29)]
snps2 <- GRanges(seqnames = dat$chr_num,ranges = IRanges(dat$start, dat$end))
hits <- findOverlaps(snps2, gff_chr)
q <- queryHits(hits)
s <- subjectHits(hits)
annots <- split(as.character(gff_chr$type[s]), q)
dat$annotation <- NA
dat$annotation[as.integer(names(annots))] <-sapply(annots, function(x) paste(unique(x), collapse=";"))
dat$annotation[is.na(dat$annotation)] <- "intergenic"


pathOutputAnnot="/scratch/project_2010912/ezel/diff_meth_Full_Annot.csv"
write.csv(dat, pathOutputAnnot, row.names = FALSE)