########################
#######ANNOTATE THE DIFFERENTIALLY METHYLATED GENES
#we get the salmon gene data from salmoBASE
#we get the gene naming formats to match with our methylation reports,
#and we add that naming convention to our DMR db 


#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
###
source("~/rSalmon/diffMeth.R")
##paths
path_annot = "/scratch/project_2010912/ezel/salmon_assembly_report.txt"
#for the gene info
path_annot_genes ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"

#download the annot file
download.file("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/237/065/GCF_905237065.1_Ssal_v3.1/GCF_905237065.1_Ssal_v3.1_assembly_report.txt",
              destfile = path_annot)
#read the annot data
annot_data = read.table(path_annot,
                        sep = "\t",
                        comment.char = "#",
                        fill = TRUE,
                        header = FALSE,
                        stringsAsFactors = FALSE)

#get the V1=chromosome names, V7= RefSeq genomic accession number
maps = annot_data[, c(7, 1)]
colnames(maps) = c("RefSeq", "ChrName")
#get the chromosomes, get rid of the "contig"s
maps = maps[grep("NC_", maps$RefSeq), ]
# lookup for gene name matching with the report
#lookup has "ssa01" "ssa02" etc
lookup = maps$ChrName
#pair them with RefSeqs, NC_059442.1 NC_059443.1 "ssa01"     "ssa02 etc
names(lookup) = maps$RefSeq
#remove ssa's
lookup = gsub("ssa", "", lookup)
names(lookup) =names(lookup)

#we turn our 'methylDiff' object of significant diff methylated regions to a df
alldif_df = as.data.frame(all.diff) 
alldif_raw = getData(alldif_df)
alldif_raw$chr_num = lookup[alldif_raw$chr]
alldif_raw$chr_num = as.character(alldif_raw$chr_num)
#we added the "ssa__" format chromosome number to our diff-methy dataset

#convert our diff-meth object to GRanges obj
sitesGr = as(all.diff, "GRanges")
#rename our ssa genes to NCBI format
#we do all these to turn ssa into simply 1,2 etc
map_names = mapSeqlevels(seqlevels(sitesGr), "NCBI")
refs = colnames(map_names)
trg = as.character(1:length(refs))
names(trg) = refs
#name_dict is still a 'GRanges' obj
name_dict =renameSeqlevels(sitesGr, trg)

#creat intron,exon etc mappings
gene=txdbmaker::makeTxDbFromGFF(path_annot_genes)
#gene=genes(gene)
gene_gr <- genes(gene)
exon_gr <- exons(gene)
intron_gr <- unlist(intronsByTranscript(gene))
#WE MIGHT WANT TO MAKE SURE THESE UPST-DOWNST VALUES ARE GOOD
gene_parts <- list(
  promoters = promoters(gene, upstream=2000, downstream=200),
  exons     = exon_gr,
  introns   = intron_gr
)

proms <- gene_parts$promoters
exs   <- gene_parts$exons
ints  <- gene_parts$introns
p_hits <- findOverlaps(name_dict, proms)
e_hits <- findOverlaps(name_dict, exs)
i_hits <- findOverlaps(name_dict, ints)

##
# convert GRanges  into a standard table
final_table <- as.data.frame(name_dict)
# initialize 
final_table$is_promoter <- 0
final_table$is_exon     <- 0
final_table$is_intron   <- 0
#fill
final_table$is_promoter[queryHits(p_hits)] <- 1
final_table$is_exon[queryHits(e_hits)]     <- 1
final_table$is_intron[queryHits(i_hits)]   <- 1
# intergenic
final_table$is_intergenic <- ifelse(
  final_table$is_promoter == 0 & final_table$is_exon == 0 & final_table$is_intron == 0,
  1, 0
)

#get the distance to the nearest promoter
near_idx = nearest(name_dict, proms)
final_table$nearest_transcript = proms$tx_name[near_idx]
final_table$dist_to_tss = distance(name_dict, proms[near_idx])
#NOTE, we could limit this to upstream TSS?
#we can get whether the sequence we are looking at is uo or downsteram from the tss
#get the index of the TSS and whether it is on the + strand or -
tss = promoters(proms, upstream = 0, downstream = 1)
#finds the closest TSS(subjectHits) to the sites in our diffmeth results
near = distanceToNearest(name_dict, tss)
#stores the ids of the TSS hits we get
idx = subjectHits(near)
site_pos= start(name_dict)
tss_pos = start(tss[idx])
strand_dir = ifelse(as.character(strand(tss[idx])) == "+", 1, -1)
signed_dist = (site_pos - tss_pos) * strand_dir
hist(signed_dist, breaks = 100)
final_table$is_upstream = signed_dist < 0
#NOTE, we might want to visualize which are downstream/upstream 
#AND check their methylation levels


#write.csv(final_table, "SalmonAnnot.csv", row.names = FALSE)

write.csv(final_table, "SalmonAnnotFull.csv", row.names = FALSE)


