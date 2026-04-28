#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)

#paths
path= "/scratch/project_2010912/ezel/methylBase_7x_40ind.txt.bgz"
#to save our diff methylation results
path_diff_meth="/scratch/project_2010912/ezel/diff_meth_results_T1_vs_T2.csv"
#to get the annotations and save them
trg_dir = "/scratch/project_2010912/ezel/"
#for raw data with chromosome numbers
path_annot = paste0(trg_dir, "salmon_assembly_report.txt")
#for the gene info
path_annot_genes ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
#
db = readMethylDB(path)
n=db@num.records
#this is just a random sampling we do to work on a smaller subset of the data first
inds = sample(1:n,1e5)
db= db[inds]
all_ids =db@sample.ids

#differential mehtylation analysis for T1 vs T2
#DIFF METHYLATION ANALYSIS

#find the datapoints from T1
is_t1 = grepl("T1", db@sample.ids)
#if the datapoint is from t1, =0, else, =1
new_treatment = ifelse(is_t1, 0, 1)
data.frame(ID = db@sample.ids, Treatment = new_treatment)
dmf=calculateDiffMeth(db)
#filter the results where the difference is < than 25 and qvalue >.01
all.diff=getMethylDiff(dmf,difference=25,qvalue=0.01,type="all")
write.csv(all.diff,
          file =path_diff_meth,
          row.names = FALSE)
########################
#######ANNOTATE THE DIFFERENTIALLY METHYLATED GENES
#we get the salmon gene data from salmoBASE
#we get the gene naming formats to match with our methylation reports,
#and we add that naming convention to our DMR db 

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
#find the missing code
gene=txdbmaker::makeTxDbFromGFF(pathAnnot)geneDf = as.data.frame(gene)
gene_parts <- list(
  promoters = promoters(gene, upstream=2000, downstream=200),
  exons     = exons(gene),
  introns   = unlist(intronsByTranscript(gene))
)
#convert our diff-meth object to GRanges obj
sitesGr = as(all.diff, "GRanges")
#rename our ssa genes to NCBI format
#we do all these to turn ssa into simply 1,2 etc
map_names = mapSeqlevels(seqlevels(sitesGr), "NCBI")
refs = colnames(map_names)
trg = as.character(1:length(refs))
names(trg) = refs
#name_dict is still a 'GRanges' obj
name_dict =renameSeqlevels(sitesGr, map2)

#creat intron,exon etc mappings

#WE MIGHT WANT TO MAKE SURE THESE UPST-DOWNST VALUES ARE GOOD
gene_parts <- list(
  promoters = promoters(gene, upstream=2000, downstream=200),
  exons     = exons(gene),
  introns   = unlist(intronsByTranscript(gene))
)

proms <- gene_parts$promoters
exs   <- gene_parts$exons
ints  <- gene_parts$introns
p_hits <- findOverlaps(name_dict, proms)
e_hits <- findOverlaps(name_dict, exs)
i_hits <- findOverlaps(name_dict, ints)

##
# convert 42 GRanges  into a standard table
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

near_idx <- nearest(name_dict, proms)
# Extract the Transcript name
final_table$nearest_transcript <- proms$tx_name[near_idx]

final_table$dist_to_tss <- distance(sitesGr2, proms[near_idx])

#########???????????????????????
geneDf$seqnames <- as.character(geneDf$seqnames)


###
allDiffFinal <- merge(alldif_raw,
                      geneDf,
                      by.x = "chr_num",
                      by.y = "seqnames")
#NOW this version will map to way too many genes, so if we have chromose 28 for example,
#it will map to all genes on that so we have to restrict that
annot <- allDiffFinal[allDiffFinal$start.x >= allDiffFinal$start.y &
                        allDiffFinal$start.x <= allDiffFinal$end.y, ]
# we do this
#but then this only gets the parts that are exactly in the gene so our
#dataset of 42 becomes a dataset of 16
# if we want to get the promoters etc we have to do sth else

