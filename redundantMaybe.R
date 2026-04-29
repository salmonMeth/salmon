#problem: we dont have geneDf
#########???????????????????????
###
#find the missing code
gene=txdbmaker::makeTxDbFromGFF(path_annot_genes)
gene=genes(gene)
gene_df = as.data.frame(gene)
gene_df$seqnames = as.character(gene_df$seqnames)

allDiffFinal <- merge(alldif_raw,
                      geneDf,
                      by.x = "chr_num",
                      by.y = "seqnames")
#NOW this version will map to way too many genes, so if we have chromose 28 for example,
#it will map to all genes on that so we have to restrict that
annot <- allDiffFinal[allDiffFinal$start.x >= allDiffFinal$start.y &
                        allDiffFinal$start.x <= allDiffFinal$end.y, ]
