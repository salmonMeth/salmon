library(GenomicRanges)



path_diff_cpg="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_updated_homeologs_8thJune.csv"

path_diff_seg="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_segments_difference_10.csv"

diff_cpg=read.csv(path_diff_cpg)

diff_seg=read.csv(path_diff_seg)

diff_cpg$bonferroni_pvalue <- p.adjust(
  diff_cpg$pvalue,
  method = "bonferroni"
)
diff_seg$bonferroni_pvalue <- p.adjust(
  diff_seg$pvalue,
  method = "bonferroni"
)


#find the matches

cpg_gr <- GRanges(
  seqnames = diff_cpg$chr,
  ranges = IRanges(
    start = diff_cpg$start,
    end = diff_cpg$end
  )
)

seg_gr <- GRanges(
  seqnames = diff_seg$chr,
  ranges = IRanges(
    start = diff_seg$start,
    end = diff_seg$end
  )
)

# Find CpGs overlapping segments
hits <- findOverlaps(cpg_gr, seg_gr)

matching_pairs <- cbind(
  diff_cpg[queryHits(hits), c("chr", "start", "end", "pvalue", "qvalue")],
  setNames(
    diff_seg[subjectHits(hits), c("chr", "start", "end", "pvalue", "qvalue")],
    paste0("seg_", c("chr", "start", "end", "pvalue", "qvalue"))
  )
)

quantile(diff_cpg$qvalue,probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9))
quantile(diff_seg$qvalue,probs=c(.1,.2,.3,.4,.5,.6,.7,.8,.9))
path_match= "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_seg_cpg_matches.csv"


write.csv(matching_pairs, 
          file = path_match, 
          row.names = FALSE)