args <- commandArgs(trailingOnly = TRUE)
target_chr <- args[1]

idx <- chr == target_chr

library(bsseq)
library(dmrseq)

bs_chr <- BSseq(
  chr = chr[idx],
  pos = pos[idx],
  M = M[idx, ],
  Cov = Cov[idx, ],
  sampleNames = colnames(M)
)

pData(bs_chr)$group <- factor(treat)

res <- dmrseq(bs_chr, testCovariate = "group")

saveRDS(res, paste0("dmr_", target_chr, ".rds"))