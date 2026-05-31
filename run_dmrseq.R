library(bsseq)
library(dmrseq)
library(BiocParallel)

bs <- readRDS("bsseq_object.rds")
bp <- MulticoreParam(workers = 4)

dmrs <- dmrseq(
  bs = bs,
  testCovariate = "group",
  BPPARAM = bp,
  maxPerms = 10   # start small for testing
)

saveRDS(dmrs, "dmrs.rds")