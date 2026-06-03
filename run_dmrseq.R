library(bsseq)
library(dmrseq)
library(BiocParallel)

bs <- readRDS("data_files/bsseq_object.rds")
bp <- MulticoreParam(workers = 4)

dmrs <- dmrseq(
  bs = bs,
  testCovariate = "group",
  BPPARAM = bp,
  maxPerms = 10   # start small
)

saveRDS(dmrs, "data_files/dmrs.rds")