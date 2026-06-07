library(bsseq)
library(changepoint)

bs_obj=readRDS( "data_files/bsseq_object.rds")

#
meth <- getMeth(bs_obj, type = "raw")
gr <- rowRanges(bs_obj)

chr <- as.character(seqnames(gr))
pos <- start(gr)
#run the changepoint function over all the samples separately

n_samples <- ncol(meth)
cpt_results <- vector("list", n_samples)
names(cpt_results) <- colnames(meth)
n_samples
for (i in seq_len(n_samples)) {
  x <- meth[, i]
  fit <- cpt.mean(x,
                  method = "PELT",
                  penalty = "MBIC")
  breaks <- cpts(fit)
  cpt_results[[i]] <- data.frame(
    sample = colnames(meth)[i],
    chr = chr[breaks],
    pos = pos[breaks],
    index = breaks
  )
}
