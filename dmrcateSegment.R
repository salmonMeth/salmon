#https://mirror.accum.se/mirror/bioconductor.org/packages/3.12/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf
#DMRcate pipeline to do segmentation
source("~/rSalmon/segmentBSSQ.R")

M <- as.matrix(df[, grep("^numCs_", colnames(df))])
Cov <- as.matrix(df[, grep("^coverage_", colnames(df))])
colnames(M) <- sample_ids
colnames(Cov) <- sample_ids
# coverage filter: at least 5 reads in at least 2 samples
keep_cov <- rowSums(Cov >= 5, na.rm = TRUE) >= 2
# missingness filter: allow up to 20% NA
keep_na <- rowSums(is.na(M)) <= (0.2 * ncol(M))
# final mask
keep <- keep_cov & keep_na
M   <- M[keep, ]
Cov <- Cov[keep, ]
df  <- df[keep, ]
#JUST FOR NOW WE CAN SET THE NAS TO 0 AND THEN MAYBE USE A BETTER METHOD:
#
M[is.na(M)] <- 0
Cov[is.na(Cov)] <- 0
beta <- (M / Cov) * 100
# avoid log(0)
beta[beta == 0] <- 0.001
beta[beta == 100] <- 99.999
Mval <- log2(beta / (100 - beta))
chr  <- df$chr
pos  <- df$start
strand <- df$strand
ord <- order(chr, pos)
Mval <- Mval[ord, ]
chr <- chr[ord]
pos <- pos[ord]
strand <- strand[ord]
stopifnot(
  nrow(Mval) == length(chr),
  nrow(Mval) == length(pos)
)
library(GenomicRanges)
sitesGr <- GRanges(
  seqnames = chr,
  ranges = IRanges(start = pos, end = pos),
  strand = strand
)
colnames(Mval)
samples <- colnames(Mval)
group <- ifelse(grepl("T1$", samples), "control",
                ifelse(grepl("T2$", samples), "treatment", NA))
table(group, useNA = "ifany")
group <- factor(group)
design <- model.matrix(~ group)
cbind(samples, group)
library(DMRcate)
test <- DMLtest.multiFactor(
  bs,
  group1 = "control",
  group2 = "treatment"
)
fit <- DMLfit.multiFactor(
  bs,
  design = data.frame(group = group),
  formula = ~ group
)
test <- DMLtest.multiFactor(fit, term = "group")
#now this part requires modification it currently doesnt find any DMRs
dmrs <- callDMR(
  test,
  p.threshold = 0.05,
  minCG = 2,
  dis.merge = 500
)
#No DMR found! Please use less stringent criteria. 
#But all of this works