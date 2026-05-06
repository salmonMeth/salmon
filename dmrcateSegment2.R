library(methylKit)
library(dmrseq)
library(bsseq)
library(GenomeInfoDb)
library(DMRcate)

source("~/rSalmon/segmentBSSQ.R")
#old stuff we need to do the segmentation, obtained from segmentBSSQ

df <- getData(db)
sample_ids <- db@sample.ids
n_samples <- length(sample_ids)
#we change the colnames of the df to include the sample.id as well
colnames(df) <- c(
  "chr", "start", "end", "strand",
  as.vector(sapply(sample_ids, function(id) {
    c(paste0("coverage_", id),
      paste0("numCs_", id),
      paste0("numTs_", id))
  }))
)
##
chr  <- df$chr
pos  <- df$start
strand <- df$strand
ord <- order(chr, pos)
chr <- chr[ord]
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
Mval <- Mval[ord, ]
treat <- db@treatment
names(treat) <- db@sample.ids
bs <- BSseq(
  chr = chr,
  pos = pos,
  M = M,
  Cov = Cov,
  sampleNames = colnames(M)
)
pData(bs)$group <- factor(treat)
########################
#new stuff begins
#this Mval is the filtered and logodds transformed methylation level matrix
#we get from the dataframe
samples <- colnames(Mval)
group <- ifelse(grepl("T1$", samples), "control",
                 ifelse(grepl("T2$", samples), "treatment", NA))
group <- factor(group, levels = c("control", "treatment"))
design <- model.matrix(~ group)
rownames(design) <- samples
methdesign <- edgeR::modelMatrixMeth(design)
colnames(methdesign) <- make.names(colnames(methdesign))
##########################
cont.mat <- limma::makeContrasts(
  T2_vs_T1 = grouptreatment,
  levels = methdesign
)

##########################

seq_annot <- sequencing.annotate(bs, methdesign, all.cov = TRUE,
                                  contrasts = TRUE,
                                  cont.matrix = cont.mat,
                                  coef = "T2_vs_T1",
                                  fdr = 0.05)