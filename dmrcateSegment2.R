library(methylKit)
library(dmrseq)
library(bsseq)
library(GenomeInfoDb)
library(DMRcate)

source("~/rSalmon/segmentBSSQ.R")
samples <- colnames(Mval)
group2 <- ifelse(grepl("T1$", samples), "control",
                 ifelse(grepl("T2$", samples), "treatment", NA))
group2 <- factor(group2, levels = c("control", "treatment"))
design2 <- model.matrix(~ group2)
rownames(design2) <- samples
methdesign2 <- edgeR::modelMatrixMeth(design2)
colnames(methdesign2) <- make.names(colnames(methdesign2))
##########################
cont.mat2 <- limma::makeContrasts(
  T2_vs_T1 = grouptreatment,
  levels = methdesign2
)

##########################
cont.mat2 <- limma::makeContrasts(
  T2_vs_T1 = group2treatment,
  levels = methdesign2
)
seq_annot2 <- sequencing.annotate(bs, methdesign2, all.cov = TRUE,
                                  contrasts = TRUE,
                                  cont.matrix = cont.mat2,
                                  coef = "T2_vs_T1",
                                  fdr = 0.05)