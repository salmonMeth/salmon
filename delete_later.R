library(bsseq)
library(GenomicRanges)
library(DSS)
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(dplyr)
#path_db_normalized="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_meth_normalized.rds"
path_db_normalized="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds"
db_norm=readRDS(path_db_normalized)
df <- getData(db_norm)
Cov <- as.matrix(df[, db_norm@coverage.index])
M <- as.matrix(df[, db_norm@numCs.index])
treatment <- db_norm@treatment

keep <- rowSums(is.na(M)) == 0 &
  rowSums(is.na(Cov)) == 0
M <- M[keep, , drop = FALSE]
Cov <- Cov[keep, , drop = FALSE]
rownames(M) <- NULL
rownames(Cov) <- NULL
colnames(M) <- db_norm@sample.ids
colnames(Cov) <- db_norm@sample.ids
str(M)
str(Cov)
####if we use the non-normalized data. just with the samples removed, we need to remove
#the NA values and then we need to adjust this gRanges object accordingly
df_keep <- df[keep, ]

gr <- GRanges(
  seqnames = df_keep$chr,
  ranges = IRanges(
    start = df_keep$start,
    end = df_keep$end
  )
)
length(gr)
dim(M)
dim(Cov)
###
#gr <- GRanges(seqnames = df$chr,ranges = IRanges(df$start, df$end))
BSobj <- BSseq(gr = gr, M = M,Cov = Cov,sampleNames = db_norm@sample.ids)
cat("Starting DMLtest 1\n")


cat("Starting DMLtest 2\n")
flush.console()
dmlTest_smoothing <- DMLtest(
  BSobj, group1 = which(treatment == 1),
  group2 = which(treatment == 0),
  smoothing=TRUE)

# Save the DML test object
#saveRDS(dmlTest_no_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/DMLtest_no_smoothing_NOTNORMALIZED.rds")
saveRDS(dmlTest_smoothing, "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/DMLtest_smoothing_checkSign.rds")


