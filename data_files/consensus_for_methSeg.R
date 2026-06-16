#we are investigating different ways of creating a consensus profile to use in methseg
#one idea is to find the CpGs that are common to most samples and then to weigh them by their coverage percentage
#or sth like that,we will see

#libs
library(GenomicRanges)
library(dplyr)
library(GenomeInfoDb)
library(methylKit)


#all this from segment_bssq.R
db=readRDS("db.rds")
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

M <- as.matrix(df[, grep("^numCs_", colnames(df))])
Cov <- as.matrix(df[, grep("^coverage_", colnames(df))])
colnames(M) <- sample_ids
colnames(Cov) <- sample_ids

#lets keep the CpGs most samples cover

present <- Cov > 0
#80% coverage?
keep80 <- rowMeans(present, na.rm = TRUE) >= 0.8
M_filt80 <- M[keep80, ]
Cov_filt80 <- Cov[keep80, ]
df_filt80 <- df[keep80, c("chr", "start", "end", "strand")]

#90% coverage?
keep90 <- rowMeans(present, na.rm = TRUE) >= 0.9
M_filt90 <- M[keep90, ]
Cov_filt90 <- Cov[keep90, ]
df_filt90 <- df[keep90, c("chr", "start", "end", "strand")]
