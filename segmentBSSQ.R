#we are trying to do identify differentially methylated regions
#using dmrseq
# permutation-based approach 
library(methylKit)
library(dmrseq)
library(bsseq)
#
source("~/rSalmon/diffMeth.R")
####
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


#extract the converage and methylatedC number matrices
M <- as.matrix(df[, grep("^numCs_", colnames(df))])
Cov <- as.matrix(df[, grep("^coverage_", colnames(df))])
colnames(M) <- sample_ids
colnames(Cov) <- sample_ids

#
#PROBLEM: we have NAs what do we do we the NAs
#we can get rid of them but there are plenty
#we do this instead
# at least 5 reads and at least 2 samples
#keep <- rowSums(Cov >= 5, na.rm = TRUE) >= 2
#we can adjust this 10
#keep3 <- keep3 & rowSums(is.na(Cov)) <= 10
#or we can allow max 20% NAs
#keep_na <- rowSums(is.na(Cov)) <= (0.2 * ncol(Cov))

#keep <- keep & keep_na
#we can also add variance filtering to get rid of uninformative rows
#var_keep <- apply(M_filt, 1, var, na.rm = TRUE) > 0
#we should try these and pick one 
#BUT THESE ALL GET STUCK BECAUSE WE HAVE NAS LEFT SPO WE GOTTA BE SUPER STRICT
chr <- df[, 1]
pos <- df[, 2]
keep <- rowSums(is.na(M)) == 0
M <- M[keep, ]
Cov <- Cov[keep, ]
chr <- chr[keep]
pos <- pos[keep]

#We can choose to smooth, we should decide?
#we need to do some manipulations for this, maybe if needed
#bs <- BSmooth(bs)
treat <- db@treatment
names(treat) <- db@sample.ids

ord <- order(chr, pos)
M <- M[ord, ]
Cov <- Cov[ord, ]
chr <- chr[ord]
pos <- pos[ord]
bs <- BSseq(
  chr = chr,
  pos = pos,
  M = M,
  Cov = Cov,
  sampleNames = colnames(M)
)
pData(bs)$group <- factor(treat)
saveRDS(bs, file = "bsseq_object.rds")
res <- dmrseq(bs, testCovariate = "group")

##################

#we will want to parallelize the chromosomes 
#to do that we will want to have the list of all chromosome names
writeLines(sort(unique(chr)), "chrom_list.txt")
#we will use the path of this chrom_list in the shell script so we need to modify that accordingly