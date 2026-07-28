library(data.table)
base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"
cs_files <- list.files(
  base_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)
all_chr <- list()
for(cs_file in cs_files){
  chr_name <- sub(
    "_segment_cs\\.csv$",
    "",
    basename(cs_file)
  )
  cov_file <- file.path(
    base_dir,
    paste0(chr_name, "_segment_coverage.csv")
  )
  cs  <- fread(cs_file)
  cov <- fread(cov_file)
  stopifnot(
    identical(cs$chr, cov$chr),
    identical(cs$start, cov$start),
    identical(cs$end, cov$end)
  )
  sample_cols <- names(cs)[4:ncol(cs)]
  out <- cs[, .(chr, start, end)]
  out$strand <- "*"
  for(s in sample_cols){
    cov_name <- paste0("coverage_", sub("^X", "", s))
    cs_name  <- paste0("numCs_", sub("^X", "", s))
    ts_name  <- paste0("numTs_", sub("^X", "", s))
    out[[cov_name]] <- cov[[s]]
    out[[cs_name]]  <- cs[[s]]
    out[[ts_name]]  <- cov[[s]] - cs[[s]]
  }
  all_chr[[chr_name]] <- out
}
meth_df <- rbindlist(all_chr)
sample_ids <- sub("^X", "", names(cs)[4:ncol(cs)])
head(sample_ids)
treatment <- ifelse(grepl("_T2$", sample_ids), 1, 0)
####################
coverage.index <- grep("^coverage_", names(meth_df))
numCs.index    <- grep("^numCs_", names(meth_df))
numTs.index    <- grep("^numTs_", names(meth_df))
library(methylKit)

mb <- new(
  "methylBase",
  meth_df,
  sample.ids = sample_ids,
  assembly = "Ssal_v3.1",
  treatment = treatment,
  context = "CpG",
  resolution = "region",
  coverage.index = coverage.index,
  numCs.index = numCs.index,
  numTs.index = numTs.index,
  destranded = FALSE
)

saveRDS(mb, file = "/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_segments.rds")

data <- getData(mb)

cs_cols <- mb@numCs.index
ts_cols <- mb@numTs.index

# extract the actual matrices
cs <- as.matrix(data[, ..cs_cols])
ts <- as.matrix(data[, ..ts_cols])

dim(cs)
dim(ts)
cntlist <- lapply(seq_len(nrow(cs)), function(i){
  
  matrix(
    c(cs[i, ], ts[i, ]),
    ncol = 2,
    dimnames = list(
      mb@sample.ids,
      c("numCs", "numTs")
    )
  )
  
})
phis <- sapply(cntlist, function(x) {
  methylKit:::estimatePhi(
    counts = x,
    modelMat = modelMat,
    treatment = mb@treatment
  )[1]
})


summary(phis)
mean(phis, na.rm=TRUE)
median(phis,na.rm=TRUE)
quantile(phis, c(0.9, 0.95, 0.99), na.rm=TRUE)
mean(phis > 1,na.rm=TRUE)     # any overdispersion
mean(phis > 1.5,na.rm=TRUE)   # moderate
mean(phis > 2,na.rm=TRUE)     # substantial
mean(phis > 5,na.rm=TRUE)    # extreme

####################
#the phi values for the segments seem to be much greater than the ones for the cpgs
#which is odd?
#we can check the distribution of variances for the CpGs vs the segments

#variance of the CpGs
path_db="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_meth_normalized.rds"
mb=readRDS(path_db)
df=getData(mb)
str(df)
n_samples <- sum(grepl("^coverage", names(df)))
n_samples
meth_mat <- sapply(1:n_samples, function(i){
  df[[paste0("numCs", i)]] /
    df[[paste0("coverage", i)]]
})
dim(meth_mat)
cpg_variance <- apply(
  meth_mat,
  1,
  var,
  na.rm = TRUE
)
summary(cpg_variance)
hist(
  cpg_variance,
  breaks = 100,
  main = "Variance of CpG methylation",
  xlab = "Variance"
)

#for the segments
#	

path_db_seg="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_segments.rds"
mb_seg=readRDS(path_db_seg)
df_seg=getData(mb_seg)
str(df_seg)
n_samples_seg <- sum(grepl("^coverage", names(df_seg)))
n_samples_seg
cs_cols <- grep("^numCs_", names(df_seg), value=TRUE)
cov_cols <- grep("^coverage_", names(df_seg), value=TRUE)
meth_mat_seg <- sapply(seq_along(cs_cols), function(i){
  
  df_seg[[cs_cols[i]]] /
    df_seg[[cov_cols[i]]]
  
})
cpg_variance_seg <- apply(
  meth_mat_seg,
  1,
  var,
  na.rm = TRUE
)
summary(cpg_variance_seg)
hist(
  cpg_variance_seg,
  breaks = 100,
  main = "Variance of segment methylation",
  xlab = "Variance"
)

quantile(
  cpg_variance,
  c(0.5,0.75,0.9,0.95,0.99),
  na.rm=TRUE
)

quantile(
  cpg_variance_seg,
  c(0.5,0.75,0.9,0.95,0.99),
  na.rm=TRUE
)
#doesnt seem to be related to the actual variance
#is phi changing with coverage?
cov_cols <- grep("^coverage_", names(df_seg), value = TRUE)
meth_mat_seg <- sapply(seq_along(cs_cols), function(i){
  df_seg[[cs_cols[i]]] / df_seg[[cov_cols[i]]]
})
dim(meth_mat_seg)
seg_var <- apply(
  meth_mat_seg,
  1,
  var,
  na.rm = TRUE
)
length(seg_var)
keep <- !is.na(phis) & !is.na(mean_cov_seg)
plot(
  log10(mean_cov_seg[keep]),
  log10(phis[keep]),
  pch = 16,
  xlab = "log10(mean coverage)",
  ylab = "log10(phi)"
)
cor(
  log10(mean_cov_seg[keep]),
  log10(phis[keep]),
  method="spearman"
)
cor(
  log10(df_seg$end - df_seg$start)[keep],
  log10(phis[keep]),
  method="spearman"
)