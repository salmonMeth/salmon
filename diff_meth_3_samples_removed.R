#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(dplyr)

path_db="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db.rds"

path_diff_meth="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_3_samples_removed_with_overdispersion_difference_10.csv"

db=readRDS(path_db)

###################################
#REMOVE THE SAMPLES and create a new rds object
samples_to_remove <- c("746_T2", "176_T1", "861_T2")
samples_to_keep <- setdiff(db@sample.ids, samples_to_remove)

db_samples_removed <- reorganize(
  db,
  sample.ids = samples_to_keep,
  treatment = db@treatment[match(samples_to_keep, db@sample.ids)]
)
###################################
###################################
n = db@num.records
# Remove specific samples
samples_to_remove <- c("746_T2", "176_T1", "861_T2")

# Indices of samples to remove
rm.idx <- which(db@sample.ids %in% samples_to_remove)

# Remove the samples
db <- reorganize(
  db,
  sample.ids = db@sample.ids[-rm.idx],
  treatment  = db@treatment[-rm.idx]
)


########################################
########################################
#we want to analyze the data before to see if we need preprocessing
dat=getData(db)
raw_list <- lapply(seq_along(db@sample.ids), function(i) {
  x <- data.frame(
    chr = dat$chr,
    start = dat$start,
    end = dat$end,
    strand = dat$strand,
    coverage = dat[[paste0("coverage", i)]],
    numCs = dat[[paste0("numCs", i)]],
    numTs = dat[[paste0("numTs", i)]]
  )
  new(
    "methylRaw",
    x,
    sample.id = db@sample.ids[i],
    assembly = db@assembly,
    context = db@context,
    resolution = "base"
  )
})

#plot the methylation distribution
# extract methylation categories for one methylRaw object
#takes in a methylraw object

get_meth_categories <- function(mr) {
  
  dat <- as.data.frame(mr)
  
  meth <- dat$numCs / dat$coverage * 100
  
  category <- cut(
    meth,
    breaks = c(-Inf, 20, 40, 60, 80, 100),
    labels = c(
      "0-20%",
      "20-40%",
      "40-60%",
      "60-80%",
      "80-100%"
    ),
    include.lowest = TRUE
  )
  
  category <- as.character(category)
  
  category[is.na(category)] <- "NA coverage"
  
  factor(
    category,
    levels = c(
      "NA coverage",
      "0-20%",
      "20-40%",
      "40-60%",
      "60-80%",
      "80-100%"
    )
  )
}
all_categories <- unlist(
  lapply(raw_list, get_meth_categories)
)

sample_ids <- db@sample.ids

t1_idx <- grepl("T1", sample_ids)
t2_idx <- grepl("T2", sample_ids)

t1_categories <- unlist(
  lapply(raw_list[t1_idx], get_meth_categories)
)

t2_categories <- unlist(
  lapply(raw_list[t2_idx], get_meth_categories)
)
make_pie <- function(x, title) {
  
  df <- as.data.frame(table(x))
  colnames(df) <- c("category", "count")
  
  df <- df[df$count > 0, ]
  
  df$percent <- df$count / sum(df$count) * 100
  
  ###labels
  # select labels: largest 2 + NA coverage
  largest2 <- order(df$percent, decreasing = TRUE)[1:2]
  
  na_idx <- which(df$category == "NA coverage")
  
  label_idx <- unique(c(largest2, na_idx))
  
  df$label <- ""
  df$label[label_idx] <- paste0(
    round(df$percent[label_idx], 1),
    "%"
  )
  
  ggplot(df, aes(x = 1, y = percent, fill = category)) +
    geom_col(
      width = 1,
      color = "black",
      linewidth = 0.5
    )+
    coord_polar(theta = "y") +
    theme_void() +
    labs(title = title, fill = "Methylation") +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 2.1,
      nudge_x = 0.15
    )
}
p_all <- make_pie(all_categories, "All samples")

p_t1 <- make_pie(t1_categories, "T1 samples")

p_t2 <- make_pie(t2_categories, "T2 samples")

p_all_t1_t2 <- p_all + p_t1 + p_t2
ggsave(
  filename = "/scratch/project_2010912/ezel/salmon/data_files/use_in_thesis/methylation_pie_all_T1_T2.png",
  plot = p_all_t1_t2,
  width = 12,
  height = 5,
  dpi = 300
)

#############################
##########################
#coverage pie charts
get_cov_categories <- function(mr) {
  
  dat <- as.data.frame(mr)
  
  cov <- dat$coverage
  
  category <- rep(NA_character_, length(cov))
  
  # NA coverage
  category[is.na(cov)] <- "NA coverage"
  
  # coverage bins for non-NA values
  category[!is.na(cov)] <- as.character(
    cut(
      cov[!is.na(cov)],
      breaks = c(-Inf, 0, 10, 50, 100, 500, Inf),
      labels = c(
        "0",
        "1-10",
        "10-50",
        "50-100",
        "100-500",
        ">500"
      ),
      include.lowest = TRUE
    )
  )
  
  factor(
    category,
    levels = c(
      "NA coverage",
      "0",
      "1-10",
      "10-50",
      "50-100",
      "100-500",
      ">500"
    )
  )
}

all_cov_categories <- unlist(
  lapply(raw_list, get_cov_categories)
)

t1_cov_categories <- unlist(
  lapply(raw_list[t1_idx], get_cov_categories)
)

t2_cov_categories <- unlist(
  lapply(raw_list[t2_idx], get_cov_categories)
)

p_cov_all <- make_pie(all_cov_categories, "Coverage - All samples")
p_cov_t1  <- make_pie(t1_cov_categories, "Coverage - T1 samples")
p_cov_t2  <- make_pie(t2_cov_categories, "Coverage - T2 samples")

p_cov_all=p_cov_all + p_cov_t1 + p_cov_t2

ggsave(
  filename = "/scratch/project_2010912/ezel/salmon/data_files/use_in_thesis/coverage_pie_all_T1_T2.png",
  plot = p_cov_all,
  width = 12,
  height = 5,
  dpi = 300
)
################################
################################
###### we can also plot the methylation&coverage histograms of specific samples
#replace 1 with somethin in range(1:92) to get any of the 92 samples

getMethylationStats(raw_list[[1]], plot=TRUE)
getCoverageStats(raw_list[[1]], plot=TRUE)

##################################
##################################
########### we can apply the builtin coverage filtering function
# but because it uses quantile internally, without the na.rm=TRUE
#argument, we need to first get rid of the na values ourselves

raw_list_obj <- new(
  "methylRawList",
  raw_list,
  treatment = db@treatment
)
raw_list_clean <- lapply(raw_list_obj, function(mr) {
  
  dat <- as.data.frame(mr)
  dat <- dat[!is.na(dat$coverage), ]
  
  new(
    "methylRaw",
    dat,
    sample.id = mr@sample.id,
    assembly = mr@assembly,
    context = mr@context
  )
})

raw_list_clean <- new(
  "methylRawList",
  raw_list_clean,
  treatment = raw_list_obj@treatment
)
filtered_obj=filterByCoverage(raw_list_clean,lo.count=NULL,lo.perc=1, hi.count=NULL,hi.perc=99.9)
normalized_obj=normalizeCoverage(filtered_obj)
#now we have to unite the samples
meth=methylKit::unite(filtered_obj, destrand=FALSE)
meth_normal=methylKit::unite(normalized_obj, destrand=FALSE)
###save
path_meth_normal="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_meth_normalized.rds"
saveRDS(meth_normal, path_meth_normal)
####analyse the sd of methylation
pm=percMethylation(meth_normal) # get percent methylation matrix
mds=matrixStats::rowSds(pm) # calculate standard deviation of CpGs
head(meth_normal[mds>20,])
hist(mds,col="cornflowerblue",xlab="Std. dev. per CpG", breaks=40, xlim=c(0,max(mds)))

####check the dendrogram 
hc <- clusterSamples(meth_normal, dist = "correlation", method = "ward", plot = FALSE)

label_cols <- ifelse(meth_normal@treatment == 0, "steelblue", "firebrick")

plot(hc, labels = FALSE, hang = -1, main = "Sample clustering")

labs <- hc$labels[hc$order]
cols <- label_cols[match(labs, meth_normal@sample.ids)]

text(x = 1:length(labs), y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]), labels = labs, srt = 90, adj = 1, xpd = TRUE, cex = 0.5, col = cols)
###denodgram for the unfilterd&unnormalized data

mb_orig <- new(
  "methylBase",
  dat,
  sample.ids = db@sample.ids,
  assembly = db@assembly,
  context = db@context,
  treatment = db@treatment,
  resolution = "base"
)
#PCA on filtered&normalized
pc=PCASamples(meth_normal,obj.return = TRUE, adj.lim=c(1,1))

# Differential methylation analysis for T1 vs T2

dmf=calculateDiffMeth(db, num.cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")), overdispersion = "MN")

all.diff = getMethylDiff(dmf, difference = 10, qvalue = 0.01, type = "all")

all.diff.df <- getData(all.diff)

write.csv(
  all.diff.df,
  file = path_diff_meth,
  row.names = FALSE
)
###########################################
###########################################
###########################################
####checking for overdispersion
#using the normalized Cpgs
path_db="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_meth_normalized.rds"
mb=readRDS(path_db)
Tcols <- mb@numTs.index
Ccols <- Tcols - 1
# extract count matrix
cntlist <- split(as.matrix(getData(mb)[, c(Ccols, Tcols)]),
                 1:nrow(getData(mb)))
vars <- data.frame(treatment = mb@treatment)
modelMat <- model.matrix(~ treatment, vars)
phis <- sapply(cntlist, function(x) {
  estimatePhi(x,
              modelMat = modelMat,
              treatment = mb@treatment)[1]
})
methylKit:::estimatePhi
phis <- sapply(cntlist, function(x)
  methylKit:::estimatePhi(
    counts = x,
    modelMat = modelMat,
    treatment = mb@treatment
  )[1]
)
summary(phis)
mean(phis)
median(phis)
quantile(phis, c(0.9, 0.95, 0.99))
mean(phis > 1)      # any overdispersion
mean(phis > 1.5)    # moderate
mean(phis > 2)      # substantial
mean(phis > 5)      # extreme


#######for the segments
####checking for overdispersion
path_db="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/segment_matrix_3_samples_removed.rds"
mb=readRDS(path_db)
Tcols <- mb@numTs.index
Ccols <- Tcols - 1
# extract count matrix
cntlist <- split(as.matrix(getData(mb)[, c(Ccols, Tcols)]),
                 1:nrow(getData(mb)))
vars <- data.frame(treatment = mb@treatment)
modelMat <- model.matrix(~ treatment, vars)


phis <- sapply(cntlist, function(x)
  methylKit:::estimatePhi(
    counts = x,
    modelMat = modelMat,
    treatment = mb@treatment
  )[1]
)
summary(phis)
mean(phis)
median(phis)
quantile(phis, c(0.9, 0.95, 0.99))
mean(phis > 1)      # any overdispersion
mean(phis > 1.5)    # moderate
mean(phis > 2)      # substantial
mean(phis > 5)      # extreme