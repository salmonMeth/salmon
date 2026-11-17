#model pedigrees
install.packages("nadiv")
library(data.table)
library(nadiv)

remove_samples <- c("746_T2", "176_T1", "861_T2")

meta <- meta[!sample %in% remove_samples]

# Check
meta[, .N, by = TimePoint]
meta[sample %in% remove_samples]

# Updated counts
meta[, .(
  n_samples = .N,
  n_animals = uniqueN(Running.ID),
  n_dams = uniqueN(ID_ma),
  n_sires = uniqueN(ID_pa)
)]

rbindlist(list(
  meta_T1[, .(
    TimePoint = "T1",
    n_samples = .N,
    n_animals = uniqueN(Running.ID)
  )],
  meta_T2[, .(
    TimePoint = "T2",
    n_samples = .N,
    n_animals = uniqueN(Running.ID)
  )]
))
ped0 <- unique(meta[, .(
  ID   = as.character(Running.ID),
  Dam  = as.character(ID_ma),
  Sire = as.character(ID_pa)
)])

# Check that every animal has only one parent combination
parent_check <- ped0[, .(
  n_combinations = .N
), by = ID]

parent_check[n_combinations > 1]

# Check missing parents
ped0[, .(
  missing_dam = sum(is.na(Dam) | Dam == ""),
  missing_sire = sum(is.na(Sire) | Sire == "")
)]

# Prepare pedigree
ped <- nadiv::prepPed(as.data.frame(ped0))

# Construct relationship matrix
A <- nadiv::makeA(ped)

A_dense <- as.matrix(A)

dim(A_dense)
count_siblings <- function(x) {
  
  offspring <- unique(x[, .(
    ID = as.character(Running.ID),
    Dam = as.character(ID_ma),
    Sire = as.character(ID_pa)
  )])
  
  pairs <- CJ(
    ID1 = offspring$ID,
    ID2 = offspring$ID,
    unique = TRUE
  )[ID1 < ID2]
  
  pairs <- merge(
    pairs, offspring,
    by.x = "ID1", by.y = "ID"
  )
  setnames(pairs, c("Dam", "Sire"), c("Dam1", "Sire1"))
  
  pairs <- merge(
    pairs, offspring,
    by.x = "ID2", by.y = "ID"
  )
  setnames(pairs, c("Dam", "Sire"), c("Dam2", "Sire2"))
  
  pairs[, relationship := fcase(
    Dam1 == Dam2 & Sire1 == Sire2, "Full siblings",
    Dam1 == Dam2 | Sire1 == Sire2, "Half siblings",
    default = "Unrelated"
  )]
  
  pairs[, .N, by = relationship]
}

count_siblings(meta_T1)
count_siblings(meta_T2)

meta_T1[, animal := factor(
  as.character(Running.ID),
  levels = rownames(A_dense)
)]

meta_T2[, animal := factor(
  as.character(Running.ID),
  levels = rownames(A_dense)
)]

# Verify that all sampled animals are represented
stopifnot(!anyNA(meta_T1$animal))
stopifnot(!anyNA(meta_T2$animal))

keep=read.csv("/scratch/project_2010912/ezel/salmon/data_files/high_methyl_variability_segments/segments_variable_T1_or_T2.csv")

library(data.table)

data_dir <- "/scratch/project_2010912/ezel/salmon/data_files/animal_model"

cs_files <- list.files(
  data_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)

cov_files <- list.files(
  data_dir,
  pattern = "_segment_coverage\\.csv$",
  full.names = TRUE
)

cs_list <- lapply(cs_files, fread)
cov_list <- lapply(cov_files, fread)

cs_all <- rbindlist(cs_list, fill = TRUE)
cov_all <- rbindlist(cov_list, fill = TRUE)
key_cols <- c("chr", "start", "end")

setkeyv(cs_all, key_cols)
setkeyv(cov_all, key_cols)

stopifnot(
  nrow(cs_all) == nrow(cov_all),
  identical(cs_all[, ..key_cols], cov_all[, ..key_cols])
)
library(data.table)

keep_dt <- as.data.table(keep)

# Keep only the segment identifiers and coordinates
keep_dt <- unique(keep_dt[, .(segment_id, chr, start, end)])

# Read the CS and coverage files
cs_files <- list.files(
  data_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)

cov_files <- list.files(
  data_dir,
  pattern = "_segment_coverage\\.csv$",
  full.names = TRUE
)

cs_all <- rbindlist(lapply(cs_files, fread), fill = TRUE)
cov_all <- rbindlist(lapply(cov_files, fread), fill = TRUE)

# Match the segments in keep to the count tables
key_cols <- c("chr", "start", "end")
library(data.table)

keep_dt <- as.data.table(keep)[, .(segment_id, chr, start, end)]

# Match the selected segments
cs_keep <- merge(
  keep_dt, cs_all,
  by = c("chr", "start", "end"),
  all = FALSE
)

cov_keep <- merge(
  keep_dt, cov_all,
  by = c("chr", "start", "end"),
  all = FALSE
)

# Check matching
cat("Segments in keep:", uniqueN(keep_dt$segment_id), "\n")
cat("Matched in CS:", uniqueN(cs_keep$segment_id), "\n")
cat("Matched in coverage:", uniqueN(cov_keep$segment_id), "\n")

# Filter to the retained T1 and T2 samples
samples_T1 <- meta_T1$sample
samples_T2 <- meta_T2$sample

cs_T1 <- cs_keep[, c("segment_id", samples_T1), with = FALSE]
cov_T1 <- cov_keep[, c("segment_id", samples_T1), with = FALSE]

cs_T2 <- cs_keep[, c("segment_id", samples_T2), with = FALSE]
cov_T2 <- cov_keep[, c("segment_id", samples_T2), with = FALSE]

library(data.table)

meta_T1 <- as.data.table(meta_T1)
meta_T2 <- as.data.table(meta_T2)

meta_T1[, animal := as.character(Running.ID)]
meta_T2[, animal := as.character(Running.ID)]

meta_T1[, sample := as.character(sample)]
meta_T2[, sample := as.character(sample)]

make_model_data <- function(cs, cov, metadata) {
  
  cs <- copy(as.data.table(cs))
  cov <- copy(as.data.table(cov))
  metadata <- copy(as.data.table(metadata))
  
  # Convert count tables to long format
  cs_long <- melt(
    cs,
    id.vars = "segment_id",
    variable.name = "sample",
    value.name = "meth"
  )
  
  cov_long <- melt(
    cov,
    id.vars = "segment_id",
    variable.name = "sample",
    value.name = "coverage"
  )
  
  # Combine methylated counts and coverage
  dat <- merge(
    cs_long,
    cov_long,
    by = c("segment_id", "sample")
  )
  
  # Add metadata
  dat <- merge(
    dat,
    metadata[, .(
      sample,
      animal,
      Weight_sc,
      ID_ma,
      ID_pa
    )],
    by = "sample"
  )
  
  # Calculate unmethylated counts
  dat[, unmeth := coverage - meth]
  
  # Remove observations with missing or invalid counts
  dat <- dat[
    !is.na(meth) &
      !is.na(coverage) &
      coverage > 0 &
      meth >= 0 &
      meth <= coverage
  ]
  
  dat[, `:=`(
    meth = as.integer(meth),
    unmeth = as.integer(unmeth),
    coverage = as.integer(coverage)
  )]
  
  return(dat[])
}
dat_T1 <- make_model_data(
  cs = cs_T1,
  cov = cov_T1,
  metadata = meta_T1
)

dat_T2 <- make_model_data(
  cs = cs_T2,
  cov = cov_T2,
  metadata = meta_T2
)
cat("T1 observations:", nrow(dat_T1), "\n")
cat("T1 segments:", uniqueN(dat_T1$segment_id), "\n")
cat("T1 animals:", uniqueN(dat_T1$animal), "\n")

cat("T2 observations:", nrow(dat_T2), "\n")
cat("T2 segments:", uniqueN(dat_T2$segment_id), "\n")
cat("T2 animals:", uniqueN(dat_T2$animal), "\n")

# Inspect
head(dat_T1)
head(dat_T2)

# Check coverage and methylated counts
summary(dat_T1$coverage)
summary(dat_T1$meth)
summary(dat_T1$unmeth)

summary(dat_T2$coverage)

saveRDS(
  dat_T1,
  file.path(data_dir, "keep_segments_model_data_T1.rds")
)

saveRDS(
  dat_T2,
  file.path(data_dir, "keep_segments_model_data_T2.rds")
)

# Animal IDs in each dataset
ids_T1 <- sort(unique(as.character(dat_T1$animal)))
ids_T2 <- sort(unique(as.character(dat_T2$animal)))

# Check that all animals are present in the pedigree matrix
stopifnot(all(ids_T1 %in% rownames(A_dense)))
stopifnot(all(ids_T2 %in% rownames(A_dense)))

# Subset A
A_T1 <- A_dense[ids_T1, ids_T1, drop = FALSE]
A_T2 <- A_dense[ids_T2, ids_T2, drop = FALSE]
# Dimensions
dim(A_T1)
dim(A_T2)

# Confirm row and column ordering
stopifnot(identical(rownames(A_T1), ids_T1))
stopifnot(identical(colnames(A_T1), ids_T1))

stopifnot(identical(rownames(A_T2), ids_T2))
stopifnot(identical(colnames(A_T2), ids_T2))

# Check positive definiteness
eigen(A_T1, symmetric = TRUE, only.values = TRUE)$values |> summary()
eigen(A_T2, symmetric = TRUE, only.values = TRUE)$values |> summary()

saveRDS(A_T1, file.path(data_dir, "A_T1.rds"))
saveRDS(A_T2, file.path(data_dir, "A_T2.rds"))
library(brms)

fit_animal_model <- function(dat, A, seg_id) {
  
  dat_seg <- copy(dat[segment_id == seg_id])
  
  # Identify animals
  ids <- sort(unique(as.character(dat_seg$animal)))
  
  # Subset relationship matrix
  A_seg <- A[ids, ids, drop = FALSE]
  
  # Match factor levels
  dat_seg[, animal := factor(as.character(animal), levels = ids)]
  
  stopifnot(identical(levels(dat_seg$animal), rownames(A_seg)))
  stopifnot(identical(rownames(A_seg), colnames(A_seg)))
  
  # Model formula
  formula <- bf(
    meth | trials(coverage) ~
      1 + Weight_sc +
      (1 | gr(animal, cov = A))
  )
  
  model <- brm(
    formula = formula,
    data = dat_seg,
    data2 = list(A = A_seg),
    family = brms::beta_binomial(link = "logit"),
    chains = 4,
    cores = 4,
    iter = 2000,
    warmup = 1000,
    seed = 123,
    control = list(adapt_delta = 0.95),
    refresh = 100
  )
  
  return(model)
}