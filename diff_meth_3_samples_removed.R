#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)

path_db="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db.rds"

path_diff_meth="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_3_samples_removed.csv"

db=readRDS(path_db)
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
# Differential methylation analysis for T1 vs T2

dmf=calculateDiffMeth(db, num.cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))

all.diff = getMethylDiff(dmf, difference = 25, qvalue = 0.01, type = "all")

write.csv(all.diff, 
          file = path_diff_meth, 
          row.names = FALSE)