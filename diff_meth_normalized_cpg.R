#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)

path_meth_normal="/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_meth_normalized.rds"
db = readRDS(path_meth_normal)

n=db@num.records
#this is just a random sampling we do to work on a smaller subset of the data first
#inds = sample(1:n,1e5)
#db= db[inds]
all_ids =db@sample.ids
#differential mehtylation analysis for T1 vs T2
#DIFF METHYLATION ANALYSIS

#find the datapoints from T1
is_t1 = grepl("T1", db@sample.ids)
#if the datapoint is from t1, =0, else, =1
new_treatment = ifelse(is_t1, 0, 1)
data.frame(ID = db@sample.ids, Treatment = new_treatment)

dmf=calculateDiffMeth(db, num.cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
#filter the results where the difference is < than 25 and qvalue >.01

all.diff = getMethylDiff(dmf, difference = 25, qvalue = 0.01, type = "all")

write.csv(getData(all.diff), 
          file = path_diff_meth, 
          row.names = FALSE)
