#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(dplyr)
db=readRDS("/scratch/project_2010912/ezel/salmon/data_files/objects_intermediary/db_3_samples_removed.rds")
n=db@num.records
all_ids =db@sample.ids
head(all_ids)
length(all_ids)
is_t1 = grepl("T1", db@sample.ids)
#if the datapoint is from t1, =0, else, =1
new_treatment = ifelse(is_t1, 0, 1)
data.frame(ID = db@sample.ids, Treatment = new_treatment)
#dmf=calculateDiffMeth(db, adjust="bonferroni",num.cores = as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")))
#filter the results where the difference is < than 25 and qvalue >.01

#all.diff = getMethylDiff(dmf, difference = 10, qvalue = 0.01, type = "all")

#all.diff <- getData(all.diff)

#write.csv(all.diff, file = "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_mk_bonf_correction.csv", row.names = FALSE)
db_mem <- unite(db)

dmf_DSS=calculateDiffMethDSS(db_mem, adjust="bonferroni")
#filter the results where the difference is < than 25 and qvalue >.01

all.diff_DSS = getMethylDiff(dmf_DSS, difference = 10, qvalue = 0.01, type = "all")
all.diff_DSS <- getData(all.diff_DSS)

write.csv(all.diff_DSS, file = "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/diff_meth_DSS_via_mk_bonf_correction.csv", row.names = FALSE)
