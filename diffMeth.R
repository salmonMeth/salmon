#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)

#paths
path= "/scratch/project_2010912/ezel/methylBase_7x_40ind.txt.bgz"
#to save our diff methylation results
#path_diff_meth="/scratch/project_2010912/ezel/diff_meth_results_T1_vs_T2.csv"
path_diff_meth="/scratch/project_2010912/ezel/diff_meth_results_T1_vs_T2_Full.csv"

#to get the annotations and save them
trg_dir = "/scratch/project_2010912/ezel/"
#for raw data with chromosome numbers
path_annot = paste0(trg_dir, "salmon_assembly_report.txt")
#for the gene info
path_annot_genes ="/scratch/project_2010912/ezel/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz"
#
db = readMethylDB(path)
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
all.diff=getMethylDiff(dmf,difference=25,qvalue=0.01,type="all")
write.csv(getData(all.diff),
          file =path_diff_meth,
          row.names = FALSE)
