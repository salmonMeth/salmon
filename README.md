So, the latest order of documents:

1. diffMeth.R : uses a script "parallelDiffMeth.sh" to use parallelization to get all the DMCs in the data(T1~T2), the output of this is at "/scratch/project_2010912/ezel/diff_meth_results_T1_vs_T2_Full.csv", and there are also some visualization scripts.
2. Islands_And_Annot.R : takes that diff_meth_results_T1_vs_T2_Full.csv and annotates them and checks if any of the DMCs are in CpG islands.
3. Homeolog_Annot.R: uses the result from step 2 and adds information about the homeologs, using the homeologs.xlsx file found in the supplement of the homeolog/duplication paper.
