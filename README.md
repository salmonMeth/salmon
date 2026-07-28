So, the latest order of documents:


# R scripts

## Initial analysis
1. diff_meth.R : uses a script "parallelDiffMeth.sh" to use parallelization to get all the DMCs in the data(T1~T2), the output of this is at "/scratch/project_2010912/ezel/diff_meth_results_T1_vs_T2_Full.csv", and there are also some visualization scripts.
2. islands_and_annot.R : takes that diff_meth_results_T1_vs_T2_Full.csv and annotates them and checks if any of the DMCs are in CpG islands.
3. homeolog_annot.R: uses the result from step 2 and adds information about the homeologs, using the homeologs.xlsx file found in the supplement of the homeolog/duplication paper.

## Differentially methylated regions

4. segment_bssq.R: used this to generate the bsseq object "bsseq_object.rds" and the plot for then mean methylation level
5. run_dmrseq.R: uses the bsseq object and parallelizes the segmentation.
6. tile.R : contains the code needed to tile the genome, and saves the resulting matrices containing the tiled CpG information as well as the coordinates of each tile in the "methylation_tiles_all.rds" file.
7. smooth_plot.R : takes the tiles and applies two different kinds of smoothing, loess and averaging and also contains some plots to compare the results.

## Segmentation 

This no longer uses the tiled matrices but goes back to using the methylation&coverage matrices obtained from "db.rds", can be found in segment_bssq.R
---
NOTE: all the results obtained here and downstream analysis is done using the default parameters
except the ones specified as different params in the data files.
---
8. methseg.R : extracts the methylation and coverage matrices and loops over all the chromosomes to find segments using the methSeg function, and saves the results as bed files that can be viewed in IGV.
-The results of these are in the segmentation_output/bed, segmentation_output/diagnostics file.
-Some stats we exctracted from these consensus segments are in summaries
9. methseg_per_sample.R : creates separate segmentations for each sample, and extracts some stats.
-The results of these are stored in the  segmentation_output/per_sample file
-Some stats regarding the per sample segmentation are in the csv files in this directory
10. annot_segments.R : we annotated the segments we got using the consensus profile. Also contains the code used to get the stats and create the plots for the group_level_stats_by_K (K being the number of segment groups)
-The results of this are in data_files/segment_annotations.
- This uses the median consensus profile
---
11. methSeg_T1_vs_T2.R : creates the segments, diagnostic plots, and summary stats for T1 and T2 samples separately, using the consensus profile based on the coverage weighted mean of the subsetted samples
12. segment_overlap.R : compares the segments obtained from T1, T2 and cumulative consensus profiles
13. segment_disagree.R : finds the chunks where there is no 3way agreement among all the sections, merges the overlapping/adjacent ones and then annotates them
14. plot_segment_overlap.R : plots the overlapping parts of the T1, T2 and cumulative segments for all chromosomes.

## Methylation of the segments

15. segment_methyl.R : calculates the methylation levels using the segments obtained from the cumulative coverage weighted consensus

## Testing the segment approach (against individual CpGs) using a model

16. construct_cs_cov_matrices.R : constructs the methylation&coverage count matrices we need to run the model.

17. model_segment_counts.R : has the code needed to model the methylated&unmethylated Cpg counts using weight, ma&pa ID has covariates.
Has code both for the pooled model and one where we fit two different models to T1 and T2 samples. Also fits two different models, one with weight scaled and one not.

## Comparing the DSS DMR/DML with the Methylkit segments

18. dss_dmr_dml.R : has the code to execute the DSS DMR/DML as well, compares the differentially methylated CpGs found by DSS to the Methylkit ones, their locations, the difference detected by each method, and compares the DSS segments to the Methylkit ones

## Annotating the differentially methylated bases/segments:

19. annot_diff_meth_DLM_DMR.R : annotates the differentially methylated bases/segments identified by DSS and Methylkit, 
also has some code to compare the distribution of the genomic features identified and has some plots.

20. gene_ontology_diff_meth.R : does the gene ontology for the annotated bases/segments
## Assesing the model performance

19. model_problem_segments.R : checks out the segments where the model didn't converge, trying to identify the coverage, coefs etc. Also has code
to compare the performance of different models (with weight scaled or unscaled) regarding how many segments throw errors.

----------------------------------------------------------------------
----------------------------------------------------------------------

# Data files/plots

## Segmentation using consensus profiles:

*** All this except the ones specified otherwise are using the default parameters.

### Consensus using the median:

a) data_files/segmentation_output/beds : contains the segmented chromosomes obtained using the median methylation value over all the samples as consensus.

b) data_files/segmentation_output/diagnostics : contains the diagnostic plots for those.

c) data_files/segmentation_output/summaries: contains the summary stats of the segments, per chromosome.

### Consensus using coverage weighted mean:

a) data_files/segmentation_output_cov_weighted_consensus/beds: contains the segmented chromosomes obtained using the coverage weighted mean methylation value over all the samples as consensus.

b) data_files/segmentation_output_cov_weighted_consensus/diagnostics: contains the diagnostic plots for those.

c) data_files/segmentation_output_cov_weighted_consensus/summaries: contains the summary stats for these, the details of the contents are listed in "data_files/segmentation_output/summaries/colum_names.Rmd"

***

## Segmentation of each sample individually

a) data_files/segmentation_output/per_sample: contains one folder for each sample, within those folders there are beds and the diagnostic plots

b) data_files/segmentation_output/per_sample_different_params: contains folders named after the different parameters tried for methSeg
NOTE: this wasnt done on all the samples, just some were chosen randomly to see the differences.

c) data_files/segmentation_output/sample_summaries: contains the summaries for each sample, default parameters

## Annotations for the segments (based on the median consensus profile)

a) data_files/segment_annotations: contains one file per chromosome.

Each file annotates the segments found based on the median consensus profile, adding new columns to the beds present in "data_files/segmentation_output/beds"
Mainly looks at the number of genes, promoters found in each segment as well as their lengths and the number of each per base pairs(can be greater than 1 due to overlap)

## Bed files of the segments suitable for UCSC browser

a) /data_files/segmentation_output/beds_fixed: contains the beds for the segmentation creates separately using the T1 coverage weighted consensus, T2 coverage weighted consensus, and the one
created from all the samples, the "fixed" refers to the chromosome number which is in these files in NCIB format. The purpose of these is to be used in the UCSC browser.

## Results from the comparison of the segments seperately obtained from T1, T2 and the cumulative consensus

a) data_files/segment_overlap_results : contains information on all the chromosomes regarding the 3way and 2way intersections

b) data_files/segment_overlap_results/disagreement_regions: contains the disagreeing segments and their annotations

##3way segmentation matches

a)data_files/segmentation_3way_overlap: the csv files contain information about the mathcing segments found for each chromosome,
the cumulative consensus is used as the reference.

## Segment methylation levels

a)data_files/segment_methylation : has one csv file per chromosome, rows are the segments
of the chromosome and there is one column per sample.

b)data_files/objects_intermediary/segment_matrix.rds : has the methylation matrix for the segments, one row per segment, one col per sample


## To run the model count model

a) data_files/segment_cs_cov_3_samples_removed: has the counts needed to run the weight model, two files per chromosome,
each row is one segment, columns have the 92 samples

b) data_files/model_results/model_results_T1_T2_separated : has the UNSCALED weight model outputs, one file per chromosome.

c)data_files/model_results/weight_rescaled_T1_T2_separated : has the model outputs for the scaled weight model.
Also contains the subfolder:

c.2) data_files/model_results/weight_rescaled_T1_T2_separated/problem_segments: which has some stats and example model outputs of the segments where the model
didn't converge.

## Differential methylation using methseg

a)  data_files/objects_intermediary/db_meth_normalized.rds: is the filtered and normalized (using the builtin functions from MethylKit) methylBase object

## Diff. methylated base/segment annotation/gene ontology

### Annotations

a)  diff_meth_results_base_segments/annotated:

a.1)  	dss_diff_meth_bases.csv : has the annotation for all the bases identified by DSS

a.2)    methylkit_diff_meth_bases.csv : has the annotation for all the bases identified by methylkit

a.3)    methylkit_diff_meth_segs.csv : has the annotation for all the segments identified by methylKit.

### Genes comparison:

b.1)  	genes/dss_diff_meth_bases_genes.csv 

b.2) 	genes/methylkit_diff_meth_bases_genes.csv

b.3 )	methylkit_diff_meth_segs_genes.csv : these have the gene names extracted from the annotation files, one name per line
### Gene ontology

c.1) genes/gene_ontology/zebrafish_orthologs_all_methods.csv : has the ontology for the genes found by all 3 methods, based on the zebrafish onthology
## Things that are in data_files/objects_intermediary


### Base based (as opposed to segments) objects 

Raw objects

a) db.rds: contains the info for all 95 samples

b) db_3_samples_removed.rds : the same info but the bad samples removed

For the methylated base and total coverage *counts*

a) cpg_coverage_total_matrix_3_samples_removed.rds

b) cpg_meth_total_matrix_3_samples_removed.rds

Percent methylation

a) 	cpg_meth_percent_matrix_3_samples_removed.rds

Filtered&normalized using Methylkit's builtin functions

a) 	db_meth_normalized.rds : has the coverage, methylated and unmethylated base counts for the bases that passed the normalization&filtering

### Using DSS

The results of differentially methylated base detection, w and w/o smoothing

a) DMLtest_no_smoothing.rds

b) DMLtest_smoothing.rds

Construcing differentially methylated regions using the results of the DMLtest

a) DSS_DMRs_no_smoothing.rds

b) DSS_DMRs_smoothing.rds

### Segments 

a) 	db_segments.rds has the coverage, methylated and unmethylated base counts accumulated over the segments

# Visualization

a) data_files/segmentation_output/T1_vs_T2_plots: contains various plots to compare the segmentations of individual samples, color coded for T1 vs T2 samples

a.1) data_files/segmentation_output/T1_vs_T2_plots/segment_group_num_per_chr: compares the total number of segment groups found for each chromosome for each sample.

Each plot contains the plots of multiple chromosomes.

a.2) data_files/segmentation_output/T1_vs_T2_plots/total_segment_num_per_chr: compares the total number of segments per chromosome for each sample.

Each plot contains the plots of multiple chromosomes.

a.3) data_files/segmentation_output/T1_vs_T2_plots/mean_methylation_per_chromosome: compares the mean methylation value of each segment group for all possible segment group numbers found in the data per chromosome.

There is one plot for every chromosome here.

## Comparison based on the number of segments, across all chromosomes

The aim here is try to see if there are any similarities between chromosomes that are found to have the same number of totssl segment groups
Looks at gene density (using the annotations), mean methylation level and average segment length.
Uses the segments obtained using the median consensus profile.

a) data_files/group_level_stats_by_K has one subdirectory for every segment group number present in the segment data.
a.1) data_files/group_level_stats_by_K/K_2
a.2) data_files/group_level_stats_by_K/K_3
a.3) data_files/group_level_stats_by_K/K_4
a.4) data_files/group_level_stats_by_K/K_8

Each subdirectory contains 4 files:

a) gene_density_boxplot: to compare the mean gene density per segment for segments belonging to a specific segment group
b) methylation_boxplot: to compare the mean methylation level
c) segment_length_boxplot: to compare the mean segment length
d) group_stats.csv : contains some stats regarding the segment groups. The empty cells are due to there being only one chromosome with the specified number of segment groups, hence standard deviation etc can't be calculated.


## 3way segment overlap to compare how the segments found using just T1 or T2 differ from the ones from all the samples

a) data_files/segmentation_3way_overlap/__plots : folders contains the visualizations to compare how the matching segments align