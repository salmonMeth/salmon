#libraries
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(ggplot2)
library(ggrepel)
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
#
cols <- c(
  "#2E004B", # Midnight Purple (Darkest)
  "#4B0082", # Indigo (Original)
  "#6A0DAD", # Deep Violet
  "#8A2BE2", # Blue Violet (Original)
  "#4169E1", # Royal Blue (Original)
  "#0000CD", # Medium Blue (Original)
  "#6495ED", # Cornflower Blue (Transition)
  "#9370DB", # Medium Purple (Light-ish)
  "#B0C4DE", # Light Steel Blue (Light)
  "#E6E6FA"  # Lavender (Lightest/End)
)#
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

all.diff = getMethylDiff(dmf, difference = 25, qvalue = 0.01, type = "all")

write.csv(all.diff, 
          file = path_diff_meth, 
          row.names = FALSE)

#visualize
#reload the csv file
diff_data <- read.csv(path_diff_meth, stringsAsFactors = FALSE)
accession_nums <- as.numeric(gsub("NC_0|\\.1", "", diff_data$chr))
diff_data$simple_id <- accession_nums - 59441
#Volcano Plot from the loaded CSV
plot(diff_data$meth.diff, -log10(diff_data$qvalue),
     pch = 20, 
     main = "Volcano Plot: T1 vs T2",
     xlab = "Methylation Difference (%)", 
     ylab = "-log10(Q-value)",
     col = ifelse(abs(diff_data$meth.diff) > 25 & diff_data$qvalue < 0.01, cols[1],cols[3]))

###
#advanced volcano plot
ggplot(diff_data, aes(x = meth.diff, y = -log10(qvalue))) +
  geom_point(alpha = 0.2, color = cols[5], size = 1) +
  geom_point(data = top_hits, color = cols[2], size = 2) +
  geom_hline(yintercept = y_cutoff, linetype = "dashed", color = cols[2], alpha = 0.5) +
  
  # Adjusting the scales
  scale_x_continuous(breaks = seq(-100, 100, by = 20)) + # Fewer vertical lines
  scale_y_continuous(n.breaks = 5) + # Tells R to only pick ~5 nice numbers for the Y-axis
  
  geom_text_repel(data = top_hits, aes(label = simple_id), size = 3, max.overlaps = 50) +
  
  theme_light() + 
  theme(
    # Clean up the grid appearance
    panel.grid.major = element_line(color = cols[9], linewidth = 0.2),
    panel.grid.minor = element_blank(),
  )
