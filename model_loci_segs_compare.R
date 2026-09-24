#comparing the DHARMA and sig(the phi parameter extracted from the model, which is the overdispersion parameter)
#for the beta binomial model

#the sigma parameter, as well as covergence, positive definite hessian , AIC BIC
#stats all come from the rds files

loci_model=readRDS("/scratch/project_2010912/ezel/salmon/data_files/model_loci_results/glmm_chr1_results.rds")

summary(loci_model$sigma_T1)
#check how many models didnt converge
sum(loci_model$convergence_T1==1,na.rm=TRUE)
#the dharma stats come from 
path_dharma_res="/scratch/project_2010912/ezel/salmon/data_files/model_loci_results/DHARMa/DHARMa_chr1_results.tsv"
dharma_res=read.delim(path_dharma_res)
#check how many have significant overdispersion error p values
sum(dharma_res$DHARMa_dispersion_pvalue_T1<.05,na.rm=TRUE)

sum(dharma_res$DHARMa_outlier_pvalue_T1<.05,na.rm = TRUE)

non_uniform_loci=(dharma_res[dharma_res$DHARMa_uniformity_pvalue_T1<.05,])$start
over_under_dispersed_loci=(dharma_res[dharma_res$DHARMa_dispersion_pvalue_T1<.05,])$start

length(intersect(non_uniform_loci,over_under_dispersed_loci))
####
#investigate the DHARMA results for the segments
dharma_segs=read.csv("/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated_beta_binom/DHARMa_results_summary/NC_059442.1_DHARMa.csv",header=TRUE)
sum(dharma_segs$DHARMa_dispersion_pvalue_T1<.05,na.rm=TRUE)

sum(dharma_segs$DHARMa_outlier_pvalue_T1<.05,na.rm = TRUE)

non_uniform_segs_start=(dharma_segs[dharma_segs$DHARMa_uniformity_pvalue_T1<.05,])$start
over_under_dispersed_segs_start=(dharma_segs[dharma_segs$DHARMa_dispersion_pvalue_T1<.05,])$start
non_uniform_segs_end=(dharma_segs[dharma_segs$DHARMa_uniformity_pvalue_T1<.05,])$end
over_under_dispersed_segs_end=(dharma_segs[dharma_segs$DHARMa_dispersion_pvalue_T1<.05,])$end
length(intersect(non_uniform_segs_start,over_under_dispersed_segs_start))
over_under_d_segs=c(over_under_dispersed_segs_start,over_under_dispersed_segs_end)

#####################################
####################################
#dharma for the segments
library(DHARMa)

base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated_beta_binom"

# Separate output directory
out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated_beta_binom/DHARMa_results_summary"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_dharma <- function(model) {
  
  if (is.null(model)) {
    return(list(
      uniformity_statistic = NA,
      uniformity_pvalue = NA,
      dispersion_statistic = NA,
      dispersion_pvalue = NA,
      outlier_statistic = NA,
      outlier_pvalue = NA,
      error = "Model is NULL"
    ))
  }
  
  tryCatch({
    
    sim <- simulateResiduals(model)
    
    u <- testUniformity(sim)
    d <- testDispersion(sim)
    o <- testOutliers(sim)
    
    list(
      uniformity_statistic = unname(u$statistic),
      uniformity_pvalue = u$p.value,
      
      dispersion_statistic = unname(d$statistic),
      dispersion_pvalue = d$p.value,
      
      outlier_statistic = unname(o$statistic),
      outlier_pvalue = o$p.value,
      
      error = NA_character_
    )
    
  }, error = function(e) {
    
    list(
      uniformity_statistic = NA,
      uniformity_pvalue = NA,
      dispersion_statistic = NA,
      dispersion_pvalue = NA,
      outlier_statistic = NA,
      outlier_pvalue = NA,
      error = conditionMessage(e)
    )
  })
}


# All chromosome model files
model_files <- list.files(
  base_dir,
  pattern = "^model_.*\\.rds$",
  full.names = TRUE
)


# Loop over chromosomes
for (model_file in model_files) {
  
  chr <- sub(
    "^model_(.*)\\.rds$",
    "\\1",
    basename(model_file)
  )
  
  cat("Processing", chr, "\n")
  
  # Load chromosome object
  check <- readRDS(model_file)
  
  chromosome_results <- list()
  
  
  # Loop over segments
  for (segment_name in names(check)) {
    
    segment <- check[[segment_name]]
    
    cat("  ", segment_name, "\n")
    
    # Coordinates
    chr_segment <- segment$location$chr
    start <- segment$location$start
    end <- segment$location$end
    
    
    # Run DHARMa for T1
    t1 <- run_dharma(segment$model_T1)
    
    # Run DHARMa for T2
    t2 <- run_dharma(segment$model_T2)
    
    
    # Store results
    chromosome_results[[length(chromosome_results) + 1]] <- data.frame(
      
      chr = chr_segment,
      start = start,
      end = end,
      
      DHARMa_uniformity_statistic_T1 =
        t1$uniformity_statistic,
      
      DHARMa_uniformity_pvalue_T1 =
        t1$uniformity_pvalue,
      
      DHARMa_dispersion_statistic_T1 =
        t1$dispersion_statistic,
      
      DHARMa_dispersion_pvalue_T1 =
        t1$dispersion_pvalue,
      
      DHARMa_outlier_statistic_T1 =
        t1$outlier_statistic,
      
      DHARMa_outlier_pvalue_T1 =
        t1$outlier_pvalue,
      
      DHARMa_error_T1 =
        t1$error,
      
      
      DHARMa_uniformity_statistic_T2 =
        t2$uniformity_statistic,
      
      DHARMa_uniformity_pvalue_T2 =
        t2$uniformity_pvalue,
      
      DHARMa_dispersion_statistic_T2 =
        t2$dispersion_statistic,
      
      DHARMa_dispersion_pvalue_T2 =
        t2$dispersion_pvalue,
      
      DHARMa_outlier_statistic_T2 =
        t2$outlier_statistic,
      
      DHARMa_outlier_pvalue_T2 =
        t2$outlier_pvalue,
      
      DHARMa_error_T2 =
        t2$error,
      
      stringsAsFactors = FALSE
    )
  }
  

  chromosome_results <- do.call(
    rbind,
    chromosome_results
  )

  output_file <- file.path(
    out_dir,
    paste0(chr, "_DHARMa.csv")
  )
  
  write.csv(
    chromosome_results,
    output_file,
    row.names = FALSE
  )
  
  cat("Saved:", output_file, "\n\n")
}