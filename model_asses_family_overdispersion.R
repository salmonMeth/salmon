library(data.table)
library(glmmTMB)
library(DHARMa)


path_base="/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated_beta_binom"


model_files <- list.files(path_base,pattern="^model_.*\\.rds$",full.names=TRUE)

path_out <- "/scratch/project_2010912/ezel/salmon/data_files/model_analysis_family_overdispersion"
dir.create(path_out,showWarnings=FALSE,recursive=TRUE)

diagnostic_results <- list()
counter <- 1

chromosome_totals <- list()
chrom_counter <- 1

for(model_file in model_files){
  
  chr_name <- sub("^model_(.*)\\.rds$","\\1",basename(model_file))
  cat("\nProcessing",chr_name,"\n")
  
  chr_out <- file.path(path_out,chr_name)
  dir.create(chr_out,showWarnings=FALSE,recursive=TRUE)
  
  models <- readRDS(model_file)
  
  total_models <- 0
  ma_problem_models <- 0
  pa_problem_models <- 0
  total_problem_models <- 0
  
  for(i in seq_along(models)){
    cat(i,"/",length(models),"\r")
    
    for(tp in c("T1","T2")){
      
      model_name <- paste0("model_",tp)
      if(!is.list(models[[i]][[model_name]])) next
      
      total_models <- total_models+1
      
      m_i <- models[[i]][[model_name]]
      
      res <- tryCatch(
        simulateResiduals(m_i),
        error=function(e){
          message("\nResidual simulation failed: ",chr_name," ",i," ",tp,": ",e$message)
          NULL
        }
      )
      
      if(is.null(res)) next
      
      mf <- model.frame(m_i)
      
      # ID_ma
      ma_levene_p <- NA_real_
      ma_max_sd <- NA_real_
      ma_min_sd <- NA_real_
      ma_sd_ratio <- NA_real_
      ma_max_sd_ID <- NA
      ma_min_sd_ID <- NA
      ma_issue <- FALSE
      
      if("ID_ma" %in% names(mf)){
        
        ma_test <- tryCatch(
          testCategorical(res,catPred=factor(mf$ID_ma),plot=FALSE),
          error=function(e){
            message("\nID_ma test failed: ",chr_name," ",i," ",tp,": ",e$message)
            NULL
          }
        )
        
        if(!is.null(ma_test) && !is.null(ma_test$homogeneity)){
          
          ma_levene_p <- ma_test$homogeneity$`Pr(>F)`[1]
          ma_issue <- !is.na(ma_levene_p) && ma_levene_p < 0.05
          
          ma_sd <- tapply(
            res$scaledResiduals,
            factor(mf$ID_ma),
            sd,
            na.rm=TRUE
          )
          
          ma_sd <- ma_sd[is.finite(ma_sd)]
          
          if(length(ma_sd)>0){
            
            ma_max_sd <- max(ma_sd)
            ma_min_sd <- min(ma_sd)
            ma_sd_ratio <- ma_max_sd/ma_min_sd
            
            ma_max_sd_ID <- names(ma_sd)[which.max(ma_sd)]
            ma_min_sd_ID <- names(ma_sd)[which.min(ma_sd)]
          }
          
          if(ma_issue){
            
            ma_problem_models <- ma_problem_models+1
            
            png(
              file.path(chr_out,paste0("model_",i,"_",tp,"_ID_ma.png")),
              width=1800,height=1400,res=200
            )
            
            plotResiduals(
              res,
              form=factor(mf$ID_ma)
            )
            
            dev.off()
          }
        }
      }
      
      # ID_pa
      pa_levene_p <- NA_real_
      pa_max_sd <- NA_real_
      pa_min_sd <- NA_real_
      pa_sd_ratio <- NA_real_
      pa_max_sd_ID <- NA
      pa_min_sd_ID <- NA
      pa_issue <- FALSE
      
      if("ID_pa" %in% names(mf)){
        
        pa_test <- tryCatch(
          testCategorical(res,catPred=factor(mf$ID_pa),plot=FALSE),
          error=function(e){
            message("\nID_pa test failed: ",chr_name," ",i," ",tp,": ",e$message)
            NULL
          }
        )
        
        if(!is.null(pa_test) && !is.null(pa_test$homogeneity)){
          
          pa_levene_p <- pa_test$homogeneity$`Pr(>F)`[1]
          pa_issue <- !is.na(pa_levene_p) && pa_levene_p < 0.05
          
          pa_sd <- tapply(
            res$scaledResiduals,
            factor(mf$ID_pa),
            sd,
            na.rm=TRUE
          )
          
          pa_sd <- pa_sd[is.finite(pa_sd)]
          
          if(length(pa_sd)>0){
            
            pa_max_sd <- max(pa_sd)
            pa_min_sd <- min(pa_sd)
            pa_sd_ratio <- pa_max_sd/pa_min_sd
            
            pa_max_sd_ID <- names(pa_sd)[which.max(pa_sd)]
            pa_min_sd_ID <- names(pa_sd)[which.min(pa_sd)]
          }
          
          if(pa_issue){
            
            pa_problem_models <- pa_problem_models+1
            
            png(
              file.path(chr_out,paste0("model_",i,"_",tp,"_ID_pa.png")),
              width=1800,height=1400,res=200
            )
            
            plotResiduals(
              res,
              form=factor(mf$ID_pa)
            )
            
            dev.off()
          }
        }
      }
      
      # Save only problematic models
      if(ma_issue || pa_issue){
        
        total_problem_models <- total_problem_models+1
        
        diagnostic_results[[length(diagnostic_results)+1]] <- data.frame(
          chromosome=chr_name,
          location_chr=models[[i]]$location$chr,
          model_index=i,
          model_type=tp,
          segment_start=models[[i]]$location$start,
          segment_end=models[[i]]$location$end,
          
          ID_ma_levene_p=ma_levene_p,
          ID_ma_max_sd=ma_max_sd,
          ID_ma_min_sd=ma_min_sd,
          ID_ma_sd_ratio=ma_sd_ratio,
          ID_ma_max_sd_ID=ma_max_sd_ID,
          ID_ma_min_sd_ID=ma_min_sd_ID,
          
          ID_pa_levene_p=pa_levene_p,
          ID_pa_max_sd=pa_max_sd,
          ID_pa_min_sd=pa_min_sd,
          ID_pa_sd_ratio=pa_sd_ratio,
          ID_pa_max_sd_ID=pa_max_sd_ID,
          ID_pa_min_sd_ID=pa_min_sd_ID
        )
      }
    }
  }
  
  cat("\n")
  
  # Chromosome-level summary
  chromosome_totals[[length(chromosome_totals)+1]] <- data.frame(
    chromosome=chr_name,
    model_type=tp,
    total_models=total_models,
    problematic_models=total_problem_models,
    proportion_problematic=total_problem_models/total_models,
    ID_ma_problematic=ma_problem_models,
    ID_ma_proportion=ma_problem_models/total_models,
    ID_pa_problematic=pa_problem_models,
    ID_pa_proportion=pa_problem_models/total_models
  )
}

# Save problematic models
if(length(diagnostic_results)>0){
  
  diagnostic_results <- do.call(rbind,diagnostic_results)
  
  write.csv(
    diagnostic_results,
    file=file.path(
      path_out,
      "DHARMa_ID_ma_ID_pa_levene_problems.csv"
    ),
    row.names=FALSE
  )
  
  print(diagnostic_results)
  
}else{
  
  cat("\nNo models showed significant Levene tests.\n")
}

# Save chromosome-level summary
chromosome_totals <- do.call(rbind,chromosome_totals)

write.csv(
  chromosome_totals,
  file=file.path(
    path_out,
    "DHARMa_levene_chromosome_summary.csv"
  ),
  row.names=FALSE
)

###get the chr stats from the big csv
problem_file <- file.path(
  path_out,
  "DHARMa_ID_ma_ID_pa_levene_problems.csv"
)

problem_models <- read.csv(problem_file)

summary_results <- list()

for(model_file in model_files){
  
  chr_name <- sub(
    "^model_(.*)\\.rds$",
    "\\1",
    basename(model_file)
  )
  
  models <- readRDS(model_file)
  
  for(tp in c("T1","T2")){
    
    model_name <- paste0("model_",tp)
    
    total_models <- sum(
      sapply(
        models,
        function(x) is.list(x[[model_name]])
      )
    )
    
    problems <- problem_models[
      problem_models$chromosome == chr_name &
        problem_models$model_type == tp,
    ]
    
    problematic_models <- nrow(problems)
    
    ma_problematic <- sum(
      !is.na(problems$ID_ma_levene_p) &
        problems$ID_ma_levene_p < 0.05
    )
    
    pa_problematic <- sum(
      !is.na(problems$ID_pa_levene_p) &
        problems$ID_pa_levene_p < 0.05
    )
    
    summary_results[[length(summary_results)+1]] <- data.frame(
      chromosome=chr_name,
      model_type=tp,
      total_models=total_models,
      problematic_models=problematic_models,
      proportion_problematic=problematic_models/total_models,
      ID_ma_problematic=ma_problematic,
      ID_ma_proportion=ma_problematic/total_models,
      ID_pa_problematic=pa_problematic,
      ID_pa_proportion=pa_problematic/total_models
    )
  }
}

summary_results <- do.call(rbind,summary_results)

write.csv(
  summary_results,
  file=file.path(
    path_out,
    "DHARMa_levene_chromosome_summary_T1_T2.csv"
  ),
  row.names=FALSE
)

print(summary_results)


#check why we have zeros for some ma_sd or pa_sd values

problem_file <- file.path(
  path_out,
  "DHARMa_ID_ma_ID_pa_levene_problems.csv"
)

problems <- read.csv(problem_file)

# ID_ma minimum SD = 0
ma_zero <- problems[
  !is.na(problems$ID_ma_min_sd) &
    problems$ID_ma_min_sd == 0,
]

# ID_pa minimum SD = 0
pa_zero <- problems[
  !is.na(problems$ID_pa_min_sd) &
    problems$ID_pa_min_sd == 0,
]

cat("ID_ma min SD = 0:",nrow(ma_zero),"\n")
cat("ID_pa min SD = 0:",nrow(pa_zero),"\n")

problems <- read.csv(
  file.path(path_out,"DHARMa_ID_ma_ID_pa_levene_problems.csv")
)


zero_sd <- problems[
  problems$ID_ma_min_sd == 0 | problems$ID_pa_min_sd == 0,
]

for(j in seq_len(nrow(zero_sd))){
  
  row <- zero_sd[j,]
  
  model_file <- model_files[
    grepl(
      paste0("^model_",row$chromosome,"\\.rds$"),
      basename(model_files)
    )
  ]
  
  models <- readRDS(model_file)
  
  m_i <- models[[row$model_index]][[
    paste0("model_",row$model_type)
  ]]
  
  cat("\n==============================\n")
  cat("Chromosome:",row$chromosome,"\n")
  cat("Model:",row$model_index,"\n")
  cat("Type:",row$model_type,"\n")
  cat("Segment:",row$segment_start,"-",row$segment_end,"\n")
  
  cat("\nRandom effects:\n")
  print(VarCorr(m_i))
  
  cat("\nFixed effects:\n")
  print(fixef(m_i))
  
  cat("\nFull model summary:\n")
  print(summary(m_i))
}

#find the %s of the over/under and OK dispersed segments


for(f in csv_files){
  
  x <- read.csv(f)
  
  x$DHARMa_dispersion_class <- toupper(
    trimws(x$DHARMa_dispersion_class)
  )
  
  chr <- unique(x$chr)
  
  for(tp in unique(x$TimePoint)){
    
    y <- x[x$TimePoint == tp,]
    
    n_total <- nrow(y)
    
    n_over <- sum(
      y$DHARMa_dispersion_class == "OVERDISPERSION",
      na.rm=TRUE
    )
    
    n_under <- sum(
      y$DHARMa_dispersion_class == "UNDERDISPERSION",
      na.rm=TRUE
    )
    
    n_ok <- sum(
      y$DHARMa_dispersion_class == "OK",
      na.rm=TRUE
    )
    
    results[[length(results)+1]] <- data.frame(
      chr=chr[1],
      TimePoint=tp,
      total_segments=n_total,
      overdispersion_n=n_over,
      overdispersion_percent=100*n_over/n_total,
      underdispersion_n=n_under,
      underdispersion_percent=100*n_under/n_total,
      ok_n=n_ok,
      ok_percent=100*n_ok/n_total
    )
  }
}

results <- do.call(rbind,results)

results <- results[order(results$chr,results$TimePoint),]

write.csv(
  results,
  file=file.path(
    path_out,
    "DHARMa_dispersion_summary_all_chromosomes.csv"
  ),
  row.names=FALSE
)
# Overdispersion segments
over_segments <- do.call(rbind,lapply(csv_files,function(f){
  
  x <- read.csv(f)
  x <- x[x$DHARMa_dispersion_class=="OVERDISPERSION",]
  
  data.frame(
    chr=x$chr,
    start=x$start,
    end=x$end,
    TimePoint=x$TimePoint
  )
}))

# Family variance heterogeneity segments
hetero_segments <- unique(data.frame(
  chr=problems$location_chr,
  start=problems$segment_start,
  end=problems$segment_end,
  TimePoint=problems$model_type
))

# Function to calculate overlap statistics
get_overlap_stats <- function(tp){
  
  over <- unique(over_segments[
    over_segments$TimePoint==tp,
  ])
  
  hetero <- unique(hetero_segments[
    hetero_segments$TimePoint==tp,
  ])
  
  over$key <- paste(over$chr,over$start,over$end,sep="_")
  hetero$key <- paste(hetero$chr,hetero$start,hetero$end,sep="_")
  
  over_keys <- unique(over$key)
  hetero_keys <- unique(hetero$key)
  
  both <- intersect(over_keys,hetero_keys)
  over_only <- setdiff(over_keys,hetero_keys)
  hetero_only <- setdiff(hetero_keys,over_keys)
  
  data.frame(
    TimePoint=tp,
    overdispersion_total=length(over_keys),
    heterogeneity_total=length(hetero_keys),
    both=length(both),
    overdispersion_only=length(over_only),
    heterogeneity_only=length(hetero_only)
  )
}

overlap_results <- rbind(
  get_overlap_stats("T1"),
  get_overlap_stats("T2")
)

print(overlap_results)

write.csv(
  overlap_results,
  file=file.path(
    path_out,
    "overdispersion_family_heterogeneity_overlap_T1_T2.csv"
  ),
  row.names=FALSE
)