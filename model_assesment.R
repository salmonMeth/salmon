library(data.table)
library(glmmTMB)
library(DHARMa)


path_base="/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated"

#for one chr

chr_name <- "NC_059442.1"

models <- readRDS(
  file.path(
    path_base,
    paste0("model_", chr_name, ".rds")
  )
)

safe_min <- function(x){
  if(length(x)==0 || all(is.na(x))) return(NA_real_)
  min(x,na.rm=TRUE)
}

safe_max <- function(x){
  if(length(x)==0 || all(is.na(x))) return(NA_real_)
  max(x,na.rm=TRUE)
}

safe_mean <- function(x){
  if(length(x)==0 || all(is.na(x))) return(NA_real_)
  mean(x,na.rm=TRUE)
}

safe_sd <- function(x){
  if(length(x)==0 || all(is.na(x))) return(NA_real_)
  sd(x,na.rm=TRUE)
}

assess_model <- function(model_object){
  
  if(is.null(model_object)) return(NULL)
  
  out <- list()
  
  out$convergence_code <- tryCatch(model_object$fit$convergence,error=function(e) NA)
  out$convergence_message <- tryCatch(model_object$fit$message,error=function(e) NA)
  out$pdHess <- if(is.null(model_object$sdr)) NA else model_object$sdr$pdHess
  
  out$AIC <- tryCatch(AIC(model_object),error=function(e) NA)
  
  res <- tryCatch(residuals(model_object,type="pearson"),error=function(e) numeric(0))
  
  out$n_obs <- length(res)
  out$pearson_mean <- safe_mean(res)
  out$pearson_sd <- safe_sd(res)
  out$pearson_min <- safe_min(res)
  out$pearson_max <- safe_max(res)
  out$n_abs2 <- sum(abs(res)>2,na.rm=TRUE)
  out$n_abs3 <- sum(abs(res)>3,na.rm=TRUE)
  
  vc <- tryCatch(VarCorr(model_object),error=function(e) NULL)
  
  out$sd_ID_ma <- NA_real_
  out$sd_ID_pa <- NA_real_
  
  if(!is.null(vc)){
    if("ID_ma" %in% names(vc$cond))
      out$sd_ID_ma <- attr(vc$cond$ID_ma,"stddev")
    if("ID_pa" %in% names(vc$cond))
      out$sd_ID_pa <- attr(vc$cond$ID_pa,"stddev")
  }
  
  out$large_random_effect <- 
    (!is.na(out$sd_ID_ma) && out$sd_ID_ma>0.01) ||
    (!is.na(out$sd_ID_pa) && out$sd_ID_pa>0.01)
  
  out$small_random_effect_1e4 <- 
    (!is.na(out$sd_ID_ma) &&
       !is.na(out$sd_ID_pa) &&
       out$sd_ID_ma<1e-4 &&
       out$sd_ID_pa<1e-4)
  
  pred <- tryCatch(predict(model_object,type="response"),error=function(e) numeric(0))
  
  out$pred_mean <- safe_mean(pred)
  out$pred_sd <- safe_sd(pred)
  out$pred_min <- safe_min(pred)
  out$pred_max <- safe_max(pred)
  
  out$DHARMa_uniformity <- NA_real_
  out$DHARMa_outliers <- NA_real_
  out$DHARMa_dispersion <- NA_real_
  out$DHARMa_dispersion_stat <- NA_real_
  out$DHARMa_dispersion_class <- NA_character_
  
  tryCatch({
    
    sim <- simulateResiduals(model_object,plot=FALSE,conditional=TRUE)
    
    out$DHARMa_uniformity <- testUniformity(sim)$p.value
    out$DHARMa_outliers <- testOutliers(sim)$p.value
    
    disp_test <- testDispersion(sim)
    
    out$DHARMa_dispersion <- disp_test$p.value
    out$DHARMa_dispersion_stat <- unname(disp_test$statistic)
    
    out$DHARMa_dispersion_class <- if(disp_test$p.value>=0.05){
      "OK"
    }else if(disp_test$statistic>1){
      "OVERDISPERSION"
    }else{
      "UNDERDISPERSION"
    }
    
  },error=function(e){
    out$DHARMa_dispersion_class <- "DHARMa_ERROR"
  })
  
  out
}
results <- vector("list", length(models)*2)

counter <- 1

for(i in seq_along(models)){
  
  cat(i,"/",length(models),"\r")
  
  for(tp in c("T1","T2")){
    
    model_name <- paste0("model_",tp)
    
    if(!is.list(models[[i]][[model_name]])) next
    
    info <- assess_model(models[[i]][[model_name]])
    
    if(is.null(info)) next
    
    results[[counter]] <- data.table(
      
      model_index=i,
      
      chr=models[[i]]$location$chr,
      start=models[[i]]$location$start,
      end=models[[i]]$location$end,
      TimePoint=tp,
      
      convergence_code=info$convergence_code,
      convergence_message=info$convergence_message,
      pdHess=info$pdHess,
      
      AIC=info$AIC,
      n_obs=info$n_obs,
      
      pearson_mean=info$pearson_mean,
      pearson_sd=info$pearson_sd,
      pearson_min=info$pearson_min,
      pearson_max=info$pearson_max,
      n_abs2=info$n_abs2,
      n_abs3=info$n_abs3,
      
      sd_ID_ma=info$sd_ID_ma,
      sd_ID_pa=info$sd_ID_pa,
      
      large_random_effect=info$large_random_effect,
      small_random_effect_1e4=info$small_random_effect_1e4,
      
      pred_mean=info$pred_mean,
      pred_sd=info$pred_sd,
      pred_min=info$pred_min,
      pred_max=info$pred_max,
      
      DHARMa_uniformity=info$DHARMa_uniformity,
      DHARMa_dispersion=info$DHARMa_dispersion,
      DHARMa_dispersion_stat=info$DHARMa_dispersion_stat,
      DHARMa_dispersion_class=info$DHARMa_dispersion_class,
      DHARMa_outliers=info$DHARMa_outliers
    )
    
    counter <- counter+1
  }
}

qc <- rbindlist(results,fill=TRUE)

write_qc_report <- function(qc,chr_name,out_dir,parent_effects=NULL){  
  sink(file.path(out_dir, paste0("QC_report_", chr_name, ".txt")))
  
  cat("MODEL QUALITY CONTROL REPORT\n")
  cat("============================\n\n")
  
  cat("Chromosome:", chr_name, "\n")
  cat("Models assessed:", nrow(qc), "\n")
  cat("Segments:", length(unique(qc$model_index)), "\n\n")
  
  
  cat("SUMMARY BY TIMEPOINT\n")
  cat("===================\n\n")
  
  print(
    qc[, .(
      n_models=.N,
      convergence_fail=sum(convergence_code!=0, na.rm=TRUE),
      non_positive_hessian=sum(!pdHess, na.rm=TRUE),
      mean_AIC=mean(AIC, na.rm=TRUE),
      median_AIC=median(AIC, na.rm=TRUE)
    ), by=TimePoint]
  )
  
  
  cat("\n\nCONVERGENCE\n")
  cat("===========\n\n")
  
  print(
    qc[, table(convergence_code, useNA="ifany"), by=TimePoint]
  )
  
  cat("\nOptimizer messages:\n")
  
  print(
    qc[, .N, by=.(TimePoint, convergence_message)][
      order(TimePoint, -N)
    ]
  )
  
  
  cat("\n\nHESSIAN STATUS\n")
  cat("==============\n\n")
  
  print(
    qc[, table(pdHess, useNA="ifany"), by=TimePoint]
  )
  
  
  cat("\n\nAIC SUMMARY\n")
  cat("===========\n\n")
  print(
    qc[, .(
      min_AIC=min(AIC, na.rm=TRUE),
      q1_AIC=quantile(AIC, 0.25, na.rm=TRUE),
      median_AIC=median(AIC, na.rm=TRUE),
      mean_AIC=mean(AIC, na.rm=TRUE),
      q3_AIC=quantile(AIC, 0.75, na.rm=TRUE),
      max_AIC=max(AIC, na.rm=TRUE)
    ), by=TimePoint]
  )
  
  cat("\n\nPEARSON RESIDUALS\n")
  cat("=================\n\n")
  
  print(
    qc[, .(
      mean_res_sd=mean(pearson_sd, na.rm=TRUE),
      median_res_sd=median(pearson_sd, na.rm=TRUE),
      mean_abs2=mean(n_abs2, na.rm=TRUE),
      mean_abs3=mean(n_abs3, na.rm=TRUE)
    ), by=TimePoint]
  )
  
  
  cat("\n\nRANDOM EFFECT SUMMARY\n")
  cat("=====================\n\n")
  
  print(
    qc[, .(
      mean_ID_ma=mean(sd_ID_ma, na.rm=TRUE),
      median_ID_ma=median(sd_ID_ma, na.rm=TRUE),
      mean_ID_pa=mean(sd_ID_pa, na.rm=TRUE),
      median_ID_pa=median(sd_ID_pa, na.rm=TRUE)
    ), by=TimePoint]
  )
  
  if(!is.null(parent_effects)){
    
    cat("\n\nPARENT EFFECT DISTRIBUTION SUMMARY\n")
    cat("==================================\n\n")
    
    print(
      parent_effects[
        ,
        .(
          n=.N,
          fraction_near_zero=mean(abs(effect)<1e-4,na.rm=TRUE),
          mean_abs_effect=mean(abs(effect),na.rm=TRUE),
          median_abs_effect=median(abs(effect),na.rm=TRUE),
          max_abs_effect=max(abs(effect),na.rm=TRUE)
        ),
        by=.(dispersion_class,parent_type)
      ]
    )
  }
  cat("\n\nPREDICTED METHYLATION BY DISPERSION CLASS\n")
  cat("========================================\n\n")
  
  print(
    qc[
      DHARMa_dispersion_class != "OK",
      .(
        n_models=.N,
        mean_prediction=mean(pred_mean,na.rm=TRUE),
        median_prediction=median(pred_mean,na.rm=TRUE),
        min_prediction=min(pred_mean,na.rm=TRUE),
        max_prediction=max(pred_mean,na.rm=TRUE),
        
        mean_pred_sd=mean(pred_sd,na.rm=TRUE),
        median_pred_sd=median(pred_sd,na.rm=TRUE)
      ),
      by=.(TimePoint,DHARMa_dispersion_class)
    ]
  )
  
  cat("\nRandom effects > 0.01\n")
  cat("======================\n\n")
  
  print(
    qc[, .(
      large_RE=sum(large_random_effect, na.rm=TRUE),
      total=.N
    ), by=TimePoint]
  )
  
  
  cat("\nBoth ID_ma and ID_pa < 1e-4\n")
  cat("============================\n\n")
  
  print(
    qc[, .(
      both_RE_near_zero=sum(small_random_effect_1e4, na.rm=TRUE),
      total=.N
    ), by=TimePoint]
  )
  
  
  cat("\n\nPREDICTIONS\n")
  cat("===========\n\n")
  
  print(
    qc[, .(
      mean_prediction=mean(pred_mean, na.rm=TRUE),
      sd_prediction=mean(pred_sd, na.rm=TRUE),
      min_prediction=min(pred_min, na.rm=TRUE),
      max_prediction=max(pred_max, na.rm=TRUE)
    ), by=TimePoint]
  )
  
  
  cat("\n\nDHARMa DIAGNOSTICS\n")
  cat("==================\n\n")
  
  print(
    qc[, .(
      uniformity_fail=sum(DHARMa_uniformity<0.05, na.rm=TRUE),
      dispersion_fail=sum(DHARMa_dispersion<0.05, na.rm=TRUE),
      outlier_fail=sum(DHARMa_outliers<0.05, na.rm=TRUE),
      total=.N
    ), by=TimePoint]
  )
  
  
  cat("\n\nUNDER/OVER DISPERSION\n")
  cat("=====================\n\n")
  
  print(
    qc[, .N, by=.(TimePoint, DHARMa_dispersion_class)]
  )
  
  under <- qc[
    DHARMa_dispersion_class=="UNDERDISPERSION",
    .(
      model_index,
      chr,
      start,
      end,
      TimePoint,
      pred_mean,
      pred_sd,
      pred_min,
      pred_max,
      sd_ID_ma,
      sd_ID_pa,
      DHARMa_dispersion,
      DHARMa_dispersion_stat
    )
  ]
  
  cat("\n\nUNDERDISPERSED SEGMENTS\n")
  cat("=======================\n\n")
  cat("Number:", nrow(under), "\n\n")
  
  if(nrow(under)>0)
    print(under[seq_len(min(50, nrow(under)))])
  over <- qc[
    DHARMa_dispersion_class=="OVERDISPERSION",
    .(
      model_index,
      chr,
      start,
      end,
      TimePoint,
      pred_mean,
      pred_sd,
      pred_min,
      pred_max,
      sd_ID_ma,
      sd_ID_pa,
      DHARMa_dispersion,
      DHARMa_dispersion_stat
    )
  ]
  
  cat("\n\nOVERDISPERSED SEGMENTS\n")
  cat("======================\n\n")
  cat("Number:", nrow(over), "\n\n")
  
  if(nrow(over)>0)
    print(over[seq_len(min(50, nrow(over)))])
  
  cat("\n\nRANDOM EFFECT SDs BY DISPERSION CLASS\n")
  cat("=====================================\n\n")
  
  print(
    qc[
      DHARMa_dispersion_class != "OK",
      .(
        n_models=.N,
        mean_ID_ma=mean(sd_ID_ma,na.rm=TRUE),
        median_ID_ma=median(sd_ID_ma,na.rm=TRUE),
        min_ID_ma=min(sd_ID_ma,na.rm=TRUE),
        max_ID_ma=max(sd_ID_ma,na.rm=TRUE),
        
        mean_ID_pa=mean(sd_ID_pa,na.rm=TRUE),
        median_ID_pa=median(sd_ID_pa,na.rm=TRUE),
        min_ID_pa=min(sd_ID_pa,na.rm=TRUE),
        max_ID_pa=max(sd_ID_pa,na.rm=TRUE)
      ),
      by=.(TimePoint,DHARMa_dispersion_class)
    ]
  )
  
  large <- qc[
    large_random_effect==TRUE,
    .(
      model_index,
      chr,
      start,
      end,
      TimePoint,
      sd_ID_ma,
      sd_ID_pa
    )
  ]
  
  cat("\n\nLARGE RANDOM EFFECT SEGMENTS\n")
  cat("============================\n\n")
  cat("Number:", nrow(large), "\n\n")
  
  if(nrow(large)>0)
    print(large[seq_len(min(50, nrow(large)))])
  
  
  sink()
}


####get the parental effects to plot or examine them
extract_parent_effects <- function(models,qc){
  
  parent_effects <- list()
  k <- 1
  
  for(i in seq_along(models)){
    
    for(tp in c("T1","T2")){
      
      model_name <- paste0("model_",tp)
      
      if(!is.list(models[[i]][[model_name]]))
        next
      
      cls <- qc[
        model_index==i & TimePoint==tp,
        DHARMa_dispersion_class
      ]
      re <- ranef(models[[i]][[model_name]],condVar=FALSE)$cond
      
      if(!is.null(re$ID_ma)){
        
        ma <- data.table(
          model_index=i,
          TimePoint=tp,
          dispersion_class=cls,
          parent=rownames(re$ID_ma),
          effect=as.numeric(re$ID_ma[,1]),
          parent_type="ID_ma"
        )
        
        parent_effects[[k]] <- ma
        k <- k+1
      }
      
      if(!is.null(re$ID_pa)){
        
        pa <- data.table(
          model_index=i,
          TimePoint=tp,
          dispersion_class=cls,
          parent=rownames(re$ID_pa),
          effect=as.numeric(re$ID_pa[,1]),
          parent_type="ID_pa"
        )
        
        parent_effects[[k]] <- pa
        k <- k+1
      }
    }
  }
  
  return(rbindlist(parent_effects,fill=TRUE))
}
parent_effects=extract_parent_effects(models,qc)

write_qc_report(qc,chr_name,out_dir,parent_effects)
####plotting

##

save_qc_plots <- function(segment_RE,qc,out_dir){
  
  dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
  
  labels <- c(
    UNDERDISPERSION="Under",
    OK="OK",
    OVERDISPERSION="Over"
  )
  
  segment_RE <- copy(segment_RE)
  
  segment_RE[,ID_ma_sd_plot:=ifelse(ID_ma_sd<=0,1e-6,ID_ma_sd)]
  segment_RE[,ID_pa_sd_plot:=ifelse(ID_pa_sd<=0,1e-6,ID_pa_sd)]
  
  
  ## maternal random effects
  
  if(sum(!is.na(segment_RE$ID_ma_sd_plot))>0){
    
    p_ma <- ggplot(
      segment_RE,
      aes(dispersion_class,ID_ma_sd_plot,fill=dispersion_class)
    )+
      geom_boxplot(outlier.alpha=.3)+
      scale_y_log10()+
      scale_x_discrete(labels=labels)+
      scale_fill_discrete(labels=labels)+
      labs(
        x="Dispersion class",
        y="SD of maternal random effects",
        fill="Dispersion"
      )+
      theme_gray()
    
    ggsave(
      file.path(out_dir,"ID_ma_random_effect_SD.png"),
      p_ma,
      width=7,
      height=5,
      dpi=300
    )
  }
  
  
  ## paternal random effects
  
  if(sum(!is.na(segment_RE$ID_pa_sd_plot))>0){
    
    p_pa <- ggplot(
      segment_RE,
      aes(dispersion_class,ID_pa_sd_plot,fill=dispersion_class)
    )+
      geom_boxplot(outlier.alpha=.3)+
      scale_y_log10()+
      scale_x_discrete(labels=labels)+
      scale_fill_discrete(labels=labels)+
      labs(
        x="Dispersion class",
        y="SD of paternal random effects",
        fill="Dispersion"
      )+
      theme_gray()
    
    ggsave(
      file.path(out_dir,"ID_pa_random_effect_SD.png"),
      p_pa,
      width=7,
      height=5,
      dpi=300
    )
  }
  
  
  ## predicted methylation
  
  if(sum(!is.na(qc$pred_mean))>0){
    
    p_meth <- ggplot(
      qc,
      aes(
        DHARMa_dispersion_class,
        pred_mean,
        fill=DHARMa_dispersion_class
      )
    )+
      geom_boxplot(outlier.alpha=.3)+
      scale_x_discrete(labels=labels)+
      scale_fill_discrete(labels=labels)+
      labs(
        x="Dispersion class",
        y="Mean predicted methylation",
        fill="Dispersion"
      )+
      theme_gray()
    
    ggsave(
      file.path(out_dir,"predicted_methylation_by_dispersion.png"),
      p_meth,
      width=7,
      height=5,
      dpi=300
    )
  }
  
  
  ## predicted methylation range
  
  if(all(c("pred_min","pred_max") %in% colnames(qc))){
    
    qc_range <- melt(
      qc,
      id.vars=c("TimePoint","DHARMa_dispersion_class"),
      measure.vars=c("pred_min","pred_max"),
      variable.name="prediction_type",
      value.name="prediction"
    )
    
    p_range <- ggplot(
      qc_range,
      aes(
        DHARMa_dispersion_class,
        prediction,
        fill=prediction_type
      )
    )+
      geom_boxplot(outlier.alpha=.2)+
      scale_x_discrete(labels=labels)+
      labs(
        x="Dispersion class",
        y="Predicted methylation",
        fill="Prediction"
      )+
      theme_gray()
    
    ggsave(
      file.path(out_dir,"prediction_range_by_dispersion.png"),
      p_range,
      width=7,
      height=5,
      dpi=300
    )
  }
  
  cat("Plots saved to:",out_dir,"\n")
}

save_qc_plots(parent_effects,qc,path_out)

#########################
#########################
#########################
#loop over all chrs
model_files <- list.files(path_base,pattern="^model_.*\\.rds$",full.names=TRUE)

path_out <- "/scratch/project_2010912/ezel/salmon/data_files/model_analysis"

dir.create(path_out,showWarnings=FALSE,recursive=TRUE)
#####
for(model_file in model_files){
  
  chr_name <- sub("^model_(.*)\\.rds$","\\1",basename(model_file))
  cat("\nProcessing",chr_name,"\n")
  
  chr_out <- file.path(path_out,chr_name)
  dir.create(chr_out,showWarnings=FALSE,recursive=TRUE)
  ###
  
  models <- readRDS(model_file)
  results <- vector("list",length(models)*2)
  counter <- 1
  
  for(i in seq_along(models)){
    cat(i,"/",length(models),"\r")
    
    for(tp in c("T1","T2")){
      model_name <- paste0("model_",tp)
      if(!is.list(models[[i]][[model_name]])) next
      
      info <- assess_model(models[[i]][[model_name]])
      if(is.null(info)) next
      
      results[[counter]] <- data.table(
        model_index=i,
        chr=models[[i]]$location$chr,
        start=models[[i]]$location$start,
        end=models[[i]]$location$end,
        TimePoint=tp,
        convergence_code=info$convergence_code,
        convergence_message=info$convergence_message,
        pdHess=info$pdHess,
        AIC=info$AIC,
        n_obs=info$n_obs,
        pearson_mean=info$pearson_mean,
        pearson_sd=info$pearson_sd,
        pearson_min=info$pearson_min,
        pearson_max=info$pearson_max,
        n_abs2=info$n_abs2,
        n_abs3=info$n_abs3,
        sd_ID_ma=info$sd_ID_ma,
        sd_ID_pa=info$sd_ID_pa,
        large_random_effect=info$large_random_effect,
        small_random_effect_1e4=info$small_random_effect_1e4,
        pred_mean=info$pred_mean,
        pred_sd=info$pred_sd,
        pred_min=info$pred_min,
        pred_max=info$pred_max,
        DHARMa_uniformity=info$DHARMa_uniformity,
        DHARMa_dispersion=info$DHARMa_dispersion,
        DHARMa_dispersion_stat=info$DHARMa_dispersion_stat,
        DHARMa_dispersion_class=info$DHARMa_dispersion_class,
        DHARMa_outliers=info$DHARMa_outliers
      )
      counter <- counter+1
    }
  }
  
  qc <- rbindlist(results,fill=TRUE)
  
  parent_effects <- extract_parent_effects(models,qc)
  
  segment_RE <- parent_effects[,.(ID_ma_sd=sd(effect[parent_type=="ID_ma"],na.rm=TRUE),
                                  ID_pa_sd=sd(effect[parent_type=="ID_pa"],na.rm=TRUE)),
                               by=.(model_index,TimePoint,dispersion_class)]
  
  fwrite(qc,file.path(chr_out,paste0("QC_table_",chr_name,".csv")))
  fwrite(parent_effects,file.path(chr_out,paste0("parent_effects_",chr_name,".csv")))
  
  write_qc_report(qc,chr_name,chr_out,parent_effects)
  save_qc_plots(segment_RE,qc,chr_out)
  
  cat("\nFinished",chr_name,"\n")
}