library ( glmmTMB )
library(data.table)

#we are modeling the methyl&unmethyl CpG count using weight, mom and dad IDs

path_ids="/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"
meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

# this model pools T1 and T2 samples
################################################
################################################
#loop over all chrs
base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"
out_dir <- file.path(base_dir, "model_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cs_files <- list.files(
  base_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)
for(cs_file in cs_files){
  
  chr_name <- sub(
    "_segment_cs\\.csv$",
    "",
    basename(cs_file)
  )
  
  cov_file <- file.path(
    base_dir,
    paste0(chr_name, "_segment_coverage.csv")
  )
  
  if(!file.exists(cov_file)){
    cat("Missing coverage file for", chr_name, "\n")
    next
  }
  
  cat("Processing", chr_name, "\n")
  
  meth <- fread(cs_file)
  cov <- fread(cov_file)
  
  meth_long <- melt(
    meth,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "coverage"
  )
  
  dat <- merge(
    meth_long,
    cov_long,
    by = c("chr", "start", "end", "sample")
  )
  
  dat <- merge(
    dat,
    meta,
    by = "sample",
    all.x = TRUE
  )
  
  dat[, unmeth_reads := coverage - meth_reads]
  
  regions <- unique(
    dat[, .(chr, start, end)]
  )
  
  fit_region <- function(i){
    
    d <- dat[
      chr == regions$chr[i] &
        start == regions$start[i] &
        end == regions$end[i]
    ]
    
    tryCatch({
      
      model <- glmmTMB(
        cbind(meth_reads, unmeth_reads) ~
          Weight.g +
          (1 | ID_ma) +
          (1 | ID_pa),
        data = d,
        family = betabinomial(link = "logit")
      )
      
      coefs <- summary(model)$coefficients$cond
      
      if(!"Weight.g" %in% rownames(coefs))
        return(NULL)
      
      data.table(
        chr = regions$chr[i],
        start = regions$start[i],
        end = regions$end[i],
        estimate = unname(coefs["Weight.g", "Estimate"]),
        SE = unname(coefs["Weight.g", "Std. Error"]),
        z = unname(coefs["Weight.g", "z value"]),
        pvalue = unname(coefs["Weight.g", "Pr(>|z|)"])
      )
      
    }, error = function(e){
      NULL
    })
  }
  
  results <- rbindlist(
    lapply(seq_len(nrow(regions)), fit_region),
    fill = TRUE
  )
  
  out_file <- file.path(
    out_dir,
    paste0("model_", chr_name, ".csv")
  )
  
  fwrite(
    results,
    out_file
  )
  
  cat(
    "Saved",
    nrow(results),
    "regions to",
    out_file,
    "\n\n"
  )
}

cat("Done.\n")

#################################################
#################################################
#################################################
#this model separates T1 and T2 samples.
base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"
out_dir <- file.path(base_dir, "model_results_T1_T2_separated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cs_files <- list.files(
  base_dir,
  pattern = "_segment_cs\\.csv$",
  full.names = TRUE
)
for(cs_file in cs_files){
  
  chr_name <- sub(
    "_segment_cs\\.csv$",
    "",
    basename(cs_file)
  )
  
  cov_file <- file.path(
    base_dir,
    paste0(chr_name, "_segment_coverage.csv")
  )
  
  if(!file.exists(cov_file)){
    cat("Missing coverage file for", chr_name, "\n")
    next
  }
  
  cat("Processing", chr_name, "\n")
  
  meth <- fread(cs_file)
  cov <- fread(cov_file)
  
  meth_long <- melt(
    meth,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = c("chr", "start", "end"),
    variable.name = "sample",
    value.name = "coverage"
  )
  
  dat <- merge(
    meth_long,
    cov_long,
    by = c("chr", "start", "end", "sample")
  )
  
  dat <- merge(
    dat,
    meta,
    by = "sample",
    all.x = TRUE
  )
  
  dat[, unmeth_reads := coverage - meth_reads]
  
  regions <- unique(
    dat[, .(chr, start, end)]
  )
  #this is for the scaled weight variable
  fit_region_sc <- function(i){
    
    d <- dat[
      chr == regions$chr[i] &
        start == regions$start[i] &
        end == regions$end[i]
    ]
    
    d_T1 <- d[TimePoint == "T1"]
    d_T2 <- d[TimePoint == "T2"]
    
    out <- data.table(
      chr = regions$chr[i],
      start = regions$start[i],
      end = regions$end[i],
      
      estimate_T1 = NA_real_,
      SE_T1 = NA_real_,
      z_T1 = NA_real_,
      pvalue_T1 = NA_real_,
      status_T1 = "NOT_RUN",
      
      estimate_T2 = NA_real_,
      SE_T2 = NA_real_,
      z_T2 = NA_real_,
      pvalue_T2 = NA_real_,
      status_T2 = "NOT_RUN"
    )
    
    ####################################################
    # T1
    
    if(
      nrow(d_T1) >= 10 &&
      length(unique(d_T1$Weight.g)) > 1
    ){
      
      tryCatch({
        
        d_T1$Weight_sc <- as.numeric(
          scale(d_T1$Weight.g)
        )
        
        model_T1 <- glmmTMB(
          cbind(
            meth_reads,
            unmeth_reads
          ) ~
            Weight_sc +
            (1 | ID_ma) +
            (1 | ID_pa),
          data = d_T1,
          family = betabinomial(link = "logit")
        )
        
        coefs <- summary(model_T1)$coefficients$cond
        
        out$estimate_T1 <- unname(
          coefs["Weight_sc", "Estimate"]
        )
        
        out$SE_T1 <- unname(
          coefs["Weight_sc", "Std. Error"]
        )
        
        out$z_T1 <- unname(
          coefs["Weight_sc", "z value"]
        )
        
        out$pvalue_T1 <- unname(
          coefs["Weight_sc", "Pr(>|z|)"]
        )
        
        if(
          !is.null(model_T1$sdr) &&
          isTRUE(model_T1$sdr$pdHess)
        ){
          out$status_T1 <- "OK"
        } else {
          out$status_T1 <- "NON_POSITIVE_HESSIAN"
        }
        
      }, error = function(e){
        
        out$status_T1 <<- paste0(
          "ERROR: ",
          conditionMessage(e)
        )
        
      })
      
    } else {
      
      out$status_T1 <- "INSUFFICIENT_DATA"
      
    }
    
    ####################################################
    # T2
    
    if(
      nrow(d_T2) >= 10 &&
      length(unique(d_T2$Weight.g)) > 1
    ){
      
      tryCatch({
        
        d_T2$Weight_sc <- as.numeric(
          scale(d_T2$Weight.g)
        )
        
        model_T2 <- glmmTMB(
          cbind(
            meth_reads,
            unmeth_reads
          ) ~
            Weight_sc +
            (1 | ID_ma) +
            (1 | ID_pa),
          data = d_T2,
          family = betabinomial(link = "logit")
        )
        
        coefs <- summary(model_T2)$coefficients$cond
        
        out$estimate_T2 <- unname(
          coefs["Weight_sc", "Estimate"]
        )
        
        out$SE_T2 <- unname(
          coefs["Weight_sc", "Std. Error"]
        )
        
        out$z_T2 <- unname(
          coefs["Weight_sc", "z value"]
        )
        
        out$pvalue_T2 <- unname(
          coefs["Weight_sc", "Pr(>|z|)"]
        )
        
        if(
          !is.null(model_T2$sdr) &&
          isTRUE(model_T2$sdr$pdHess)
        ){
          out$status_T2 <- "OK"
        } else {
          out$status_T2 <- "NON_POSITIVE_HESSIAN"
        }
        
      }, error = function(e){
        
        out$status_T2 <<- paste0(
          "ERROR: ",
          conditionMessage(e)
        )
        
      })
      
    } else {
      
      out$status_T2 <- "INSUFFICIENT_DATA"
      
    }
    
    return(out)
  }
  
  fit_region <- function(i){
    
    d <- dat[
      chr == regions$chr[i] &
        start == regions$start[i] &
        end == regions$end[i]
    ]
    
    d_T1 <- d[TimePoint == "T1"]
    d_T2 <- d[TimePoint == "T2"]
    
    out <- data.table(
      chr = regions$chr[i],
      start = regions$start[i],
      end = regions$end[i],
      
      estimate_T1 = NA_real_,
      SE_T1 = NA_real_,
      z_T1 = NA_real_,
      pvalue_T1 = NA_real_,
      status_T1 = NA_character_,
      
      estimate_T2 = NA_real_,
      SE_T2 = NA_real_,
      z_T2 = NA_real_,
      pvalue_T2 = NA_real_,
      status_T2 = NA_character_
    )
    
    ## ---------- T1 ----------
    tryCatch({
      
      model_T1 <- glmmTMB(
        cbind(meth_reads, unmeth_reads) ~
          Weight.g +
          (1 | ID_ma) +
          (1 | ID_pa),
        data = d_T1,
        family = betabinomial(link = "logit")
      )
      
      coefs <- summary(model_T1)$coefficients$cond
      
      out$estimate_T1 <- unname(coefs["Weight.g","Estimate"])
      out$SE_T1 <- unname(coefs["Weight.g","Std. Error"])
      out$z_T1 <- unname(coefs["Weight.g","z value"])
      out$pvalue_T1 <- unname(coefs["Weight.g","Pr(>|z|)"])
      
      if(isTRUE(model_T1$sdr$pdHess)){
        out$status_T1 <- "OK"
      } else {
        out$status_T1 <- "NON_POSITIVE_HESSIAN"
      }
      
    }, error = function(e){
      
      out$status_T1 <- paste0(
        "ERROR: ",
        conditionMessage(e)
      )
      
    })
    
    ## ---------- T2 ----------
    tryCatch({
      
      model_T2 <- glmmTMB(
        cbind(meth_reads, unmeth_reads) ~
          Weight.g +
          (1 | ID_ma) +
          (1 | ID_pa),
        data = d_T2,
        family = betabinomial(link = "logit")
      )
      
      coefs <- summary(model_T2)$coefficients$cond
      
      out$estimate_T2 <- unname(coefs["Weight.g","Estimate"])
      out$SE_T2 <- unname(coefs["Weight.g","Std. Error"])
      out$z_T2 <- unname(coefs["Weight.g","z value"])
      out$pvalue_T2 <- unname(coefs["Weight.g","Pr(>|z|)"])
      
      if(isTRUE(model_T2$sdr$pdHess)){
        out$status_T2 <- "OK"
      } else {
        out$status_T2 <- "NON_POSITIVE_HESSIAN"
      }
      
    }, error = function(e){
      
      out$status_T2 <- paste0(
        "ERROR: ",
        conditionMessage(e)
      )
      
    })
    
    out
  }
  
  results <- rbindlist(
    lapply(seq_len(nrow(regions)), fit_region_sc),
    fill = TRUE
  )
  
  out_file <- file.path(
    out_dir,
    paste0("model_", chr_name, ".csv")
  )
  
  fwrite(
    results,
    out_file
  )
  
  cat(
    "Saved",
    nrow(results),
    "regions to",
    out_file,
    "\n\n"
  )
}

cat("Done.\n")


#################################################
#################################################
#examine the non converging results
result_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed/model_results_T1_T2_separated"

files <- list.files(
  result_dir,
  pattern = "^model_.*\\.csv$",
  full.names = TRUE
)
bad_T1 <- rbindlist(
  lapply(files, function(f){
    x <- fread(f)
    x[
      status_T1 != "OK"
    ][, file := basename(f)]
  }),
  fill = TRUE
)
bad_T2 <- rbindlist(
  lapply(files, function(f){
    x <- fread(f)
    x[
      status_T2 != "OK"
    ][, file := basename(f)]
  }),
  fill = TRUE
)
str(bad_T1)
str(bad_T2)
bad_T1_T2 <- rbindlist(
  lapply(files, function(f){
    x <- fread(f)
    x[
      status_T2 != "OK" & (status_T1 != "OK" )
    ][, file := basename(f)]
  }),
  fill = TRUE
)
str(bad_T1_T2)
segment_summary=fread("/scratch/project_2010912/ezel/salmon/data_files/segmentation_output_3_samples_removed/summaries/chromosome_summary.csv")
tot_segs=sum(segment_summary$n_segments)
nrow(bad_T1)/tot_segs
nrow(bad_T2)/tot_segs

