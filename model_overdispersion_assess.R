#we want to check the estimated overdispersion parameters of the betabinomial model
#to see if we can instead use a binomial
library(data.table)
library(glmmTMB)
base_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segment_cs_cov_3_samples_removed"

path_ids <- "/scratch/project_2010912/ezel/hannu_methylation.selectedIDs_10Dec24.txt"

out_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated"

out_dir_2="/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated/overdispersion"

dir.create(out_dir,recursive = TRUE,showWarnings = FALSE)
dir.create(out_dir_2,recursive = TRUE,showWarnings = FALSE)


#########################################################

meta <- fread(path_ids)
meta[, sample := paste(Running.ID, TimePoint, sep = "_")]

#########################################################
#model

fit_region_phi <- function(i){
  
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
    
    phi_T1 = NA_real_,
    status_T1 = "NOT_RUN",
    
    phi_T2 = NA_real_,
    status_T2 = "NOT_RUN"
  )
  
  #####################################################
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
      
      out$phi_T1 <- sigma(model_T1)
      
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
  
  #####################################################
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
      
      out$phi_T2 <- sigma(model_T2)
      
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
  
  out
}

#########################################################
#loop to fit the model

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
    paste0(
      chr_name,
      "_segment_coverage.csv"
    )
  )
  
  if(!file.exists(cov_file)){
    
    cat(
      "Missing coverage file for",
      chr_name,
      "\n"
    )
    
    next
  }
  
  cat(
    "\nProcessing",
    chr_name,
    "\n"
  )
  meth <- fread(cs_file)
  
  cov <- fread(cov_file)
  
  meth_long <- melt(
    meth,
    id.vars = c(
      "chr",
      "start",
      "end"
    ),
    variable.name = "sample",
    value.name = "meth_reads"
  )
  
  cov_long <- melt(
    cov,
    id.vars = c(
      "chr",
      "start",
      "end"
    ),
    variable.name = "sample",
    value.name = "coverage"
  )
  
  dat <- merge(
    meth_long,
    cov_long,
    by = c(
      "chr",
      "start",
      "end",
      "sample"
    )
  )
  
  dat <- merge(
    dat,
    meta,
    by = "sample",
    all.x = TRUE
  )
  
  dat[
    ,
    unmeth_reads := coverage - meth_reads
  ]
  
  regions <- unique(
    dat[
      ,
      .(
        chr,
        start,
        end
      )
    ]
  )
    # Fit the model
  
  phi_results <- rbindlist(
    lapply(
      seq_len(nrow(regions)),
      fit_region_phi
    ),
    fill = TRUE
  )
  
  #####################################################
  # SAVE
  fwrite(
    phi_results,
    file.path(
      out_dir_2,
      paste0(
        chr_name,
        "_overdispersion.csv"
      )
    )
  )
  
  cat(
    "Saved",
    nrow(phi_results),
    "segments\n"
  )
}

#########################################################
# get the estimated overdispersion parameter

phi_files <- list.files(
  out_dir_2,
  pattern = "_overdispersion\\.csv$",
  full.names = TRUE
)

all_phi <- rbindlist(
  lapply(
    phi_files,
    fread
  ),
  fill = TRUE
)
all_phi[, phi_T1 := as.numeric(phi_T1)]
all_phi[, phi_T2 := as.numeric(phi_T2)]

fwrite(
  all_phi,
  file.path(
    out_dir_2,
    "all_overdispersion.csv"
  )
)
########################################################
#get rho and stats related to rho
##############
##############
#t1
t1_ok <- all_phi[
  status_T1 == "OK" &
    !is.na(phi_T1)
]

t1_bad <- all_phi[
  status_T1 != "OK" &
    !is.na(phi_T1)
]

t1_ok[, rho := 1/(phi_T1 + 1)]
t1_bad[, rho := 1/(phi_T1 + 1)]
##############
##############
#t1
t2_ok <- all_phi[
  status_T2 == "OK" &
    !is.na(phi_T2)
]

t2_bad <- all_phi[
  status_T2 != "OK" &
    !is.na(phi_T2)
]

t2_ok[, rho := 1/(phi_T2 + 1)]
t2_bad[, rho := 1/(phi_T2 + 1)]

p <- c(0,.01,.05,.25,.5,.75,.95,.99,1)

data.table(
  Quantile = names(quantile(t1_ok$rho, p, na.rm=TRUE)),
  T1_OK = signif(quantile(t1_ok$rho, p, na.rm=TRUE), 3),
  T1_BAD = signif(quantile(t1_bad$rho, p, na.rm=TRUE), 3)
)
data.table(
  Quantile = names(quantile(t2_ok$rho, p, na.rm=TRUE)),
  T2_OK = signif(quantile(t2_ok$rho, p, na.rm=TRUE), 3),
  T2_BAD = signif(quantile(t2_bad$rho, p, na.rm=TRUE), 3)
)
##
rbind(
  data.table(
    group = "T1_OK",
    n = nrow(t1_ok),
    median_rho = median(t1_ok$rho, na.rm=TRUE),
    mean_rho = mean(t1_ok$rho, na.rm=TRUE)
  ),
  data.table(
    group = "T1_BAD",
    n = nrow(t1_bad),
    median_rho = median(t1_bad$rho, na.rm=TRUE),
    mean_rho = mean(t1_bad$rho, na.rm=TRUE)
  ),
  data.table(
    group = "T2_OK",
    n = nrow(t2_ok),
    median_rho = median(t2_ok$rho, na.rm=TRUE),
    mean_rho = mean(t2_ok$rho, na.rm=TRUE)
  ),
  data.table(
    group = "T2_BAD",
    n = nrow(t2_bad),
    median_rho = median(t2_bad$rho, na.rm=TRUE),
    mean_rho = mean(t2_bad$rho, na.rm=TRUE)
  )
)

p <- c(0,.01,.05,.25,.5,.75,.95,.99,1)

make_stats <- function(x, group){
  
  q <- round(
    quantile(x, probs = p, na.rm = TRUE),
    3
  )
  
  data.table(
    group = group,
    n = sum(!is.na(x)),
    mean = round(mean(x, na.rm = TRUE), 3),
    median = round(median(x, na.rm = TRUE), 3),
    sd = round(sd(x, na.rm = TRUE), 3),
    min = q[1],
    q01 = q[2],
    q05 = q[3],
    q25 = q[4],
    q50 = q[5],
    q75 = q[6],
    q95 = q[7],
    q99 = q[8],
    max = q[9]
  )
}

rho_phi_summary <- rbindlist(list(
  
  make_stats(t1_ok$phi_T1,  "PHI_T1_OK"),
  make_stats(t1_bad$phi_T1, "PHI_T1_BAD"),
  
  make_stats(t2_ok$phi_T2,  "PHI_T2_OK"),
  make_stats(t2_bad$phi_T2, "PHI_T2_BAD"),
  
  make_stats(t1_ok$rho,  "RHO_T1_OK"),
  make_stats(t1_bad$rho, "RHO_T1_BAD"),
  
  make_stats(t2_ok$rho,  "RHO_T2_OK"),
  make_stats(t2_bad$rho, "RHO_T2_BAD")
  
), fill = TRUE)

fwrite(
  rho_phi_summary,
  "/scratch/project_2010912/ezel/salmon/data_files/model_results/weight_rescaled_T1_T2_separated/overdispersion/rho_phi_summary.csv"
)

rho_phi_summary