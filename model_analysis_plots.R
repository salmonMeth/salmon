library(data.table)
library(ggplot2)

path_out <- "/scratch/project_2010912/ezel/salmon/data_files/model_analysis"

qc_files <- list.files(
  path_out,
  pattern="^QC_table_.*\\.csv$",
  recursive=TRUE,
  full.names=TRUE
)

qc_all <- rbindlist(
  lapply(qc_files,fread),
  fill=TRUE
)

qc_all <- qc_all[
  DHARMa_dispersion_class %in% c(
    "OK",
    "OVERDISPERSION",
    "UNDERDISPERSION"
  )
]

qc_all[,.(n=.N,
          mean_ID_ma=mean(sd_ID_ma,na.rm=TRUE),
          median_ID_ma=median(sd_ID_ma,na.rm=TRUE),
          mean_ID_pa=mean(sd_ID_pa,na.rm=TRUE),
          median_ID_pa=median(sd_ID_pa,na.rm=TRUE)),
       by=DHARMa_dispersion_class]


ggsave(
  file.path(path_out,"maternal_RE_by_dispersion.png"),
  p_ma,
  width=7,
  height=5,
  dpi=300
)

qc_all[,.(n=.N,
          frac_ID_ma_zero=mean(sd_ID_ma<1e-4,na.rm=TRUE),
          frac_ID_pa_zero=mean(sd_ID_pa<1e-4,na.rm=TRUE)),
       by=DHARMa_dispersion_class]

qc_all[,.(median_ma=median(sd_ID_ma,na.rm=TRUE),
          median_pa=median(sd_ID_pa,na.rm=TRUE),
          mean_ma=mean(sd_ID_ma,na.rm=TRUE),
          mean_pa=mean(sd_ID_pa,na.rm=TRUE)),
       by=DHARMa_dispersion_class]

qc_all[,segment_length:=end-start]

qc_all[,.(n=.N,
          median_length=median(segment_length,na.rm=TRUE),
          mean_length=mean(segment_length,na.rm=TRUE),
          min_length=min(segment_length,na.rm=TRUE),
          max_length=max(segment_length,na.rm=TRUE)),
       by=DHARMa_dispersion_class]



library(data.table)
library(ggplot2)

# segment lengths by dispersion class
len_all <- qc_all[
  is.finite(start) & is.finite(end),
  .(
    Dispersion = factor(
      DHARMa_dispersion_class,
      levels=c("UNDERDISPERSION","OK","OVERDISPERSION"),
      labels=c("Under","OK","Over")
    ),
    Length = end - start
  )
]


# save raw segment lengths
fwrite(
  len_all,
  file.path(path_out,"segment_lengths_by_dispersion.csv")
)


# summary statistics
len_stats <- len_all[
  ,
  .(
    n=.N,
    mean_length=mean(Length, na.rm=TRUE),
    median_length=median(Length, na.rm=TRUE),
    sd_length=sd(Length, na.rm=TRUE),
    min_length=min(Length, na.rm=TRUE),
    max_length=max(Length, na.rm=TRUE),
    q25=quantile(Length,0.25,na.rm=TRUE),
    q75=quantile(Length,0.75,na.rm=TRUE)
  ),
  by=Dispersion
]


fwrite(
  len_stats,
  file.path(path_out,"segment_length_stats_by_dispersion.csv")
)


# plot (no filtering)
p_len <- ggplot(
  len_all,
  aes(
    x=Dispersion,
    y=Length,
    fill=Dispersion
  )
) +
  geom_boxplot(
    outlier.shape=NA
  ) +
  scale_y_log10() +
  scale_fill_manual(
    values=c(
      "Under"="#F8766D",
      "OK"="#00BA38",
      "Over"="#619CFF"
    )
  ) +
  theme_gray() +
  labs(
    x="Dispersion class",
    y="Segment length (bp)",
    fill="Dispersion"
  )


ggsave(
  file.path(path_out,"segment_length_by_dispersion.png"),
  p_len,
  width=7,
  height=5,
  dpi=300
)


###sd ma pa
library(data.table)
library(ggplot2)

# Function to prepare SD data
make_sd_plot_data <- function(dt, col) {
  x <- dt[
    is.finite(get(col)),
    .(
      Dispersion = factor(
        DHARMa_dispersion_class,
        levels=c("UNDERDISPERSION","OK","OVERDISPERSION"),
        labels=c("Under","OK","Over")
      ),
      SD = get(col)
    )
  ]
  
  # remove boxplot outliers within each class
  x_clean <- x[
    ,
    {
      q1 <- quantile(SD, 0.25)
      q3 <- quantile(SD, 0.75)
      iqr <- q3 - q1
      
      .SD[SD >= q1 - 1.5*iqr &
            SD <= q3 + 1.5*iqr]
    },
    by=Dispersion
  ]
  
  return(x_clean)
}


# Prepare maternal and paternal
ma_clean <- make_sd_plot_data(qc_all, "sd_ID_ma")
pa_clean <- make_sd_plot_data(qc_all, "sd_ID_pa")


# use same y limits for comparison
all_sd <- c(ma_clean$SD, pa_clean$SD)
ylims <- quantile(all_sd, c(0.01, 0.99), na.rm=TRUE)


# plotting function
make_boxplot <- function(df, ylabel) {
  ggplot(
    df,
    aes(
      x=Dispersion,
      y=SD,
      fill=Dispersion
    )
  ) +
    geom_boxplot(outlier.shape=NA) +
    scale_y_log10() +
    coord_cartesian(
      ylim=ylims
    ) +
    scale_fill_manual(
      values=c(
        "Under"="#F8766D",
        "OK"="#00BA38",
        "Over"="#619CFF"
      )
    ) +
    theme_gray() +
    labs(
      x="Dispersion class",
      y=ylabel,
      fill="Dispersion"
    )
}


p_ma <- make_boxplot(
  ma_clean,
  "Maternal random effect SD"
)

p_pa <- make_boxplot(
  pa_clean,
  "Paternal random effect SD"
)


# save
ggsave(
  file.path(path_out,"maternal_SD_by_dispersion_clean.png"),
  p_ma,
  width=7,
  height=5,
  dpi=300
)

ggsave(
  file.path(path_out,"paternal_SD_by_dispersion_clean.png"),
  p_pa,
  width=7,
  height=5,
  dpi=300
)

#########
# create segment length and methylation bins
seg_meth <- qc_all[
  is.finite(start) &
    is.finite(end) &
    is.finite(methylation_level),
  .(
    Length = end - start,
    Methylation = methylation_level
  )
]


# split methylation into 3 equal-sized groups
seg_meth[
  ,
  Meth_group := cut(
    Methylation,
    breaks = quantile(
      Methylation,
      probs=c(0, 1/3, 2/3, 1),
      na.rm=TRUE
    ),
    include.lowest=TRUE,
    labels=c("Low","Medium","High")
  )
]


# save the assigned groups
fwrite(
  seg_meth,
  file.path(path_out,"segment_length_by_methylation_group.csv")
)


# summary statistics
seg_meth_stats <- seg_meth[
  ,
  .(
    n=.N,
    mean_length=mean(Length),
    median_length=median(Length),
    sd_length=sd(Length),
    min_length=min(Length),
    max_length=max(Length)
  ),
  by=Meth_group
]


fwrite(
  seg_meth_stats,
  file.path(path_out,"segment_length_stats_by_methylation_group.csv")
)


# boxplot
p_len_meth <- ggplot(
  seg_meth,
  aes(
    x=Meth_group,
    y=Length,
    fill=Meth_group
  )
) +
  geom_boxplot(
    outlier.shape=NA
  ) +
  scale_y_log10() +
  scale_fill_manual(
    values=c(
      "Low"="#F8766D",
      "Medium"="#00BA38",
      "High"="#619CFF"
    )
  ) +
  theme_gray() +
  labs(
    x="Methylation level group",
    y="Segment length (bp)",
    fill="Methylation"
  )


ggsave(
  file.path(path_out,"segment_length_by_methylation_group.png"),
  p_len_meth,
  width=7,
  height=5,
  dpi=300
)

cov_dir <- "/scratch/project_2010912/ezel/salmon/data_files/segmentation/segment_cs_cov_3_samples_removed"



p_med <- ggplot(
  cov_disp,
  aes(x=Dispersion,y=median_coverage,fill=Dispersion)
) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_log10() +
  scale_fill_manual(values=c(
    "Under"="#F8766D",
    "OK"="#00BA38",
    "Over"="#619CFF"
  )) +
  theme_gray() +
  labs(
    x="Dispersion class",
    y="Median segment coverage",
    fill="Dispersion"
  )

ggsave(
  file.path(path_out,"segment_median_coverage_by_dispersion.png"),
  p_med,
  width=7,
  height=5,
  dpi=300
)



cov_files <- list.files(
  cov_dir,
  pattern="*_segment_coverage\\.csv$",
  full.names=TRUE
)

cov_all <- rbindlist(lapply(cov_files,fread),fill=TRUE)

cov_cols <- setdiff(colnames(cov_all),c("chr","start","end"))

cov_all[,`:=`(
  mean_coverage=rowMeans(.SD,na.rm=TRUE),
  median_coverage=apply(.SD,1,median,na.rm=TRUE),
  q25_coverage=apply(.SD,1,quantile,probs=0.25,na.rm=TRUE),
  q75_coverage=apply(.SD,1,quantile,probs=0.75,na.rm=TRUE),
  sd_coverage=apply(.SD,1,sd,na.rm=TRUE),
  min_coverage=apply(.SD,1,min,na.rm=TRUE),
  max_coverage=apply(.SD,1,max,na.rm=TRUE)
),.SDcols=cov_cols]

cov_all[,cv_coverage:=sd_coverage/mean_coverage]

cov_seg <- cov_all[,.(chr,start,end,
                      mean_coverage,
                      median_coverage,
                      q25_coverage,
                      q75_coverage,
                      sd_coverage,
                      min_coverage,
                      max_coverage,
                      cv_coverage)]

cov_disp <- merge(
  cov_seg,
  qc_all[,.(chr,start,end,DHARMa_dispersion_class)],
  by=c("chr","start","end")
)

cov_disp[,Dispersion:=factor(
  DHARMa_dispersion_class,
  levels=c("UNDERDISPERSION","OK","OVERDISPERSION"),
  labels=c("Under","OK","Over")
)]

fwrite(
  cov_disp,
  file.path(path_out,"segment_coverage_by_dispersion_full.csv")
)

cov_stats <- cov_disp[,.(n=.N,
                         mean_coverage=mean(mean_coverage,na.rm=TRUE),
                         median_coverage=median(median_coverage,na.rm=TRUE),
                         mean_q25=mean(q25_coverage,na.rm=TRUE),
                         mean_q75=mean(q75_coverage,na.rm=TRUE),
                         mean_sd=mean(sd_coverage,na.rm=TRUE),
                         mean_cv=mean(cv_coverage,na.rm=TRUE)),
                      by=Dispersion]

fwrite(
  cov_stats,
  file.path(path_out,"segment_coverage_stats_by_dispersion.csv")
)

p_cov <- ggplot(
  cov_disp,
  aes(x=Dispersion,y=median_coverage,fill=Dispersion)
) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_log10() +
  scale_fill_manual(values=c(
    "Under"="#F8766D",
    "OK"="#00BA38",
    "Over"="#619CFF"
  )) +
  theme_gray() +
  labs(
    x="Dispersion class",
    y="Median segment coverage",
    fill="Dispersion"
  )

ggsave(
  file.path(path_out,"segment_median_coverage_by_dispersion.png"),
  p_cov,
  width=7,
  height=5,
  dpi=300
)

p_cv <- ggplot(
  cov_disp,
  aes(x=Dispersion,y=cv_coverage,fill=Dispersion)
) +
  geom_boxplot(outlier.shape=NA) +
  scale_y_log10() +
  scale_fill_manual(values=c(
    "Under"="#F8766D",
    "OK"="#00BA38",
    "Over"="#619CFF"
  )) +
  theme_gray() +
  labs(
    x="Dispersion class",
    y="Coverage CV",
    fill="Dispersion"
  )

ggsave(
  file.path(path_out,"segment_coverage_CV_by_dispersion.png"),
  p_cv,
  width=7,
  height=5,
  dpi=300
)