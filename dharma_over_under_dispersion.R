library(DHARMa)

library(DHARMa)
model_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated_beta_binom"
dharm_dir <- file.path(model_dir, "DHARMa_results")
dir.create(dharm_dir, showWarnings = FALSE)
# Find all model files
model_files <- list.files(
  model_dir,
  pattern = "^model_.*\\.rds$",
  full.names = TRUE
)
# Store results here
summary_results <- list()
for (f in model_files) {
  file_name <- tools::file_path_sans_ext(basename(f))
  cat("\n========================================\n")
  cat("FILE:", file_name, "\n")
  cat("========================================\n")
  # Load the container
  x <- readRDS(f)
  # Output directory for this chromosome/contig
  file_dir <- file.path(dharm_dir, file_name)
  dir.create(file_dir, showWarnings = FALSE)
  # Loop over genomic regions
  for (region_name in names(x)) {
    cat("  Region:", region_name, "\n")
    region <- x[[region_name]]
    # Create directory for this region
    region_dir <- file.path(file_dir, region_name)
    dir.create(region_dir, showWarnings = FALSE)
    # Loop over T1 and T2
    for (timepoint in c("model_T1", "model_T2")) {
      cat("    ", timepoint, "\n")
      # Check that model exists
      if (!timepoint %in% names(region)) {
        cat("      Model missing -- skipping\n")
        next
      }
      mod <- region[[timepoint]]
      # Check that it is actually a glmmTMB model
      if (!inherits(mod, "glmmTMB")) {
        cat("      Not a glmmTMB model -- skipping\n")
        next
      }
      # Run DHARMa safely
      dharm <- tryCatch(
        {
          simulateResiduals(
            fittedModel = mod
          )
        },
        error = function(e) {
          cat("      ERROR:", conditionMessage(e), "\n")
          NULL
        }
      )
      # If DHARMa failed, record failure and continue
      if (is.null(dharm)) {
        summary_results[[length(summary_results) + 1]] <- data.frame(
          file = file_name,
          region = region_name,
          timepoint = timepoint,
          status = "DHARMa_FAILED",
          uniformity_p = NA_real_,
          dispersion_p = NA_real_,
          outlier_p = NA_real_
        )
        next
      }
      # Save DHARMa object
      saveRDS(
        dharm,
        file = file.path(
          region_dir,
          paste0(timepoint, "_DHARMa.rds")
        )
      )
      # Save plot
      pdf(
        file = file.path(
          region_dir,
          paste0(timepoint, "_DHARMa.pdf")
        ),
        width = 8,
        height = 8
      )
      plot(dharm)
      dev.off()
      # Run DHARMa tests
      uniformity <- tryCatch(
        testUniformity(dharm),
        error = function(e) NULL
      )
      dispersion <- tryCatch(
        testDispersion(dharm),
        error = function(e) NULL
      )
      outliers <- tryCatch(
        testOutliers(dharm),
        error = function(e) NULL
      )
      # Save summary
      summary_results[[length(summary_results) + 1]] <- data.frame(
        file = file_name,
        region = region_name,
        timepoint = timepoint,
        status = "OK",
        uniformity_p = if (!is.null(uniformity))
          uniformity$p.value else NA_real_,
        dispersion_p = if (!is.null(dispersion))
          dispersion$p.value else NA_real_,
        outlier_p = if (!is.null(outliers))
          outliers$p.value else NA_real_
      )
      # Free memory
      rm(mod, dharm)
      gc()
    }
  }
}
# Combine summary results
summary_results <- do.call(
  rbind,
  summary_results
)
# Save summary table
write.csv(
  summary_results,
  file = file.path(
    dharm_dir,
    "DHARMa_summary.csv"
  ),
  row.names = FALSE
)
cat("\nFinished!\n")

##############################
##############################
library(DHARMa)

dharm_dir <- "/scratch/project_2010912/ezel/salmon/data_files/model_results_T1_T2_separated_beta_binom/DHARMa_results"

dharm_files <- list.files(
  dharm_dir,
  pattern = "_DHARMa\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

dispersion_results <- list()

for (f in dharm_files) {
  
  dharm <- readRDS(f)
  
  timepoint <- sub("_DHARMa\\.rds$", "", basename(f))
  region <- basename(dirname(f))
  model_file <- basename(dirname(dirname(f)))
  
  result <- tryCatch(
    {
      test <- testDispersion(dharm)
      
      data.frame(
        model_file = model_file,
        region = region,
        timepoint = timepoint,
        p_value = test$p.value,
        statistic = unname(test$statistic),
        status = "tested",
        error = NA_character_
      )
    },
    
    error = function(e) {
      
      data.frame(
        model_file = model_file,
        region = region,
        timepoint = timepoint,
        p_value = NA_real_,
        statistic = NA_real_,
        status = "failed",
        error = conditionMessage(e)
      )
    }
  )
  
  dispersion_results[[length(dispersion_results) + 1]] <- result
}

dispersion_results <- do.call(
  rbind,
  dispersion_results
)

#get the problematic dispersions
dispersion_results$dispersion_status <- ifelse(
  dispersion_results$status == "failed",
  "failed",
  ifelse(
    dispersion_results$p_value < 0.05 &
      dispersion_results$statistic > 1,
    "overdispersed",
    ifelse(
      dispersion_results$p_value < 0.05 &
        dispersion_results$statistic < 1,
      "underdispersed",
      "OK"
    )
  )
)
overdispersed <- subset(
  dispersion_results,
  dispersion_status == "overdispersed"
)

underdispersed <- subset(
  dispersion_results,
  dispersion_status == "underdispersed"
)

###########get the plots for the over/underdispersed
# Recreate output directories
over_dir <- file.path(dharm_dir, "overdispersed")
under_dir <- file.path(dharm_dir, "underdispersed")

dir.create(over_dir, showWarnings = FALSE)
dir.create(under_dir, showWarnings = FALSE)


save_selected_plots <- function(results, status, output_dir) {
  
  selected <- results[
    results$dispersion_status == status &
      results$status == "tested",
    ,
    drop = FALSE
  ]
  
  cat("\n", status, ":", nrow(selected), "models\n")
  
  saved <- 0
  missing <- 0
  
  for (i in seq_len(nrow(selected))) {
    
    model_file <- selected$model_file[i]
    region <- selected$region[i]
    timepoint <- selected$timepoint[i]
    
    # Exact path to the saved DHARMa object
    f <- file.path(
      dharm_dir,
      model_file,
      region,
      paste0(timepoint, "_DHARMa.rds")
    )
    
    if (!file.exists(f)) {
      cat("MISSING:", f, "\n")
      missing <- missing + 1
      next
    }
    
    # Load DHARMa object
    dharm <- readRDS(f)
    
    # Unique PDF name
    pdf_file <- file.path(
      output_dir,
      paste0(
        model_file, "_",
        region, "_",
        timepoint, "_DHARMa.pdf"
      )
    )
    
    # Save plot
    pdf(
      pdf_file,
      width = 8,
      height = 8
    )
    
    plot(dharm)
    
    dev.off()
    
    saved <- saved + 1
    
    if (saved %% 25 == 0) {
      cat("  saved:", saved, "\n")
    }
  }
  
  cat(
    "\nFinished", status,
    "| saved:", saved,
    "| missing:", missing, "\n"
  )
}


# Overdispersed
save_selected_plots(
  dispersion_results,
  "overdispersed",
  over_dir
)

# Underdispersed
save_selected_plots(
  dispersion_results,
  "underdispersed",
  under_dir
)


#find the ones with the highest/lowest over/under dispersion
# Read the dispersion results
dispersion_results <- read.csv(
  file.path(dharm_dir, "dispersion_all_models.csv"),
  stringsAsFactors = FALSE
)

# Top 20 overdispersed
top20_over <- dispersion_results[
  dispersion_results$dispersion_status == "overdispersed",
]

top20_over <- top20_over[
  order(-top20_over$statistic),
]

top20_over <- head(top20_over, 20)


# Top 20 underdispersed
top20_under <- dispersion_results[
  dispersion_results$dispersion_status == "underdispersed",
]

top20_under <- top20_under[
  order(top20_under$statistic),
]

top20_under <- head(top20_under, 20)

write.csv(
  top20_over,
  file.path(dharm_dir, "top20_overdispersed_models.csv"),
  row.names = FALSE
)

write.csv(
  top20_under,
  file.path(dharm_dir, "top20_underdispersed_models.csv"),
  row.names = FALSE
)

#####################
#####################
#plotting dispersion status vs segment lengths
library(ggplot2)
library(dplyr)
library(stringr)

dispersion_results <- read.csv(
  file.path(dharm_dir, "dispersion_all_models.csv"),
  stringsAsFactors = FALSE
)

plot_data <- dispersion_results %>%
  filter(
    status == "tested",
    !is.na(statistic)
  ) %>%
  mutate(
    start = as.numeric(str_extract(region, "(?<=_)[0-9]+(?=_)")),
    end = as.numeric(str_extract(region, "(?<=_)[0-9]+$")),
    segment_length = end - start + 1
  )

p <- ggplot(
  plot_data,
  aes(
    x = segment_length,
    y = statistic,
    color = dispersion_status,
    alpha = dispersion_status
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  scale_x_log10() +
  scale_color_manual(
    values = c(
      "OK" = "grey40",
      "overdispersed" = "#D73027",
      "underdispersed" = "#4575B4"
    )
  ) +
  scale_alpha_manual(
    values = c(
      "OK" = 0.12,
      "overdispersed" = 0.9,
      "underdispersed" = 0.9
    ),
    guide = "none"
  ) +
  facet_wrap(
    ~ timepoint,
    ncol = 2
  ) +
  labs(
    x = "Segment length (bp)",
    y = "DHARMa dispersion statistic",
    color = "Dispersion status"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    panel.spacing = unit(1.5, "lines"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

ggsave(
  filename = file.path(dharm_dir, "segment_length_vs_dispersion_T1_T2.png"),
  plot = p,
  width = 12,
  height = 6,
  dpi = 300
)