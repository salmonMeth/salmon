#we want to visualize the matching segments(3way matches) to see how much overlap there is,
# how they differ etc


#######looping over all the chromosomes
library(data.table)
library(ggplot2)


match_dir <-"/scratch/project_2010912/ezel/salmon/data_files/segmentation_3way_overlap"
relative_tolerance <- 0.2
matches_per_plot <- 20



##################################################
# Get all the chromosome match files

match_files <- list.files(
  match_dir,
  pattern = paste0(
    "_matches_relTol_",
    relative_tolerance,
    "\\.csv$"
  ),
  full.names = TRUE
)


cat(
  "Found",
  length(match_files),
  "chromosome files\n"
)



##################################################
# Loop =

for(match_file in match_files){
  
  chr_to_plot <- sub(
    paste0(
      "_matches_relTol_",
      relative_tolerance,
      "\\.csv$"
    ),
    "",
    basename(match_file)
  )
  
  
  cat(
    "Plotting",
    chr_to_plot,
    "\n"
  )
  
  matches <- fread(match_file)
  
  if(nrow(matches) == 0){
    
    cat(
      "Skipping",
      chr_to_plot,
      "- no matches\n"
    )
    
    next
    
  }
  #split the chromosomes into multiple groups
  #for ease of viewing
  matches[, plot_group :=
            ceiling(seq_len(.N) / matches_per_plot)]
  
  
  
  ##
  
  plot_dir <- file.path(
    match_dir,
    paste0(
      chr_to_plot,
      "_match_plots"
    )
  )
  
  
  dir.create(
    plot_dir,
    showWarnings = FALSE
  )
  
  
  
  ############################################
  # plot
  
  for(g in unique(matches$plot_group)){
    
    
    plot_matches <- matches[
      plot_group == g
    ]
    
    plot_segments <- rbindlist(list(
      
      data.table(
        sample = "Weighted",
        start = plot_matches$Weighted_start,
        end = plot_matches$Weighted_end,
        match_id = seq_len(nrow(plot_matches)),
        score = plot_matches$score
      ),
      
      data.table(
        sample = "T1",
        start = plot_matches$T1_start,
        end = plot_matches$T1_end,
        match_id = seq_len(nrow(plot_matches)),
        score = plot_matches$score
      ),
      
      data.table(
        sample = "T2",
        start = plot_matches$T2_start,
        end = plot_matches$T2_end,
        match_id = seq_len(nrow(plot_matches)),
        score = plot_matches$score
      )
      
    ))
    
    
    
    plot_segments[, sample := factor(
      sample,
      levels = c(
        "Weighted",
        "T1",
        "T2"
      )
    )]
    
    
    p <- ggplot(
      plot_segments,
      aes(
        x = start,
        xend = end,
        y = as.numeric(sample),
        yend = as.numeric(sample),
        color = score
      )
    ) +
      
      geom_segment(
        linewidth = 2.5
      ) +
      
      geom_segment(
        aes(
          x = start,
          xend = start,
          y = as.numeric(sample) - 0.25,
          yend = as.numeric(sample) + 0.25
        ),
        linewidth = 1
      ) +
      
      geom_segment(
        aes(
          x = end,
          xend = end,
          y = as.numeric(sample) - 0.25,
          yend = as.numeric(sample) + 0.25
        ),
        linewidth = 1
      ) +
      
      scale_y_continuous(
        breaks = 1:3,
        labels = c(
          "Weighted",
          "T1",
          "T2"
        ),
        limits = c(
          0.5,
          3.5
        )
      ) +
      
      scale_color_gradientn(
        colours = c(
          "#08306B",
          "#00CED1",
          "#7FFFD4",
          "#FF69B4",
          "#C71585"
        ),
        values = scales::rescale(
          c(
            0,
            0.02,
            0.08,
            0.15,
            0.2
          )
        ),
        limits = c(
          0,
          0.2
        ),
        oob = scales::squish,
        name = "Similarity score"
      ) +
      
      theme_classic() +
      
      theme(
        axis.text.y = element_text(
          size = 12
        ),
        plot.title = element_text(
          size = 14
        )
      ) +
      
      labs(
        title = paste(
          chr_to_plot,
          "ABC matches",
          g
        ),
        x = "Genomic position",
        y = ""
      )
    
    
    
    ggsave(
      filename = file.path(
        plot_dir,
        paste0(
          "matches_group_",
          g,
          ".png"
        )
      ),
      plot = p,
      width = 14,
      height = 6,
      dpi = 300
    )
    
  }
  
}

cat("Finished all chromosomes\n")