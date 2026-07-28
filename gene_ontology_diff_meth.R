#############gene ontology for the diff methylated bases/segments
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(methylKit)
library(genomation)
library(GenomeInfoDb)
library(dplyr)
library(gprofiler2)
library(biomaRt)

############

path_dml_dss="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/dss_diff_meth_bases_genes.csv"
path_dml_mk="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/methylkit_diff_meth_bases_genes.csv"
path_dmr_mk="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/methylkit_diff_meth_segs_genes.csv"

methylkit_dml <- read.csv(path_dml_mk)
dss_dml <- read.csv(path_dml_dss)
dmr_seg <- read.csv(path_dmr_mk)


#############################################
#############################################
ensembl <- useEnsembl(
  biomart = "genes",
  version = NULL
)

####to find the datasets we need we can grep them as such
datasets[grep("lucius|esox|pike", datasets$dataset, ignore.case = TRUE), ]
# dataset : elucius_gene_ensembl 
salmon <- useEnsembl(
  biomart = "genes",
  dataset = "ssalar_gene_ensembl"
)
zebrafish <- useEnsembl(
  biomart = "genes",
  dataset = "drerio_gene_ensembl"
)

pike <- useEnsembl(
  biomart = "genes",
  dataset = "elucius_gene_ensembl"
)

###we want to find  homologs
attrs <- listAttributes(salmon)

#zebrafish
attrs[grep("drerio", attrs$name), ]


#pike
attrs[grep("elucius", attrs$name), ]


#organism=the ontology we want to use, zebrafish, pike etc
get_orthologs <- function(genes, mart, organism) {
  
  attrs <- c(
    "ensembl_gene_id",
    "external_gene_name",
    paste0(organism, "_homolog_ensembl_gene"),
    paste0(organism, "_homolog_associated_gene_name"),
    paste0(organism, "_homolog_orthology_type"),
    paste0(organism, "_homolog_orthology_confidence")
  )
  
  res <- getBM(
    attributes = attrs,
    filters = "ensembl_gene_id",
    values = unique(genes),
    mart = mart
  )
  
  ortholog_col <- paste0(organism, "_homolog_ensembl_gene")
  confidence_col <- paste0(organism, "_homolog_orthology_confidence")
  ## remove genes with no ortholog
  res <- res[res[[ortholog_col]] != "", ]
  
  ## keep only high-confidence orthologs
  res <- res[res[[confidence_col]] == 1, ]
  
  res
}


get_gene_list <- function(x){
  
  genes <- unique(
    na.omit(x$gene_id)
  )
  
  # remove multiple genes if separated by ;
  genes <- unlist(
    strsplit(genes, ";")
  )
  
  unique(genes)
}


genes_methylkit <- get_gene_list(methylkit_dml)

genes_dss <- get_gene_list(dss_dml)

genes_dmr <- get_gene_list(dmr_seg)

org="elucius"
mk_zfish  <- get_orthologs(genes_methylkit, salmon, org)
dss_zfish <- get_orthologs(genes_dss, salmon, org)
dmr_zfish <- get_orthologs(genes_dmr, salmon,org)

# Add method labels
mk_zfish$method  <- "MethylKit"
dss_zfish$method <- "DSS"
dmr_zfish$method <- "Segment"

write.csv(
  mk_zfish,
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_mk_bases.csv",
  row.names = FALSE
)
write.csv(
  dss_zfish,
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_dss_bases.csv",
  row.names = FALSE
)
write.csv(
  dmr_zfish,
  "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_mk_segs.csv",
  row.names = FALSE
)


#############################
##
#####################################
#WORD CLOUD V2
path_orth_dss="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_dss_bases.csv"
path_orth_mk="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_mk_bases.csv"
path_orth_mk_segs="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_mk_segs.csv"
path_word_cloud="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orth_word_cloud.png"



#we can modify this function to plot the gene names of salmon 
#instead of the pike orthologues
#by changing the line
# ortholog_col <- paste0(org,"_homolog_associated_gene_name") with 
#ortholog_col <- "external_gene_name"
make_methylation_triangle_plot <- function(
    mk_file,
    dss_file,
    seg_file,
    org = "elucius",
    save_path = "methylation_triangle_plot.png",
    title = "Genes associated with differential methylation"
) {
  library(grid)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  
  #####################################
  # Read ortholog files
  #####################################
  
  mk_orth <- read.csv(
    mk_file,
    stringsAsFactors = FALSE
  )
  
  dss_orth <- read.csv(
    dss_file,
    stringsAsFactors = FALSE
  )
  
  seg_orth <- read.csv(
    seg_file,
    stringsAsFactors = FALSE
  )
  
  
  #####################################
  # Extract ortholog gene names
  #####################################
  
  #ortholog_col <- paste0(org,"_homolog_associated_gene_name")
  ortholog_col <- "external_gene_name"
  
  for(df in list(mk_orth, dss_orth, seg_orth)){
    if(!ortholog_col %in% colnames(df)){
      stop(
        paste(
          "Column",
          ortholog_col,
          "not found"
        )
      )
    }
  }
  
  
  mk <- mk_orth %>%
    filter(
      !is.na(.data[[ortholog_col]]),
      .data[[ortholog_col]] != ""
    ) %>%
    pull(all_of(ortholog_col)) %>%
    unique()
  
  
  dss <- dss_orth %>%
    filter(
      !is.na(.data[[ortholog_col]]),
      .data[[ortholog_col]] != ""
    ) %>%
    pull(all_of(ortholog_col)) %>%
    unique()
  
  
  seg <- seg_orth %>% 
    filter(
      !is.na(.data[[ortholog_col]]),
      .data[[ortholog_col]] != ""
    ) %>%
    pull(all_of(ortholog_col)) %>%
    unique()
  
  
  
  #####################################
  # Create gene table
  #####################################
  
  gene_table <- data.frame(
    gene = unique(c(mk, dss, seg)),
    stringsAsFactors = FALSE
  )
  
  
  gene_table$MethylKit <- gene_table$gene %in% mk
  gene_table$DSS        <- gene_table$gene %in% dss
  gene_table$Segment    <- gene_table$gene %in% seg
  
  
  gene_table$group <- apply(
    gene_table[,c(
      "MethylKit",
      "DSS",
      "Segment"
    )],
    1,
    function(x){
      
      paste(
        names(x)[x],
        collapse = "+"
      )
      
    }
  )
  
  
  gene_table$group <- recode(
    gene_table$group,
    
    "MethylKit" =
      "MethylKit only",
    
    "DSS" =
      "DSS only",
    
    "Segment" =
      "Segment only",
    
    "MethylKit+DSS" =
      "MethylKit + DSS",
    
    "MethylKit+Segment" =
      "MethylKit + Segment",
    
    "DSS+Segment" =
      "DSS + Segment",
    
    "MethylKit+DSS+Segment" =
      "All three"
  )
  
  
  # number of methods
  gene_table$n_methods <-
    rowSums(
      gene_table[,c(
        "MethylKit",
        "DSS",
        "Segment"
      )]
    )
  
  
  #####################################
  # Font size and transparency
  #####################################
  
  gene_table <- gene_table %>%
    mutate(
      
      plot_size = case_when(
        
        n_methods == 1 ~ 3.5,
        
        n_methods == 2 ~ 6,
        
        n_methods == 3 ~ 10
        
      ),
      
      plot_alpha = case_when(
        
        n_methods == 1 ~ 1,
        
        n_methods == 2 ~ 0.95,
        
        n_methods == 3 ~ 0.75
        
      )
      
    )
  #####################################
  # Split groups
  #####################################
  
  mk_only <- gene_table %>%
    filter(group == "MethylKit only")
  
  dss_only <- gene_table %>%
    filter(group == "DSS only")
  
  seg_only <- gene_table %>%
    filter(group == "Segment only")
  
  mk_dss <- gene_table %>%
    filter(group == "MethylKit + DSS")
  
  mk_seg <- gene_table %>%
    filter(group == "MethylKit + Segment")
  
  dss_seg <- gene_table %>%
    filter(group == "DSS + Segment")
  
  all_three <- gene_table %>%
    filter(group == "All three")
  
  
  set.seed(123)
  
  
  place_group <- function(df, xmin, xmax, ymin, ymax){
    
    if(nrow(df)>0){
      
      df$x <- runif(
        nrow(df),
        xmin,
        xmax
      )
      
      df$y <- runif(
        nrow(df),
        ymin,
        ymax
      )
      
    }
    
    df
  }
  
  
  mk_only <- place_group(
    mk_only,
    -1.5,-0.3,
    -1,0
  )
  
  dss_only <- place_group(
    dss_only,
    -0.5,0.5,
    0.5,1.4
  )
  
  seg_only <- place_group(
    seg_only,
    0.3,1.5,
    -1,0
  )
  
  
  # intersections closer to centre
  
  mk_dss <- place_group(
    mk_dss,
    -0.5,0,
    0.2,0.7
  )
  
  mk_seg <- place_group(
    mk_seg,
    -0.4,0,
    -0.6,-0.1
  )
  
  dss_seg <- place_group(
    dss_seg,
    0,0.5,
    -0.6,-0.1
  )
  
  all_three <- place_group(
    all_three,
    -0.15,0.15,
    -0.35,0.1
  )
  
  
  #####################################
  plot_df <- bind_rows(
    mk_only,
    dss_only,
    seg_only,
    mk_dss,
    mk_seg,
    dss_seg,
    all_three
  )
  
  
  plot_df$plot_size <- case_when(
    plot_df$n_methods == 1 ~ 3.5,
    plot_df$n_methods == 2 ~ 6,
    plot_df$n_methods == 3 ~ 10
  )
  
  
  plot_df$plot_alpha <- case_when(
    plot_df$n_methods == 1 ~ 1,
    plot_df$n_methods == 2 ~ 0.95,
    plot_df$n_methods == 3 ~ 0.75
  )
  
  cat(
    "Total unique genes plotted:",
    nrow(plot_df),
    "\n"
  )
  
  cat(
    "Genes missing coordinates:",
    sum(is.na(plot_df$x) | is.na(plot_df$y)),
    "\n"
  )
  #####################################
  # Plot
  #####################################
  
  p <- ggplot(
    plot_df,
    aes(
      x = x,
      y = y,
      label = gene,
      color = group
    )
  ) +
    
    geom_text_repel(
      aes(
        size = plot_size,
        alpha = plot_alpha,
        color = group
      ),
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 0.05,
      point.padding = 0.02,
      force = 0.15,
      min.segment.length = Inf,
      seed = 123
    ) +
    
    scale_color_manual(
      values=c(
        "MethylKit only" = "#A020F0",
        "DSS only" = "#FF7F00",
        "Segment only" = "#0072B2",
        "MethylKit + DSS" = "#009E73",
        "MethylKit + Segment" = "#CC79A7",
        "DSS + Segment" = "#F0E442",
        "All three" = "#444444"
      ),
      labels=c(
        "MethylKit only" = "MK",
        "DSS only" = "DSS",
        "Segment only" = "DMR",
        "MethylKit + DSS" = "MK + DSS",
        "MethylKit + Segment" = "MK + DMR",
        "DSS + Segment" = "DSS + DMR",
        "All three" = "All methods"
      )
    )+
    
    scale_size_identity() +
    scale_alpha_identity()+
    
    guides(
      color = guide_legend(
        title = "Overlap group",
        nrow = 2,
        byrow = TRUE,
        override.aes = list(
          size = 5,
          alpha = 1
        )
      ),
      size = "none",
      alpha = "none"
    ) +
    
    theme_void() +
    
    theme(
      legend.position = "bottom",
      
      legend.text = element_text(
        size = 10
      ),
      
      legend.title = element_text(
        size = 12,
        face = "bold"
      ),
      
      legend.key.width = unit(
        1.2,
        "cm"
      ),
      
      legend.spacing.x = unit(
        0.3,
        "cm"
      ),
      
      plot.title = element_text(
        hjust = 0.5,
        size = 18,
        face = "bold"
      )
    ) +
    
    labs(
      title = title
    )
  
  ggsave(
    filename = save_path,
    plot = p,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  return(
    list(
      plot = p,
      gene_table = gene_table
    )
  )
}


result <- make_methylation_triangle_plot(path_orth_mk, path_orth_dss,path_orth_mk_segs,
                                         org = "elucius",path_word_cloud
)

result$plot

