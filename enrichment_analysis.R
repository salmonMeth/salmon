#########################
########################
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
library(tidyr)

##############################
##############################
#some extra code
###########
##########we want to check the integrity of the annotations
extract_gene_ids <- function(files, out_dir = ".") {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  output_files <- character(length(files))
  
  for (i in seq_along(files)) {
    
    message("Processing: ", files[i])
    
    dat <- read.csv(
      files[i],
      stringsAsFactors = FALSE
    )
    
    if (!"gene_id" %in% names(dat)) {
      warning("Skipping ", files[i], ": no gene_id column.")
      next
    }
    
    # Split entries like "gene1;gene2"
    genes <- unlist(strsplit(dat$gene_id, ";", fixed = TRUE))
    
    genes <- sort(unique(
      genes[
        !is.na(genes) &
          genes != ""
      ]
    ))
    
    genes <- data.frame(
      gene_id = genes,
      stringsAsFactors = FALSE
    )
    
    out_file <- file.path(
      out_dir,
      paste0(
        tools::file_path_sans_ext(basename(files[i])),
        "_genes.csv"
      )
    )
    
    write.csv(
      genes,
      out_file,
      row.names = FALSE,
      quote = FALSE
    )
    
    output_files[i] <- out_file
  }
  
  invisible(output_files)
}
name1="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/methylkit_diff_meth_bases.csv"
name2="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/dss_diff_meth_bases.csv"
name3="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/methylkit_diff_meth_segs.csv"
names=c(name1,name2,name3)
out_dir="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes"

#################################
################################
################################
#now the actual enrichment analysis code


path_segs ="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_mk_segs.csv"
path_mk= "/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_mk_bases.csv"
path_dss="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/gene_ontology/pike_orthologs_dss_bases.csv"

mk_orth=read.csv(path_mk)
dss_orth=read.csv(path_dss)
seg_orth=read.csv(path_segs)
extract_gene_sets <- function(df, org) {
  
  original <- df %>%
    dplyr::select(
      ensembl_gene_id,
      external_gene_name
    ) %>%
    dplyr::filter(!is.na(ensembl_gene_id)) %>%
    dplyr::distinct()
  
  orthologue <- df %>%
    dplyr::select(
      homolog_ensembl_gene = all_of(paste0(org, "_homolog_ensembl_gene")),
      homolog_gene_name = all_of(paste0(org, "_homolog_associated_gene_name"))
    ) %>%
    dplyr::filter(!is.na(homolog_ensembl_gene)) %>%
    dplyr::distinct()
  
  list(
    original = original,
    orthologue = orthologue
  )
}
org="elucius"

mk_sets  <- extract_gene_sets(mk_orth, org)
dss_sets <- extract_gene_sets(dss_orth, org)
seg_sets <- extract_gene_sets(seg_orth, org)

library(gprofiler2)

salmon_GO <- gost(
  query = mk_sets$original$ensembl_gene_id,
  organism = "ssalar"
)

pike_GO <- gost(
  query = mk_sets$orthologue$homolog_ensembl_gene,
  organism = "elucius"
)

out_gosttable_pike="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/enrichment/pike_enrichment_table.png"
out_gosttable_salmon="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/enrichment/salmon_enrichment_table.png"


sal_tab <- publish_gosttable(
  salmon_GO,
  use_colors = TRUE,
  ggplot = TRUE
)

ggsave(
  filename = out_gosttable_salmon,
  plot = sal_tab,
  width = 10,
  height = 3,
  units = "in",
  dpi = 300
)

library(gprofiler2)
library(ggplot2)
library(dplyr)

save_enrichment_tables <- function(gene_sets, method_name,
                                   base_dir="/scratch/project_2010912/ezel/salmon/data_files/diff_meth_results_base_segments/annotated/genes/enrichment",
                                   salmon_org="ssalar", pike_org="elucius") {
  
  out_dir <- file.path(base_dir, method_name)
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  
  salmon_GO <- gost(gene_sets$original$ensembl_gene_id, organism=salmon_org, evcodes=TRUE)
  pike_GO <- gost(gene_sets$orthologue$homolog_ensembl_gene, organism=pike_org, evcodes=TRUE)
  
  save_plot <- function(go, file, w, h) {
    if(!is.null(go$result)) {
      p <- publish_gosttable(go, use_colors=TRUE, ggplot=TRUE) +
        ggplot2::theme(text=ggplot2::element_text(size=12))
      ggsave(file, p, width=w, height=h, units="in", dpi=300)
    }
  }
  
  save_plot(salmon_GO, file.path(out_dir,"salmon_enrichment_table.png"), 10, 4)
  save_plot(pike_GO, file.path(out_dir,"pike_enrichment_table.png"), 10, 6)
  
  save_results <- function(go, map, id_col, name_col, file) {
    if(!is.null(go$result)) {
      res <- go$result
      lookup <- setNames(map[[name_col]], map[[id_col]])
      
      res$intersection_gene_names <- sapply(res$intersection, function(x) {
        genes <- unlist(strsplit(x, ","))
        names <- lookup[genes]
        names[is.na(names) | names==""] <- genes[is.na(names) | names==""]
        paste(names, collapse=";")
      })
      
      res <- res %>%
        dplyr::mutate(dplyr::across(where(is.list), ~sapply(., paste, collapse=";")))
      
      write.csv(res, file, row.names=FALSE)
    }
  }
  
  salmon_map <- gene_sets$original %>%
    dplyr::select(ensembl_gene_id, external_gene_name)
  
  pike_map <- gene_sets$orthologue %>%
    dplyr::select(homolog_ensembl_gene, homolog_gene_name)
  
  save_results(salmon_GO, salmon_map,
               "ensembl_gene_id", "external_gene_name",
               file.path(out_dir,"salmon_enrichment_results.csv"))
  
  save_results(pike_GO, pike_map,
               "homolog_ensembl_gene", "homolog_gene_name",
               file.path(out_dir,"pike_enrichment_results.csv"))
  
  invisible(list(salmon=salmon_GO, pike=pike_GO, output_dir=out_dir))
}


save_enrichment_tables(mk_sets, "MethylKit")
save_enrichment_tables(dss_sets, "DSS")
save_enrichment_tables(seg_sets, "Segments")


################################
################################
