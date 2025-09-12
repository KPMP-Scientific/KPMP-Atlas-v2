# loading libraries -------------------------------------------------------

#for loading data
library(Seurat)
library(BPCells)

#for converting seurat objects to h5ad
library(hdf5r)
library(dplyr)


# helper functions --------------------------------------------------------

seurat_write_h5 <- function(seurat_obj = NULL, file = NULL, assays = NULL, reductions_map= list('UMAP'='UMAP')){
  
  if(is.null(file)){
    stop('No such file or directory')
  }
    
    seurat_to_h5(seurat_obj = seurat_obj, file = file, h5 = h5, assays = assays, reductions_map = reductions_map)
  
}


seurat_to_h5 <- function(seurat_obj, h5, file, assays, reductions_map){
  
  #save the relevant assay data

  lapply(assays, function(assay){

    slot <- 'counts'
    if(assay == 'SCT') slot <- 'data'

    mat <- GetAssayData(seurat_obj, slot = slot, assay = assay)

    write_matrix_anndata_hdf5(
    mat,
    file,
    group = paste('assays', assay, sep = '/')
    )
  })

  h5 <- h5file(file, mode = "a")

  tryCatch({
  
  # --- save the cell annotations
  df_to_h5(df = slot(seurat_obj, name = 'meta.data'), h5 = h5, gr_name = 'obs')
  
  # --- save gene information
  if(!is.null(rownames(seurat_obj))) {
    df_to_h5(df = data.frame(row.names = rownames(seurat_obj)), h5 = h5, gr_name = 'var')
  }
  
  reductions_to_h5(seurat_obj, h5, reductions_map)
  
  print('INFO: Successfully converted Seurat object to h5')
  print('INFO: final h5 file is:\n')
  print(h5)

  },
  error = function(e) print(e),
  finally = {
    h5close(h5)
  })
  
}


df_to_h5 <- function(df, h5, gr_name=NULL){
  if(!any(grepl(gr_name, names(h5)))){
    h5df <- h5$create_group(gr_name)
  }else{
    h5df <- h5[[gr_name]]
  }
  
  h5df[['index']] = rownames(df)
  if(ncol(df)>0){
    h5df[['colnames']] = colnames(df)
  }
  # factor to levels,character to levels,logical to levels
  for(k in names(df)){
    if(is.factor(df[[k]])){
      h5df[[k]]<- as.integer(df[[k]]) - 1L # for 0 begin
      h5df[[paste0(k,'_levels')]]<- levels(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'category'
    }
    if(is.character(df[[k]])){
      str_to_lvl <- factor(df[[k]])
      h5df[[k]]<- as.integer(str_to_lvl) - 1L
      h5df[[paste0(k,'_levels')]]<- levels(str_to_lvl)
      h5attr(h5df[[k]], 'origin_dtype') = 'string'
    }
    if(is.logical(df[[k]])){
      h5df[[k]] <- as.integer(df[[k]])
      h5attr(h5df[[k]], 'origin_dtype') = 'bool'
    }
    if(any(is.numeric(df[[k]]),is.integer(df[[k]]))){
      h5df[[k]] <- df[[k]]
      h5attr(h5df[[k]], 'origin_dtype') = 'number'
    }
  }
}

reductions_to_h5 <- function(seurat_obj, h5, reductions_map){
  
  h5red <- h5$create_group('reductions')
  
  lapply(names(reductions_map), function(reduction){
    if(reduction %in% names(slot(seurat_obj, name = 'reductions'))){
      
      mat <- slot(slot(seurat_obj, name = 'reductions')[[reduction]], 'cell.embeddings')
      
      h5mat <- h5red$create_group(reductions_map[[reduction]])
      
      for(k in 1:(length(ncol(mat)) + 1)){
        if(any(is.numeric(mat[,k]), is.integer(mat[,k]))){
          
          h5mat[[colnames(mat)[k]]] <- mat[,k]
          h5attr(h5mat[[colnames(mat)[k]]], 'origin_dtype') = 'number'
          
        }
      }
    }
    
  })
}


# define input and output paths -------------------------------------------
if(exists("snakemake")){
  seurat_fp <- snakemake@input$data
  data_dir <- snakemake@input$data_dir
  output_fp <- snakemake@output[[1]]
  
  assays <- snakemake@params[[1]]
  reductions  <- snakemake@params[[2]]
  
}else{
  seurat_fp <- 'sds/data/KPMPv2/snRNA/Kidney_AtlasV2_Seurat_05162024.rds'
  data_dir <- 'sds/data/KPMPv2/snRNA/full_kidney_count_set_0424'

  output_fp <- "test.h5ad"
  
  reductions  <- list( 'full.umap'= 'umap', 'pca.full'= 'pca')
  
  assays <- c('RNA', 'SCT')
}


# load data ---------------------------------------------------------------

if(grepl('.RDS$', seurat_fp,  ignore.case = TRUE)) seurat_obj <- readRDS(seurat_fp)


# reformat metadata --------------------------------------------------------
# library(readxl)
# metadata <- read_excel(clusters_fp, sheet = "Table S - sn annot", skip = 4) %>% select(-contains('...')) %>% 
#   mutate_at(vars(starts_with('v2.subclass')), ~stringr::str_replace_all(., '/', '_')) %>% # replace / with _ in subclass names
#   mutate_at(vars(starts_with('v2.subclass')), ~stringr::str_replace_all(., '\\+', ''))  %>%  # remove special chars from subclass names
#   mutate_at(vars(starts_with('v2.subclass')), ~stringr::str_replace_all(., 'ï', 'i'))

# #remove all columns from seurat_obj that are in metadata (except for v2.clusters)
# metadata_cols <- intersect(colnames(metadata), colnames(seurat_obj@meta.data))
# metadata_cols <- metadata_cols[!metadata_cols %in% 'v2.clusters']

# #remove metadata_cols from seurat_obj
# seurat_obj@meta.data <- seurat_obj@meta.data %>% select(-all_of(metadata_cols))

# #left_join metadata with seurat_obj on v2.clusters
# new_annotations <- left_join(seurat_obj@meta.data %>% tibble::rownames_to_column('cell_id'), metadata, by = c('v2.clusters' = 'v2.clusters')) %>% 
#   tibble::column_to_rownames('cell_id')

technical = c('percent.er', 'percent.mt', 'source', 'assay', 'percent_cortex', 'percent_medulla', 'region_level1', 
              'region_level2', 'location', 'laterality', 'protocol', 'tissue_type_full', 'tissue_type', 'atlas_version')
clinical = c('age_binned', 'sex', 'race', 'KDIGO_stage', 'baseline_eGFR_binned', 'proteinuria_binned', 'A1c_binned', 'albuminuria_binned',
             'diabetes_history', 'diabetes_duration', 'hypertension_history', 'hypertension_duration', 'on_RAAS_blockade')

condition = c('condition', 'condition_level1', 'condition_level2', 'condition_level3')
celltype = c('v2.clusters', 'v2.subclass.full', 'v2.subclass.sp', 'v2.subclass.l3', 'v2.subclass.l2', 'v2.subclass.l1', 'v2.state.l2', 'v2.state.l1', 'v2.class', 'v2.structure', 'integration_stable')

# #remove trailing whitespace from all character columns
# new_annotations <- new_annotations %>% mutate_if(is.character, ~stringr::str_trim(.))

# #select columns in technical, clinical, uid, condition, celltype
# new_annotations <- new_annotations %>% select(all_of(technical), all_of(clinical), all_of(uid), all_of(condition), all_of(celltype))

# #rename columns in technical, clinical, uid to start with tech_, clin_, uid_
# new_annotations <- new_annotations %>% rename_at(vars(technical), ~paste0('tech_', .)) %>% 
#   rename_at(vars(clinical), ~paste0('clin_', .)) %>% 
#   rename_at(vars(uid), ~paste0('uid_', .)) %>%
#   rename('uid_donor_id' = 'uid_patient')

#rename celltype columns


# new_annotations <- new_annotations %>% rename('condition_1' = 'condition_level1',
#                                             'condition_2' = 'condition_level2',
#                                             'condition_3' = 'condition_level3',
#                                             'condition_full' = 'condition',
#                                             'celltype_integration_stable' = 'integration_stable')


#replace NA with 'NA' in all character columns
uid = c('patient', 'specimen', 'library', 'v2.clusters')
seurat_obj@meta.data <- seurat_obj@meta.data %>% select(all_of(uid)) %>%
  mutate_if(is.character, ~replace(., is.na(.), 'NA')) %>%
  mutate_if(is.character, ~stringr::str_trim(.))

print('INFO: Added metadata to seurat object')


# add counts matrix -------------------------------------------------------

print('INFO: Adding counts matrix to seurat object')
# load count matrix from disk
counts.mat <- open_matrix_dir(dir = data_dir)
seurat_obj$RNA$counts <- counts.mat


# dev ---------------------------------------------------------------------

assays <- intersect(assays, Assays(seurat_obj))

if(length(assays) == 0) {
  stop('None of the assays provided can be found in the Seurat object')
}

seurat_write_h5(seurat_obj, output_fp, assays = assays, reductions_map = reductions)



