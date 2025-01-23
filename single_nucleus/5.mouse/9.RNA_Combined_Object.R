library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(BPCells)
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")


#Load and combine individual objects
KB.all <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_C.Rds")
KB.all
#70956 features across 330539 samples

#Load Individual objects
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
KB.str <- KB@meta.data

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_C.Rds")
KB.imm <- KB@meta.data

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_C.Rds")
KB.ec <- KB@meta.data

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")
KB.pt <- KB@meta.data

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")
KB.dt <- KB@meta.data

col.to.use <- intersect(colnames(KB.str), c(colnames(KB.imm), colnames(KB.ec), colnames(KB.pt), colnames(KB.dt)))
meta <- rbind(KB.str[,col.to.use], KB.imm[,col.to.use], KB.ec[,col.to.use], KB.pt[,col.to.use], KB.dt[,col.to.use])

KB.all <- subset(KB.all, cells = rownames(meta))
KB.all@meta.data <- meta[colnames(KB.all),]
DefaultAssay(KB.all) <- "RNA"
KB.all <- JoinLayers(object = KB.all)
KB.all[["sketch"]] <- NULL

# Write the matrix to a directory
write_matrix_dir(
  mat = KB.all[["RNA"]]$counts,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/kidney_counts_112024",
  overwrite = TRUE
)

kidney.mat <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/kidney_counts_112024")
KB.all[["RNA"]]$counts <- kidney.mat
KB.all <- NormalizeData(KB.all)
KB.all[["pca"]] <- NULL
KB.all[["integrated.rpca"]] <- KB.all[["integrated.rpca.full"]]
KB.all[["integrated.rpca.full"]] <- NULL
KB.all[["umap"]] <- KB.all[["umap.full"]]
KB.all[["umap.full"]] <- NULL

KB.all <- RunUMAP(KB.all, reduction = "integrated.rpca", dims = 1:50)
DimPlot(KB.all, reduction = "umap", pt.size = 0.5, label = TRUE, 
        group.by = "v2.clusters", repel = TRUE) + NoLegend()

KB.all
#An object of class Seurat 
#35478 features across 315764 samples within 1 assay 
#Active assay: RNA (35478 features, 2000 variable features)
#2 layers present: data, counts
#2 dimensional reductions calculated: integrated.rpca, umap


#update cluster metadata
meta <- KB.all@meta.data
cl.meta <- read.delim("mouse_IRI/Mouse_Kidney_AtlasV2_Cluster_Metadata_11152024.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.all@meta.data <- meta

DimPlot(KB.all, label = T, label.size = 3, reduction = "umap", 
        group.by = "v2.subclass.l1", alpha = 0.1) + NoLegend()
DimPlot(KB.all, label = T, label.size = 3, reduction = "umap", 
        group.by = "v2.subclass.l2", alpha = 0.1, repel = TRUE) + NoLegend()
DimPlot(KB.all, label = T, label.size = 3, reduction = "umap", 
        group.by = "v2.subclass.l3", alpha = 0.1, repel = TRUE) + NoLegend()

saveRDS(KB.all, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
