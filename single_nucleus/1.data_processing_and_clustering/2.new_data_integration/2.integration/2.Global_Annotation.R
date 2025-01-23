library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

###Prepare Reference Data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_Atlas_V2_07-2023_Object_Downsampled_SCT_Reference_B.Rda")
KB

###Prepare Query Data 
nKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_Atlas_V2_04-2024-subset_Object.Rds")
nKB[["RNA"]]$data <- NormalizeData(nKB[["RNA"]]$counts)

###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KB,
  query = nKB,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

nKB <- MapQuery(
  anchorset = anchors,
  query = nKB,
  reference = KB,
  refdata = list(
    v2.subclass.l1 = "v2.subclass.l1",
    v2.subclass.l3 = "v2.subclass.l3",
    v2.clusters = "v2.clusters"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(nKB, reduction = "ref.umap", group.by = "predicted.v2.subclass.l1", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(nKB, reduction = "ref.umap", group.by = "predicted.v2.clusters", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()

hist(nKB$predicted.v2.subclass.l1.score)
hist(nKB$predicted.v2.subclass.l3.score)
hist(nKB$predicted.v2.clusters.score)

#remove prediction assays
nKB[["prediction.score.v2.subclass.l1"]] <- NULL
nKB[["prediction.score.v2.subclass.l3"]] <- NULL
nKB[["prediction.score.v2.clusters"]] <- NULL

saveRDS(nKB, file = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_Atlas_V2_04-2024-subset_Object_global-predictions.Rds")
