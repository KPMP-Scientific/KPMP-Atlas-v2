library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)
library(tidyverse)
library(patchwork)
library(pagoda2)
library(BPCells)
options(future.globals.maxSize = 1e10)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")



###Global Object
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/v2_object_filtered/Kidney_Atlas_V2_07-2023_Object_Filtered_B.Rds")

#update cluster metadata
meta <- KB@meta.data
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB@meta.data <- meta

Idents(KB)
KB <- subset(KB, downsample = 1000)
KB[["RNA"]]$counts <- as(KB[["RNA"]]$counts, "dgCMatrix")

countMatrix <- KB[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 4983, maxit=1000)
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:4983]

KB <- SCTransform(KB, verbose = TRUE)
VariableFeatures(KB) <- sn.od.genes
KB <- RunPCA(KB, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("pca_", 1:50, sep = "")
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KB, reduction = "pca")),]
colnames(cell.embeddings) <- paste("pca_", 1:50, sep = "")
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 30L,
              min.dist = 0.3, umap.method = "uwot-learn")
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, repel = TRUE) + NoLegend()
KB <- FindNeighbors(
  object = KB,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])

KB[["sketch"]] <- NULL
KB[["v1.pca"]] <- NULL
KB[["v1.umap"]] <- NULL
KB[["pca.full"]] <- NULL
KB[["full.umap"]] <- NULL

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_Atlas_V2_07-2023_Object_Downsampled_SCT_Reference_B.Rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_Atlas_V2_07-2023_Object_Downsampled_SCT_Reference_B.Rda")




###Stromal object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k200.rda")
sn.od.genes <- VariableFeatures(KB.STR)
table(Idents(KB.STR))

#update cluster metadata
meta <- KB.STR@meta.data
meta$str.v2.clusters <- paste0("S_",meta$str.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$str.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.STR@meta.data <- meta

Idents(KB.STR)
KB.STR <- subset(KB.STR, downsample = 1000)
KB <- SCTransform(KB.STR, verbose = TRUE)
VariableFeatures(KB) <- sn.od.genes
KB <- RunPCA(KB, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("pca_", 1:50, sep = "")
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KB, reduction = "pca")),]
colnames(cell.embeddings) <- paste("pca_", 1:50, sep = "")
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 30L,
              min.dist = 0.3, umap.method = 'uwot', return.model = TRUE)
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, repel = TRUE) + NoLegend()
KB <- FindNeighbors(
  object = KB,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_SCT_Reference_full.Rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_SCT_Reference_full.Rda")




###Immune object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k200.rda")
sn.od.genes <- VariableFeatures(KB.IMM)
table(Idents(KB.IMM))

#update cluster metadata
meta <- KB.IMM@meta.data
meta$imm.v2.clusters <- paste0("I_",meta$imm.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$imm.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.IMM@meta.data <- meta

Idents(KB.IMM)
KB.IMM <- subset(KB.IMM, downsample = 1000)
KB <- SCTransform(KB.IMM, verbose = TRUE)
VariableFeatures(KB) <- sn.od.genes
KB <- RunPCA(KB, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("pca_", 1:50, sep = "")
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KB, reduction = "pca")),]
colnames(cell.embeddings) <- paste("pca_", 1:50, sep = "")
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 30L,
              min.dist = 0.3, umap.method = 'uwot', return.model = TRUE)
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, repel = TRUE) + NoLegend()
KB <- FindNeighbors(
  object = KB,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_SCT_Reference_full.Rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_SCT_Reference_full.Rda")



###Vascular object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_vasculature_k200.rda")
sn.od.genes <- VariableFeatures(KB.EC)
table(Idents(KB.EC))

#update cluster metadata
meta <- KB.EC@meta.data
meta$ec.v2.clusters <- paste0("E_",meta$ec.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$ec.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.EC@meta.data <- meta

Idents(KB.EC)
KB.EC <- subset(KB.EC, downsample = 1000)
KB <- SCTransform(KB.EC, verbose = TRUE)
VariableFeatures(KB) <- sn.od.genes
KB <- RunPCA(KB, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("pca_", 1:50, sep = "")
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KB, reduction = "pca")),]
colnames(cell.embeddings) <- paste("pca_", 1:50, sep = "")
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 30L,
              min.dist = 0.3, umap.method = 'uwot', return.model = TRUE)
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, repel = TRUE) + NoLegend()
KB <- FindNeighbors(
  object = KB,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_SCT_Reference.Rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_SCT_Reference.Rda")




###PT object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200.rda")
sn.od.genes <- VariableFeatures(KB.PT)
table(Idents(KB.PT))

#update cluster metadata
meta <- KB.PT@meta.data
meta$pt.v2.clusters <- paste0("P_",meta$pt.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$pt.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.PT@meta.data <- meta

Idents(KB.PT)
KB.PT <- subset(KB.PT, downsample = 1000)
KB <- SCTransform(KB.PT, verbose = TRUE)
VariableFeatures(KB) <- sn.od.genes
KB <- RunPCA(KB, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v
colnames(pca.loadings) <- paste("pca_", 1:50, sep = "")
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KB, reduction = "pca")),]
colnames(cell.embeddings) <- paste("pca_", 1:50, sep = "")
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 30L,
              min.dist = 0.3, umap.method = 'uwot', return.model = TRUE)
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, repel = TRUE) + NoLegend()
KB <- FindNeighbors(
  object = KB,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Proximal-tubules_subset_filtered_SCT_Reference.Rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Proximal-tubules_subset_filtered_SCT_Reference.Rda")




###DT object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200.rda")
sn.od.genes <- VariableFeatures(KB.DT)
table(Idents(KB.DT))

#update cluster metadata
meta <- KB.DT@meta.data
meta$dt.v2.clusters <- paste0("D_",meta$dt.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$dt.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.DT@meta.data <- meta

Idents(KB.DT)
KB.DT <- subset(KB.DT, downsample = 1000)
KB <- SCTransform(KB.DT, verbose = TRUE)
VariableFeatures(KB) <- sn.od.genes
KB <- RunPCA(KB, verbose = FALSE)
pca.loadings <- p2$misc$PCA$v
colnames(pca.loadings) <- paste("pca_", 1:50, sep = "")
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])
cell.embeddings <- p2$reductions$PCA[rownames(Embeddings(KB, reduction = "pca")),]
colnames(cell.embeddings) <- paste("pca_", 1:50, sep = "")
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 30L,
              min.dist = 0.3, umap.method = 'uwot', return.model = TRUE)
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, repel = TRUE) + NoLegend()
KB <- FindNeighbors(
  object = KB,
  reduction = "pca",
  dims = 1:50,
  graph.name = "pca.annoy.neighbors", 
  k.param = 100,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

#FindNeighbors wipes the pca loadings, so re-add it
KB@reductions$pca@feature.loadings <- pca.loadings
Loadings(object = KB[["pca"]])

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Distal-tubules_subset_filtered_SCT_Reference.Rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Distal-tubules_subset_filtered_SCT_Reference.Rda")


