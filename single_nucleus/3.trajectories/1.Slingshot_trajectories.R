library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
library(anndata)
library(SeuratWrappers)
library(magrittr)
library(cowplot)
library(slingshot)
library(Polychrome)
library(dendextend)
library("RColorBrewer")
library(gplots)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")
source("misc/utils.R")


###Fibroblast Trajectories
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_04192024.rds")
meta <- KB@meta.data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
load("color_factors_v2-clusters.robj")
KB.STR <- KB
KB.STR@meta.data <- meta[rownames(KB.STR@meta.data),]

###Subset to cortical set 1 - clusters 16-19 (perivascular fibroblasts)
KB.STR <- subset(KB.STR, v2.clusters %in% paste0("S_",c(16:19)))
countMatrix <- KB.STR[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#1144 overdispersed genes
p2$calculatePcaReduction(nPcs = 30, n.odgenes = 800, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:800]
VariableFeatures(KB.STR) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.STR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.STR))
KB.STR <- RunUMAP(object = KB.STR, reduction = "pca", dims = 1:30, n.neighbors = 20L,
                  min.dist = 0.01)

v2.fn.cols <- setNames(c("#C8FF00","#FFD300",
                         "#02AD24","#00FF00"),
                       c("S_16","S_17",
                         "S_18","S_19"))
Idents(KB.STR) <- "v2.clusters"
DimPlot(KB.STR, reduction = "umap", pt.size = 0.7, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.STR))],
) + NoLegend()

save(KB.STR, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_16-19_0424-newData.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_16-19_0424-newData.rda")

#Slingshot
sce <- as.SingleCellExperiment(KB.STR)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus="S_16")
SlingshotDataSet(sce)
#lineages: 1 
#Lineage1: S_16  S_17  S_18  S_19  

#curves: 1 
#Curve1: Length: 16.465	Samples: 9271

# Plot the lines and lineages
pdf(file='trajectories/S16-19-Lineages_Slingshot.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()

pdf(file='trajectories/S16-19-Lineages_Slingshot_2.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')
dev.off()

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_16-19_slingshot_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_16-19_slingshot_0424-newData.rda")

KB.STR <- AddMetaData(KB.STR, metadata = slingPseudotime(sce), col.name = c("pseudotime.lineage1"))
save(KB.STR, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_16-19_0424-newData.rda")


###Subset to cortical set 2 (interstitial fibroblasts)
KB.STR <- KB
KB.STR@meta.data <- meta[rownames(KB.STR@meta.data),]
KB.STR <- subset(KB.STR, v2.clusters %in% paste0("S_",c(9:14)))

countMatrix <- KB.STR[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#1429 overdispersed genes
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 600, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1429]
VariableFeatures(KB.STR) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.STR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.STR))
KB.STR <- RunUMAP(object = KB.STR, reduction = "pca", dims = 1:30, n.neighbors = 20L,
                  min.dist = 0.01)

v2.fn.cols <- setNames(adjustcolor(c("#9A4D42","#02AD24","#0000FF","#009FFF","#00479E","#FF0000"), alpha.f = 0.3),
                       c("S_9","S_10","S_11","S_12","S_13","S_14"))

Idents(KB.STR) <- "v2.clusters"
DimPlot(KB.STR, reduction = "umap", pt.size = 0.7, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.STR))],
) + NoLegend()

save(KB.STR, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_0424-newData.rda")

#Slingshot
sce <- as.SingleCellExperiment(KB.STR)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("S_9"), end.clus = c("S_12","S_14"))
SlingshotDataSet(sce)
#lineages: 4 
#Lineage1: S_9  S_13  S_12  
#Lineage2: S_9  S_11  
#Lineage3: S_9  S_10  
#Lineage4: S_9  S_14  

#curves: 4 
#Curve1: Length: 8.0891	Samples: 11982.29
#Curve2: Length: 9.6694	Samples: 14996.18
#Curve3: Length: 6.1791	Samples: 9811.18
#Curve4: Length: 6.8633	Samples: 9849.27

# Plot the lines and lineages
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(2,3))
lines(SlingshotDataSet(sce), lwd = 3, col = "darkred", type = 'curves', linInd = c(1,4))

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "gray30", type = 'curves', linInd = c(2,3))
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1,4))


save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_slingshot_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_slingshot_0424-newData.rda")

KB.STR <- AddMetaData(KB.STR, metadata = slingPseudotime(sce), col.name = c("pseudotime.lineage1","pseudotime.lineage2","pseudotime.lineage3","pseudotime.lineage4"))
KB.STR <- AddMetaData(KB.STR, metadata = slingAvgPseudotime(sce), col.name = c("avg.pseudotime"))
save(KB.STR, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_0424-newData.rda")


###Subset Cortical set 2 Lineage 1 (S9, S12, S13)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_0424-newData.rda")
l1 <- rownames(KB.STR@meta.data[!is.na(KB.STR@meta.data$pseudotime.lineage1), ])
KB.STR <- subset(KB.STR, cells = c(l1))
table(KB.STR$v2.clusters)
KB.STR <- subset(KB.STR, v2.clusters %in% paste0("S_", c(9,12,13)))

KB.STR@meta.data
countMatrix <- KB.STR[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#794 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 794, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:794]
VariableFeatures(KB.STR) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.STR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.STR))
KB.STR <- RunUMAP(object = KB.STR, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                  min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(adjustcolor(c("#9A4D42","#02AD24","#0000FF","#009FFF","#00479E","#FF0000"), alpha.f = 0.3),
                       c("S_9","S_10","S_11","S_12","S_13","S_14"))

DimPlot(KB.STR, reduction = "umap", pt.size = 0.7, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.STR))],
) + NoLegend()
FeaturePlot(KB.STR, reduction = "umap", features = c("pseudotime.lineage1"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.STR, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_Lineage1_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_Lineage1_0424-newData.rda")

                                             
###Subset Cortical set 2 Lineage 4
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_0424-newData.rda")
l4 <- rownames(KB.STR@meta.data[!is.na(KB.STR@meta.data$pseudotime.lineage4), ])
KB.STR <- subset(KB.STR, cells = c(l4))
table(KB.STR$v2.clusters)
KB.STR <- subset(KB.STR, v2.clusters %in% paste0("S_", c(9,11,14)))

KB.STR@meta.data
countMatrix <- KB.STR[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#900 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 900, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:900]
VariableFeatures(KB.STR) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.STR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.STR))
KB.STR <- RunUMAP(object = KB.STR, reduction = "pca", dims = 1:20, n.neighbors = 50L,
                  min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(adjustcolor(c("#9A4D42","#02AD24","#0000FF","#009FFF","#00479E","#FF0000"), alpha.f = 0.3),
                       c("S_9","S_10","S_11","S_12","S_13","S_14"))

DimPlot(KB.STR, reduction = "umap", pt.size = 0.7, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.STR))],
) + NoLegend()
FeaturePlot(KB.STR, reduction = "umap", features = c("pseudotime.lineage4"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()


save(KB.STR, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_Lineage4_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_Lineage4_0424-newData.rda")





                                             
###PT Trajectories
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
meta <- KB@meta.data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")
load("color_factors_v2-clusters.robj")
KB.PT <- KB
KB.PT@meta.data <- meta[rownames(KB.PT@meta.data),]

###Subset to PT lineages
KB.PT <- subset(KB.PT, v2.clusters %in% paste0("P_",c(4,8,10,12:22)))

countMatrix <- KB.PT[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#5559 overdispersed genes
p2$calculatePcaReduction(nPcs = 30, n.odgenes = 2000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2000]
VariableFeatures(KB.PT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.PT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.PT))
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:20, n.neighbors = 20L, spread = 3,
                 min.dist = 0.01)

v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#FE8F42","#DC5E93","#A10300","#858567","#783FC1","#FFD300","#9A4D42"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

Idents(KB.PT) <- "v2.clusters"
DimPlot(KB.PT, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.PT))],
) + NoLegend()

save(KB.PT, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")

#Slingshot
sce <- as.SingleCellExperiment(KB.PT)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("P_15"))
SlingshotDataSet(sce)
#lineages: 5 
#Lineage1: P_15  P_19  P_21  P_13  P_4  P_8  P_10  
#Lineage2: P_15  P_19  P_21  P_13  P_16  P_17  
#Lineage3: P_15  P_19  P_14  P_18  P_22  
#Lineage4: P_15  P_19  P_21  P_12  
#Lineage5: P_15  P_19  P_20  

#curves: 5 
#Curve1: Length: 50.628	Samples: 250514.4
#Curve2: Length: 34.397	Samples: 67338.86
#Curve3: Length: 21.246	Samples: 47450.98
#Curve4: Length: 31.695	Samples: 62681.92
#Curve5: Length: 19.852	Samples: 40413.34

plot(as.SlingshotDataSet(sce), col=c("black", "red", "green", "blue", "cyan"))

# Plot the lines and lineages
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1,4))
lines(SlingshotDataSet(sce), lwd = 3, col = "darkred", type = 'curves', linInd = c(2,5))
lines(SlingshotDataSet(sce), lwd = 3, col = "gray30", type = 'curves', linInd = c(3))

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_slingshot_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_slingshot_0424-newData.rda")

KB.PT <- AddMetaData(KB.PT, metadata = slingPseudotime(sce), col.name = c("pseudotime.lineage1","pseudotime.lineage2","pseudotime.lineage3","pseudotime.lineage4","pseudotime.lineage5"))
KB.PT <- AddMetaData(KB.PT, metadata = slingAvgPseudotime(sce), col.name = c("avg.pseudotime"))
save(KB.PT, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")


###Subset Lineage 2 (P_15/14/19/20/21/16/17)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
l2 <- rownames(KB.PT@meta.data[!is.na(KB.PT@meta.data$pseudotime.lineage2), ])
KB.PT <- subset(KB.PT, cells = c(l2))
table(KB.PT$v2.clusters)
KB.PT <- subset(KB.PT, v2.clusters %in% paste0("P_", c(15,14,19,20,21,16,17)))

KB.PT@meta.data
countMatrix <- KB.PT[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#3651 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 1000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1000]
VariableFeatures(KB.PT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.PT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.PT))
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                  min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#FE8F42","#DC5E93","#A10300","#858567","#783FC1","#FFD300","#9A4D42"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

DimPlot(KB.PT, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.PT))],
) + NoLegend()
FeaturePlot(KB.PT, reduction = "umap", features = c("pseudotime.lineage2"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.PT, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage2_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage2_0424-newData.rda")

###Subset Lineage 1 (P_15/14/19/20/21/13/4)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
l1 <- rownames(KB.PT@meta.data[!is.na(KB.PT@meta.data$pseudotime.lineage1), ])
KB.PT <- subset(KB.PT, cells = c(l1))
table(KB.PT$v2.clusters)
KB.PT <- subset(KB.PT, v2.clusters %in% paste0("P_", c(15,14,19,20,21,13,4)))
KB.PT <- subset(KB.PT, downsample = 10000) #Cap the number of P_4 cells

KB.PT@meta.data
countMatrix <- KB.PT[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#3267 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 1000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1000]
VariableFeatures(KB.PT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.PT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.PT))
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                 min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#FE8F42","#DC5E93","#A10300","#858567","#783FC1","#FFD300","#9A4D42"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

DimPlot(KB.PT, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.PT))],
) + NoLegend()
FeaturePlot(KB.PT, reduction = "umap", features = c("pseudotime.lineage1"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.PT, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage1_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage1_0424-newData.rda")

###Subset Lineage 5 (P_15/14/19/20/21/18/22)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
l5 <- rownames(KB.PT@meta.data[!is.na(KB.PT@meta.data$pseudotime.lineage5), ])
KB.PT <- subset(KB.PT, cells = c(l5))
table(KB.PT$v2.clusters)
KB.PT <- subset(KB.PT, v2.clusters %in% paste0("P_", c(15,14,19,21,18,22)))

KB.PT@meta.data
countMatrix <- KB.PT[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#2706 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 1000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1000]
VariableFeatures(KB.PT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.PT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.PT))
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                 min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#FE8F42","#DC5E93","#A10300","#858567","#783FC1","#FFD300","#9A4D42"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

DimPlot(KB.PT, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.PT))],
) + NoLegend()
FeaturePlot(KB.PT, reduction = "umap", features = c("pseudotime.lineage5"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.PT, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage5_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage5_0424-newData.rda")

###Subset Lineage 4 (P_15/14/19/21/12/10)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
l4 <- rownames(KB.PT@meta.data[!is.na(KB.PT@meta.data$pseudotime.lineage4), ])
KB.PT <- subset(KB.PT, cells = c(l4))
table(KB.PT$v2.clusters)
KB.PT <- subset(KB.PT, v2.clusters %in% paste0("P_", c(15,14,19,21,12,10)))
KB.PT <- subset(KB.PT, downsample = 10000) #to limit number of PT-S3 (P_10)

KB.PT@meta.data
countMatrix <- KB.PT[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#3291 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 1000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1000]
VariableFeatures(KB.PT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.PT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.PT))
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                 min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#FE8F42","#DC5E93","#A10300","#858567","#783FC1","#FFD300","#9A4D42"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

DimPlot(KB.PT, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.PT))],
) + NoLegend()
FeaturePlot(KB.PT, reduction = "umap", features = c("pseudotime.lineage4"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.PT, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage4_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage4_0424-newData.rda")






###TAL Trajectories
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_04192024.rds")
meta <- KB@meta.data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
load("color_factors_v2-clusters.robj")
KB.DT <- KB
KB.DT@meta.data <- meta[rownames(KB.DT@meta.data),]

###Subset lineages
KB.TAL <- subset(KB.DT, v2.clusters %in% paste0("D_",c(5:8,11:13,15:17)))
rm(KB.DT)
gc(reset=TRUE)

countMatrix <- KB.TAL[["RNA"]]$counts
countMatrix <- as(countMatrix, "dgCMatrix")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#4319 overdispersed genes
p2$calculatePcaReduction(nPcs = 30, n.odgenes = 3000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3000]
VariableFeatures(KB.TAL) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.TAL[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.TAL))
KB.TAL <- RunUMAP(object = KB.TAL, reduction = "pca", dims = 1:20, n.neighbors = 20L, spread = 3,
                  min.dist = 0.01)

v2.fn.cols <- setNames(c("#005300","#B1CC71","#886C00","#02AD24","#009FFF","#14F9FF","#1F9698",
                         "#DC5E93","#FFD300","#FE8F42"),
                       c("D_5","D_6","D_7","D_8","D_11","D_12","D_13",
                         "D_15","D_16","D_17"))


Idents(KB.TAL) <- "v2.clusters"
DimPlot(KB.TAL, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.TAL))],
) + NoLegend()

save(KB.TAL, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_0424-newData_TALA.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_0424-newData_TALA.rda")

#Slingshot
sce <- as.SingleCellExperiment(KB.TAL)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("D_16"))
SlingshotDataSet(sce)
#lineages: 3 
#Lineage1: D_16  D_17  D_6  D_7  D_5  
#Lineage2: D_16  D_17  D_6  D_7  D_8  
#Lineage3: D_16  D_15  

#curves: 3 
#Curve1: Length: 46.859	Samples: 75258.17
#Curve2: Length: 40.274	Samples: 78778.65
#Curve3: Length: 27.935	Samples: 42101

# Plot the lines and lineages
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1,2))
lines(SlingshotDataSet(sce), lwd = 3, col = "darkred", type = 'curves', linInd = c(3))

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_slingshot_0424-newData_TALA.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_slingshot_0424-newData_TALA.rda")

KB.TAL <- AddMetaData(KB.TAL, metadata = slingPseudotime(sce), col.name = c("pseudotime.lineage1","pseudotime.lineage2","pseudotime.lineage3"))
KB.TAL <- AddMetaData(KB.TAL, metadata = slingAvgPseudotime(sce), col.name = c("avg.pseudotime"))
save(KB.TAL, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_0424-newData_TALA.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_0424-newData_TALA.rda")

###Subset Lineage 3 (D_17,D_16,D_15)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_0424-newData_TALA.rda")
l3 <- rownames(KB.TAL@meta.data[!is.na(KB.TAL@meta.data$pseudotime.lineage3), ])
KB.TAL <- subset(KB.TAL, cells = c(l3))
table(KB.TAL$v2.clusters)
KB.TAL <- subset(KB.TAL, v2.clusters %in% paste0("D_", c(15:17)))

KB.TAL@meta.data
KB.TAL[["RNA"]]$counts <- as(KB.TAL[["RNA"]]$counts, "dgCMatrix")
countMatrix <- KB.TAL[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#2834 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 1000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1000]
VariableFeatures(KB.TAL) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.TAL[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.TAL))
KB.TAL <- RunUMAP(object = KB.TAL, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                  min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#005300","#B1CC71","#886C00","#02AD24","#009FFF","#14F9FF","#1F9698",
                         "#DC5E93","#FFD300","#FE8F42"),
                       c("D_5","D_6","D_7","D_8","D_11","D_12","D_13",
                         "D_15","D_16","D_17"))

DimPlot(KB.TAL, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.TAL))],
) + NoLegend()
FeaturePlot(KB.TAL, reduction = "umap", features = c("pseudotime.lineage3"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.TAL, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_Lineage3_0424-newData_TALA.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_Lineage3_0424-newData.rda_TALA")

###Subset Lineage 1 (D_16, D_17,D_7,D_5)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_0424-newData_TALA.rda")
l1 <- rownames(KB.TAL@meta.data[!is.na(KB.TAL@meta.data$pseudotime.lineage1), ])
KB.TAL <- subset(KB.TAL, cells = c(l1))
table(KB.TAL$v2.clusters)
KB.TAL <- subset(KB.TAL, v2.clusters %in% paste0("D_", c(5,17,16)))

KB.TAL@meta.data
KB.TAL[["RNA"]]$counts <- as(KB.TAL[["RNA"]]$counts, "dgCMatrix")
countMatrix <- KB.TAL[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#3094 overdispersed genes
p2$calculatePcaReduction(nPcs = 40, n.odgenes = 1000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1000]
VariableFeatures(KB.TAL) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.TAL[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.TAL))
KB.TAL <- RunUMAP(object = KB.TAL, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                  min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#005300","#B1CC71","#886C00","#02AD24","#009FFF","#14F9FF","#1F9698",
                         "#DC5E93","#FFD300","#FE8F42"),
                       c("D_5","D_6","D_7","D_8","D_11","D_12","D_13",
                         "D_15","D_16","D_17"))

DimPlot(KB.TAL, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.TAL))],
) + NoLegend()
FeaturePlot(KB.TAL, reduction = "umap", features = c("pseudotime.lineage1"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.TAL, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_Lineage1_0424-newData_TALA.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_Lineage1_0424-newData_TALA.rda")







###Myeloid Trajectories
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
meta <- KB@meta.data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")
KB.IMM <- KB
KB.IMM@meta.data <- meta[rownames(KB.IMM@meta.data),]
load("color_factors_v2-clusters.robj")
load("color_factors.robj")

###Subset to myeloid lineages - circulating monocyte lineage
KB.mye <- subset(KB.IMM, v2.clusters %in% paste0("I_",c(16:19,22)))
countMatrix <- KB.mye[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#430 overdispersed genes
p2$calculatePcaReduction(nPcs = 30, n.odgenes = 430, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:430]
VariableFeatures(KB.mye) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.mye[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.mye))
KB.mye <- RunUMAP(object = KB.mye, reduction = "pca", dims = 1:30, n.neighbors = 20L, spread = 3,
                  min.dist = 0.01)

v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("I_14","I_15","I_16","I_17",
                         "I_18","I_19","I_20","I_21","I_22"))

DimPlot(KB.mye, reduction = "umap", pt.size = 1, label = TRUE, alpha = 0.4,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.mye))],
) + NoLegend()

save(KB.mye, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_0424-newData.rda")

#Slingshot
sce <- as.SingleCellExperiment(KB.mye)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("I_22"))
SlingshotDataSet(sce)
#lineages: 3 
#Lineage1: I_22  I_16  I_18  
#Lineage2: I_22  I_16  I_19  
#Lineage3: I_22  I_16  I_17  

#curves: 3 
#Curve1: Length: 34.842	Samples: 3510.67
#Curve2: Length: 32.552	Samples: 3625.06
#Curve3: Length: 30.045	Samples: 2946.29

# Plot the lines and lineages
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_slingshot_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_slingshot_0424-newData.rda")

KB.mye <- AddMetaData(KB.mye, metadata = slingPseudotime(sce), col.name = c("pseudotime.lineage1","pseudotime.lineage2","pseudotime.lineage3"))
KB.mye <- AddMetaData(KB.mye, metadata = slingAvgPseudotime(sce), col.name = c("avg.pseudotime"))
save(KB.mye, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_0424-newData.rda")

###Subset to circulating monocyte lineage 1
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_0424-newData.rda")
l1 <- rownames(KB.mye@meta.data[!is.na(KB.mye@meta.data$pseudotime.lineage1), ])
KB.mye <- subset(KB.mye, cells = c(l1))
table(KB.mye$v2.clusters)
KB.mye <- subset(KB.mye, v2.clusters %in% paste0("I_", c(22,16,18)))

KB.mye@meta.data
countMatrix <- KB.mye[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#287 overdispersed genes
p2$calculatePcaReduction(nPcs = 20, n.odgenes = 287, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:287]
VariableFeatures(KB.mye) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.mye[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.mye))
KB.mye <- RunUMAP(object = KB.mye, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                  min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("I_14","I_15","I_16","I_17",
                         "I_18","I_19","I_20","I_21","I_22"))

DimPlot(KB.mye, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.mye))],
) + NoLegend()
FeaturePlot(KB.mye, reduction = "umap", features = c("pseudotime.lineage1"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.mye, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_Lineage1_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_Lineage1_0424-newData.rda")

###Subset to circulating monocyte lineage 3
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_0424-newData.rda")
l3 <- rownames(KB.mye@meta.data[!is.na(KB.mye@meta.data$pseudotime.lineage3), ])
KB.mye <- subset(KB.mye, cells = c(l3))
table(KB.mye$v2.clusters)
KB.mye <- subset(KB.mye, v2.clusters %in% paste0("I_", c(22,17)))

KB.mye@meta.data
countMatrix <- KB.mye[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#136 overdispersed genes
p2$calculatePcaReduction(nPcs = 20, n.odgenes = 136, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:136]
VariableFeatures(KB.mye) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.mye[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.mye))
KB.mye <- RunUMAP(object = KB.mye, reduction = "pca", dims = 1:20, n.neighbors = 20L,
                  min.dist = 0.01)

#set color palette
v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("I_14","I_15","I_16","I_17",
                         "I_18","I_19","I_20","I_21","I_22"))

DimPlot(KB.mye, reduction = "umap", pt.size = 0.7, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.mye))],
) + NoLegend()
FeaturePlot(KB.mye, reduction = "umap", features = c("pseudotime.lineage3"), pt.size = 0.7, label = FALSE, alpha = 0.6,
            repel = TRUE,
) + NoLegend()

save(KB.mye, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_Lineage3_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_Lineage3_0424-newData.rda")

###Subset to resident MAC lineage
KB.mye <- subset(KB.IMM, v2.clusters %in% paste0("I_",c(14,15)))
table(Idents(KB.mye))
KB.mye <- subset(KB.mye, downsample = 10000) #limit res macs to 10k cells

countMatrix <- KB.mye[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#463 overdispersed genes
p2$calculatePcaReduction(nPcs = 30, n.odgenes = 463, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:463]
VariableFeatures(KB.mye) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.mye[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.mye))
KB.mye <- RunUMAP(object = KB.mye, reduction = "pca", dims = 1:30, n.neighbors = 20L, spread = 3,
                  min.dist = 0.01)

v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("I_14","I_15","I_16","I_17",
                         "I_18","I_19","I_20","I_21","I_22"))

DimPlot(KB.mye, reduction = "umap", pt.size = 1, label = TRUE, alpha = 0.4,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.mye))],
) + NoLegend()

save(KB.mye, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I14-15_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I14-15_0424-newData.rda")


#Slingshot
sce <- as.SingleCellExperiment(KB.mye)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("I_14"))
SlingshotDataSet(sce)
#lineages: 1 
#Lineage1: I_14  I_15  

#curves: 1 
#Curve1: Length: 33.445	Samples: 11142

# Plot the lines and lineages
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')

plot(reducedDims(sce)$UMAP, col = v2.fn.cols[sce$v2.clusters], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I14-15_slingshot_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I14-15_slingshot_0424-newData.rda")

KB.mye <- AddMetaData(KB.mye, metadata = slingPseudotime(sce), col.name = c("pseudotime.lineage1"))
save(KB.mye, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I14-15_0424-newData.rda")
