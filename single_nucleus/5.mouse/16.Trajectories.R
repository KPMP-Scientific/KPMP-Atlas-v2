library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(pagoda2)
library(BPCells)
library(slingshot)

options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

###Stromal Trajectories
##Perivascular Set
#Load Combined Data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
table(KB$v2.clusters)
KB <- subset(KB, v2.clusters %in% paste0("S_",c(14:19)))

##Integrate and Cluster
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- ScaleData(KB)
KB <- RunPCA(KB)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)
KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:20)

DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) 
table(Idents(KB))
v2.fn.cols <- setNames(c("#C8FF00","#FFD300",
                         "#02AD24","#02AD24","#02AD24","#00FF00"),
                       c("S_14","S_15",
                         "S_16","S_17","S_18","S_19"))


pdf(file='mouse_IRI/trajectories/pvFIB-Lineages_UMAP.pdf',width=8,height=6)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_pvFIB-Subset.Rds")

#Trajectory Analysis
sce <- as.SingleCellExperiment(KB)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("S_14"))

SlingshotDataSet(sce)
#lineages: 5 
#Lineage1: P_14  P_16  P_10  P_7  P_8  P_9  
#Lineage2: P_14  P_15  P_4  P_5  
#Lineage3: P_14  P_15  P_3  
#Lineage4: P_14  P_16  P_17  
#Lineage5: P_14  P_15  P_13  

# Plot the lines and lineages
pdf(file='mouse_IRI/trajectories/pvFIB-Lineages_Slingshot.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_pvFIB-Subset_slingshot.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_pvFIB-Subset_slingshot.rda")



##Interstitial Set
#Load Combined Data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
table(KB$v2.clusters)
KB <- subset(KB, v2.clusters %in% paste0("S_",c(5:12)))

##Integrate and Cluster
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- ScaleData(KB)
KB <- RunPCA(KB)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)
KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:30)

DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) 
table(Idents(KB))
v2.fn.cols <- setNames(adjustcolor(c("#9A4D42","#0000FF","#0000FF","#00ffb7","#00479E","#009FFF","#FF0000","#FF0000"), alpha.f = 0.3),
                       c("S_5","S_6","S_7","S_8","S_9","S_10","S_11","S_12"))


pdf(file='mouse_IRI/trajectories/intFIB-Lineages_UMAP.pdf',width=8,height=6)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_intFIB-Subset.Rds")


#Trajectory Analysis
KB[["RNA"]] <- as(KB[["RNA"]], Class = "Assay")
sce <- as.SingleCellExperiment(KB)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("S_5"), end.clus=c("S_10","S_12"))

SlingshotDataSet(sce)
#lineages: 3 
#Lineage1: S_5  S_9  S_7  S_8  S_11  S_12  
#Lineage2: S_5  S_9  S_7  S_10  
#Lineage3: S_5  S_6  


# Plot the lines and lineages
pdf(file='mouse_IRI/trajectories/intFIB-Lineages_Slingshot.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_intFIB-Subset_slingshot.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_intFIB-Subset_slingshot.rda")




### aPT Trajectory
#Load Combined Data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")
table(KB$v2.clusters)
KB <- subset(KB, v2.clusters %in% paste0("P_",c(3,4,5,7,8,9,10,13:17)))

##Integrate and Cluster
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- ScaleData(KB)
KB <- RunPCA(KB)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)
KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:30)

DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) 

v2.fn.cols <- setNames(c("#00479E","#009FFF","#009FFF","#0000FF","#0000FF","#0000FF","#0000FF",
                         "#FE8F42","#B1CC71","#00FFBE","#14F9FF","#14F9FF"),
                       c("P_3","P_4","P_5","P_7","P_8","P_9","P_10",
                         "P_13","P_14","P_15","P_16","P_17"))

pdf(file='mouse_IRI/trajectories/PT-Lineages_UMAP.pdf',width=8,height=6)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-aPT-Subset.Rds")


#Trajectory Analysis
sce <- as.SingleCellExperiment(KB)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("P_14"))

SlingshotDataSet(sce)
#lineages: 5 
#Lineage1: P_14  P_16  P_10  P_7  P_8  P_9  
#Lineage2: P_14  P_15  P_4  P_5  
#Lineage3: P_14  P_15  P_3  
#Lineage4: P_14  P_16  P_17  
#Lineage5: P_14  P_15  P_13  

# Plot the lines and lineages
pdf(file='mouse_IRI/trajectories/PT-Lineages_Slingshot.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()

pdf(file='mouse_IRI/trajectories/PT-Lineages_Slingshot_2.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')
dev.off()

pdf(file='mouse_IRI/trajectories/PT-Lineages_Slingshot_3.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1,2,3,5))
lines(SlingshotDataSet(sce), lwd = 3, col = "darkred", type = 'curves', linInd = c(4))
dev.off()

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-aPT-Subset_slingshot.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-aPT-Subset_slingshot.rda")






### aTAL Trajectory  -------------------------------------------------------
#Load Combined Data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")
table(KB$v2.clusters)
KB <- subset(KB, v2.clusters %in% paste0("D_",c(6,7,14:16)))

##Integrate and Cluster
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB, nfeatures = 500)
KB <- ScaleData(KB)
KB <- RunPCA(KB, npcs = 20)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)
KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:20)

DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) 
DimPlot(KB, group.by = "source", reduction = "umap", label = TRUE) 


###cluster using Pagoda2
KB[["RNA"]] <- JoinLayers(KB[["RNA"]])

countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k200infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k200infomap, col.name = "rpca_k200_infomap_atal")
DimPlot(KB, group.by = "rpca_k200_infomap_atal", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()

Idents(KB) <- "rpca_k200_infomap_atal"
cl15.markers <- FindMarkers(KB, ident.1 = "15", min.pct = 0.25, logfc.threshold = 0.25)

#remove 2 small outlier clusters from Kirita et al
KB <- subset(KB, rpca_k200_infomap_atal %in% c(12,15), invert = TRUE)

KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:20)

DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) 
DimPlot(KB, group.by = "source", reduction = "umap", label = TRUE) 

Idents(KB) <- "v2.clusters"
v2.fn.cols <- setNames(c("#005300","#B1CC71",
                         "#DC5E93","#FFD300","#FE8F42"),
                       c("D_6","D_7",
                         "D_16","D_14","D_15"))

pdf(file='mouse_IRI/trajectories/TAL-Lineages_UMAP.pdf',width=8,height=6)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_TALA-aTAL-Subset.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_TALA-aTAL-Subset.Rds")

#Trajectory Analysis
sce <- as.SingleCellExperiment(KB)
sce <- slingshot(sce, clusterLabels = 'v2.clusters', reducedDim = 'UMAP', start.clus=c("D_14"))

SlingshotDataSet(sce)
#lineages: 2 
#Lineage1: D_14  D_15  D_6  D_7  
#Lineage2: D_14  D_15  D_16  

plot(as.SlingshotDataSet(sce), col=c("black", "red", "green", "blue", "cyan"))

# Plot the lines and lineages
pdf(file='mouse_IRI/trajectories/TALA-Lineages_Slingshot.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()

pdf(file='mouse_IRI/trajectories/TALA-Lineages_Slingshot_2.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'lineages')
dev.off()

pdf(file='mouse_IRI/trajectories/TALA-Lineages_Slingshot_3.pdf',width=8,height=6)
plot(reducedDims(sce)$UMAP, col = v2.fn.cols[as.character(sce$v2.clusters)], pch=16, cex = 0.6)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1))
lines(SlingshotDataSet(sce), lwd = 3, col = "darkred", type = 'curves', linInd = c(2))
dev.off()

save(sce, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_TALA-aTAL-Subset_slingshot.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_TALA-aTAL-Subset_slingshot.rda")
