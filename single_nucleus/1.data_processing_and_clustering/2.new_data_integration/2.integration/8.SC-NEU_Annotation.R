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

nKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_Atlas_V2_04-2024-subset_Object_global-predictions.Rds")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/v2_object_filtered/Kidney_Atlas_V2_012024_Object_Filtered.Rds")

nKB.sc <- subset(nKB, predicted.v2.subclass.l1 %in% c("NEU")) #subset to SC/NEU cells
nKB.sc <- subset(nKB.sc, predicted.v2.subclass.l1.score > 0.8) #subset to high confidence cells


##Merge with full object
KB.sc <- subset(KB, v2.subclass.l1 %in% c("NEU")) #subset to SC/NEU cells
KB.sc$group <- "ref"

nKB.sc$group <- "query"
KB <- merge(x = KB.sc, y = nKB.sc, merge.dr = FALSE)
KB <- JoinLayers(KB)
KB[["sketch"]] <- NULL

##Clustering in Pagoda2 using integrated PCs
#Filtered Count Matrix from Seurat
#Filtered Count Matrix from Seurat
countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)
KB[["RNA"]]$counts <- countMatrix

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#1 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 10, n.odgenes = 100, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:100]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:10, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_sc <- k100infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.8,
        group.by = "pagoda_k100_infomap_sc", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.8,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()

#Query cells are well overlapped with reference cells

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_SC-NEU_k100_round1_newData-0424.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_SC-NEU_k100_round1_newData-0424.rda")
save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_SC-NEU_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_SC-NEU_subset_filtered_0424-newData.rda")

