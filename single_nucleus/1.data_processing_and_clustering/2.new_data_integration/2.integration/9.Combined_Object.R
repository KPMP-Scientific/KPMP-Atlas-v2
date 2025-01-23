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

nKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_Atlas_V2_04-2024-subset_Object.Rds")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/v2_object_filtered/Kidney_Atlas_V2_012024_Object_Filtered.Rds")

KB.all <- merge(x = KB, y = nKB, merge.dr = FALSE)
rm(nKB)
KB.all
#73176 features across 1477957 samples

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
KB.str <- KB@meta.data

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")
KB.imm <- KB@meta.data

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")
KB.ec <- KB@meta.data

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")
KB.pt <- KB@meta.data

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
KB.dt <- KB@meta.data

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_SC-NEU_subset_filtered_0424-newData.rda")
KB.sc <- KB@meta.data

col.to.use <- intersect(colnames(KB.str), c(colnames(KB.imm), colnames(KB.ec), colnames(KB.pt), colnames(KB.dt), colnames(KB.sc)))
col.to.use <- col.to.use[!col.to.use %in% grep("predicted", col.to.use, value = T)]
col.to.use <- col.to.use[!col.to.use %in% "set"]
meta <- rbind(KB.str[,col.to.use], KB.imm[,col.to.use], KB.ec[,col.to.use], KB.pt[,col.to.use], KB.dt[,col.to.use], KB.sc[,col.to.use])

KB.all <- subset(KB.all, cells = rownames(meta))
KB.all <- JoinLayers(object = KB.all)
KB.all[["sketch"]] <- NULL

# Write the matrix to a directory
write_matrix_dir(
  mat = KB.all[["RNA"]]$counts,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/full_kidney_count_set_0424",
  overwrite = TRUE
)

kidney.mat <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/full_kidney_count_set_0424")
KB.all[["RNA"]]$counts <- kidney.mat
KB.all[["RNA"]]$data.1 <- NULL
KB.all[["RNA"]]$data <- NULL
KB.all <- NormalizeData(KB.all)
KB.all
#An object of class Seurat 
#36588 features across 1388643 samples within 1 assay 
#Active assay: RNA (36588 features, 0 variable features)
#2 layers present: counts, data

save(KB.all, file = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_04192024.rda")


###Global UMAP visualization
##Create a sketch data set
KB.all <- FindVariableFeatures(KB.all)
KB.all <- SketchData(
  object = KB.all,
  ncells = 200000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
KB.all

DefaultAssay(KB.all) <- "sketch"


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.all[["sketch"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#5241 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 5241, maxit=1000)

#pagoda2 PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:5241]
VariableFeatures(KB.all) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.all[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.all))

KB.all <- RunUMAP(object = KB.all, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.2, return.model = T)

DimPlot(KB.all, reduction = "umap", pt.size = 0.5) 
DimPlot(KB.all, reduction = "umap", pt.size = 0.5, label = TRUE, 
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()


#Add pc loadings
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("pca_", 1:50, sep = "")
KB.all@reductions$pca@feature.loadings <- pca.loadings
Embeddings(KB.all, reduction = "pca")


#Extend results to the full datasets
KB.all <- ProjectData(
  object = KB.all,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = NULL
)

DefaultAssay(KB.all) <- "RNA"

#Update metadata
meta <- rbind(KB.str[,col.to.use], KB.imm[,col.to.use], KB.ec[,col.to.use], KB.pt[,col.to.use], KB.dt[,col.to.use], KB.sc[,col.to.use])
KB.all@meta.data <- meta[rownames(KB.all@meta.data),]

table(KB.all$v2.subclass.l1)
KB.all@meta.data[rownames(KB.sc),]$v2.clusters <- "N_1"

#update cluster metadata
meta <- KB.all@meta.data
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_04192024.txt")
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

DefaultAssay(KB.all) <- "RNA"
KB.all[["pca"]] <- KB.all[["pca.full"]]
KB.all[["umap"]] <- KB.all[["full.umap"]]
VariableFeatures(KB.all) <- sn.od.genes
KB.all[["sketch"]] <- NULL
KB.all[["pca.full"]] <- NULL
KB.all[["full.umap"]] <- NULL
saveRDS(KB.all, "~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_04192024.rds")

KB.all
#An object of class Seurat 
#36588 features across 1388643 samples within 1 assay 
#Active assay: RNA (36588 features, 5241 variable features)
#2 layers present: counts, data
#2 dimensional reductions calculated: pca, umap
