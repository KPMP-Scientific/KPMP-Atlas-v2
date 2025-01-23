library(Seurat)
library(ggplot2)
library(pagoda2)
library(spacexr)
library(future)
plan("multisession", workers = 10)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")
load("color_factors_v2-clusters.robj")

nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")

####RCTD
###Proximal Epithelial Sub-clusters
##Load Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
table(ref$v2.subclass.l1)
ref <- subset(ref, v2.subclass.l1 %in% c("POD","PEC","PT","DTL","ATL"))
ref <- subset(ref, v2.state.l2 %in% c("degenerative","cycling"), invert = TRUE)
ref <- subset(ref, v2.subclass.sp %in% c("-"), invert = TRUE)
Idents(ref) <- "v2.subclass.sp"
ref <- subset(ref, downsample = 5000)
table(ref$v2.subclass.sp)
ref <- subset(ref, v2.subclass.sp %in% "aDTL", invert = TRUE) #overpredicts

#extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
ref$v2.subclass.sp <- gsub("/","",ref$v2.subclass.sp)
cluster <- as.factor(ref$v2.subclass.sp)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

#RCTD method with full mode
nano.obj.sub <- subset(nano.obj, v2.subclass.l1 %in% c("POD","PEC","PT","TL"))
table(nano.obj.sub$v2.subclass.l1)
counts <- nano.obj.sub[["Nanostring"]]$counts
coords <- data.frame(x = rep(1, length(colnames(counts))),
                     y = rep(2, length(colnames(counts))))
rownames(coords) <- colnames(counts)
query <- SpatialRNA(coords, counts, colSums(counts), use_fake_coords = TRUE)

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
save(RCTD, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Proximal-Epithelial_Combined_RCTD_08052024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Proximal-Epithelial_Combined_RCTD_08052024.rda")

#add l3 prediction weights
rctd.cls <- RCTD@results$weights
rctd.cls <- 100 * sweep(rctd.cls, 1, rowSums(rctd.cls), "/")
rctd.cls <- as.data.frame(rctd.cls)
maxWeight.l3 <- as.numeric(apply(rctd.cls, 1, max))
names(maxWeight.l3) <- rownames(rctd.cls)
maxWeight.l3.ct <- colnames(rctd.cls)[apply(rctd.cls, 1, which.max)]
names(maxWeight.l3.ct) <- rownames(rctd.cls)

rctd.cls <- as.data.frame(rctd.cls)[colnames(nano.obj.sub),]
rctd.cls[is.na(rctd.cls)] <- 0
rownames(rctd.cls) <- colnames(nano.obj.sub)

nano.obj.sub[["l3.predictions"]] <- CreateAssayObject(counts = t(as.matrix(rctd.cls)))
nano.obj.sub$maxWeight.l3 <- as.numeric(maxWeight.l3[rownames(nano.obj.sub@meta.data)])
nano.obj.sub$maxCellType.l3 <- as.character(maxWeight.l3.ct[rownames(nano.obj.sub@meta.data)])
table(nano.obj.sub$maxCellType.l3)

saveRDS(nano.obj.sub, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_pEpi_Sub.RDS")
nano.obj.sub <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_pEpi_Sub.RDS")

###Distal Tubule Sub-clusters
##Load Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
table(ref$v2.subclass.l1)
ref <- subset(ref, v2.subclass.l1 %in% c("TAL","DCT","CNT"))
ref <- subset(ref, v2.state.l2 %in% c("degenerative","cycling"), invert = TRUE)
ref <- subset(ref, v2.subclass.sp %in% c("-"), invert = TRUE)
Idents(ref) <- "v2.subclass.sp"
ref <- subset(ref, downsample = 5000)
table(ref$v2.subclass.sp)

#extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
ref$v2.subclass.sp <- gsub("/","",ref$v2.subclass.sp)
cluster <- as.factor(ref$v2.subclass.sp)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

#RCTD method with full mode
nano.obj.sub <- subset(nano.obj, v2.subclass.l1 %in% c("TAL","DCT/CNT"))
table(nano.obj.sub$v2.subclass.l1)
counts <- nano.obj.sub[["Nanostring"]]$counts
coords <- data.frame(x = rep(1, length(colnames(counts))),
                     y = rep(2, length(colnames(counts))))
rownames(coords) <- colnames(counts)
query <- SpatialRNA(coords, counts, colSums(counts), use_fake_coords = TRUE)

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
save(RCTD, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Distal-Epithelial_Combined_RCTD_08052024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Distal-Epithelial_Combined_RCTD_08052024.rda")

#add l3 prediction weights
rctd.cls <- RCTD@results$weights
rctd.cls <- 100 * sweep(rctd.cls, 1, rowSums(rctd.cls), "/")
rctd.cls <- as.data.frame(rctd.cls)
maxWeight.l3 <- as.numeric(apply(rctd.cls, 1, max))
names(maxWeight.l3) <- rownames(rctd.cls)
maxWeight.l3.ct <- colnames(rctd.cls)[apply(rctd.cls, 1, which.max)]
names(maxWeight.l3.ct) <- rownames(rctd.cls)

rctd.cls <- as.data.frame(rctd.cls)[colnames(nano.obj.sub),]
rctd.cls[is.na(rctd.cls)] <- 0
rownames(rctd.cls) <- colnames(nano.obj.sub)

nano.obj.sub[["l3.predictions"]] <- CreateAssayObject(counts = t(as.matrix(rctd.cls)))
nano.obj.sub$maxWeight.l3 <- as.numeric(maxWeight.l3[rownames(nano.obj.sub@meta.data)])
nano.obj.sub$maxCellType.l3 <- as.character(maxWeight.l3.ct[rownames(nano.obj.sub@meta.data)])
table(nano.obj.sub$maxCellType.l3)

saveRDS(nano.obj.sub, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_dEpi_Sub.RDS")
nano.obj.sub <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_dEpi_Sub.RDS")

###Collecting Duct Sub-clusters
##Load Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
table(ref$v2.subclass.l1)
ref <- subset(ref, v2.subclass.l1 %in% c("PC","IC"))
ref <- subset(ref, v2.state.l2 %in% c("degenerative","cycling"), invert = TRUE)
ref <- subset(ref, v2.subclass.sp %in% c("-"), invert = TRUE)
Idents(ref) <- "v2.subclass.sp"
ref <- subset(ref, downsample = 5000)
table(ref$v2.subclass.sp)

#extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
ref$v2.subclass.sp <- gsub("/","",ref$v2.subclass.sp)
cluster <- as.factor(ref$v2.subclass.sp)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

#RCTD method with full mode
nano.obj.sub <- subset(nano.obj, v2.subclass.l1 %in% c("PC","IC"))
table(nano.obj.sub$v2.subclass.l1)
counts <- nano.obj.sub[["Nanostring"]]$counts
coords <- data.frame(x = rep(1, length(colnames(counts))),
                     y = rep(2, length(colnames(counts))))
rownames(coords) <- colnames(counts)
query <- SpatialRNA(coords, counts, colSums(counts), use_fake_coords = TRUE)

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
save(RCTD, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Connecting-Epithelial_Combined_RCTD_08052024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Connecting-Epithelial_Combined_RCTD_08052024.rda")

#add l3 prediction weights
rctd.cls <- RCTD@results$weights
rctd.cls <- 100 * sweep(rctd.cls, 1, rowSums(rctd.cls), "/")
rctd.cls <- as.data.frame(rctd.cls)
maxWeight.l3 <- as.numeric(apply(rctd.cls, 1, max))
names(maxWeight.l3) <- rownames(rctd.cls)
maxWeight.l3.ct <- colnames(rctd.cls)[apply(rctd.cls, 1, which.max)]
names(maxWeight.l3.ct) <- rownames(rctd.cls)

rctd.cls <- as.data.frame(rctd.cls)[colnames(nano.obj.sub),]
rctd.cls[is.na(rctd.cls)] <- 0
rownames(rctd.cls) <- colnames(nano.obj.sub)

nano.obj.sub[["l3.predictions"]] <- CreateAssayObject(counts = t(as.matrix(rctd.cls)))
nano.obj.sub$maxWeight.l3 <- as.numeric(maxWeight.l3[rownames(nano.obj.sub@meta.data)])
nano.obj.sub$maxCellType.l3 <- as.character(maxWeight.l3.ct[rownames(nano.obj.sub@meta.data)])
table(nano.obj.sub$maxCellType.l3)

saveRDS(nano.obj.sub, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_cEpi_Sub.RDS")
nano.obj.sub <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_cEpi_Sub.RDS")

###EC-IMM Sub-clusters
##Load Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
table(ref$v2.subclass.l1)
ref <- subset(ref, v2.subclass.l1 %in% c("EC","FIB","VSM/P","Lymphoid","Myeloid"))
ref <- subset(ref, v2.state.l2 %in% c("degenerative","cycling"), invert = TRUE)
ref <- subset(ref, v2.subclass.sp %in% c("-"), invert = TRUE)
Idents(ref) <- "v2.subclass.sp"
ref <- subset(ref, downsample = 5000)
table(ref$v2.subclass.sp)

#extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
ref$v2.subclass.sp <- gsub("/","",ref$v2.subclass.sp)
cluster <- as.factor(ref$v2.subclass.sp)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

#RCTD method with full mode
nano.obj.sub <- subset(nano.obj, v2.subclass.l1 %in% c("EC","FIB","VSMP","Lymphoid","Myeloid","ERY","cycling"))
table(nano.obj.sub$v2.subclass.l1)
counts <- nano.obj.sub[["Nanostring"]]$counts
coords <- data.frame(x = rep(1, length(colnames(counts))),
                     y = rep(2, length(colnames(counts))))
rownames(coords) <- colnames(counts)
query <- SpatialRNA(coords, counts, colSums(counts), use_fake_coords = TRUE)

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
save(RCTD, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/EC-IMM_Combined_RCTD_08052024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/EC-IMM_Combined_RCTD_08052024.rda")

#add l3 prediction weights
rctd.cls <- RCTD@results$weights
rctd.cls <- 100 * sweep(rctd.cls, 1, rowSums(rctd.cls), "/")
rctd.cls <- as.data.frame(rctd.cls)
maxWeight.l3 <- as.numeric(apply(rctd.cls, 1, max))
names(maxWeight.l3) <- rownames(rctd.cls)
maxWeight.l3.ct <- colnames(rctd.cls)[apply(rctd.cls, 1, which.max)]
names(maxWeight.l3.ct) <- rownames(rctd.cls)

rctd.cls <- as.data.frame(rctd.cls)[colnames(nano.obj.sub),]
rctd.cls[is.na(rctd.cls)] <- 0
rownames(rctd.cls) <- colnames(nano.obj.sub)

nano.obj.sub[["l3.predictions"]] <- CreateAssayObject(counts = t(as.matrix(rctd.cls)))
nano.obj.sub$maxWeight.l3 <- as.numeric(maxWeight.l3[rownames(nano.obj.sub@meta.data)])
nano.obj.sub$maxCellType.l3 <- as.character(maxWeight.l3.ct[rownames(nano.obj.sub@meta.data)])
table(nano.obj.sub$maxCellType.l3)

saveRDS(nano.obj.sub, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_EC-IMM_Sub.RDS")
nano.obj.sub <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_EC-IMM_Sub.RDS")

###Combined Object
nano.p.epi <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_pEpi_Sub.RDS")
nano.non.epi <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_EC-IMM_Sub.RDS")
nano.c.epi <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_cEpi_Sub.RDS")
nano.d.epi <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_dEpi_Sub.RDS")

nano.obj <- merge(x = nano.p.epi, y = c(nano.d.epi,nano.c.epi,nano.non.epi))
table(nano.obj$maxCellType.l3)
table(nano.obj$v2.subclass.l1)

nano.obj.comb <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")
nano.obj.comb <- subset(nano.obj.comb, cells = colnames(nano.obj))
nano.obj.comb[["l3.predictions"]] <- nano.obj[["l3.predictions"]]

nano.obj.comb <- AddMetaData(nano.obj.comb, metadata = nano.obj@meta.data[,34:37])

saveRDS(nano.obj.comb, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")






#### Clustering
library(Seurat)
library(ggplot2)
library(pagoda2)
library(spacexr)
library("corrplot")

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")
load("color_factors_v2-clusters.robj")

###V2 atlas variable genes
load("AtlasV2_VariableFeatures.rda")


###Global Clustering 
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")

###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#90 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference atlas overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
od.genes <- unique(c(pt.od.genes,dt.od.genes,imm.od.genes,str.od.genes,ec.od.genes,neu.od.genes))
sn.od.genes <- od.genes[od.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done separately from K100

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap

saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Pagoda_K100.RDS")

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap <- k100infomap[rownames(nano.obj@meta.data)]
nano.obj@meta.data$pagoda_k50_infomap <- k50infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k50_infomap") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxCellType.l3") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")



###Proximal Epithelial Sub-clusters
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")
nano.obj <- subset(nano.obj, v2.subclass.l1 %in% c("POD","PEC","PT","TL"))

###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#19 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference atlas overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- pt.od.genes[pt.od.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap

saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Pagoda_K100_PT.RDS")

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap.pt <- k100infomap[rownames(nano.obj@meta.data)]
nano.obj@meta.data$pagoda_k50_infomap.pt <- k50infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap.pt") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k50_infomap.pt") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxCellType.l3") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_PT.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_PT.RDS")

###Annotate to subclass.l3
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_PT.RDS")
nano.obj@meta.data

##Generate table of cluster cell types
Idents(nano.obj) <- "pagoda_k100_infomap.pt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxCellType.l3))), 2)
    max.len = 2
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(query = prop.table(table(seurat.obj.sub$maxCellType.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxCellType.l3)))
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(nano.obj, celltype)
max.idents

##Check marker genes
pt.markers <- c(
  "NPHS1","NPHS2","ST6GALNAC3","PODXL",                         #POD
  "PTGDS","BST2","TPPP3","CHI3L1",                              #dPOD
  "ALDH1A2","CFH","FAM155A","CLDN1",                            #PEC
  
  "LRP2","CUBN",                                                #PT
  "SLC5A12","SLC22A6","HNF4A",                                  #S1/S2
  "RALYL","PCDH15","PRODH2","SLC22A8",                          #S1
  "SLC34A1","ANKS1B","SLC5A10","SLC5A11",                       #S2                                  
  "SLC5A8","GPM6A","SLC22A24","SLC7A13",                        #PT-S3
  
  "IGFBP7","SPP1","ITGB8","CDH6","TMEM178B","ALPK2","HAVCR1",
  "ITGB3",
  "CST3","CLU","VIM","PIGR",#"APOE",                            #aPT2
  "IL32","SOX4","VCAM1","MMP7","SOX9","CCL2",                   #aPT2
  "DCC",                                                        #aPT1
  "GDA","GLIS1",                                                #aPT-S1/S2
  "DLGAP1","PROM1",                                             #frPTS1/S2
  "APBB1IP","ROBO2","COL23A1","MEG3",                           #frPTS1/S2
  "LSAMP","KCNIP1","NRXN3","WNT2B",                             #frPTS1/S2
  "KCTD16","SPON1",                                             #aPT-S3
  "NRG3","FAM189A1","DTNA","KITLG","GRM8",                      #frPT-S3
  "TOP2A","MKI67",                                              #cycling
  "EGR1","FOS",                                                 #dPT
  
  "PAX8-AS1","SLC44A5","CRYAB","TACSTD2",                       #TL
  "ABCA13",                                                     #DTL
  "AQP1", "UNC5D","LRRC4C","DENND2A","SLCO1A2",                 #DTL2
  "IRX3", "SERPINA1",                                           #aDTL2
  "SIM2",                                                       #DTL1/3/ATL
  "ADGRL3","JAG1","SMAD9","ID1",                                #DTL1
  "SLC14A2","FAM169A","SMOC2",                                  #DTL3
  "ABCA4","BCL6","AKR1B1","SH3GL3",                             #DTL3/ATL
  "CLCNKA","PROX1","CACNA1C","COLEC10",                         #ATL
  "SOD3","HSPA2","CSRP2",                                       #dATL
  "CST3","APOE","GATM","ALDOB","CLU","DEFB1",                   #Injury
  "SPP1","IGFBP7","CLU"
)

Idents(nano.obj) <- "pagoda_k100_infomap.pt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))

DotPlot(nano.obj, features = unique(pt.markers[pt.markers %in% rownames(nano.obj)]), dot.scale = 8) + RotatedAxis()


##Add ref atlas correlations
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
var.genes <- VariableFeatures(ref)
Idents(object = ref) <- "v2.subclass.sp"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","aDTL","DTL1","DTL3","ATL",
  "M-TAL","C-TAL","MD","frTAL","aTAL1","aTAL2","DCT","aDCT","CNT","aCNT","PC",
  "IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
  "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM","M-FIB","C/M-FIB",
  "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
  "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
  "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
  "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N","-"))

ref <- subset(ref, idents = c("POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","DTL1","DTL3","ATL"))
ref.ds <- subset(ref, downsample = 1000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")
ref.ds <- NormalizeData(ref.ds)

var.genes <- pt.od.genes[pt.od.genes %in% rownames(nano.obj)]
ref.ds <- ScaleData(ref.ds, features = var.genes, assay = "RNA")

nano.obj <- NormalizeData(nano.obj)
nano.obj <- ScaleData(nano.obj, features = var.genes)

Idents(nano.obj) <- "pagoda_k100_infomap.pt"
levels(Idents(object = nano.obj)) <- paste("CL", levels(Idents(object = nano.obj)), sep = "")
ave.nano.obj <- AverageExpression(nano.obj, features = var.genes, assay = "Nanostring", layer = "scale.data")
ave.ref <- AverageExpression(ref.ds, features = var.genes, assays = "RNA", layer = "scale.data")
library("corrplot")
ave.cor <- cor(cbind(as.data.frame(ave.ref),as.data.frame(ave.nano.obj)))
ave.cor
rownames(ave.cor) <- gsub("RNA.","",rownames(ave.cor))
colnames(ave.cor) <- gsub("Nanostring.","",colnames(ave.cor))
ave.cor <- ave.cor[1:10,11:21]

#Generate table
Idents(nano.obj) <- "pagoda_k100_infomap.pt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxCor <- function(ave.cor, celltype) {
  max.cor <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    if (ct %in% colnames(ave.cor) == FALSE) {
      max.cor <- data.frame(ref.cor.name = c(NA,NA),
                            ref.cor.value = c(NA,NA))
    } else {
      max.cor <- data.frame(ref.cor.name = tail(names(sort(ave.cor[,ct])), 2),
                            ref.cor.value = tail(sort(ave.cor[,ct]), 2))
    } 
    
  }))
  return(max.cor)
}

max.cor <- getMaxCor(ave.cor, celltype)
max.idents$cor.subclass <- max.cor$ref.cor.name
max.idents$cor.value <- max.cor$ref.cor.value

#Write out table
write.table(max.idents, file="nanostring/PT_Kidney_atlasV2_subclass_level3_overlaps_08-07-2024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

###Cluster annotations
#1 - PT-S1/S2
#2 - aPT2
#3 - PT-S3
#4 - aPT1
#5 - frPT
#6 - POD
#7 - aPT1
#8 - aPT1
#9 - PEC
#10 - TL
#11 - frPT

###Update experiment metadata - check region localization
nano.obj@meta.data
colnames(nano.obj@meta.data)[colnames(nano.obj@meta.data) == "orig.ident"] <- "library"
nano.obj@meta.data$library <- gsub("[.]", "_", nano.obj@meta.data$Run_Tissue_name)
nano.obj@meta.data$fov <- paste0(nano.obj@meta.data$library, "_", nano.obj@meta.data$fov)
meta <- nano.obj@meta.data
exp.meta <- read.delim("nanostring/Kidney_AtlasV2_Nanostring_Experiment_Metadata.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level1","region_level2",
         "age_binned","sex","race","tissue_preservation","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$fov, exp.meta$fov)]
}

nano.obj@meta.data <- meta
table(nano.obj$region_level1,nano.obj$pagoda_k100_infomap.pt)

#Add annotations to nano.obj
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_PT.RDS")
Idents(nano.obj) <- "pagoda_k100_infomap.pt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
cls <- c(1:11)
new.cls <- c("PT-S1/S2","aPT2","PT-S3","aPT1","frPT","POD","aPT1","aPT1","PEC","TL","frPT")

Idents(nano.obj) <- plyr::mapvalues(Idents(nano.obj), from = cls, to = new.cls)
Idents(nano.obj) <- factor(Idents(nano.obj), levels = c("POD","PEC","PT-S1/S2","PT-S3","aPT2","aPT1","frPT",
                                                        "TL"))

table(Idents(nano.obj))
nano.obj$v2.subclass.l3 <- Idents(nano.obj)

v2.scl3.cols.mod <- v2.scl2.cols
names(v2.scl3.cols.mod) <- gsub("PT-S1","PT-S1/S2",names(v2.scl3.cols.mod))
names(v2.scl3.cols.mod) <- gsub("DTL","TL",names(v2.scl3.cols.mod))

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,raster=FALSE,
        group.by = "v2.subclass.l3", repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl3.cols.mod[levels(Idents(object = nano.obj))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_PT_anno.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_PT_anno.RDS")





###Distal Epithelial Sub-clusters
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")
table(nano.obj$v2.subclass.l1)
nano.obj <- subset(nano.obj, v2.subclass.l1 %in% c("TAL","DCT/CNT","PC","IC"))

###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#32 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference atlas overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- dt.od.genes[dt.od.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap

saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Pagoda_K100_DT.RDS")

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap.dt <- k100infomap[rownames(nano.obj@meta.data)]
nano.obj@meta.data$pagoda_k50_infomap.dt <- k50infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap.dt") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k50_infomap.dt") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxCellType.l3") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_DT.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_DT.RDS")



###Annotate to subclass.l3
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_DT.RDS")
nano.obj@meta.data

##Generate table of cluster cell types
Idents(nano.obj) <- "pagoda_k100_infomap.dt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxCellType.l3))), 2)
    max.len = 2
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(query = prop.table(table(seurat.obj.sub$maxCellType.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxCellType.l3)))
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(nano.obj, celltype)
max.idents

##Check marker genes
dt.markers <- c(
  "CASR","SLC12A1","UMOD","EGF","ESRRB",                  #TAL
  "ANK2","CLCNKA","KCTD16","BMP6","RIPOR3","CLDN14",      #M-TAL
  "PHACTR1","SLCO3A1","CXCL12","CNTN1","CABP1",           #TAL-A
  "KCNMB2","RGS6",                                        #C/M-TAL-A
  "ENOX1","CALCR","RBM20","PDE3A",                        #C-TAL-A
  "DACH1","LRMDA",                                        #TAL-B
  
  "TENM4","FGF14","PXDNL","GRM8",                         #C/M-TAL-B
  "KCNJ10","TMEM52B","CLDN16","TMEM207","JAG1",           #TAL-B
  "SCN7A","COL8A1","LINGO2",                              #C-TAL-B                          
  "BBOX1","NOS1","ROBO2","TMPRSS4",                       #MD
  
  "CREB5","COL4A1","ITGA3","HIVEP3","CYTOR",              #aTAL1
  "CD44","DIAPH3","SYT16","HAVCR1",                       #aTAL1                                      
  
  "ITGB6","NRP1","TFPI",                                  #aTAL
  
  "HIF1A","ADAMTS1","DNAH5",                              #aTAL2
  
  "ITGB8","PROM1","ARHGAP26","RNF144B","TMPRSS4","RHEX",  #frTAL
  "EZH2","CENPP","MKI67",#"TOP2A",                        #cycling
  
  "SLC12A3","CNNM2","KLHL3","TRPM6",                       #DCT
  "GRAMD1B","ADAMTS17","ZNF385D",                          #DCT1
  
  "UNC5C","HMCN1","FSTL4","TRPV5",                         #DCT2
  "CACNA1D", "ACSL4",                                      #frDCT
  "FGF13","IGF2BP2","FAM155A","NRG1",                      #aDCT 
  
  "SLC8A1",                                                #DCT2 / CNT
  "HSD11B2","CALB1","ANGPT1","SNTG1",                      #CNT
  "CTSD",                                                  #dCNT
  "RAPGEF5","DLGAP1","BIRC3",                              #aCNT
  
  "GATA3","AQP2","PAPPA",                                  #PC
  "SCNN1G","SCNN1B",
  "SGCD","STC1",                                           #CNT-PC
  "FGF12","PTCSC3","CNTNAP3B",                             #CCD/OMCD PC
  "MCTP1","CPAMD8","GABRA2",                               #OMCD-PC
  "GREB1L",                                                #OMCD-PC/IMCD
  "SLC14A2","HS3ST5","TENM3","TGFBR3","AKR1C1","FAT3",    #IMCD
  "AOC1","TFPI2","LCN2",                                   #dIMCD
  "ADIRF","SNCG","DHRS2","GPX2","TP63",#"FXYD3",           #PapE                              #PapE
  
  "ATP6V0D2", "ATP6V1C2", "CLNK",                          #IC
  "SLC26A7", "SLC4A1",                                     #IC-A                                   
  "HS6ST3","NRXN3", "NXPH2", "LEF1",                       #CCD-IC-A
  
  "FAM184B","ADTRP","AQP6","STAP1",                        #OMCD-IC-A
  
  "SLC4A9", "SLC26A4", "INSRR", "TLDC2",                #IC-B
  
  "CKB","COX8A",#"PEBP1","UQCRB",                         #Injury
  "CST3","DEFB1","SPP1","IGFBP7","CLU"
)


Idents(nano.obj) <- "pagoda_k100_infomap.dt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))

DotPlot(nano.obj, features = unique(dt.markers[dt.markers %in% rownames(nano.obj)]), dot.scale = 8) + RotatedAxis()


##Add ref atlas correlations
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
var.genes <- VariableFeatures(ref)
Idents(object = ref) <- "v2.subclass.sp"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","aDTL","DTL1","DTL3","ATL",
  "M-TAL","C-TAL","MD","frTAL","aTAL1","aTAL2","DCT","aDCT","CNT","aCNT","PC",
  "IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
  "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM","M-FIB","C/M-FIB",
  "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
  "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
  "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
  "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N","-"))

ref <- subset(ref, idents = c("M-TAL","C-TAL","MD","frTAL","aTAL1","aTAL2","DCT","aDCT","CNT","aCNT","PC",
                              "IMCD","IC-A","tPC-IC","IC-B"))
ref.ds <- subset(ref, downsample = 1000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")
ref.ds <- NormalizeData(ref.ds)

var.genes <- dt.od.genes[dt.od.genes %in% rownames(nano.obj)]
ref.ds <- ScaleData(ref.ds, features = var.genes, assay = "RNA")

nano.obj <- NormalizeData(nano.obj)
nano.obj <- ScaleData(nano.obj, features = var.genes)

Idents(nano.obj) <- "pagoda_k100_infomap.dt"
levels(Idents(object = nano.obj)) <- paste("CL", levels(Idents(object = nano.obj)), sep = "")
ave.nano.obj <- AverageExpression(nano.obj, features = var.genes, assay = "Nanostring", layer = "scale.data")
ave.ref <- AverageExpression(ref.ds, features = var.genes, assays = "RNA", layer = "scale.data")
library("corrplot")
ave.cor <- cor(cbind(as.data.frame(ave.ref),as.data.frame(ave.nano.obj)))
ave.cor
rownames(ave.cor) <- gsub("RNA.","",rownames(ave.cor))
colnames(ave.cor) <- gsub("Nanostring.","",colnames(ave.cor))
ave.cor <- ave.cor[1:15,16:31]

#Generate table
Idents(nano.obj) <- "pagoda_k100_infomap.dt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxCor <- function(ave.cor, celltype) {
  max.cor <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    if (ct %in% colnames(ave.cor) == FALSE) {
      max.cor <- data.frame(ref.cor.name = c(NA,NA),
                            ref.cor.value = c(NA,NA))
    } else {
      max.cor <- data.frame(ref.cor.name = tail(names(sort(ave.cor[,ct])), 2),
                            ref.cor.value = tail(sort(ave.cor[,ct]), 2))
    } 
    
  }))
  return(max.cor)
}

max.cor <- getMaxCor(ave.cor, celltype)
max.idents$cor.subclass <- max.cor$ref.cor.name
max.idents$cor.value <- max.cor$ref.cor.value

#Write out table
write.table(max.idents, file="nanostring/DT_Kidney_atlasV2_subclass_level3_overlaps_08-07-2024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

###Cluster annotations (check region info)
#1 - aTAL
#2 - aCNT
#3 - DCT
#4 - ambiguous (only 2 cycling cells)
#5 - CNT
#6 - dTAL
#7 - M-PC
#8 - IC-A
#9 - ambiguous
#10 - C-PC
#11 - frTAL
#12 - IC-B
#13 - frTAL
#14 - C-TAL
#15 - M-TAL
#16 - MD

###Update experiment metadata
nano.obj@meta.data
colnames(nano.obj@meta.data)[colnames(nano.obj@meta.data) == "orig.ident"] <- "library"
nano.obj@meta.data$library <- gsub("[.]", "_", nano.obj@meta.data$Run_Tissue_name)
nano.obj@meta.data$fov <- paste0(nano.obj@meta.data$library, "_", nano.obj@meta.data$fov)
meta <- nano.obj@meta.data
exp.meta <- read.delim("nanostring/Kidney_AtlasV2_Nanostring_Experiment_Metadata.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level1","region_level2",
         "age_binned","sex","race","tissue_preservation","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$fov, exp.meta$fov)]
}

nano.obj@meta.data <- meta
table(nano.obj$region_level1,nano.obj$pagoda_k100_infomap.dt)

#Add annotations to nano.obj
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_DT.RDS")
Idents(nano.obj) <- "pagoda_k100_infomap.dt"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
cls <- c(1:16)
new.cls <- c("aTAL","aCNT","DCT","ambiguous","CNT","dTAL","M-PC","IC-A","ambiguous",
             "C-PC","frTAL","IC-B","frTAL","C-TAL","M-TAL","MD")

Idents(nano.obj) <- plyr::mapvalues(Idents(nano.obj), from = cls, to = new.cls)
Idents(nano.obj) <- factor(Idents(nano.obj), levels = c("M-TAL","C-TAL","MD","aTAL","frTAL","dTAL","DCT",
                                                        "CNT","aCNT","C-PC","M-PC",
                                                        "IC-A","IC-B","ambiguous"))

table(Idents(nano.obj))
nano.obj$v2.subclass.l3 <- Idents(nano.obj)

v2.scl3.cols.mod <- v2.scl2.cols
names(v2.scl3.cols.mod) <- gsub("aTAL1","aTAL",names(v2.scl3.cols.mod))

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,raster=FALSE,
        group.by = "v2.subclass.l3", repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl3.cols.mod[levels(Idents(object = nano.obj))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_DT_anno.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_DT_anno.RDS")








###Endothelial Sub-clusters
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")
table(nano.obj$v2.subclass.l1)
nano.obj <- subset(nano.obj, v2.subclass.l1 %in% c("EC"))

###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#10 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference atlas overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- ec.od.genes[ec.od.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap

saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Pagoda_K100_EC.RDS")

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap.ec <- k100infomap[rownames(nano.obj@meta.data)]
nano.obj@meta.data$pagoda_k50_infomap.ec <- k50infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap.ec") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k50_infomap.ec") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxCellType.l3") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_EC.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_EC.RDS")

###Annotate to subclass.l3
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_EC.RDS")
nano.obj@meta.data

##Generate table of cluster cell types
Idents(nano.obj) <- "pagoda_k100_infomap.ec"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxCellType.l3))), 2)
    max.len = 2
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(query = prop.table(table(seurat.obj.sub$maxCellType.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxCellType.l3)))
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(nano.obj, celltype)
max.idents

##Check marker genes
ec.markers <- c("PECAM1","PTPRB","FLT1",                                   #Broad EC
                "EMCN","HECW2","ITGA8",                                    #EC-GC
                "EHD3","RXFP1", "MGP","SOST",                              #EC-GC      
                'PTCHD4',"ZMAT3",                                          #aEC-GC
                "AQP1","FILIP1","H19","ESM1","SLC45A4",                    #EC-GC-FILIP1+
                
                "PDE3A","SULF1","NKAIN2","NOS1",                           #EC-AA
                "ADAMTS6","MCTP1","PALMD","SLC14A1","ITIH5",               #EC-DVR  
                #"LYPD6B","EPHA3","STXBP6","CP",
                
                "CEACAM1","PLVAP","DNASE1L3",                              #PTC/AVR
                "COL15A1","FABP5","ALPL", #"CD36",                         #M-EC-PTC
                "GPM6A","NR2F2",                                           #PTC/AVR
                "ZNF385D","RANBP3L","EDIL3","TEK","GAS6","CD9",            #EC-AVR
                "MX2","RSAD2","ISG15","IFIT1",                             #iaEC-AVR
                
                "SLCO2A1",                                                
                "VWF","RYR3","ADGRG6","CPE","TRABD2B",                     #EC-V
                "ADAMTSL1","CMKLR1",                                       #EC-V/EC-PCV
                "DOK6",                                                    #EC-PCV
                "KDR","FLT1",                                              #EC-PTC
                "NAV3","OSMR","SYCP2L",                                    #C-EC-PTC
                "AFAP1L1","USP31","MYO1B","LAMA4","NETO2",                 #angEC-PTC
                "SLC6A6","FUT8","ATP13A3","AFF3",                          #EC-EA
                "IFITM3","HLA-DRA","CAVIN2","CCL14","CA4",                 #dEC-PTC
                
                'ICAM1',"TNFAIP3",'CCL2','SELE',"CXCL2",                   #infEC-PTC
                "VCAM1",                                                   #inf/iaEC-PTC
                "CXCL10","GBP1","CXCL11","CTSS",                           #iaEC-PTC
                "MMRN1","CD36","TBX1","PROX1",                             #EC-LYM
                "TOP2A","MKI67","CENPF"                                    #cycling
                )


Idents(nano.obj) <- "pagoda_k100_infomap.ec"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))

DotPlot(nano.obj, features = unique(ec.markers[ec.markers %in% rownames(nano.obj)]), dot.scale = 8) + RotatedAxis()


##Add ref atlas correlations
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
var.genes <- VariableFeatures(ref)
Idents(object = ref) <- "v2.subclass.sp"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","aDTL","DTL1","DTL3","ATL",
  "M-TAL","C-TAL","MD","frTAL","aTAL1","aTAL2","DCT","aDCT","CNT","aCNT","PC",
  "IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
  "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM","M-FIB","C/M-FIB",
  "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
  "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
  "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
  "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N","-"))

ref <- subset(ref, idents = c("EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
                              "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM"))
ref.ds <- subset(ref, downsample = 1000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")
ref.ds <- NormalizeData(ref.ds)

var.genes <- ec.od.genes[ec.od.genes %in% rownames(nano.obj)]
ref.ds <- ScaleData(ref.ds, features = var.genes, assay = "RNA")

nano.obj <- NormalizeData(nano.obj)
nano.obj <- ScaleData(nano.obj, features = var.genes)

Idents(nano.obj) <- "pagoda_k100_infomap.ec"
levels(Idents(object = nano.obj)) <- paste("CL", levels(Idents(object = nano.obj)), sep = "")
ave.nano.obj <- AverageExpression(nano.obj, features = var.genes, assay = "Nanostring", layer = "scale.data")
ave.ref <- AverageExpression(ref.ds, features = var.genes, assays = "RNA", layer = "scale.data")
library("corrplot")
ave.cor <- cor(cbind(as.data.frame(ave.ref),as.data.frame(ave.nano.obj)))
ave.cor
rownames(ave.cor) <- gsub("RNA.","",rownames(ave.cor))
colnames(ave.cor) <- gsub("Nanostring.","",colnames(ave.cor))
ave.cor <- ave.cor[1:13,14:20]

#Generate table
Idents(nano.obj) <- "pagoda_k100_infomap.ec"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxCor <- function(ave.cor, celltype) {
  max.cor <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    if (ct %in% colnames(ave.cor) == FALSE) {
      max.cor <- data.frame(ref.cor.name = c(NA,NA),
                            ref.cor.value = c(NA,NA))
    } else {
      max.cor <- data.frame(ref.cor.name = tail(names(sort(ave.cor[,ct])), 2),
                            ref.cor.value = tail(sort(ave.cor[,ct]), 2))
    } 
    
  }))
  return(max.cor)
}

max.cor <- getMaxCor(ave.cor, celltype)
max.idents$cor.subclass <- max.cor$ref.cor.name
max.idents$cor.value <- max.cor$ref.cor.value

#Write out table
write.table(max.idents, file="nanostring/EC_Kidney_atlasV2_subclass_level3_overlaps_08-07-2024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

###Cluster annotations (check region info)
#1 - EC
#2 - EC-PTC
#3 - EC-GC
#4 - EC-LYM
#5 - EC-AVR/EC-V
#6 - EC-AEA/DVR
#7 - M-EC-PTC

###Update experiment metadata
nano.obj@meta.data
colnames(nano.obj@meta.data)[colnames(nano.obj@meta.data) == "orig.ident"] <- "library"
nano.obj@meta.data$library <- gsub("[.]", "_", nano.obj@meta.data$Run_Tissue_name)
nano.obj@meta.data$fov <- paste0(nano.obj@meta.data$library, "_", nano.obj@meta.data$fov)
meta <- nano.obj@meta.data
exp.meta <- read.delim("nanostring/Kidney_AtlasV2_Nanostring_Experiment_Metadata.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level1","region_level2",
         "age_binned","sex","race","tissue_preservation","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$fov, exp.meta$fov)]
}

nano.obj@meta.data <- meta
table(nano.obj$region_level1,nano.obj$pagoda_k100_infomap.ec)

#Add annotations to nano.obj
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_EC.RDS")
Idents(nano.obj) <- "pagoda_k100_infomap.ec"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
cls <- c(1:7)
new.cls <- c("EC","EC-PTC","EC-GC","EC-LYM","EC-AVR/EC-V","EC-AEA/DVR","M-EC-PTC"
)

Idents(nano.obj) <- plyr::mapvalues(Idents(nano.obj), from = cls, to = new.cls)
Idents(nano.obj) <- factor(Idents(nano.obj), levels = c("EC","EC-GC","EC-AEA/DVR","M-EC-PTC","EC-AVR/EC-V","EC-PTC","EC-LYM"))

table(Idents(nano.obj))
nano.obj$v2.subclass.l3 <- Idents(nano.obj)

v2.scl3.cols.mod <- v2.scl2.cols
names(v2.scl3.cols.mod) <- gsub("EC-AA","EC-AEA/DVR",names(v2.scl3.cols.mod))
names(v2.scl3.cols.mod) <- gsub("EC-AVR","EC-AVR/EC-V",names(v2.scl3.cols.mod))
names(v2.scl3.cols.mod) <- gsub("C-EC-PTC","EC-PTC",names(v2.scl3.cols.mod))

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.4,raster=FALSE,
        group.by = "v2.subclass.l3", repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl3.cols.mod[levels(Idents(object = nano.obj))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_EC_anno.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_EC_anno.RDS")





###Stromal Sub-clusters
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")
table(nano.obj$v2.subclass.l1)
nano.obj <- subset(nano.obj, v2.subclass.l1 %in% c("FIB", "VSMP"))

###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#9 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference atlas overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- str.od.genes[str.od.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap

saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Pagoda_K100_STR.RDS")

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap.str <- k100infomap[rownames(nano.obj@meta.data)]
nano.obj@meta.data$pagoda_k50_infomap.str <- k50infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap.str") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k50_infomap.str") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxCellType.l3") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_STR.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_STR.RDS")

###Annotate to subclass.l3
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_STR.RDS")
nano.obj@meta.data

##Generate table of cluster cell types
Idents(nano.obj) <- "pagoda_k50_infomap.str"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxCellType.l3))), 2)
    max.len = 2
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(query = prop.table(table(seurat.obj.sub$maxCellType.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxCellType.l3)))
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(nano.obj, celltype)
max.idents

##Check marker genes
str.markers <- c(
  "DCN","C7","PDGFRA",                                     #Pan FIB
  "SYT1","TNC","IGFBP2","RARRES2",                         #Pan Medullary FIB
  "CA8","HPSE2","GABRG3",                                  #IM/OM-FIB
  "MCTP2","SNAP25","BDKRB2",                               #IM-FIB
  "KIF26B","FREM1","KCNQ3","LOXHD1",                       #OM & C/M-FIB
  "SPP1","TIMP1",                                          #dM-FIB
  "COL1A2","COL3A1","COL1A1",                              #dFIB & MYOF
  "ADAMTSL1","KCNIP1","ADARB2",                            #C/M-FIB
  "RYR2","ZNF536","SEMA3D","ACTG2",# "PAMR1",              #IM-pvMYOF
  
  "NEGR1","LAMA2","ABCA8","MEG3",                           #Pan cortical FIB
  "CCN1","CCN2","ELL2","SAMHD1","SLC2A3",                   #C-FIB (interstitial fib)
  "GRIN2A","EMID1",                                         #C-FIB-PATH
  "SELENOP","LUM","CXCL12",'GGT5',"ECRG4",                  #C-FIB-OSMRlo
  "OSMR","SOD2","UGCG","IL1R1","CCL2",                      #C-FIB-OSMRhi
  "RELB","CXCL10","CCL19",                                  #C-FIB-OSMRhi
  "SULF1","GLI2","NTM","INHBA","FAP","POSTN",               #C-MYOF
  "SPARC","BGN","VIM",                                      #dFIB
  
  "FLRT2","COL12A1","FGF14",                                #Pan pvFIB
  "PDZRN4",'IGF1','ADAMTS3',"RSPO3","WNT5B",                #pvFIB-RSPO3+
  "C3","EBF2","SFRP2","CD34","PI16",                        #pvFIB-PI16+
  "ITGBL1","PLD5","CNTN5",                                  #pvFIB & pvMYOF
  "MGAT4C","EPHA3",                                         #pvFIB
  "ITGA8",                                                  #pvMYOF & VSMC
  "ADGRL3","MYH11","ACTA2",'KCNMA1',"PCDH7",                #pvMYOF
  "PRUNE2","MYOCD",#"SYNPO2",'MACROD2',                     #pvMYOF
  
  
  "PDGFRB","SLCO3A1",                                       #Pan VSM markers
  "SAMD12","ROBO1","PIEZO2",                                #MC & REN
  "DAAM2","GATA3","DCC","POSTN","IL1RL1",                   #MC
  "ROBO2","REN","KCNQ5","SLCO2A1","SPOCK3",                 #REN
  "MYH11","NTRK3",'RGS6','KCNMA1',"ADRA1A","MCAM",          #VSMC
  "NOTCH3",                                                 #VSMC &VSMC/P
  "RGS5",'PLCB4',"FRMD3",                                   #VSMC/P
  "ADGRB3","SLC38A11","C2CD2","DNAH11","SLC6A1",            #M-VSMC/P
  'PDE1C',"STEAP4","RXFP1",                                 #VSMC/P (cortical)
  'FLNA', 'TAGLN',"ACTA2","VIM","ACTB",                     #dVSMC
  "CD36","PLIN1","LPL","ADIPOQ","FABP4"                     #Ad
)


Idents(nano.obj) <- "pagoda_k50_infomap.str"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))

DotPlot(nano.obj, features = unique(str.markers[str.markers %in% rownames(nano.obj)]), dot.scale = 8) + RotatedAxis()


##Add ref atlas correlations
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
var.genes <- VariableFeatures(ref)
Idents(object = ref) <- "v2.subclass.sp"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","aDTL","DTL1","DTL3","ATL",
  "M-TAL","C-TAL","MD","frTAL","aTAL1","aTAL2","DCT","aDCT","CNT","aCNT","PC",
  "IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
  "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM","M-FIB","C/M-FIB",
  "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
  "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
  "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
  "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N","-"))

ref <- subset(ref, idents = c("M-FIB","C/M-FIB",
                              "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
                              "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P"))
ref.ds <- subset(ref, downsample = 1000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")
ref.ds <- NormalizeData(ref.ds)

var.genes <- str.od.genes[str.od.genes %in% rownames(nano.obj)]
ref.ds <- ScaleData(ref.ds, features = var.genes, assay = "RNA")

nano.obj <- NormalizeData(nano.obj)
nano.obj <- ScaleData(nano.obj, features = var.genes)

Idents(nano.obj) <- "pagoda_k50_infomap.str"
levels(Idents(object = nano.obj)) <- paste("CL", levels(Idents(object = nano.obj)), sep = "")
ave.nano.obj <- AverageExpression(nano.obj, features = var.genes, assay = "Nanostring", layer = "scale.data")
ave.ref <- AverageExpression(ref.ds, features = var.genes, assays = "RNA", layer = "scale.data")
library("corrplot")
ave.cor <- cor(cbind(as.data.frame(ave.ref),as.data.frame(ave.nano.obj)))
ave.cor
rownames(ave.cor) <- gsub("RNA.","",rownames(ave.cor))
colnames(ave.cor) <- gsub("Nanostring.","",colnames(ave.cor))
ave.cor <- ave.cor[1:16,17:57]

#Generate table
Idents(nano.obj) <- "pagoda_k50_infomap.str"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxCor <- function(ave.cor, celltype) {
  max.cor <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    if (ct %in% colnames(ave.cor) == FALSE) {
      max.cor <- data.frame(ref.cor.name = c(NA,NA),
                            ref.cor.value = c(NA,NA))
    } else {
      max.cor <- data.frame(ref.cor.name = tail(names(sort(ave.cor[,ct])), 2),
                            ref.cor.value = tail(sort(ave.cor[,ct]), 2))
    } 
    
  }))
  return(max.cor)
}

max.cor <- getMaxCor(ave.cor, celltype)
max.idents$cor.subclass <- max.cor$ref.cor.name
max.idents$cor.value <- max.cor$ref.cor.value

#Write out table
write.table(max.idents, file="nanostring/STR_Kidney_atlasV2_subclass_level3_overlaps_08-07-2024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

###Cluster annotations (check region info)
#1 - VSMC/P
#2 - FIB
#3 - FIB
#4 - C-FIB-OSMRlo
#5 - C-FIB-OSMRhi
#6 - pvFIB-PI16+
#7 - FIB
#8 - C-FIB-OSMRlo
#9 - C-FIB-OSMRlo
#10 - C-MYOF
#11 - FIB
#12 - FIB
#13 - FIB
#14 - FIB
#15 - FIB
#16 - FIB
#17 - FIB
#18 - FIB
#19 - pvFIB-RSPO3+
#20 - FIB
#21 - pvFIB-RSPO3+
#22 - FIB
#23 - FIB
#24 - FIB
#25 - pvMYOF
#26 - MC
#27 - FIB
#28 - FIB
#29 - FIB
#30 - FIB
#31 - FIB
#32 - FIB
#33 - REN
#34 - FIB
#35 - VSMC
#36 - FIB
#37 - FIB
#38 - FIB
#39 - FIB
#40 - FIB
#41 - FIB

###Update experiment metadata
nano.obj@meta.data
colnames(nano.obj@meta.data)[colnames(nano.obj@meta.data) == "orig.ident"] <- "library"
nano.obj@meta.data$library <- gsub("[.]", "_", nano.obj@meta.data$Run_Tissue_name)
nano.obj@meta.data$fov <- paste0(nano.obj@meta.data$library, "_", nano.obj@meta.data$fov)
meta <- nano.obj@meta.data
exp.meta <- read.delim("nanostring/Kidney_AtlasV2_Nanostring_Experiment_Metadata.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level1","region_level2",
         "age_binned","sex","race","tissue_preservation","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$fov, exp.meta$fov)]
}

nano.obj@meta.data <- meta
table(nano.obj$region_level1,nano.obj$pagoda_k50_infomap.str)

Idents(nano.obj) <- "pagoda_k50_infomap.str"
CL16_markers <- FindMarkers(nano.obj, ident.1 = "16", min.pct = 0.25, logfc.threshold = 0.25)
DotPlot(ref.ds, features = c("LCN2","MT1X"), dot.scale = 8) + RotatedAxis()



#Add annotations to nano.obj
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_STR.RDS")
Idents(nano.obj) <- "pagoda_k50_infomap.str"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
cls <- c(1:41)
new.cls <- c("VSMC/P","FIB","FIB","C-FIB-OSMRlo","C-FIB-OSMRhi","pvFIB-PI16+",
  "FIB","C-FIB-OSMRlo","C-FIB-OSMRlo","C-MYOF","FIB","FIB","FIB","FIB","FIB",
  "FIB","FIB","FIB","pvFIB-RSPO3+","FIB","pvFIB-RSPO3+","FIB","FIB","FIB",
  "pvMYOF","MC","FIB","FIB","FIB","FIB","FIB","FIB","REN","FIB","VSMC","FIB",
  "FIB","FIB","FIB","FIB","FIB"
)

Idents(nano.obj) <- plyr::mapvalues(Idents(nano.obj), from = cls, to = new.cls)
Idents(nano.obj) <- factor(Idents(nano.obj), levels = c("FIB","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","pvFIB-RSPO3+",
                                                        "pvFIB-PI16+","pvMYOF","MC","REN","VSMC","VSMC/P"))

table(Idents(nano.obj))
nano.obj$v2.subclass.l3 <- Idents(nano.obj)

v2.scl3.cols.mod <- v2.scl3.cols
names(v2.scl3.cols.mod) <- gsub("C-FIB","FIB",names(v2.scl3.cols.mod))

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.4,raster=FALSE,
        group.by = "v2.subclass.l3", repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl3.cols.mod[levels(Idents(object = nano.obj))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_STR_anno.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_STR_anno.RDS")






###Immune Sub-clusters
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")
table(nano.obj$v2.subclass.l1)
nano.obj <- subset(nano.obj, v2.subclass.l1 %in% c("Lymphoid", "Myeloid"))

###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#13 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference atlas overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- imm.od.genes[imm.od.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap

saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Pagoda_K100_IMM.RDS")

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap.imm <- k100infomap[rownames(nano.obj@meta.data)]
nano.obj@meta.data$pagoda_k50_infomap.imm <- k50infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap.imm") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k50_infomap.imm") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxCellType.l3") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_IMM.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_IMM.RDS")

###Annotate to subclass.l3
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_IMM.RDS")
nano.obj@meta.data

##Generate table of cluster cell types
Idents(nano.obj) <- "pagoda_k100_infomap.imm"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxCellType.l3))), 2)
    max.len = 2
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(query = prop.table(table(seurat.obj.sub$maxCellType.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxCellType.l3)))
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(nano.obj, celltype)
max.idents

##Check marker genes
imm.markers <- c("PTPRC",                                             #Broad Immune
                 "BANK1","MS4A1","CD37","CD79A",                      #B Cells
                 "IGKC","XBP1","MZB1","JCHAIN",#"SDC1",               #PL Cells
                 "CD96","CD247","BCL11B","THEMIS",                    #T
                 "INPP4B","TRAC","CD3D",                              #T
                 "IL7R",                                              #T/MAIT/ILC3
                 "LEF1","CD4","SELL",                                 #Naïve Th
                 "SLC4A10","KLRB1","CCR6",                            #MAIT
                 "PCDH9","TOX2","KIT","RORC",                         #ILC3
                 "IKZF2","RTKN2","IL2RA","CTLA4",#"FOXP3",            #T-REG
                 "CD8A",                                              #CD8+
                 "GZMK",                                              #CD8+ TEM/TRM
                 "CCL5","SAMD3","GZMA","CCL4",                        #CD8+ & NK 
                 "NKG7","KLRD1","GNLY","GZMB","CX3CR1",               #CD8+ TEM/TEMRA & NK
                 "GZMH",                                              #CD8+ TEM/TEMRA
                 "TXK","KLRF1","NCAM1",#"PDGFD",                      #NK
                 
                 
                 "HBB","HBA2","HBA1",                                #Erythrocyte
                 "CPA3","IL18R1","TPSB2","TPSAB1",                   #MAST
                 "CD163","MSR1","CSF1R","CD14",                      #MAC
                 "MRC1","F13A1","STAB1","CD163L1","LYVE1",           #resMAC-LYVE1+
                 "HLA-DPA1","C1QA","C1QB",                           #resMAC-HLAIIhi
                 "HIF1A","NAMPT","PLAUR","ITGAX","HBEGF","OSM",      #moMAC-HBEGF+
                 "PSTPIP2","CXCL10","CXCL9","CCL2","CCL3","IL1B",    #moMAC-CXCL10+
                 "GPNMB","SPP1","APOC1","PLA2G7","CD68","CAPG",      #moFAM
                 "HMOX1","TREM2",                                    #moFAM
                 "C3","KCNQ3","ADGRB3","VASH1","CX3CR1",             #moMAC-C3+
                 "CLEC10A","FCER1A","CD1C",                          #cDC2
                 "TCF7L2","COTL1","FCGR3A",                          #ncMON
                 "FCN1",                                             #MON/ncMON
                 "VCAN","LYZ","CD36",                                #MON
                 "LAMP3","SLCO5A1","CCR7",#"EBI3","CCL19",           #mDC
                 "WDFY4","CADM1","CLEC9A","BATF3",                   #cDC1
                 "BCL11A","CLEC4C","IL3RA","PLD4",#"LILRA4",         #pDC
                 "S100A9","FCGR3B","S100A8","IFITM2",                #N
                 "TOP2A","MKI67",                                    #cycling
                 "NRXN1","GRIK2","CDH19","NCAM2"                     #SC/NEU
                 
)


Idents(nano.obj) <- "pagoda_k100_infomap.imm"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))

DotPlot(nano.obj, features = unique(imm.markers[imm.markers %in% rownames(nano.obj)]), dot.scale = 8) + RotatedAxis()

table(Idents(nano.obj))

##Add ref atlas correlations
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
var.genes <- VariableFeatures(ref)
Idents(object = ref) <- "v2.subclass.sp"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","aDTL","DTL1","DTL3","ATL",
  "M-TAL","C-TAL","MD","frTAL","aTAL1","aTAL2","DCT","aDCT","CNT","aCNT","PC",
  "IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
  "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM","M-FIB","C/M-FIB",
  "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
  "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
  "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
  "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N","-"))

ref <- subset(ref, idents = c("B","PL","T","MAIT","ILC3",
                              "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
                              "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N"))
ref.ds <- subset(ref, downsample = 1000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")
ref.ds <- NormalizeData(ref.ds)

var.genes <- imm.od.genes[imm.od.genes %in% rownames(nano.obj)]
ref.ds <- ScaleData(ref.ds, features = var.genes, assay = "RNA")

nano.obj <- NormalizeData(nano.obj)
nano.obj <- ScaleData(nano.obj, features = var.genes)

Idents(nano.obj) <- "pagoda_k100_infomap.imm"
levels(Idents(object = nano.obj)) <- paste("CL", levels(Idents(object = nano.obj)), sep = "")
ave.nano.obj <- AverageExpression(nano.obj, features = var.genes, assay = "Nanostring", layer = "scale.data")
ave.ref <- AverageExpression(ref.ds, features = var.genes, assays = "RNA", layer = "scale.data")
library("corrplot")
ave.cor <- cor(cbind(as.data.frame(ave.ref),as.data.frame(ave.nano.obj)))
ave.cor
rownames(ave.cor) <- gsub("RNA.","",rownames(ave.cor))
colnames(ave.cor) <- gsub("Nanostring.","",colnames(ave.cor))
ave.cor <- ave.cor[1:23,24:53]

#Generate table
Idents(nano.obj) <- "pagoda_k50_infomap.imm"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxCor <- function(ave.cor, celltype) {
  max.cor <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    if (ct %in% colnames(ave.cor) == FALSE) {
      max.cor <- data.frame(ref.cor.name = c(NA,NA),
                            ref.cor.value = c(NA,NA))
    } else {
      max.cor <- data.frame(ref.cor.name = tail(names(sort(ave.cor[,ct])), 2),
                            ref.cor.value = tail(sort(ave.cor[,ct]), 2))
    } 
    
  }))
  return(max.cor)
}

max.cor <- getMaxCor(ave.cor, celltype)
max.idents$cor.subclass <- max.cor$ref.cor.name
max.idents$cor.value <- max.cor$ref.cor.value

#Write out table
write.table(max.idents, file="nanostring/IMM_Kidney_atlasV2_subclass_level3_overlaps_08-07-2024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

###Cluster annotations (check region info)
#1 - resMAC
#2 - T
#3 - moMAC/DC
#4 - MAST
#5 - ncMON/N
#6 - MAST
#7 - CD8+ TEM/TEMRA/NK
#8 - B
#9 - PL
#10 - moMAC/DC
#11 - PL
#12 - moMAC/DC
#13 - ambiguous (1 cell)
#14 - ambiguous (1 cell)

#K50 clusters
#7 - mDC
#19 - infMAC
#25 - infMAC

###Update experiment metadata
nano.obj@meta.data
colnames(nano.obj@meta.data)[colnames(nano.obj@meta.data) == "orig.ident"] <- "library"
nano.obj@meta.data$library <- gsub("[.]", "_", nano.obj@meta.data$Run_Tissue_name)
nano.obj@meta.data$fov <- paste0(nano.obj@meta.data$library, "_", nano.obj@meta.data$fov)
meta <- nano.obj@meta.data
exp.meta <- read.delim("nanostring/Kidney_AtlasV2_Nanostring_Experiment_Metadata.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level1","region_level2",
         "age_binned","sex","race","tissue_preservation","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$fov, exp.meta$fov)]
}

nano.obj@meta.data <- meta
table(nano.obj$region_level1,nano.obj$pagoda_k100_infomap.imm)

#Add annotations to nano.obj
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_IMM.RDS")
Idents(nano.obj) <- "pagoda_k100_infomap.imm"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
cls <- c(1:14)
new.cls <- c("resMAC","T","moMAC/DC","MAST","ncMON/N","MAST","CD8+ TEM/TEMRA/NK",
  "B","PL","moMAC/DC","PL","moMAC/DC","ambiguous","ambiguous")

Idents(nano.obj) <- plyr::mapvalues(Idents(nano.obj), from = cls, to = new.cls)
Idents(nano.obj) <- factor(Idents(nano.obj), levels = c("B","PL","T","CD8+ TEM/TEMRA/NK","MAST","resMAC","moMAC/DC","ncMON/N",
                                                        "ambiguous"))

table(Idents(nano.obj))
nano.obj$v2.subclass.l3 <- Idents(nano.obj)

v2.scl3.cols.mod <- v2.scl3.cols
names(v2.scl3.cols.mod) <- gsub("resMAC-LYVE1+","resMAC",names(v2.scl3.cols.mod))
names(v2.scl3.cols.mod) <- gsub("MON","moMAC/DC",names(v2.scl3.cols.mod))
names(v2.scl3.cols.mod) <- gsub("T-REG","T",names(v2.scl3.cols.mod))
names(v2.scl3.cols.mod) <- gsub("NK","CD8+ TEM/TEMRA/NK",names(v2.scl3.cols.mod))

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.4,raster=FALSE,
        group.by = "v2.subclass.l3", repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl3.cols.mod[levels(Idents(object = nano.obj))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))


#K50 clusters
#7 - mDC
#19 - moMAC-INF
#25 - moMAC-INF


moMAC.INF <- rownames(nano.obj@meta.data[nano.obj@meta.data$pagoda_k50_infomap.imm %in% c("19","25"),])
mDC <- rownames(nano.obj@meta.data[nano.obj@meta.data$pagoda_k50_infomap.imm %in% c("7"),])

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.4,raster=FALSE,
        cells.highlight = moMAC.INF) 
DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.4,raster=FALSE,
        cells.highlight = mDC) 

nano.obj <- SetIdent(nano.obj, cells = moMAC.INF, value = "moMAC-INF")
nano.obj <- SetIdent(nano.obj, cells = mDC, value = "mDC")
nano.obj$v2.subclass.l3 <- Idents(nano.obj)

Idents(nano.obj) <- factor(Idents(nano.obj), levels = c("B","PL","T","CD8+ TEM/TEMRA/NK","MAST","resMAC","moMAC/DC","moMAC-INF","ncMON/N","mDC",
                                                        "ambiguous"))

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.4,raster=FALSE,
        group.by = "v2.subclass.l3", repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl3.cols.mod[levels(Idents(object = nano.obj))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_IMM_anno.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_IMM_anno.RDS")
