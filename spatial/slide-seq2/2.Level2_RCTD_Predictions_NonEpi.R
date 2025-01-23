# Slide-Seq - Combined Analyses: nonEpi Level 3 RCTD Predictions --------------------------
library(Seurat)
library(SeuratData)
library(Matrix)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(foreach)
library(BPCells)
library(spacexr)
library(pagoda2)
library("corrplot")


options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")
load("color_factors_v2-clusters.robj")

kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object_filtered.Rds")
DefaultAssay(kss) <- "Spatial"
counts <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_counts_09-02-2023")
kss[["Spatial"]]$counts <- counts[,colnames(kss)]
kss <- NormalizeData(kss)

###Non-Epithelial Sub-clusters
##Load Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
ref <- subset(ref, v2.subclass.l1 %in% c("EC","FIB","VSM/P","Lymphoid","Myeloid"))
ref <- subset(ref, v2.state.l2 %in% c("degenerative","cycling"), invert = TRUE)
ref <- subset(ref, v2.subclass.sp %in% c("-"), invert = TRUE)
Idents(ref) <- "v2.subclass.sp"
table(Idents(ref))
ref <- subset(ref, downsample = 5000)

#extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]$counts
ref$v2.subclass.sp <- gsub("/","",ref$v2.subclass.sp)
cluster <- as.factor(ref$v2.subclass.sp)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

#RCTD method with full mode
kss.sub <- subset(kss, v2.subclass.l1 %in% c("EC","FIB","VSMP","Lymphoid","Myeloid"))
counts <- kss.sub[["Spatial"]]$counts
counts <- as(object = counts, Class = "dgCMatrix")
coords <- data.frame(x = rep(1, length(colnames(counts))),
                     y = rep(2, length(colnames(counts))))
rownames(coords) <- colnames(counts)
query <- SpatialRNA(coords, counts, colSums(counts), use_fake_coords = TRUE)

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
save(RCTD, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Intermediate_Objects/Non-Epithelial_Combined_RCTD_07152024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Intermediate_Objects/Non-Epithelial_Combined_RCTD_07152024.rda")


#add l3 prediction weights
rctd.cls <- RCTD@results$weights
rctd.cls <- 100 * sweep(rctd.cls, 1, rowSums(rctd.cls), "/")
rctd.cls <- as.data.frame(rctd.cls)
maxWeight.l3 <- as.numeric(apply(rctd.cls, 1, max))
names(maxWeight.l3) <- rownames(rctd.cls)
maxWeight.l3.ct <- colnames(rctd.cls)[apply(rctd.cls, 1, which.max)]
names(maxWeight.l3.ct) <- rownames(rctd.cls)

rctd.cls <- as.data.frame(rctd.cls)[colnames(kss.sub),]
rctd.cls[is.na(rctd.cls)] <- 0
rownames(rctd.cls) <- colnames(kss.sub)

kss.sub[["l3.predictions"]] <- CreateAssayObject(counts = t(as.matrix(rctd.cls)))
kss.sub$maxWeight.l3 <- as.numeric(maxWeight.l3[rownames(kss.sub@meta.data)])
kss.sub$maxCellType.l3 <- as.character(maxWeight.l3.ct[rownames(kss.sub@meta.data)])
table(kss.sub$maxCellType.l3)
save(kss.sub, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_nonEpi_subset_07212024.Rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_nonEpi_subset_07212024.Rda")
