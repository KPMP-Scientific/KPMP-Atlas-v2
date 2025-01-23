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

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

# Generate a combined seurat object for all pucks
dir <- "~/hsKidAt/spatial_slideseq/processed/HuBMAP_batch3_standard_deconv/"
pucks <- list.files(dir)
pucks <- pucks[8:45]
puck.path <- paste0(dir, pucks)
dir2 <- "~/hsKidAt/spatial_slideseq/processed/HuBMAP_batch1_and_batch2/"
pucks.2 <- list.files(dir2)
puck.path.2 <- paste0(dir2, pucks.2)

pucks <- c(pucks, pucks.2)
puck.path <- c(puck.path, puck.path.2)
names(puck.path) <- pucks

ssg.list <- lapply(pucks, function(x) {
  print(paste("Running for Puck:", x))
  
  path <- as.character(puck.path[x])
  
  expr <- readMM(paste(path,"raw_counts.mtx", sep = "/"))
  positions <- data.frame(read.delim(paste(path,"spatial_locs.tsv", sep = "/"), row.names = 3))
  cell.meta <- data.table(read.delim(paste(path,"cell_metadata.tsv", sep = "/")))
  gene.meta <- data.table(read.delim(paste(path,"gene_metadata.tsv", sep = "/")))
  cell.meta <- as.data.frame(cell.meta)
  rownames(cell.meta) <- cell.meta$cell_ID
  rownames(x = expr) <- toupper(gene.meta$gene_ID)
  colnames(x = expr) <- cell.meta$cell_ID
  rownames(x = expr) <- gsub("_","-", rownames(x = expr))
  expr <- expr[!rownames(expr) %in% rownames(x = expr)[duplicated(rownames(x = expr))],]
  
  expr <- as(expr, "CsparseMatrix")
  expr <- convert_matrix_type(expr, type = "uint32_t")
  write_matrix_dir(
    mat = expr, 
    dir = paste0("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/bpcells_counts_",x),
    overwrite = TRUE)
  expr <- open_matrix_dir(dir = paste0("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/bpcells_counts_",x))                 
  
  obj <- CreateSeuratObject(
    counts = expr,
    project = 'SlideSeq',
    assay = 'Spatial',
    meta.data = cell.meta
  )
  positions <- na.omit(positions)
  obj <- subset(obj, cells = rownames(positions))
  obj[["image"]] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = positions
  )
  obj <- FilterSlideSeq(object = obj, radius = 2450, do.plot = FALSE)
  #obj <- subset(x = obj, subset = nCount_Spatial > 50)
  
  return(obj) 
})

kss <- merge(x = ssg.list[[1]],
             y = c(ssg.list[[2]],ssg.list[[3]],ssg.list[[4]],ssg.list[[5]],ssg.list[[6]],ssg.list[[7]],ssg.list[[8]],ssg.list[[9]],
                   ssg.list[[10]],ssg.list[[11]],ssg.list[[12]],ssg.list[[13]],ssg.list[[14]],ssg.list[[15]],ssg.list[[16]],ssg.list[[17]],ssg.list[[18]],ssg.list[[19]],
                   ssg.list[[20]],ssg.list[[21]],ssg.list[[22]],ssg.list[[23]],ssg.list[[24]],ssg.list[[25]],ssg.list[[26]],ssg.list[[27]],ssg.list[[28]],ssg.list[[29]],
                   ssg.list[[30]],ssg.list[[31]],ssg.list[[32]],ssg.list[[33]],ssg.list[[34]],ssg.list[[35]],ssg.list[[36]],ssg.list[[37]],ssg.list[[38]],ssg.list[[39]],
                   ssg.list[[40]],ssg.list[[41]],ssg.list[[42]],ssg.list[[43]],ssg.list[[44]],ssg.list[[45]],ssg.list[[46]],ssg.list[[47]],ssg.list[[48]],ssg.list[[49]],
                   ssg.list[[50]],ssg.list[[51]],ssg.list[[52]],ssg.list[[53]],ssg.list[[54]],ssg.list[[55]],ssg.list[[56]],ssg.list[[57]],ssg.list[[58]],ssg.list[[59]],
                   ssg.list[[60]],ssg.list[[61]],ssg.list[[62]],ssg.list[[63]],ssg.list[[64]],ssg.list[[65]],ssg.list[[66]],ssg.list[[67]],ssg.list[[68]],ssg.list[[69]],
                   ssg.list[[70]],ssg.list[[71]]))  
kss
names(kss@images) <- pucks
kss[["Spatial"]] <- JoinLayers(kss[["Spatial"]])

###Save object
counts <- kss[["Spatial"]]$counts
write_matrix_dir(
  mat = counts,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_counts_08-30-2023",
  overwrite = TRUE
)
counts <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_counts_08-30-2023")
kss[["Spatial"]]$counts <- counts
kss
#An object of class Seurat 
#51531 features across 2754481 samples within 1 assay 
#Active assay: Spatial (51531 features, 0 variable features)
#1 layer present: counts
#71 images present:

saveRDS(
  object = kss,
  file = "Kidney_Slide-seq_Spatial_Atlas_V2_08-2023_Object.Rds",
  destdir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_08-11-2023"
)



###Add in experiment metadata
kss@meta.data
colnames(kss@meta.data)[colnames(kss@meta.data) == "orig.ident"] <- "library"
kss@meta.data$library <- sub("\\-.*", "", kss@meta.data$cell_ID)
kss@meta.data <- kss@meta.data[,c(1,2,79:86)]

meta <- kss@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Slide-seq_Experiment_Metadata.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level1",
         "age","age_binned","sex","race","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
kss@meta.data <- meta
saveRDS(kss,file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_08-11-2023/Kidney_Slide-seq_Spatial_Atlas_V2_08-2023_Object.Rds")
#kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_08-11-2023/Kidney_Slide-seq_Spatial_Atlas_V2_08-2023_Object.Rds")


#Load Reference data
ref <- readRDS("~/hsKidAt/blake_LTS/scratch_backup/tmp_scratch/v2_object/Kidney_Atlas_V2_07-2023_Object_Filtered.Rds")
DefaultAssay(ref) <- "sketch"
var.genes <- VariableFeatures(ref)
DefaultAssay(ref) <- "RNA"



###Spatial Deconvolution using RCTD: subclass.l1
#prepare a list of puck objects
ssg.list <- lapply(pucks, function(x) {
  print(paste("Running for Puck:", x))
  
  path <- as.character(puck.path[x])
  
  positions <- data.frame(read.delim(paste(path,"spatial_locs.tsv", sep = "/"), row.names = 3))
  cell.meta <- data.table(read.delim(paste(path,"cell_metadata.tsv", sep = "/")))
  gene.meta <- data.table(read.delim(paste(path,"gene_metadata.tsv", sep = "/")))
  cell.meta <- as.data.frame(cell.meta)
  rownames(cell.meta) <- cell.meta$cell_ID
  expr <- open_matrix_dir(dir = paste0("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/bpcells_counts_",x))                 
  
  obj <- CreateSeuratObject(
    counts = expr,
    project = 'SlideSeq',
    assay = 'Spatial',
    meta.data = cell.meta
  )
  positions <- na.omit(positions)
  obj <- subset(obj, cells = rownames(positions))
  obj[["image"]] <- new(
    Class = 'SlideSeq',
    assay = "Spatial",
    coordinates = positions
  )
  obj <- FilterSlideSeq(object = obj, radius = 2450, do.plot = FALSE)
  #obj <- subset(x = obj, subset = nCount_Spatial > 50)
  names(obj@images) <- x
  return(obj) 
})
names(ssg.list) <- pucks

##Spatial Deconvolution using RCTD
#Reference data - update ref annotations
meta <- ref@meta.data
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
ref@meta.data <- meta

ref$v2.subclass.l1 <- gsub("/","",ref$v2.subclass.l1) 
ref$v2.subclass.l1 <- gsub("DTL","TL",ref$v2.subclass.l1) 
ref$v2.subclass.l1 <- gsub("ATL","TL",ref$v2.subclass.l1) 

Idents(object = ref) <- "v2.subclass.l1"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT","TL","TAL","DCT","CNT","PC","PapE","IC","EC","FIB",
  "VSMP","Ad","Lymphoid","Myeloid","NEU"))
ref <- subset(ref, idents = c("PapE","Ad","NEU"), invert = TRUE)
ref.ds <- subset(ref, downsample = 10000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")

# extract information to pass to the RCTD Reference function
counts <- ref.ds[["RNA"]]$counts
cluster <- as.factor(ref.ds$v2.subclass.l1)
names(cluster) <- colnames(ref.ds)
nUMI <- ref.ds$nCount_RNA
names(nUMI) <- colnames(ref.ds)
reference <- Reference(counts, cluster, nUMI)

# run RCTD on each puck separately
ssg.list <- lapply(pucks, function(x) {
  print(paste("Running for Puck:", x))
  
  #RCTD method with full mode
  slide.seq <- ssg.list[[x]]
  counts <- slide.seq[["Spatial"]]$counts
  coords <- GetTissueCoordinates(slide.seq)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  RCTD <- create.RCTD(query, reference, max_cores = 60)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  
  #add l1 prediction weights
  rctd.cls <- RCTD@results$weights
  rctd.cls <- 100 * sweep(rctd.cls, 1, rowSums(rctd.cls), "/")
  rctd.cls <- as.data.frame(rctd.cls)
  maxWeight.l1 <- as.numeric(apply(rctd.cls, 1, max))
  names(maxWeight.l1) <- rownames(rctd.cls)
  maxWeight.l1.ct <- colnames(rctd.cls)[apply(rctd.cls, 1, which.max)]
  names(maxWeight.l1.ct) <- rownames(rctd.cls)
  
  rctd.cls <- as.data.frame(rctd.cls)[colnames(slide.seq),]
  rctd.cls[is.na(rctd.cls)] <- 0
  rownames(rctd.cls) <- colnames(slide.seq)
  
  slide.seq[["l1.predictions"]] <- CreateAssayObject(counts = t(as.matrix(rctd.cls)))
  slide.seq$maxWeight.l1 <- as.numeric(maxWeight.l1[rownames(slide.seq@meta.data)])
  slide.seq$maxWeight.l1.ct <- as.character(maxWeight.l1.ct[rownames(slide.seq@meta.data)])
  
  #subset to confidently predicted l1
  slide.seq <- subset(slide.seq, maxWeight.l1 > 30)
  return(slide.seq) 
})

names(ssg.list) <- pucks

save(ssg.list, file = "slide-seq/Intermediate_Objects/slide-seq_puck_full_l1_predictions_Seurat_List_all-pucks_09012023.rda")
#load("slide-seq/Intermediate_Objects/slide-seq_puck_full_l1_predictions_Seurat_List_all-pucks_09012023.rda")

kss <- merge(x = ssg.list[[1]],
             y = c(ssg.list[[2]],ssg.list[[3]],ssg.list[[4]],ssg.list[[5]],ssg.list[[6]],ssg.list[[7]],ssg.list[[8]],ssg.list[[9]],
                   ssg.list[[10]],ssg.list[[11]],ssg.list[[12]],ssg.list[[13]],ssg.list[[14]],ssg.list[[15]],ssg.list[[16]],ssg.list[[17]],ssg.list[[18]],ssg.list[[19]],
                   ssg.list[[20]],ssg.list[[21]],ssg.list[[22]],ssg.list[[23]],ssg.list[[24]],ssg.list[[25]],ssg.list[[26]],ssg.list[[27]],ssg.list[[28]],ssg.list[[29]],
                   ssg.list[[30]],ssg.list[[31]],ssg.list[[32]],ssg.list[[33]],ssg.list[[34]],ssg.list[[35]],ssg.list[[36]],ssg.list[[37]],ssg.list[[38]],ssg.list[[39]],
                   ssg.list[[40]],ssg.list[[41]],ssg.list[[42]],ssg.list[[43]],ssg.list[[44]],ssg.list[[45]],ssg.list[[46]],ssg.list[[47]],ssg.list[[48]],ssg.list[[49]],
                   ssg.list[[50]],ssg.list[[51]],ssg.list[[52]],ssg.list[[53]],ssg.list[[54]],ssg.list[[55]],ssg.list[[56]],ssg.list[[57]],ssg.list[[58]],ssg.list[[59]],
                   ssg.list[[60]],ssg.list[[61]],ssg.list[[62]],ssg.list[[63]],ssg.list[[64]],ssg.list[[65]],ssg.list[[66]],ssg.list[[67]],ssg.list[[68]],ssg.list[[69]],
                   ssg.list[[70]],ssg.list[[71]]))  
kss
names(kss@images) <- pucks
kss[["Spatial"]] <- JoinLayers(kss[["Spatial"]])

###Save object
counts <- kss[["Spatial"]]$counts
write_matrix_dir(
  mat = counts,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_counts_09-02-2023",
  overwrite = TRUE
)
counts <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_counts_09-02-2023")
kss[["Spatial"]]$counts <- counts
kss
#An object of class Seurat 
#51545 features across 1668343 samples within 2 assays 
#Active assay: Spatial (51531 features, 0 variable features)
#1 layer present: counts
#1 other assay present: l1.predictions
#71 images present: 
saveRDS(
  object = kss,
  file = "Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object.Rds",
  destdir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-02-2023"
)


###Clustering on sketched dataset using Pagoda2
kss <- NormalizeData(kss)
kss <- FindVariableFeatures(kss)
kss <- SketchData(
  object = kss,
  ncells = 100000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
kss

countMatrix <- kss[["sketch"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")

# Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[2])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#513 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- var.genes[var.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
kss[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(kss))
kss <- RunUMAP(object = kss, reduction = "pca", dims = 1:50, n.neighbors = 20L,
               min.dist = 0.2, return.model = T)
kss@meta.data$pagoda_k100_infomap <- k100infomap[rownames(kss@meta.data)]

DimPlot(kss, label = T, label.size = 3, reduction = "umap", group.by = "pagoda_k100_infomap") + NoLegend()
DimPlot(kss, label = T, label.size = 3, reduction = "umap", group.by = "maxCelltype.l1") + NoLegend()


###Extend results to full data sets
#Add pc loadings
VariableFeatures(kss) <- sn.od.genes
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("PC_", 1:50, sep = "")
kss@reductions$pca@feature.loadings <- pca.loadings
Embeddings(kss, reduction = "pca")

#Extend results to the full datasets
kss <- ProjectData(
  object = kss,
  assay = "Spatial",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(pagoda_k100_infomap_full = "pagoda_k100_infomap")
)

DefaultAssay(kss) <- "Spatial"
DimPlot(kss, label = T, label.size = 3, reduction = "full.umap", raster=FALSE,
        group.by = "pagoda_k100_infomap_full", alpha = 0.01) + NoLegend()
DimPlot(kss, label = T, label.size = 3, reduction = "full.umap", raster=FALSE,
        group.by = "maxWeight.l1.ct", alpha = 0.01) + NoLegend()

saveRDS(kss, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object.Rds")


###Annotate to subclass.l1
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object.Rds")
DefaultAssay(kss) <- "Spatial"

##Generate table of top cluster cell types
Idents(kss) <- "pagoda_k100_infomap_full"
Idents(kss) <- factor(Idents(kss), levels = 1:length(levels(Idents(kss))))
levels(Idents(kss)) <- paste("CL", levels(Idents(kss)), sep = "")
celltype <- Idents(kss)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$pagoda_k100_infomap))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxWeight.l1.ct))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$pagoda_k100_infomap))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$pagoda_k100_infomap)),
                                query = prop.table(table(seurat.obj.sub$maxWeight.l1.ct))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxWeight.l1.ct)))
    if (is.null(Idents.called$ref.Var1) == TRUE) {
      Idents.called$ref.Var1 <- rep(NA, 2)
      Idents.called$ref.Freq <- rep(NA, 2)
      Idents.called$ref.Total <- rep(NA, 2)
      Idents.called$query.Var1 <- Idents.called$query.Var1
      Idents.called$query.Freq <- Idents.called$query.Freq
      Idents.called$query.Total <- Idents.called$query.Total
      Idents.called <- Idents.called[,c("ref.Var1","ref.Freq","ref.Total","query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("ref cell types called for:", ct))
    }
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$ref.Var1 <- Idents.called$ref.Var1
      Idents.called$ref.Freq <- Idents.called$ref.Freq
      Idents.called$ref.Total <- Idents.called$ref.Total
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("ref.Var1","ref.Freq","ref.Total","query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("ref.cluster", "ref.Freq","ref.Total", "query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(kss, celltype)
max.idents

##Add ref atlas correlations
#Correlate query clusters to reference cell types using variable genes
ref <- readRDS("~/hsKidAt/blake_LTS/scratch_backup/tmp_scratch/v2_object/Kidney_Atlas_V2_07-2023_Object_Filtered.Rds")
DefaultAssay(ref) <- "sketch"
var.genes <- VariableFeatures(ref)
DefaultAssay(ref) <- "RNA"

meta <- ref@meta.data
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
ref@meta.data <- meta

ref$v2.subclass.l1 <- gsub("/","",ref$v2.subclass.l1) 
ref$v2.subclass.l1 <- gsub("DTL","TL",ref$v2.subclass.l1) 
ref$v2.subclass.l1 <- gsub("ATL","TL",ref$v2.subclass.l1) 

Idents(object = ref) <- "v2.subclass.l1"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT","TL","TAL","DCT","CNT","PC","PapE","IC","EC","FIB",
  "VSMP","Ad","Lymphoid","Myeloid","NEU"))
ref <- subset(ref, idents = c("PapE","Ad","NEU"), invert = TRUE)
ref.ds <- subset(ref, downsample = 10000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")
ref.ds <- NormalizeData(ref.ds)
ref.ds <- ScaleData(ref.ds, features = var.genes, assay = "RNA")

DefaultAssay(kss) <- "sketch"
kss[["sketch"]] <- as(object = kss[["sketch"]], Class = "Assay")
kss[["sketch"]] <- NormalizeData(kss[["sketch"]])
kss[["sketch"]] <- ScaleData(kss[["sketch"]], features = sn.od.genes)

Idents(kss) <- "pagoda_k100_infomap"
sn.od.genes <- VariableFeatures(kss)
levels(Idents(object = kss)) <- paste("CL", levels(Idents(object = kss)), sep = "")
ave.kss <- AverageExpression(kss, features = sn.od.genes, assay = "sketch", slot = "scale.data")
ave.ref <- AverageExpression(ref.ds, features = sn.od.genes, assays = "RNA", slot = "scale.data")
library("corrplot")
ave.cor <- cor(cbind(as.data.frame(ave.ref),as.data.frame(ave.kss)))
ave.cor
rownames(ave.cor) <- gsub("RNA.","",rownames(ave.cor))
colnames(ave.cor) <- gsub("sketch.","",colnames(ave.cor))
ave.cor <- ave.cor[1:14,15:66]

#Generate table
Idents(kss) <- "pagoda_k100_infomap"
Idents(kss) <- factor(Idents(kss), levels = 1:length(levels(Idents(kss))))
levels(Idents(kss)) <- paste("CL", levels(Idents(kss)), sep = "")
celltype <- Idents(kss)

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
write.table(max.idents, file="slide-seq/Intermediate_Objects/Kidney_atlasV2_sketch_cluster_level1_overlaps_09-05-2023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


##Check marker genes
epi.markers <- c("PTPRQ","WT1","NPHS2",                                #normal POD
                 "CDKN1C", "SPOCK2",                                   #injured POD
                 "IGFBP2",                                             #injured POD (conserved)
                 "CLDN1", "CFH","ALDH1A2",                             #PEC
                 "LRP2",                                               #PT
                 "SLC5A12","SLC22A6",                                  #S1/S2
                 "PRODH2","SLC5A2","SLC22A8","SLC7A8",                 #S1
                 "SLC34A1","SLC5A10",                                  #S2                                  
                 "SLC5A11","SLC7A13",                                  #S3
                 
                 "ITGB8","CDH6","HAVCR1","VCAM1",                      #AS
                 "CRYAB","TACSTD2","SLC44A5",                          #TL
                 "AQP1", "UNC5D",                                     #DTL1
                 "ADGRL3","ID1",                                      #DTL2
                 #DTL3
                 "AKR1B1","SH3GL3",                                   #DTL3/ATL
                 "PROX1",                                             #ATL
                 
                 "CASR","SLC12A1","UMOD",                              #TAL
                 "PROM1",                                             #ATL-TAL
                 "ITGB6","FGF13",                                      #ATL-TAL
                 "CLDN14", "KCTD16",                        #M-TAL
                 "ANK2",                                    #M-TAL
                 "ESRRB","EGF",                                        #TAL
                 "ENOX1","TMEM207",                              #C-TAL
                 "CLDN16",                                    #C-TAL
                 "NOS1","ROBO2",                                  #MD
                 
                 "SLC12A3","TRPM6",                                #DCT
                 #DCT1
                 "SLC8A1","SCN2A","HSD11B2","CALB1",                                  #CNT
                 "TRPV5",                                                             #DCT2
                 "SCNN1G","SCNN1B",                                #CNT                   
                 
                 "GATA3","AQP2","AQP3",                            #PC
                 "PDE10A",            #C-PC
                 "KCNK13",                  #M-PC
                 "FXYD4",                                   #PC
                 "PHACTR1","SLC14A2","SOX5","PCDH7",  #IMCD
                 "TP63","GPX2", "FXYD3", "KRT5",           #pape
                 "ATP6V0D2",           #IC
                 "SLC4A1","SLC26A7",                               #IC-A
                 #CNT-IC-A
                 #C-IC-A
                 "KIT","AQP6","CALCA",           #M-IC-A
                 "SLC4A9","INSRR")     #IC-B


int.markers <- c("PECAM1","EMCN",               #EC
                 "HECW2","PLAT","ITGA8",         #EC-GC
                 "BTNL9","PALMD","AQP1","TM4SF1",                   #EC-AEA-DVR
                 "SERPINE2",  #EC-AEA
                 "SLC14A1","ENPP2",                          #EC-DVR
                 "DNASE1L3","CEACAM1",                                        #EC-AVR/PTC
                 "PITPNC1","SLCO2A1",                       #EC-PTC
                 "PLVAP","TLL1",   #EC-AVR
                 "MMRN1","PROX1",                     #EC-LYM
                 
                 "PDGFRB",                                          #VSMC/P
                 
                 "ROBO1",                                                  #MC/REN
                 "PIEZO2","POSTN",                           #MC
                 "REN","GRID2",                           #REN
                 
                 "MYH11","RGS6","MCAM",                      #C-VSMC/P               
                 "RGS5","ADGRB3",                         #VSMC/P
                 
                 "COL1A1","COL1A2","C7","DCN",               #FIB/MyoF
                 
                 "SYNPO2","PCDH7",         #MyoF
                 "MEG3","LAMA2",                                                     #FIB/MyoF
                 "COL6A3","GLI2","COL5A1",   #MyoF
                 "SLC24A3","CDH13","SMOC2","ANO3","RXFP1",                           #MyoF
                 
                 "PDGFRA",               #FIB
                 
                 "SYT1",                                                       #M-FIB
                 "PLCXD3","GREB1L",                                 #M-FIB                          
                 "ADAMTSL1","FREM1",                                  #M-FIB
                 "MGP","BGN",                                      #M-FIB
                 
                 
                 "ABI3BP","FLRT2", "FGF14","IGF1",                    #aFIB
                 
                 "TAGLN","ACTA2",                               #iVSMC/P
                 "PTPRC",                                                        #ec
                 "BANK1","MS4A1",                                  #B
                 "IGKC","MZB1",                        #pB
                 
                 "THEMIS","IL7R","CD96","CD247",                #T
                 "GNLY","NKG7","GZMA",             #NKT
                 
                 "KIT", "MS4A2",                                         #Mast
                 
                 "MRC1","CD163","CD14",                 #MAC
                 "DIAPH3","CENPF","MKI67",                      #cycMAC
                 
                 "MSR1",
                 "ITGAX","HLA-DQA1","CSF2RA",                  #cDC
                 "FLT3","CLEC9A",                                       #cDC
                 
                 "IL3RA","CLEC4C",                        #pDC
                 
                 "CTSS","FCN1","FCGR3A",   #MON
                 "S100A9","S100A8","FCGR3B"                     #NC
                 
)


Idents(kss) <- "pagoda_k100_infomap"
Idents(kss) <- factor(Idents(kss), levels = 1:length(levels(Idents(kss))))

kss.sub <- subset(kss, subset = pagoda_k100_infomap %in% c(16,38,43,46,48,50))
DotPlot(kss.sub, features = c("NPHS1","NPHS2","PODXL",                                #POD
                              "CDKN1C","SPOCK2",                                          #dPOD
                              "CLDN1", "CFH","ALDH1A2",                             #PEC
                              
                              "SLC5A12","SLC22A6",                                  #S1/S2
                              "PRODH2","SLC5A2","SLC22A8","SLC7A8",                 #S1
                              "SLC34A1","SLC5A10",                                  #S2                                  
                              "SLC5A11","SLC7A13",                                  #S3
                              
                              "ITGB8","CDH6","HAVCR1","VCAM1",                      #AS
                              
                              "CRYAB","TACSTD2","SLC44A5", "SH3GL3",                #TL
                              "EMCN","PECAM1","PDGFRB","ACTA2","C7","DCN","COL1A1", #EC/stromal
                              "PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
                              "CD96","NKG7",       
                              "MS4A2","CPA3", "KIT","S100A9",
                              "APOE","DEFB1","CST3","GATM","ALDOB",                 #Injury
                              "UMOD","SLC12A1","SLC12A3","SLC8A1","AQP2",            #Distal tubules
                              "HSD11B2","SLC26A7",
                              "MKI67","TOP2A"), dot.scale = 8) + RotatedAxis()


kss.sub <- subset(kss, subset = pagoda_k100_infomap %in% c(1:25))
DotPlot(kss.sub, features = epi.markers, dot.scale = 8) + RotatedAxis()
DotPlot(kss.sub, features = int.markers, dot.scale = 8) + RotatedAxis()

CL16_markers <- FindMarkers(object = kss, ident.1 = 16, 
                            max.cells.per.ident = 500, only.pos = TRUE)

kss.sub <- subset(kss, subset = pagoda_k100_infomap %in% c(26:52))
DotPlot(kss.sub, features = epi.markers, dot.scale = 8) + RotatedAxis()
DotPlot(kss.sub, features = int.markers, dot.scale = 8) + RotatedAxis()

CL35_markers <- FindMarkers(object = kss, ident.1 = 35, 
                            max.cells.per.ident = 500, only.pos = TRUE)

#Checked ambigous clusters showing poor prediction assignments and poor correlations
#Those showing very low DEGs were labeled as ambiguous
#Those showing DEGs were assigned to the closest subclass.l1


###Generate filtered kss object
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object.Rds")
DefaultAssay(kss) <- "Spatial"

Idents(kss) <- "pagoda_k100_infomap_full"
cls <- c(1:52)
new.cls <- c("TL","PC","FIB","TL","Myeloid","EC","EC","Lymphoid","PC","FIB","VSMP",
             "FIB","IC","PC","PC","NA","PC","FIB","NA","NA","NA","TL","NA","EC","NA",
             "NA","PC","NA","NA","NA","Myeloid","NA","NA","NA","NA","DCT","VSMP",
             "POD","TAL","TAL","TAL","FIB","PT","TAL","PapE","PT","Lymphoid","PT",
             "Lymphoid","PT","CNT","Lymphoid")
Idents(kss) <- plyr::mapvalues(Idents(kss), from = cls, to = new.cls)
Idents(kss) <- factor(Idents(kss), levels = c("POD","PEC","PT","TL","TAL","DCT","CNT",
                                              "PC","PapE","IC","EC","FIB","VSMP",
                                              "Lymphoid","Myeloid","NA"))

table(Idents(kss))
kss$v2.subclass.l1 <- Idents(kss)

v2.scl1.cols.mod <- v2.scl1.cols
names(v2.scl1.cols.mod) <- gsub("/","",names(v2.scl1.cols.mod))
names(v2.scl1.cols.mod) <- gsub("DTL","TL",names(v2.scl1.cols.mod))

DimPlot(kss, reduction = "full.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l1", repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl1.cols.mod[levels(Idents(object = kss))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))


###Add in experiment metadata
kss@meta.data
colnames(kss@meta.data)[colnames(kss@meta.data) == "orig.ident"] <- "library"
kss@meta.data$library <- sub("\\-.*", "", kss@meta.data$cell_ID)
kss@meta.data <- kss@meta.data[,c(1,2,79:94)]

meta <- kss@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Slide-seq_Experiment_Metadata.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level1",
         "age","age_binned","sex","race","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
kss@meta.data <- meta


DimPlot(kss, reduction = "full.umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "patient", repel = TRUE) #+ NoLegend()
DimPlot(kss, reduction = "full.umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level1", repel = TRUE) #+ NoLegend()
DimPlot(kss, reduction = "full.umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level1", repel = TRUE) #+ NoLegend()

saveRDS(kss, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object.Rds")



###Re-run global clustering for final UMAP visualization
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object.Rds")
kss <- subset(kss, v2.subclass.l1 %in% "NA", invert = TRUE)

DefaultAssay(kss) <- "sketch"
countMatrix <- kss[["sketch"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[2])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
library(pagoda2)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#466 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- var.genes[var.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

#pagoda2 PCA values for umap
cell.embeddings <- p2$reductions$PCA
kss[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(kss))

kss <- RunUMAP(object = kss, reduction = "pca", dims = 1:50, n.neighbors = 20L,
               min.dist = 0.2, return.model = T)

DimPlot(kss, label = T, label.size = 3, reduction = "umap", group.by = "v2.subclass.l1") + NoLegend()


###Extend results to full data sets
#Add pc loadings
VariableFeatures(kss) <- sn.od.genes
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("PC_", 1:50, sep = "")
kss@reductions$pca@feature.loadings <- pca.loadings
Embeddings(kss, reduction = "pca")


#Extend results to the full datasets
kss <- ProjectData(
  object = kss,
  assay = "Spatial",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(pagoda_k100_infomap_full = "pagoda_k100_infomap",
                 pagoda_k50_infomap_full = "pagoda_k50_infomap")
)

DefaultAssay(kss) <- "Spatial"
DimPlot(kss, label = T, label.size = 3, reduction = "full.umap.1", raster=FALSE,
        group.by = "pagoda_k100_infomap_full", alpha = 0.01) + NoLegend()
DimPlot(kss, label = T, label.size = 3, reduction = "full.umap.1", raster=FALSE,
        group.by = "v2.subclass.l1", alpha = 0.01) + NoLegend()

saveRDS(object = kss,
        file = "Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object_filtered.Rds", 
        destdir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023")
