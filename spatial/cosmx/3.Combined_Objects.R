library(Seurat)
library(SeuratData)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(BPCells)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")
load("color_factors_v2-clusters.robj")


###Merge Metadata
nano.p.epi <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_PT_anno.RDS")
nano.d.epi <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_DT_anno.RDS")
nano.ec <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_EC_anno.RDS")
nano.str <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_STR_anno.RDS")
nano.imm <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_IMM_anno.RDS")

nano.obj <- merge(x = nano.p.epi, y = c(nano.d.epi,nano.ec,nano.str,nano.imm))
table(nano.obj$v2.subclass.l3)

meta <- nano.obj@meta.data

nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined.RDS")
table(nano.obj$v2.subclass.l1)
nano.obj$v2.subclass.l3 <- as.character(nano.obj$v2.subclass.l1)
nano.obj@meta.data[rownames(meta),]$v2.subclass.l3 <- meta$v2.subclass.l3



DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "v2.subclass.l1") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxCellType.l3") + NoLegend()

#remove ambiguous (50 cells)
nano.obj <- subset(nano.obj, v2.subclass.l3 %in% "ambiguous", invert = TRUE)
nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

#Update annotations
meta <- nano.obj@meta.data
meta$l3 <- meta$v2.subclass.l3
cl.meta <- read.delim("nanostring/nanostring_V2_annotations.txt")

emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$l3, cl.meta$l3)]
}
colnames(meta)
table(meta$v2.subclass.l3)
nano.obj@meta.data <- meta




###Update experiment metadata
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

order <- c("library","nCount_Nanostring","nFeature_Nanostring",    
           "maxWeight.l1","maxWeight.l1.ct","barcode",                
           "cell_ID","fov","Area","AspectRatio","Width","Height",
           "source","assay","experiment","patient","specimen","condition_level2",
           "condition_level1","condition","region_level1","region_level2",
           "age_binned","sex","race","tissue_preservation","tissue_type","atlas_version",
           "Mean.CD298","Max.CD298","Mean.PanCK","Max.PanCK","Mean.CD45","Max.CD45",
           "Mean.CD3","Max.CD3","Mean.DAPI","Max.DAPI",
           "maxWeight.l3","maxCellType.l3","v2.subclass.full",
           "v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2",
           "v2.state.l1","v2.class","v2.structure"
)

nano.obj@meta.data <- meta[,order]

##Save Filtered subset
saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated.RDS")
nano.obj
#An object of class Seurat 
#1098 features across 68715 samples within 3 assays 
#Active assay: Nanostring (1007 features, 0 variable features)
#2 layers present: counts, data
#2 other assays present: l1.predictions, l3.predictions
#2 dimensional reductions calculated: pca, umap
#2 spatial fields of view present: R5341.S1 R5446.S2



###Re-create Seurat Object with independent FOVs
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated.RDS")
fovs <- unique(nano.obj$fov)
nano.list <- lapply(fovs, function(x) {
  print(paste("Running for FOV:", x))
  
  Idents(nano.obj) <- "fov"
  cells <- WhichCells(nano.obj, idents = x)
  counts <- nano.obj[["Nanostring"]]$counts[,cells]
  meta <- nano.obj@meta.data[cells,]
  fov <- gsub("_",".",unique(meta$library))
  coords <- GetTissueCoordinates(nano.obj[[fov]], which = "centroids")
  coords <- coords[coords$cell %in% cells,]
  seg <- GetTissueCoordinates(nano.obj[[fov]], which = "segmentation")
  seg <- seg[seg$cell %in% cells,]
  mol <- GetTissueCoordinates(nano.obj[[fov]], which = "molecules")
  mol <- mol[mol$x < max(seg$x) & mol$x > min(seg$x) & mol$y < max(seg$y) & mol$y > min(seg$y), ]
  colnames(mol) <- c("x","y","gene")
  mol <- CreateMolecules(mol)
  
  fov_data <- list(segmentation = CreateSegmentation(seg),
                   centroids = CreateCentroids(coords))
  coords <- CreateFOV(
    coords = fov_data,
    type = c("segmentation", "centroids"),
    molecules = mol,
    assay = "Nanostring")
  obj <- CreateSeuratObject(counts = counts, assay = "Nanostring", meta = meta)
  fov <- gsub("_",".",x)
  obj@images[[fov]] <- coords
  
  #Add prediction weights
  l1.predictions <- nano.obj[["l1.predictions"]]$counts[,cells]
  l3.predictions <- nano.obj[["l3.predictions"]]$counts[,cells]
  obj[["l1.predictions"]] <- CreateAssayObject(counts = l1.predictions)
  obj[["l3.predictions"]] <- CreateAssayObject(counts = l3.predictions)
  return(obj)
})

names(nano.list) <- fovs
new.nano.obj <- merge(x = nano.list[[1]],
                      y = c(nano.list[[2]],nano.list[[3]],nano.list[[4]],nano.list[[5]],nano.list[[6]],nano.list[[7]],nano.list[[8]],nano.list[[9]],
                            nano.list[[10]],nano.list[[11]],nano.list[[12]],nano.list[[13]],nano.list[[14]],nano.list[[15]],nano.list[[16]],nano.list[[17]],nano.list[[18]],nano.list[[19]],
                            nano.list[[20]],nano.list[[21]],nano.list[[22]],nano.list[[23]],nano.list[[24]],nano.list[[25]],nano.list[[26]],nano.list[[27]],nano.list[[28]],nano.list[[29]],  
                            nano.list[[30]],nano.list[[31]],nano.list[[32]],nano.list[[33]],nano.list[[34]],nano.list[[35]],nano.list[[36]],nano.list[[37]],nano.list[[38]],nano.list[[39]],
                            nano.list[[40]],nano.list[[41]]))

new.nano.obj[["Nanostring"]] <- JoinLayers(new.nano.obj[["Nanostring"]])

saveRDS(new.nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated_FOVs.RDS")
new.nano.obj
#An object of class Seurat 
#1098 features across 68715 samples within 3 assays 
#Active assay: Nanostring (1007 features, 0 variable features)
#1 layer present: counts
#2 other assays present: l1.predictions, l3.predictions
#41 spatial fields of view present:
