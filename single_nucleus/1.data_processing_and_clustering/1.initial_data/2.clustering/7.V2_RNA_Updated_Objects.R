###Update metadata for full set (including ambiquous low quality)

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

options(future.globals.maxSize = 1e9)
setwd("~/data/net/home/Human_Kidney/Atlas_V2")
source("misc/utils.R")

options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")

load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")

###Update V2 annotations
post.qc <- c(rownames(KB.IMM@meta.data),
             rownames(KB.STR@meta.data),
             rownames(KB.EC@meta.data),
             rownames(KB.PT@meta.data),
             rownames(KB.DT@meta.data))

KB.IMM@meta.data$v2.clusters <- paste0("I_",KB.IMM@meta.data$imm.v2.clusters)
KB.STR@meta.data$v2.clusters <- paste0("S_",KB.STR@meta.data$str.v2.clusters)
KB.EC@meta.data$v2.clusters <- paste0("E_",KB.EC@meta.data$ec.v2.clusters)
KB.PT@meta.data$v2.clusters <- paste0("P_",KB.PT@meta.data$pt.v2.clusters)
KB.DT@meta.data$v2.clusters <- paste0("D_",KB.DT@meta.data$dt.v2.clusters)

colnames.use <- c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt",
                  "set","v1_barcode","source","assay","experiment_long","experiment",
                  "patient","specimen","condition_level3","condition_level2",           
                  "condition_level1","condition","percent_cortex","percent_medulla",
                  "region_level2","region_level1","age_binned","sex","race","location",
                  "laterality","protocol","tissue_type_full","tissue_type","atlas_version",
                  "predicted.class.score","predicted.class","predicted.subclass.l1.score",
                  "predicted.subclass.l1","predicted.subclass.l2.score","predicted.subclass.l2",
                  "predicted.subclass.l3.score","predicted.subclass.l3",
                  "v1clusters","v1.subclass.l3","v1.subclass.l2","v1.subclass.l1","v1.state.l1",
                  "v1.state.l2","v1.class","v1.structure","v2.clusters")

###Create combined metadata tables
neu <- WhichCells(KB, idents = 100)
low.qc <- rownames(KB@meta.data)[!rownames(KB@meta.data) %in% c(post.qc,neu)]
KB@meta.data$v2.clusters <- "alq"
neu.meta <- KB@meta.data[neu,colnames.use]
neu.meta$v2.clusters <- "N_1"
low.qc.meta <- KB@meta.data[low.qc,colnames.use]
combined.meta <- rbind(KB.IMM@meta.data[,colnames.use],KB.STR@meta.data[,colnames.use],KB.EC@meta.data[,colnames.use],
                       KB.PT@meta.data[,colnames.use],KB.DT@meta.data[,colnames.use],
                       neu.meta,low.qc.meta)
dim(combined.meta)
dim(KB@meta.data)

KB@meta.data <- combined.meta[rownames(KB@meta.data),]


table(KB$v2.clusters)

Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = c(paste0("P_",1:33),
                                            paste0("D_",1:48),
                                            paste0("E_",1:22),
                                            paste0("S_",1:27),
                                            paste0("I_",1:28),
                                            paste0("N_",1),
                                            "alq"))

###Remove failed libraries
KB <- subset(KB, library %in% c("KM150", "KB54"), invert = TRUE)
KB
#36588 features across 1475982 samples

table(Idents(KB))


###update metadata
meta <- KB@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
colnames(exp.meta)
meta$KDIGO_stage <- "NA"
meta$baseline_eGFR_binned <- "NA"
meta$proteinuria._binned <- "NA"
meta$A1c_binned <- "NA"
meta$albuminuria_binned <- "NA"
meta$diabetes_history <- "NA"
meta$diabetes_duration <- "NA"
meta$hypertension_history <- "NA"
meta$hypertension_duration <- "NA"
meta$on_RAAS_blockade <- "NA"
meta$publication_status <- "NA"

emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race",
         "KDIGO_stage","baseline_eGFR_binned","proteinuria._binned","A1c_binned","albuminuria_binned",   
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration","on_RAAS_blockade",
         "location","laterality","protocol",
         "tissue_type_full","tissue_type","publication_status","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

colnames(meta)
order <- c(
          "library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","v1_barcode",
          "source","assay","experiment_long","experiment","patient","specimen",
           "condition_level3","condition_level2","condition_level1","condition",
           "percent_cortex","percent_medulla","region_level2","region_level1",
           "age","age_binned","sex","race",
           "KDIGO_stage","baseline_eGFR_binned","proteinuria._binned","A1c_binned","albuminuria_binned",   
           "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration","on_RAAS_blockade",
           "location","laterality","protocol",
           "tissue_type_full","tissue_type","publication_status","atlas_version","v1clusters","v1.subclass.l3",             
          "v1.subclass.l2","v1.subclass.l1","v1.state.l1","v1.state.l2","v1.class","v1.structure","v2.clusters")

meta <- meta[,order]
KB@meta.data <- meta


###Add v2 subclasses
###update metadata
meta <- KB@meta.data
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata.txt")

emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}

colnames(meta)
KB@meta.data <- meta

saveRDS(meta, file = "Kidney_AtlasV2_filtered_subset_07132023_metadata.RDS")


###Save object
counts <- KB[["RNA"]]$counts
write_matrix_10x_hdf5(
  counts,
  path = "~/tmp_scratch/Kidney_Atlas_V2_Combined_Counts_07-2023.h5",
  barcodes = colnames(counts),
  feature_ids = rownames(counts),
  feature_names = rownames(counts),
  feature_types = "Gene Expression"
)

###Prepare simplified kidney data objects in scratch directory for fast loading
v2.data <- open_matrix_10x_hdf5(
  path = "~/tmp_scratch/Kidney_Atlas_V2_Combined_Counts_07-2023.h5"
)

# Write the matrix to a directory
write_matrix_dir(
  mat = v2.data,
  dir = "~/tmp_scratch/v2_counts",
  overwrite = TRUE
)

v2.mat <- open_matrix_dir(dir = "~/tmp_scratch/v2_counts")

ref.pca <- Embeddings(object = KB[["ref.pca"]])
ref.umap <- Embeddings(object = KB[["ref.umap"]])
umap <- Embeddings(object = KB[["umap"]])

KB <- CreateSeuratObject(counts = v2.mat, meta.data = meta)

KB[["v1.pca"]] <- CreateDimReducObject(embeddings = ref.pca, key = "v1pca_", assay = DefaultAssay(KB))
KB[["v1.umap"]] <- CreateDimReducObject(embeddings = ref.umap, key = "v1umap_", assay = DefaultAssay(KB))
KB[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "umap_", assay = DefaultAssay(KB))
KB
#An object of class Seurat 
#36588 features across 1475982 samples within 1 assay 
#Active assay: RNA (36588 features, 0 variable features)
#1 layer present: counts

saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_07-2023_Object.Rds",
  destdir = "~/tmp_scratch/v2_object"
)





###Prepare Objects and UMAP Visualization for Annotated Subset (excluding ambiguous low quality)

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

options(future.globals.maxSize = 1e9)

setwd("~/data/net/home/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/tmp_scratch/v2_object/Kidney_Atlas_V2_07-2023_Object.Rds")
load("color_factors.robj")

Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = c(paste0("P_",1:33),
                                            paste0("D_",1:48),
                                            paste0("E_",1:22),
                                            paste0("S_",1:27),
                                            paste0("I_",1:28),
                                            paste0("N_",1),
                                            "alq"))


KB.sub <- subset(KB, v2.clusters %in% "alq", invert = TRUE)
KB.sub
#36588 features across 1058829 samples within 1 assay



##Create a sketch data set
KB.sub[["RNA"]]$data <- NormalizeData(KB.sub[["RNA"]]$counts)
KB.sub <- FindVariableFeatures(KB.sub)
KB.sub <- SketchData(
  object = KB.sub,
  ncells = 100000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
KB.sub

DefaultAssay(KB.sub) <- "sketch"


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.sub[["sketch"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#4376 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 4376, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:4376]
VariableFeatures(KB.sub) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.sub[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.sub))
KB.sub <- RunUMAP(object = KB.sub, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.2, return.model = T)

DimPlot(KB.sub, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l1", repel = TRUE) + NoLegend()
DimPlot(KB.sub, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()

#Add pc loadings
pca.loadings <- p2$misc$PCA$v[sn.od.genes,]
colnames(pca.loadings) <- paste("PC_", 1:50, sep = "")
KB.sub@reductions$pca@feature.loadings <- pca.loadings
Embeddings(KB.sub, reduction = "pca")

#Extend results to the full datasets
KB.sub <- ProjectData(
  object = KB.sub,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50)

DefaultAssay(KB.sub) <- "RNA"
DimPlot(KB.sub, reduction = "full.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l1", repel = TRUE) + NoLegend()
DimPlot(KB.sub, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()
DimPlot(KB.sub, reduction = "v1.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()


###Save to new filtered subset folder
# Write the counts layer to a directory
write_matrix_dir(mat = KB.sub[["RNA"]]$counts, dir = "~/tmp_scratch/v2_counts_filtered")
counts.mat <- open_matrix_dir(dir = "~/tmp_scratch/v2_counts_filtered")

KB.sub[["RNA"]]$counts <- counts.mat

saveRDS(
  object = KB.sub,
  file = "Kidney_Atlas_V2_07-2023_Object_Filtered.Rds",
  destdir = "~/tmp_scratch/v2_object_filtered"
)



###Make V1 Published set for sharing across Altos
KB.sub <- readRDS("~/tmp_scratch/v2_object_filtered/Kidney_Atlas_V2_07-2023_Object_Filtered.Rds")
table(KB.sub$library)
table(KB.sub$assay)
table(KB.sub$atlas_version)
table(KB.sub$publication_status)

KB.pub <- subset(KB.sub, publication_status %in% c("altos", "published"))
KB.pub
#An object of class Seurat 
#73176 features across 482319 samples within 2 assays 
#Active assay: RNA (36588 features, 2000 variable features)
#2 layers present: counts, data
#1 other assay present: sketch
#6 dimensional reductions calculated: v1.pca, umap, v1.umap, pca, pca.full, full.umap

KB.pub[["sketch"]] <- NULL
KB.pub[["umap"]] <- KB.pub[["v1.umap"]]
KB.pub[["pca"]] <- KB.pub[["v1.pca"]]
KB.pub[["v1.umap"]] <- NULL
KB.pub[["v1.pca"]] <- NULL
KB.pub[["pca.full"]] <- NULL
KB.pub[["full.umap"]] <- NULL
KB.pub

DimPlot(KB.pub, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()

to.use <- c("library","nCount_RNA","nFeature_RNA","percent.er",           
  "percent.mt","v1_barcode","source","assay","experiment_long",      
  "experiment","patient","specimen","condition_level3","condition_level2",     
  "condition_level1","condition","percent_cortex","percent_medulla","region_level2",        
  "region_level1","age_binned","sex","race",                 
  "KDIGO_stage","baseline_eGFR_binned","proteinuria._binned","A1c_binned","albuminuria_binned",   
  "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration","on_RAAS_blockade",     
  "location","laterality","protocol","tissue_type_full","tissue_type",          
  "publication_status","atlas_version","v1clusters","v1.subclass.l3","v1.subclass.l2",       
  "v1.subclass.l1","v1.state.l1","v1.state.l2","v1.class","v1.structure")
KB.pub@meta.data <- KB.pub@meta.data[,to.use]

#Convert to anndata
KB.pub[["RNA"]] <- as(object = KB.pub[["RNA"]], Class = "Assay")
expression_matrix <- KB.pub[["RNA"]]$counts
umap <- Embeddings(KB.pub, reduction = "umap")
dfobs <- KB.pub@meta.data
dfvar <- KB.pub@assays$RNA@meta.features
adata <- AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  obsm=list('X_umap'=umap)
)
write_h5ad(adata, "Human_Kidney_Atlas_V1_Expanded_snRNA-AnnData_07-2023.h5ad") #move to ~/data/pod/blake_LTS/Atlas_V2/share/
saveRDS(KB.pub, file = "~/data/pod/blake_LTS/Atlas_V2/share/V1_Data_Browser/Human_Kidney_Atlas_V1_Expanded_snRNA-Seurat_07-2023.Rds")



###Update SC/NEU in object (Needs updating!!!)
setwd("~/data/net/home/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")
source("misc/utils.R")

KB.sub <- readRDS("~/tmp_scratch/v2_object_filtered/Kidney_Atlas_V2_07-2023_Object_Filtered.Rds")
KB.sub
#73176 features across 1058829 samples

nonSC.NEU.cells <- readRDS("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_nonSC-NEU_Cell_subset.rda")
KB.sub <- subset(KB.sub, cells = nonSC.NEU.cells, invert = TRUE)
#73176 features across 1055983 samples within 2 assays

saveRDS(KB.sub, file = "~/tmp_scratch/v2_object_filtered/Kidney_Atlas_V2_07-2023_Object_Filtered_B.Rds")
saveRDS(KB.sub, file = "~/data/pod/blake_LTS/scratch_backup/tmp_scratch/v2_object/Kidney_Atlas_V2_07-2023_Object_Filtered_B.Rds")

                                             
