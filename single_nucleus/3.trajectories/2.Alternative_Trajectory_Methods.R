# PT-S1 Trajectories -----------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(patchwork)
library(pagoda2)
library(SeuratExtend)
library(monocle3)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Load PT Data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
KB <- KB.PT
rm(KB.PT)
KB <- UpdateSeuratObject(KB)

#update metadata
meta <- KB@meta.data
head(meta)

exp.meta <- read.delim("HKAv2_Experimental_Metadata_12032025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level3","region_level2",
         "region_level1","age_binned","sex","race","KDIGO_stage",
         "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
         "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
         "tissue_type_full","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

KB@meta.data <- meta


#Update annotations
meta <- KB@meta.data
cl.meta <- read.csv("HKAv2_Cluster_Metadata_09232025.csv")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)

meta <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","source","assay","experiment","patient","specimen",
                "condition_level3","condition_level2","condition_level1","condition",
                "percent_cortex","percent_medulla","region_level3","region_level2",
                "region_level1","age_binned","sex","race","KDIGO_stage",
                "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
                "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
                "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
                "tissue_type_full","tissue_type","v2.clusters","v2.subclass.full","v2.subclass.l3",
                "v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1",
                "v2.class","v2.structure")]

KB@meta.data <- meta

#Fix problem characters
unique(KB@meta.data$region_level3)
KB@meta.data[KB@meta.data$region_level3 %in% '\xca',]$region_level3 <- ""
unique(KB@meta.data$region_level3)


#Trim to only biopsy and AKI/CKD samples
KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("AKI","CKD"))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:20, n.neighbors = 20L, spread = 3,
                 min.dist = 0.01)


v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#DC5E93","#A10300"),
                       c("PT-S1","PT-S2","PT-S3",
                         "aPT-S3","aPT-S1/S2","aPT1","aPT2","frPT-S1/S2","frPT-S3"))

Idents(KB) <- "v2.subclass.l3"
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()

saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT_SeuratObject_12172025.rds")

##Save Anndata
library(anndataR)
adata <- as_AnnData(
  KB,
  assay_name = "RNA",
  x_mapping = "counts",
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")
)

adata$write_h5ad(paste0("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT_Subset_12172025.h5ad"), compression="gzip")




###Full PT Trajectory analyses
#Monocle3 Trajectory Analysis
cds <- as.cell_data_set(KB)
cds <- cluster_cells(cds, reduction_method = "UMAP")
my_clusters <- colData(cds)$v2.subclass.l3
names(my_clusters) <- rownames(colData(cds))
# Assign to the specific Monocle3 slot
cds@clusters$UMAP$clusters <- my_clusters

cds <- learn_graph(cds, use_partition = TRUE)
aPT2 <- WhichCells(KB, idents = "aPT2")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = aPT2)


p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)

plot_cells(cds,
           color_cells_by = "v2.subclass.l3",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

save_monocle_objects(cds, "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-full_Monocle3Object_12172025")

#cds <- load_monocle_objects(directory_path='/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-full_Monocle3Object_12172025/')





###PT S1/S2 Subset
S1 <- subset(KB, idents = c("PT-S1","aPT-S1/S2","aPT1","aPT2","frPT-S1/S2"))

#To ensure proper regional states, subset to only cortex samples
table(S1$region_level2)
S1 <- subset(S1, region_level2 %in% "C")

#re-compute embeddings
countMatrix <- S1[["RNA"]]$counts
countMatrix <- as(countMatrix, "dgCMatrix")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#3123 overdispersed genes
p2$calculatePcaReduction(nPcs = 20, n.odgenes = 1500, maxit=1000)
cell.embeddings <- p2$reductions$PCA
S1[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(S1))
S1 <- RunUMAP(object = S1, reduction = "pca", dims = 1:10, n.neighbors = 20L, spread = 3,
                 min.dist = 0.01)
pdf(file='trajectories/PT-S1-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_colors.pdf',width=6,height=5)
DimPlot(S1, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(S1))],
) + NoLegend()
DimPlot(S1, reduction = "umap", pt.size = 0.8, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(S1))],
) + NoLegend()
dev.off()

saveRDS(S1, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S1-2_SeuratObject_12172025.rds")
S1 <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S1-2_SeuratObject_12172025.rds")

##Save Anndata
library(anndataR)
adata <- as_AnnData(
  S1,
  assay_name = "RNA",
  x_mapping = "data",
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")
)

adata$write_h5ad("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S1-2_Subset_12172025.h5ad", compression="gzip")







###PT-S1 Subset Trajectory analyses
#Monocle3 Trajectory Analysis
library(monocle3)
S1.cds <- as.cell_data_set(S1)
S1.cds <- cluster_cells(S1.cds, reduction_method = "UMAP")
my_clusters <- colData(S1.cds)$v2.subclass.l3
names(my_clusters) <- rownames(colData(S1.cds))
# Assign to the specific Monocle3 slot
S1.cds@clusters$UMAP$clusters <- v2.fn.cols

S1.cds <- learn_graph(S1.cds, use_partition = TRUE)
aPT2 <- WhichCells(S1, idents = "aPT2")
S1.cds <- order_cells(S1.cds, reduction_method = "UMAP", root_cells = aPT2)


p1 <- plot_cells(S1.cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(S1.cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)

S1.cds@clusters$UMAP$partitions <- colors

pdf(file='trajectories/PT-S1-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Monocle3_colors.pdf',width=6,height=4)
plot_cells(S1.cds,
           color_cells_by = "v2.subclass.l3",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           alpha = 0.6)
S1.cds <- order_cells(S1.cds)
plot_cells(S1.cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.3,
           alpha = 0.15)
dev.off()
save_monocle_objects(S1.cds, "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S1-2_Monocle3Object_12172025")

#S1.cds <- load_monocle_objects(directory_path='/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S1-2_Monocle3Object_12172025/')


###Running Diffusion Map
S1 <- Palantir.RunDM(S1)
S1@reductions
DimPlot2(S1, reduction = "ms", group.by = "v2.subclass.l3", label = TRUE)


p <- DimPlot(S1, reduction = "ms", group.by = "v2.subclass.l3")
cells <- CellSelector(p)

S1 <- Palantir.Pseudotime(S1, start_cell = "KB5_TGTTGAGGTGCTCTTC-1")
ps <- S1@misc$Palantir$Pseudotime
head(ps)
# Visualize cell fate on UMAP
colnames(ps)[3] <- c("fate1")
S1@meta.data[,colnames(ps)] <- ps

pdf(file='trajectories/PT-S1-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors.pdf',width=5,height=4)
DimPlot2(S1, reduction = "ms", group.by = "v2.subclass.l3", 
         label = FALSE, theme = NoAxes(),
         cols = v2.fn.cols[levels(Idents(S1))]) + NoLegend()
DimPlot2(S1, features = "Pseudotime", reduction = "ms", 
         cols = list(continuous = "A"), theme = NoAxes())
dev.off()

saveRDS(S1, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S1-2_SeuratObject_withDM_12172025.rds")
#S1 <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S1-2_SeuratObject_withDM_12172025.rds")










# PT-S3 Trajectories -----------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(patchwork)
library(pagoda2)
library(SeuratExtend)
library(monocle3)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Load PT Data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
KB <- KB.PT
rm(KB.PT)
KB <- UpdateSeuratObject(KB)

#update metadata
meta <- KB@meta.data
head(meta)

exp.meta <- read.delim("HKAv2_Experimental_Metadata_12032025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level3","region_level2",
         "region_level1","age_binned","sex","race","KDIGO_stage",
         "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
         "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
         "tissue_type_full","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

KB@meta.data <- meta


#Update annotations
meta <- KB@meta.data
cl.meta <- read.csv("HKAv2_Cluster_Metadata_09232025.csv")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)

meta <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","source","assay","experiment","patient","specimen",
                "condition_level3","condition_level2","condition_level1","condition",
                "percent_cortex","percent_medulla","region_level3","region_level2",
                "region_level1","age_binned","sex","race","KDIGO_stage",
                "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
                "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
                "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
                "tissue_type_full","tissue_type","v2.clusters","v2.subclass.full","v2.subclass.l3",
                "v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1",
                "v2.class","v2.structure")]

KB@meta.data <- meta

#Fix problem characters
unique(KB@meta.data$region_level3)
KB@meta.data[KB@meta.data$region_level3 %in% '\xca',]$region_level3 <- ""
unique(KB@meta.data$region_level3)


#Trim to only biopsy and AKI/CKD samples
KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("AKI","CKD"))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:20, n.neighbors = 20L, spread = 3,
              min.dist = 0.01)


v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#DC5E93","#A10300"),
                       c("PT-S1","PT-S2","PT-S3",
                         "aPT-S3","aPT-S1/S2","aPT1","aPT2","frPT-S1/S2","frPT-S3"))

Idents(KB) <- "v2.subclass.l3"
pdf(file='trajectories/PT-full-S1-S3-AKI-CKD-Biopsy-Subset_subclassl3_colors.pdf',width=6,height=4)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()
#saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT_SeuratObject_12172025.rds")



###PT S3 Subset
S3 <- subset(KB, idents = c("PT-S3","aPT-S3","aPT1","aPT2","frPT-S3"))

#remoove cluster 20, apparent outlier for the S3 trajectory
S3 <- subset(S3, v2.clusters %in% "P_20", invert = TRUE)

#To ensure proper regional states, subset to only cortex samples
table(S3$region_level2)
#S3 <- subset(S3, region_level2 %in% "C")

#re-compute embeddings
countMatrix <- S3[["RNA"]]$counts
countMatrix <- as(countMatrix, "dgCMatrix")
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#2575 overdispersed genes
p2$calculatePcaReduction(nPcs = 20, n.odgenes = 800, maxit=1000)
cell.embeddings <- p2$reductions$PCA
S3[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(S3))
S3 <- RunUMAP(object = S3, reduction = "pca", dims = 1:8, n.neighbors = 20L, spread = 3,
              min.dist = 0.01)
pdf(file='trajectories/PT-S3-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_colors.pdf',width=6,height=4)
DimPlot(S3, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(S3))],
) + NoLegend()
DimPlot(S3, reduction = "umap", pt.size = 0.8, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(S3))],
) + NoLegend()
dev.off()

saveRDS(S3, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S3_SeuratObject_12172025.rds")
S3 <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S3_SeuratObject_12172025.rds")

##Save Anndata
library(anndataR)
adata <- as_AnnData(
  S3,
  assay_name = "RNA",
  x_mapping = "data",
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")
)

adata$write_h5ad("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S3_Subset_12172025.h5ad", compression="gzip")



###Running Diffusion Map
S3 <- Palantir.RunDM(S3, n_components = 8)
S3@reductions
DimPlot2(S3, reduction = "ms", group.by = "v2.subclass.l3", label = TRUE)


p <- DimPlot(S3, reduction = "ms", group.by = "v2.subclass.l3")
cells <- CellSelector(p)

S3 <- Palantir.Pseudotime(S3, start_cell = "KBCVD3_TGGGATTTCTCCTGTG-1")
ps <- S3@misc$Palantir$Pseudotime
head(ps)
# Visualize cell fate on UMAP
colnames(ps)[3] <- c("fate1")
S3@meta.data[,colnames(ps)] <- ps

pdf(file='trajectories/PT-S3-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors.pdf',width=5,height=4)
DimPlot2(S3, reduction = "ms", group.by = "v2.subclass.l3", 
         label = FALSE, theme = NoAxes(),
         cols = v2.fn.cols[levels(Idents(S3))]) + NoLegend()
DimPlot2(S3, features = "Pseudotime", reduction = "ms", 
         cols = list(continuous = "A"), theme = NoAxes())
dev.off()

saveRDS(S3, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S3_SeuratObject_withDM_12172025.rds")
#S3 <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aPT-PT-S3_SeuratObject_withDM_12172025.rds")







# TAL Trajectories -----------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(pagoda2)
library(SeuratExtend)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Load TAL Data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_0424-newData_TALA.rda")
KB <- KB.TAL
rm(KB.TAL)
KB <- UpdateSeuratObject(KB)

#update metadata
meta <- KB@meta.data
head(meta)

exp.meta <- read.delim("HKAv2_Experimental_Metadata_12032025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level3","region_level2",
         "region_level1","age_binned","sex","race","KDIGO_stage",
         "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
         "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
         "tissue_type_full","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

KB@meta.data <- meta


#Update annotations
meta <- KB@meta.data
cl.meta <- read.csv("HKAv2_Cluster_Metadata_09232025.csv")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)

meta <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","source","assay","experiment","patient","specimen",
                "condition_level3","condition_level2","condition_level1","condition",
                "percent_cortex","percent_medulla","region_level3","region_level2",
                "region_level1","age_binned","sex","race","KDIGO_stage",
                "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
                "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
                "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
                "tissue_type_full","tissue_type","v2.clusters","v2.subclass.full","v2.subclass.l3",
                "v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1",
                "v2.class","v2.structure")]

KB@meta.data <- meta

#Fix problem characters
unique(KB@meta.data$region_level3)
KB@meta.data[KB@meta.data$region_level3 %in% '\xca',]$region_level3 <- ""
unique(KB@meta.data$region_level3)


#Trim to only biopsy and AKI/CKD samples
KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("AKI","CKD"))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:10, n.neighbors = 20L, spread = 3,
              min.dist = 0.01)


v2.fn.cols <- setNames(c("#005300","#02AD24","#005300","#02AD24","#009FFF","#14F9FF","#14F9FF",
                         "#DC5E93","#FFD300","#FE8F42"),
                       c("C/M-TAL-A","C-TAL-A","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","C-TAL-B",
                         "frTAL","aTAL1","aTAL2"))

Idents(KB) <- "v2.subclass.l3"
pdf(file='trajectories/C-TAL-A-Lineage-AKI-CKD-Biopsy-Subset_subclassl3_colors.pdf',width=6,height=5)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aTAL-TALA_SeuratObject_12172025.rds")
KB <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aTAL-TALA_SeuratObject_12172025.rds")

##Save Anndata
library(anndataR)
adata <- as_AnnData(
  KB,
  assay_name = "RNA",
  x_mapping = "counts",
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")
)

adata$write_h5ad(paste0("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aTAL-TALA_Subset_12172025.h5ad"), compression="gzip")




###TAL Trajectory analyses
#Monocle3 Trajectory Analysis
cds <- as.cell_data_set(KB)
cds <- cluster_cells(cds, reduction_method = "UMAP")
my_clusters <- colData(cds)$v2.subclass.l3
names(my_clusters) <- rownames(colData(cds))
# Assign to the specific Monocle3 slot
cds@clusters$UMAP$clusters <- my_clusters

cds <- learn_graph(cds, use_partition = TRUE)
aTAL1 <- WhichCells(KB, idents = "aTAL1")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = aTAL1)


p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)

plot_cells(cds,
           color_cells_by = "v2.subclass.l3",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


pdf(file='trajectories/aTAL-TALA-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Monocle3_colors.pdf',width=6,height=4)
plot_cells(cds,
           color_cells_by = "v2.subclass.l3",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           alpha = 0.6)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.3,
           alpha = 0.15)
dev.off()



save_monocle_objects(cds, "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aTAL-TALA-full_Monocle3Object_12172025")

#cds <- load_monocle_objects(directory_path='/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aTAL-TALA-full_Monocle3Object_12172025/')







###Running Diffusion Map
KB <- Palantir.RunDM(KB)
KB@reductions
DimPlot2(KB, reduction = "ms", group.by = "v2.subclass.l3", label = TRUE)


p <- DimPlot(KB, reduction = "ms", group.by = "v2.subclass.l3")
cells <- CellSelector(p)

KB <- Palantir.Pseudotime(KB, start_cell = "KB8_ACGGTTATCTTTCTAG-1")
ps <- KB@misc$Palantir$Pseudotime
head(ps)
# Visualize cell fate on UMAP
colnames(ps)[3:4] <- c("fate1","fate2")
KB@meta.data[,colnames(ps)] <- ps
DimPlot2(KB, features = colnames(ps), reduction = "ms", 
         cols = list(continuous = "A", Entropy = "D"), theme = NoAxes())

pdf(file='trajectories/C-TAL-A-Lineage-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors.pdf',width=5,height=4)
DimPlot2(KB, reduction = "ms", group.by = "v2.subclass.l3", 
         label = FALSE, theme = NoAxes(),
         cols = v2.fn.cols[levels(Idents(KB))]) + NoLegend()
DimPlot2(KB, features = "Pseudotime", reduction = "ms", 
         cols = list(continuous = "A"), theme = NoAxes())
dev.off()


saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aTAL-TALA_SeuratObject_withDM_12172025.rds")
KB <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_aTAL-TALA_SeuratObject_withDM_12172025.rds")








# Interstitial Fibroblast Trajectories -----------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(pagoda2)
library(SeuratExtend)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Load Interstitial Fibroblast Data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_0424-newData.rda")
KB <- KB.STR
rm(KB.STR)
KB <- UpdateSeuratObject(KB)

#update metadata
meta <- KB@meta.data
head(meta)

exp.meta <- read.delim("HKAv2_Experimental_Metadata_12032025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level3","region_level2",
         "region_level1","age_binned","sex","race","KDIGO_stage",
         "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
         "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
         "tissue_type_full","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

KB@meta.data <- meta


#Update annotations
meta <- KB@meta.data
cl.meta <- read.csv("HKAv2_Cluster_Metadata_09232025.csv")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)

meta <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","source","assay","experiment","patient","specimen",
                "condition_level3","condition_level2","condition_level1","condition",
                "percent_cortex","percent_medulla","region_level3","region_level2",
                "region_level1","age_binned","sex","race","KDIGO_stage",
                "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
                "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
                "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
                "tissue_type_full","tissue_type","v2.clusters","v2.subclass.full","v2.subclass.l3",
                "v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1",
                "v2.class","v2.structure")]

KB@meta.data <- meta

#Fix problem characters
unique(KB@meta.data$region_level3)
KB@meta.data[KB@meta.data$region_level3 %in% '\xca',]$region_level3 <- ""
unique(KB@meta.data$region_level3)

KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("AKI","CKD"))


KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:20, n.neighbors = 20L, spread = 3,
              min.dist = 0.01)


v2.fn.cols <- setNames(adjustcolor(c("#02AD24","#009FFF","#00479E","#FF0000"), alpha.f = 0.3),
                       c("C-FIB","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF"))

KB$v2.subclass.l3[KB$v2.subclass.l3 %in% "C-FIB-PATH"] <- "C-FIB"
Idents(KB) <- "v2.subclass.l3"

pdf(file='trajectories/Int-FIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_colors.pdf',width=6,height=4)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_intFIB_SeuratObject_12172025.rds")

##Save Anndata
library(anndataR)
adata <- as_AnnData(
  KB,
  assay_name = "RNA",
  x_mapping = "counts",
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")
)

adata$write_h5ad(paste0("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_intFIB_Subset_12172025.h5ad"), compression="gzip")




###Full FIB Trajectory analyses
#Monocle3 Trajectory Analysis
cds <- as.cell_data_set(KB)
cds <- cluster_cells(cds, reduction_method = "UMAP")
my_clusters <- colData(cds)$v2.subclass.l3
names(my_clusters) <- rownames(colData(cds))
# Assign to the specific Monocle3 slot
cds@clusters$UMAP$clusters <- my_clusters

cds <- learn_graph(cds, use_partition = TRUE)
CFIB <- WhichCells(KB, idents = "C-FIB")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = CFIB)


p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)

plot_cells(cds,
           color_cells_by = "v2.subclass.l3",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           cell_size = 0.6,
           alpha = 0.3)
cds <- order_cells(cds)
pdf(file='trajectories/Int-FIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Monocle3.pdf',width=6.5,height=4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.4,
           alpha = 0.3)
dev.off()
save_monocle_objects(cds, "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_intFIB_Monocle3Object_12172025")

#cds <- load_monocle_objects(directory_path='/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_intFIB_Monocle3Object_12172025/')







###Running Diffusion Map
KB <- Palantir.RunDM(KB, n_components = 20)
KB@reductions
DimPlot2(KB, reduction = "dm", group.by = "v2.subclass.l3", label = TRUE, dims = c(2, 3))


p <- DimPlot(KB, reduction = "dm", dims = c(2, 3),  group.by = "v2.subclass.l3")
cells <- CellSelector(p)

KB <- Palantir.Pseudotime(KB, start_cell = "AHK23_GTCACGGGTGACTGTT-1")
ps <- KB@misc$Palantir$Pseudotime
head(ps)
# Visualize cell fate on UMAP
colnames(ps)[3] <- c("fate1")
KB@meta.data[,colnames(ps)] <- ps
DimPlot2(KB, features = colnames(ps), reduction = "dm", dims = c(2, 3), 
         cols = list(continuous = "A", Entropy = "D"), theme = NoAxes())

pdf(file='trajectories/Int-FIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors_1.pdf',width=6,height=4)
DimPlot2(KB, reduction = "dm", dims = c(2, 3), group.by = "v2.subclass.l3", 
         label = FALSE, theme = NoAxes(),
         cols = v2.fn.cols[levels(Idents(KB))]) + NoLegend()
dev.off()

pdf(file='trajectories/Int-FIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors_2.pdf',width=5,height=4)
DimPlot2(KB, features = "Pseudotime", reduction = "dm", dims = c(2, 3),  
         cols = list(continuous = "A"), theme = NoAxes())
dev.off()


saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_intFIB_SeuratObject_withDM_12172025.rds")

KB <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_intFIB_SeuratObject_withDM_12172025.rds")







# Perivascular Fibroblast Trajectories -----------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(pagoda2)
library(SeuratExtend)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Load Interstitial Fibroblast Data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_16-19_0424-newData.rda")
KB <- KB.STR
rm(KB.STR)
KB <- UpdateSeuratObject(KB)

#update metadata
meta <- KB@meta.data
head(meta)

exp.meta <- read.delim("HKAv2_Experimental_Metadata_12032025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level3","region_level2",
         "region_level1","age_binned","sex","race","KDIGO_stage",
         "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
         "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
         "tissue_type_full","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

KB@meta.data <- meta


#Update annotations
meta <- KB@meta.data
cl.meta <- read.csv("HKAv2_Cluster_Metadata_09232025.csv")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)

meta <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","source","assay","experiment","patient","specimen",
                "condition_level3","condition_level2","condition_level1","condition",
                "percent_cortex","percent_medulla","region_level3","region_level2",
                "region_level1","age_binned","sex","race","KDIGO_stage",
                "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
                "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
                "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
                "tissue_type_full","tissue_type","v2.clusters","v2.subclass.full","v2.subclass.l3",
                "v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1",
                "v2.class","v2.structure")]

KB@meta.data <- meta

#Fix problem characters
unique(KB@meta.data$region_level3)
KB@meta.data[KB@meta.data$region_level3 %in% '\xca',]$region_level3 <- ""
unique(KB@meta.data$region_level3)

KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("AKI","CKD"))


KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:10, n.neighbors = 20L, spread = 3,
              min.dist = 0.01)


v2.fn.cols <- setNames(c("#C8FF00","#FFD300",
                         "#02AD24","#00FF00"),
                       c("pvFIB-RSPO3+","pvFIB-PI16+",
                         "pvFIB","pvMYOF"))

Idents(KB) <- "v2.subclass.l3"
pdf(file='trajectories/pvFIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_colors.pdf',width=6,height=4)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_pvFIB_SeuratObject_12172025.rds")
KB <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_pvFIB_SeuratObject_12172025.rds")

##Save Anndata
library(anndataR)
adata <- as_AnnData(
  KB,
  assay_name = "RNA",
  x_mapping = "counts",
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")
)

adata$write_h5ad(paste0("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_pvFIB_Subset_12172025.h5ad"), compression="gzip")




###Full FIB Trajectory analyses
#Monocle3 Trajectory Analysis
cds <- as.cell_data_set(KB)
cds <- cluster_cells(cds, reduction_method = "UMAP")
my_clusters <- colData(cds)$v2.subclass.l3
names(my_clusters) <- rownames(colData(cds))
# Assign to the specific Monocle3 slot
cds@clusters$UMAP$clusters <- my_clusters

cds <- learn_graph(cds, use_partition = TRUE)
pvFIB <- WhichCells(KB, idents = "pvFIB-PI16+")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = pvFIB)


p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)

plot_cells(cds,
           color_cells_by = "v2.subclass.l3",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
pdf(file='trajectories/pvFIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Monocle3.pdf',width=6.5,height=4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.4,
           alpha = 0.3)
dev.off()

save_monocle_objects(cds, "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_pvFIB_Monocle3Object_01062026")

#cds <- load_monocle_objects(directory_path='/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_pvFIB_Monocle3Object_01062026/')







###Running Diffusion Map
KB <- Palantir.RunDM(KB, n_components = 10)
KB@reductions
DimPlot2(KB, reduction = "ms", group.by = "v2.subclass.l3", label = TRUE, dims = c(1, 2))


p <- DimPlot(KB, reduction = "ms", dims = c(1, 2),  group.by = "v2.subclass.l3")
cells <- CellSelector(p)

KB <- Palantir.Pseudotime(KB, start_cell = "KB29_GGAGGTAGTCATCAGT-1")
ps <- KB@misc$Palantir$Pseudotime
head(ps)
# Visualize cell fate on UMAP
colnames(ps)[3:4] <- c("fate1","fate2")
KB@meta.data[,colnames(ps)] <- ps
DimPlot2(KB, features = colnames(ps), reduction = "ms", dims = c(1, 2), 
         cols = list(continuous = "A", Entropy = "D"), theme = NoAxes())

pdf(file='trajectories/pvFIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors_1.pdf',width=6,height=4)
DimPlot2(KB, reduction = "ms", dims = c(1, 2), group.by = "v2.subclass.l3", 
         label = FALSE, theme = NoAxes(),
         cols = v2.fn.cols[levels(Idents(KB))]) + NoLegend()
dev.off()

pdf(file='trajectories/pvFIB-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors_2.pdf',width=5,height=4)
DimPlot2(KB, features = "Pseudotime", reduction = "ms", dims = c(1, 2),  
         cols = list(continuous = "A"), theme = NoAxes())
dev.off()


saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_pvFIB_SeuratObject_withDM_12172025.rds")
KB <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_pvFIB_SeuratObject_withDM_12172025.rds")








# Myeloid Trajectories -----------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(pagoda2)
library(SeuratExtend)

options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Load Interstitial Fibroblast Data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_0424-newData.rda")
KB <- KB.mye
rm(KB.mye)
KB <- UpdateSeuratObject(KB)

#update metadata
meta <- KB@meta.data
head(meta)

exp.meta <- read.delim("HKAv2_Experimental_Metadata_12032025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level3","region_level2",
         "region_level1","age_binned","sex","race","KDIGO_stage",
         "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
         "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
         "tissue_type_full","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

KB@meta.data <- meta


#Update annotations
meta <- KB@meta.data
cl.meta <- read.csv("HKAv2_Cluster_Metadata_09232025.csv")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)

meta <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","source","assay","experiment","patient","specimen",
                "condition_level3","condition_level2","condition_level1","condition",
                "percent_cortex","percent_medulla","region_level3","region_level2",
                "region_level1","age_binned","sex","race","KDIGO_stage",
                "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
                "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
                "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
                "tissue_type_full","tissue_type","v2.clusters","v2.subclass.full","v2.subclass.l3",
                "v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1",
                "v2.class","v2.structure")]

KB@meta.data <- meta

#Fix problem characters
unique(KB@meta.data$region_level3)
KB@meta.data[KB@meta.data$region_level3 %in% '\xca',]$region_level3 <- ""
unique(KB@meta.data$region_level3)

KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("AKI","CKD"))
KB <- subset(KB, v2.subclass.l3 %in% "moMAC-C3+", invert = TRUE)

KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:30, n.neighbors = 20L, spread = 3,
              min.dist = 0.01)


v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("resMAC-LYVE1+","resMAC-HLAIIhi ","moMAC-HBEGF+","moMAC-CXCL10+",
                         "moFAM","moMAC-C3+","cDC2","ncMON","MON"))

Idents(KB) <- "v2.subclass.l3"
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
pdf(file='trajectories/Myeloid-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_colors.pdf',width=6,height=4)
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = FALSE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()
dev.off()

saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_Myeloid_SeuratObject_12172025.rds")
KB <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_Myeloid_SeuratObject_12172025.rds")

##Save Anndata
library(anndataR)
adata <- as_AnnData(
  KB,
  assay_name = "RNA",
  x_mapping = "counts",
  obsm_mapping = list(X_pca = "pca", X_umap = "umap")
)

adata$write_h5ad(paste0("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_Myeloid_Subset_12172025.h5ad"), compression="gzip")




###Full FIB Trajectory analyses
#Monocle3 Trajectory Analysis
cds <- as.cell_data_set(KB)
cds <- cluster_cells(cds, reduction_method = "UMAP")
my_clusters <- colData(cds)$v2.subclass.l3
names(my_clusters) <- rownames(colData(cds))
# Assign to the specific Monocle3 slot
cds@clusters$UMAP$clusters <- my_clusters

cds <- learn_graph(cds, use_partition = TRUE)
MON <- WhichCells(KB, idents = "MON")
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = MON)


p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
wrap_plots(p1, p2)

plot_cells(cds,
           color_cells_by = "v2.subclass.l3",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
pdf(file='trajectories/Myeloid-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Monocle3.pdf',width=6.5,height=4)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.4,
           alpha = 0.3)
dev.off()

save_monocle_objects(cds, "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_Myeloid_Monocle3Object_01062026")

#cds <- load_monocle_objects(directory_path='/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_Myeloid_Monocle3Object_01062026/')







###Running Diffusion Map (Done with C3)
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_0424-newData.rda")
KB <- KB.mye
rm(KB.mye)
KB <- UpdateSeuratObject(KB)

#update metadata
meta <- KB@meta.data
head(meta)

exp.meta <- read.delim("HKAv2_Experimental_Metadata_12032025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level3","region_level2",
         "region_level1","age_binned","sex","race","KDIGO_stage",
         "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
         "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
         "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
         "tissue_type_full","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}

KB@meta.data <- meta


#Update annotations
meta <- KB@meta.data
cl.meta <- read.csv("HKAv2_Cluster_Metadata_09232025.csv")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)

meta <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","source","assay","experiment","patient","specimen",
                "condition_level3","condition_level2","condition_level1","condition",
                "percent_cortex","percent_medulla","region_level3","region_level2",
                "region_level1","age_binned","sex","race","KDIGO_stage",
                "baseline_eGFR_binned","proteinuria_binned","A1c_binned","albuminuria_binned",
                "diabetes_history","diabetes_duration","hypertension_history","hypertension_duration",
                "on_RAAS_blockade","ckd_stageC","location","laterality","protocol",
                "tissue_type_full","tissue_type","v2.clusters","v2.subclass.full","v2.subclass.l3",
                "v2.subclass.l2","v2.subclass.l1","v2.state.l2","v2.state.l1",
                "v2.class","v2.structure")]

KB@meta.data <- meta

#Fix problem characters
unique(KB@meta.data$region_level3)
KB@meta.data[KB@meta.data$region_level3 %in% '\xca',]$region_level3 <- ""
unique(KB@meta.data$region_level3)

KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 %in% c("AKI","CKD"))
KB <- subset(KB, v2.subclass.l3 %in% "moMAC-C3+", invert = TRUE)

KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:30, n.neighbors = 20L, spread = 3,
              min.dist = 0.01)


v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("resMAC-LYVE1+","resMAC-HLAIIhi ","moMAC-HBEGF+","moMAC-CXCL10+",
                         "moFAM","moMAC-C3+","cDC2","ncMON","MON"))

Idents(KB) <- "v2.subclass.l3"
DimPlot(KB, reduction = "umap", pt.size = 0.8, label = TRUE, alpha = 0.6,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB))],
) + NoLegend()

KB <- Palantir.RunDM(KB, n_components = 30)
KB@reductions
DimPlot2(KB, reduction = "ms", group.by = "v2.subclass.l3", label = TRUE, dims = c(1, 2))


p <- DimPlot(KB, reduction = "ms", dims = c(1, 2),  group.by = "v2.subclass.l3")
cells <- CellSelector(p)

KB <- Palantir.Pseudotime(KB, start_cell = "KB111_CATGCGGGTACCCACC-1")
ps <- KB@misc$Palantir$Pseudotime
head(ps)
# Visualize cell fate on UMAP
colnames(ps)[3:4] <- c("fate1","fate2")
KB@meta.data[,colnames(ps)] <- ps
DimPlot2(KB, features = colnames(ps), reduction = "ms", dims = c(1, 2), 
         cols = list(continuous = "A", Entropy = "D"), theme = NoAxes())

pdf(file='trajectories/Myeloid-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors_1.pdf',width=6,height=4)
DimPlot2(KB, reduction = "ms", dims = c(1, 2), group.by = "v2.subclass.l3", 
         label = FALSE, theme = NoAxes(),
         cols = v2.fn.cols[levels(Idents(KB))]) + NoLegend()
dev.off()

pdf(file='trajectories/Myeloid-Lineages-AKI-CKD-Biopsy-Subset_subclassl3_Palantir_colors_2.pdf',width=5,height=4)
DimPlot2(KB, features = "Pseudotime", reduction = "ms", dims = c(1, 2),  
         cols = list(continuous = "A"), theme = NoAxes())
dev.off()


saveRDS(KB, file = "/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_Myeloid_SeuratObject_withDM_12172025.rds")
KB <- readRDS("/weka/labs/kzhang/kzhang-lab-scratch/blake/HKAv2/HKAv2_Myeloid_SeuratObject_withDM_12172025.rds")







