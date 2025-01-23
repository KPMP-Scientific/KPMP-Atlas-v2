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

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")
source("misc/utils.R")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
KB
#35478 features across 315764 samples

Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = c(paste0("P_",1:31),
                                            paste0("D_",1:47),
                                            paste0("E_",1:21),
                                            paste0("S_",1:25),
                                            paste0("I_",1:27),
                                            paste0("N_",1)))
table(Idents(KB))


###Update metadata
meta <- KB@meta.data
exp.meta <- read.delim("mouse_IRI/Mouse_Experiment_Metadata_11152024.txt")
emc <- c("library","source","assay","experiment","patient",       
         "specimen","injury_condition","condition_level3","condition_level2","condition_level1","condition",      
         "age_months","age_group","sex","genotype","strain","protocol","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB@meta.data <- meta

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")


#update object
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
umap <- Embeddings(KB, reduction = "umap")
rpca <- Embeddings(KB, reduction = "integrated.rpca")
counts <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/kidney_counts_112024")
meta <- KB@meta.data

KB <- CreateSeuratObject(counts = counts, project = "Mouse Kidney Atlas V2", min.cells = 3, min.features = 200, 
                         meta.data = meta)

KB[["integrated.rpca"]] <- CreateDimReducObject(embeddings = rpca, key = "integratedrpca_", assay = DefaultAssay(KB))
KB[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "umap_", assay = DefaultAssay(KB))

KB <- NormalizeData(KB)
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")




###Final Post-V2 annotation/alignment Stats
meta <- KB@meta.data
stats <- cbind(
  data.frame(meta %>%
               group_by(library) %>%
               summarise_at(vars(nCount_RNA,nFeature_RNA), list(mean))),
  data.frame(meta %>%
               group_by(library) %>%
               tally())
)
rownames(stats) <- stats$library
stats <- stats[,-c(1,4)]
stats  

#Determine number of clusters
df <- t(table(KB$library, KB$v2.clusters))
cols <- vector()
exps <- colnames(df)

for(i in exps){
  cols[i] <- length(df[which(df[,i] > 5),i]) 
}
stats$n_clusters <- cols[rownames(stats)]

write.table(stats, file="mouse_IRI/QC_Plots/Mouse_10X_RNA_post-v2-clustering_Stats_11152024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



#cluster Stats
meta <- KB@meta.data
stats <- cbind(
  data.frame(meta %>%
               group_by(v2.clusters) %>%
               summarise_at(vars(nCount_RNA,nFeature_RNA), list(mean))),
  data.frame(meta %>%
               group_by(v2.clusters) %>%
               tally())
)
rownames(stats) <- stats$v2.clusters
stats <- stats[,-c(1,4)]
stats  


#Determine number per patient
stats2 <- cbind(
  data.frame(meta %>%
               group_by(v2.clusters, patient) %>%
               summarise(n = n()))
)
library(tidyverse)
stats2 <- stats2 %>%
  spread(patient, n, fill = 0)
rownames(stats2) <- stats2$v2.clusters

stats <- cbind(stats, stats2[rownames(stats),])

write.table(stats, file="mouse_IRI/QC_Plots/Mouse_10X_RNA_post-v2-clustering_Stats_11152024_Cluster_Stats.txt", sep = "\t", row.names=TRUE, col.names=TRUE)





###UMAP Plots
##Generate plots for all clusters and annotated clusters
source("misc/utils.R")
##Assign new exp/patient colors
load("color_factors_v2-clusters.robj")

meta <- KB@meta.data
reorder <- sample(rownames(meta))
meta <- meta[reorder,]
meta <- meta[!duplicated(meta$v2.clusters),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "v2.clusters", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "v2.subclass.l3", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "v2.subclass.l2", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "v2.subclass.l1", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "v2.subclass.sp", sort_label = F, colorset = "varibow")


v2.cl.cols <- meta$v2.clusters_color; names(v2.cl.cols) <- meta$v2.clusters_label
v2.scl3.cols <- meta$v2.subclass.l3_color; names(v2.scl3.cols) <- meta$v2.subclass.l3_label
v2.scl2.cols <- meta$v2.subclass.l2_color; names(v2.scl2.cols) <- meta$v2.subclass.l2_label
v2.scl1.cols <- meta$v2.subclass.l1_color; names(v2.scl1.cols) <- meta$v2.subclass.l1_label
v2.scsp.cols <- meta$v2.subclass.sp_color; names(v2.scsp.cols) <- meta$v2.subclass.sp_label

meta <- KB@meta.data
meta <- meta[,c("source","patient", "condition_level1", "condition_level2", "condition_level3")]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "source", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "patient", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level1", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level2", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level3", sort_label = F, colorset = "varibow")

source.cols <- meta$source_color; names(source.cols) <- meta$source_label; source.cols <- factor(na.omit(source.cols))
source.cols <- factor(source.cols[!duplicated(source.cols)])
patient.cols <- meta$patient_color; names(patient.cols) <- meta$patient_label; patient.cols <- factor(na.omit(patient.cols))
patient.cols <- factor(patient.cols[!duplicated(patient.cols)])
condl1.cols <- meta$condition_level1_color; names(condl1.cols) <- meta$condition_level1_label; condl1.cols <- factor(na.omit(condl1.cols))
condl1.cols <- factor(condl1.cols[!duplicated(condl1.cols)])
condl2.cols <- meta$condition_level2_color; names(condl2.cols) <- meta$condition_level2_label; condl2.cols <- factor(na.omit(condl2.cols))
condl2.cols <- factor(condl2.cols[!duplicated(condl2.cols)])
condl3.cols <- meta$condition_level3_color; names(condl3.cols) <- meta$condition_level3_label; condl3.cols <- factor(na.omit(condl3.cols))
condl3.cols <- factor(condl3.cols[!duplicated(condl3.cols)])

condl3.cols <- setNames(c("#E6ACAC","#806075","#145ECC","#4D6E80","#688039", "#1F9950","#45E6AF","#8936B3"), 
                        c("Sham","Ref","4hrs","12hrs","2d","2wk","4-6wk","20-24wk"))

condl2.cols <- setNames(c("#E6ACAC","#806075","#8936B3","#4D6E80","#145ECC"), 
                        c("Sham","Ref","AKI-noncycled","AKI-cycled","AKI"))
condl1.cols <- setNames(c("#E6ACAC","#145ECC"), 
                        c("Ref","AKI"))
state.l2.cols <- setNames(c("#8eab9b","#c76f67","#4747B3","#14CCCC","#B3B336","#807226","#806075","#8936B3"),
                          c("reference","degenerative","cycling","transitioning","adaptive - epi","adaptive - str","adaptive - imm","adaptive - ec"))
sex.cols <- setNames(c("#36B3B3", "#CC5252"), c("Male","Female"))
save(v2.cl.cols,v2.scl3.cols,v2.scl2.cols,v2.scl1.cols,v2.scsp.cols,condl3.cols,
     condl2.cols,condl1.cols,source.cols,patient.cols,sex.cols,
     state.l2.cols,
     file = "mouse_color_factors_v2-clusters.robj")
load("mouse_color_factors_v2-clusters.robj")


Idents(object = KB) <- "patient"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Patient_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(patient.cols[levels(Idents(KB))], 0.1), name = "Patient"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
dev.off()



#Condition (level 1)
Idents(object = KB) <- "condition_level1"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Condition_l1_umap.pdf',width=12,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1, 
        label.size = 4, repel = TRUE) + ggtitle("Condition"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB))], 0.1), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Condition (level 2)
Idents(object = KB) <- "condition_level3"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Condition_l3_umap.pdf',width=12,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Condition"
        ) + scale_color_manual(values = alpha(condl3.cols[levels(Idents(KB))], 0.1), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#State level 2
Idents(object = KB) <- "v2.state.l2"
levels(Idents(object = KB))
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_State.l2_umap.pdf',width=12,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("State.l2"
        ) + scale_color_manual(values = alpha(state.l2.cols[levels(Idents(KB))], 0.2), name = "State.l2"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Clusters
Idents(object = KB) <- "v2.clusters"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Clusters_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 2, repel = TRUE) + ggtitle("Cluster"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.cl.cols[levels(Idents(object = KB))], 0.1), name = "Cluster"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()



#Subclass level 3
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Subclass.l3_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Subclass level 1
#use human colors
load("color_factors_v2-clusters.robj")
Idents(object = KB) <- "v2.subclass.l1"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Subclass.l1_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl1.cols[levels(Idents(KB))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Subclass.l1_umap_unlabeled.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl1.cols[levels(Idents(KB))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

load("mouse_color_factors_v2-clusters.robj")

#Sex
Idents(object = KB) <- "sex"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Sex_umap.pdf',width=12,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(KB))], 0.1), name = "Sex"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Source
Idents(object = KB) <- "source"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_v2_Source_umap.pdf',width=12,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Source"
        ) + scale_color_manual(values = alpha(source.cols[levels(Idents(KB))], 0.1), name = "Sex"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()



###Plot individual v2 labels
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
meta <- KB@meta.data
load("mouse_color_factors_v2-clusters.robj")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_C.Rds")
KB@meta.data <- meta[rownames(KB@meta.data),]


Idents(KB) <- "v2.subclass.l3"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_IMM_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_IMM_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_STR_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_STR_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_C.Rds")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_EC_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_EC_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_PT_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_PT_v2_Subclass.l3_umap_Unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_DT_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='mouse_IRI/UMAP_Plots/10X_RNA_DT_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()






###Injury Timepoint Barplot V2
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
pdf(file='mouse_IRI/QC_Plots/Binned_Timepoints_Counts_Barplot.pdf',width=6,height=3)
meta <- KB@meta.data
meta <- meta[!duplicated(meta$patient),]
stats <- cbind(
  data.frame(meta %>%
               group_by(condition_level3) %>%
               tally())
)
rownames(stats) <- stats$condition_level3
barplot(stats[c("Ref","Sham","4hrs","12hrs","2d","2wk","4-6wk","20-24wk"),]$n, main = "Condition", cex.names = 1, 
        names.arg = rownames(stats[c("Ref","Sham","4hrs","12hrs","2d","2wk","4-6wk","20-24wk"),]),las=2)
dev.off()




###Cell State proportions by condition
load("color_factors.robj")
load("mouse_color_factors_v2-clusters.robj")

Idents(KB) <- "condition_level3"
order <- c("Sham","4hrs","12hrs","2d","2wk","4-6wk","20-24wk")
Idents(KB) <- factor(Idents(KB), levels = order)


pdf(file='mouse_IRI/QC_Plots/Condition_Altered-state_Barplot_Injury_Timecourse.pdf',width=5,height=6)
col.order <- order
row.order <- c("degenerative","cycling","adaptive - epi","adaptive - str","adaptive - ec","adaptive - imm")
prop1 <- prop.table(table(KB$v2.state.l2, KB$condition_level3), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(state.l2.cols[rownames(prop1)]))
dev.off()



pdf(file='mouse_IRI/QC_Plots/Sample_Altered-state_Barplot.pdf',width=26,height=6)
table(KB$condition_level3)
ref <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "Ref", ]$patient)
sham <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "Sham", ]$patient)
iri.4h <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "4hrs", ]$patient)
iri.12h <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "12hrs", ]$patient)
iri.2d <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "2d", ]$patient)
iri.2wk <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "2wk", ]$patient)
iri.4wk <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "4-6wk", ]$patient)
iri.20wk <- unique(KB@meta.data[KB@meta.data$condition_level3 %in% "20-24wk", ]$patient)


col.order <- c(ref,sham,iri.4h,iri.12h,iri.2d,iri.2wk,iri.4wk,iri.20wk)
unique(KB$patient)[!unique(KB$patient) %in% col.order]
row.order <- c("degenerative","cycling","adaptive - epi","adaptive - str","adaptive - ec","adaptive - imm")
prop1 <- prop.table(table(KB$v2.state.l2, KB$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(state.l2.cols[rownames(prop1)]))
dev.off()

write.table(prop1, file="mouse_IRI/QC_Plots/Sample_Altered-state_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


pdf(file='mouse_IRI/QC_Plots/Sample_Altered-state_Barplot_noDegen.pdf',width=26,height=6)
row.order <- c("cycling","adaptive - epi","adaptive - str","adaptive - ec","adaptive - imm")
prop1 <- prop.table(table(KB$v2.state.l2, KB$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(state.l2.cols[rownames(prop1)]))
dev.off()


pdf(file='mouse_IRI/QC_Plots/Sample_Immune_Prop_Barplot.pdf',width=26,height=6)
row.order <- c("immune cells")
prop1 <- prop.table(table(KB$v2.class, KB$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = "#E68AC9")
dev.off()


###Cluster level QC
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
order <- c(paste0("P_",1:31),
           paste0("D_",1:47),
           paste0("E_",1:21),
           paste0("S_",1:25),
           paste0("I_",1:27),
           paste0("N_",1))

Idents(object = KB) <- "v2.clusters"
Idents(object = KB) <- factor(Idents(object = KB), levels = order)

#Violin plots
pdf(file='mouse_IRI/QC_Plots/AtlasV2_Cluster_Violin_plots.pdf',width=16,height=12)
VlnPlot(object = KB, features = c("nCount_RNA", "nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, cols = v2.cl.cols[levels(Idents(KB))])
dev.off()

#Barplots 1
pdf(file='mouse_IRI/QC_Plots/AtlasV2_Cluster_Barplots.pdf',width=8,height=12)
layout(matrix(c(1,2,3,4), nrow = 4, ncol = 1, byrow = TRUE))
prop1 <- prop.table(table(KB$condition_level3, KB$v2.clusters), margin = 2)[,order]
barplot(prop1, main = "Condition Proportions", names.arg = colnames(prop1), las=2, border = NA, 
        cex.names = 0.5, col = as.character(condl2.cols[rownames(prop1)]))
prop2 <- prop.table(table(KB$patient, KB$v2.clusters), margin = 2)[,order]
barplot(prop2, main = "Patient Proportions", names.arg = colnames(prop2), border = NA,
        cex.names = 0.5, las=2, col = as.character(patient.cols[rownames(prop2)]))
prop4 <- prop.table(table(KB$sex, KB$v2.clusters)[c("Male","Female"),], margin = 2)[,order]
barplot(prop4,main = "Sex Proportions", names.arg = colnames(prop4), las=2, border = NA,
        cex.names = 0.5, col = as.character(sex.cols[rownames(prop4)]))

tab1 <- table(KB$v2.clusters)[order]
barplot(tab1, col = "gray", main = "Cluster Size",  border = NA,
        names.arg = order,las=2, cex.names = 0.5)
dev.off()

stats <- rbind(prop1,prop4,tab1)[,order]

write.table(stats, file="mouse_IRI/QC_Plots/AtlasV2_Cluster_condition-sex_Stats_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
