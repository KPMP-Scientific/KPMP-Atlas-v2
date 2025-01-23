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

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

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

write.table(stats, file="QC_Plots/10X_RNA_post-v2-clustering_Stats_12192024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



#cluster Stats
meta <- KB@meta.data
stats <- cbind(
  data.frame(meta %>%
               group_by(v2.clusters) %>%
               summarise_at(vars(nCount_RNA,nFeature_RNA,percent.er,percent.mt), list(mean))),
  data.frame(meta %>%
               group_by(v2.clusters) %>%
               tally())
)
rownames(stats) <- stats$v2.clusters
stats <- stats[,-c(1,6)]
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

write.table(stats, file="QC_Plots/10X_RNA_post-v2-clustering_Stats_12192024_Cluster_Stats.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



###UMAP Plots
load("color_factors_v2-clusters.robj")

Idents(object = KB) <- "patient"
pdf(file='UMAP_Plots/10X_snCv3_v2_Patient_umap_sub.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Patient"
        ) + scale_color_manual(values = alpha(patient.cols[levels(Idents(KB))], 0.1), name = "Patient"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + NoLegend()
dev.off()

#Condition (level 1)
Idents(object = KB) <- "condition_level1"
cells <- rownames(KB@meta.data[KB@meta.data$condition_level1 %in% c("HRT","CKD","AKI"),])
KB.sub <- subset(KB, cells = cells)
pdf(file='UMAP_Plots/10X_snCv3_v2_Condition_l1_umap_sub.pdf',width=12,height=8)
DimPlot(KB.sub, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1, 
        label.size = 4, repel = TRUE) + ggtitle("Condition"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.sub))], 0.1), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#State level 2
Idents(object = KB) <- "v2.state.l2"
levels(Idents(object = KB))
pdf(file='UMAP_Plots/10X_snCv3_v2_State.l2_umap_sub.pdf',width=12,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("State.l2"
        ) + scale_color_manual(values = alpha(state.l2.cols[levels(Idents(KB))], 0.2), name = "State.l2"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Clusters
Idents(object = KB) <- "v2.clusters"
pdf(file='UMAP_Plots/10X_snCv3_v2_Clusters_umap_sub.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 2, repel = TRUE) + ggtitle("Cluster"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.cl.cols[levels(Idents(object = KB))], 0.1), name = "Cluster"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Subclass level 3
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='UMAP_Plots/10X_snCv3_v2_Subclass.l3_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

pdf(file='UMAP_Plots/10X_snCv3_v1_Subclass.l3_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = sc.l3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

Idents(object = KB) <- "v1.subclass.l3"
pdf(file='UMAP_Plots/10X_snCv3_v1_Subclass.l3_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = sc.l3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Subclass level 1
Idents(object = KB) <- "v2.subclass.l1"
pdf(file='UMAP_Plots/10X_snCv3_v2_Subclass.l1_umap.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl1.cols[levels(Idents(KB))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='UMAP_Plots/10X_snCv3_v2_Subclass.l1_umap_unlabeled.pdf',width=10,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl1.cols[levels(Idents(KB))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#single cell subclass level 1
load("~/hsKidAt/blake_LTS/Atlas_V2/sc_objects/KPMP_PREMIERE_sc_downsampled_Jun2024.Robj")
kpmpv2small$celltype[kpmpv2small$celltype == "VSMC/MC/FIB"] <- "FIB"
kpmpv2small$celltype[kpmpv2small$celltype == "MYL"] <- "Myeloid"
kpmpv2small$celltype[kpmpv2small$celltype == "DCT/CNT"] <- "DCT"
kpmpv2small$celltype[kpmpv2small$celltype == "POD/PEC"] <- "POD"
kpmpv2small$celltype[kpmpv2small$celltype == "T"] <- "Lymphoid"
kpmpv2small$celltype[kpmpv2small$celltype == "B"] <- "Lymphoid"
kpmpv2small$celltype[kpmpv2small$celltype == "SchwannCells"] <- "NEU"

Idents(object = kpmpv2small) <- "celltype"
pdf(file='UMAP_Plots/10X_scCv3_v2_Subclass.l1_umap_unlabeled.pdf',width=10,height=8)
DimPlot(kpmpv2small, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.5,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl1.cols[levels(Idents(kpmpv2small))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Sex
Idents(object = KB) <- "sex"
pdf(file='UMAP_Plots/10X_snCv3_v2_Sex_umap_sub.pdf',width=12,height=8)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Sex"
        ) + scale_color_manual(values = alpha(sex.cols[levels(Idents(KB))], 0.1), name = "Sex"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()



###Plot clustering groups separately
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
meta <- KB@meta.data
load("color_factors_v2-clusters.robj")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")
KB@meta.data <- meta[rownames(KB@meta.data),]


Idents(KB) <- "v2.subclass.l3"
pdf(file='UMAP_Plots/10X_snCv3_IMM_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='UMAP_Plots/10X_snCv3_IMM_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
Idents(KB) <- "v2.subclass.int"
pdf(file='UMAP_Plots/10X_snCv3_IMM_v2_Subclass.integrated_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3.int"
        ) + NoLegend() + scale_color_manual(values = v2.scint.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

Idents(KB) <- "v2.clusters"
pdf(file='UMAP_Plots/10X_snCv3_IMM_v2_Clusters_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.2,
        label.size = 4, repel = TRUE) + ggtitle("Clusters"
        ) + NoLegend() + scale_color_manual(values = v2.cl.cols[levels(Idents(KB))], name = "Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='UMAP_Plots/10X_snCv3_STR_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='UMAP_Plots/10X_snCv3_STR_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
Idents(object = KB) <- "v2.clusters"
pdf(file='UMAP_Plots/10X_snCv3_STR_v2_Clusters_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Cluster"
        ) + NoLegend() + scale_color_manual(values = v2.cl.cols[levels(Idents(KB))], name = "Cluster"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='UMAP_Plots/10X_snCv3_EC_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='UMAP_Plots/10X_snCv3_EC_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
Idents(object = KB) <- "v2.clusters"
pdf(file='UMAP_Plots/10X_snCv3_EC_v2_Clusters_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Cluster"
        ) + NoLegend() + scale_color_manual(values = v2.cl.cols[levels(Idents(KB))], name = "Cluster"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='UMAP_Plots/10X_snCv3_PT_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='UMAP_Plots/10X_snCv3_PT_v2_Subclass.l3_umap_Unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
Idents(object = KB) <- "v2.clusters"
pdf(file='UMAP_Plots/10X_snCv3_PT_v2_Clusters_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Cluster"
        ) + NoLegend() + scale_color_manual(values = v2.cl.cols[levels(Idents(KB))], name = "Cluster"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
KB@meta.data <- meta[rownames(KB@meta.data),]
Idents(object = KB) <- "v2.subclass.l3"
pdf(file='UMAP_Plots/10X_snCv3_DT_v2_Subclass.l3_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='UMAP_Plots/10X_snCv3_DT_v2_Subclass.l3_umap_unlabeled.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = FALSE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KB))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
Idents(object = KB) <- "v2.clusters"
pdf(file='UMAP_Plots/10X_snCv3_DT_v2_Clusters_umap.pdf',width=7,height=6)
DimPlot(KB, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Cluster"
        ) + NoLegend() + scale_color_manual(values = v2.cl.cols[levels(Idents(KB))], name = "Cluster"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()





###Age barplot
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
pdf(file='QC_Plots/Binned_Age_Patient_Counts_Barplot.pdf',width=6,height=3)
meta <- KB@meta.data
meta <- meta[!duplicated(meta$patient),]
stats <- cbind(
  data.frame(meta %>%
               group_by(age_binned) %>%
               tally())
)
rownames(stats) <- stats$age_binned
barplot(stats$n, main = "Age Groups V2", cex.names = 1, 
        names.arg = rownames(stats),las=2)
dev.off()

#compare with age group distribution from V1
meta <- KB@meta.data[KB@meta.data$atlas_version == "V1",]
meta <- meta[!duplicated(meta$patient),]
stats <- cbind(
  data.frame(meta %>%
               group_by(age_binned) %>%
               tally())
)
rownames(stats) <- stats$age_binned
barplot(stats$n, main = "Age Groups V1", cex.names = 1, 
        names.arg = rownames(stats),las=2)




###Cell State proportions by condition
load("color_factors.robj")
load("color_factors_v2-clusters.robj")

Idents(KB) <- "condition_level1"
order <- c("HRT","RT-UCS","DM-R","NHT","AKI","CKD")
Idents(KB) <- factor(Idents(KB), levels = order)


pdf(file='QC_Plots/Condition_Altered-state_Barplot.pdf',width=5,height=6)
col.order <- order
row.order <- c("degenerative","cycling","adaptive - epi","adaptive - str","adaptive - ec","adaptive - imm","transitioning")
prop1 <- prop.table(table(KB$v2.state.l2, KB$condition_level1), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(state.l2.cols[rownames(prop1)]))
dev.off()


pdf(file='QC_Plots/Patient_Altered-state_Barplot.pdf',width=26,height=6)
table(KB$tissue_type)
dd <- unique(KB@meta.data[KB@meta.data$tissue_type %in% "Deceased Donor" &
                            KB@meta.data$condition_level1 %in% "HRT", ]$patient)
nx <- unique(KB@meta.data[KB@meta.data$tissue_type %in% "Nephrectomy" &
                            KB@meta.data$condition_level1 %in% "HRT", ]$patient)
bx <- unique(KB@meta.data[KB@meta.data$tissue_type %in% "Biopsy" &
                            KB@meta.data$condition_level1 %in% "HRT", ]$patient)
RT.UCS <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "RT-UCS", ]$patient)
NHT <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "NHT", ]$patient)
DMR <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "DM-R", ]$patient)

aki <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "AKI" &
                             KB@meta.data$tissue_type %in% "Biopsy" , ]$patient)
aki2 <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "AKI" &
                              KB@meta.data$tissue_type %in% "Deceased Donor" , ]$patient)
ckd <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "CKD" &
                             KB@meta.data$tissue_type %in% "Biopsy" , ]$patient)
ckd2 <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "CKD" &
                              KB@meta.data$tissue_type %in% "Nephrectomy" , ]$patient)
ckd3 <- unique(KB@meta.data[KB@meta.data$condition_level1 %in% "CKD" &
                              KB@meta.data$tissue_type %in% "Deceased Donor" , ]$patient)
col.order <- c(bx,dd,nx,RT.UCS,DMR,NHT,aki,aki2,ckd,ckd2,ckd3)
unique(KB$patient)[!unique(KB$patient) %in% col.order]
row.order <- c("degenerative","cycling","adaptive - epi","adaptive - str","adaptive - ec","adaptive - imm","transitioning")
prop1 <- prop.table(table(KB$v2.state.l2, KB$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(state.l2.cols[rownames(prop1)]))
dev.off()

write.table(prop1, file="QC_Plots/Patient_Altered-state_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


pdf(file='QC_Plots/Patient_Altered-state_Barplot_noDegen.pdf',width=26,height=6)
row.order <- c("cycling","adaptive - epi","adaptive - str","adaptive - ec","adaptive - imm","transitioning")
prop1 <- prop.table(table(KB$v2.state.l2, KB$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(state.l2.cols[rownames(prop1)]))
dev.off()


pdf(file='QC_Plots/Patient_Immune_Prop_Barplot.pdf',width=26,height=6)
row.order <- c("immune cells")
prop1 <- prop.table(table(KB$v2.class, KB$patient), margin = 2)[row.order,col.order]
barplot(prop1,main = "State Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = "#E68AC9")
dev.off()



###Build dendrogram
library(dendextend)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
load("color_factors_v2-clusters.robj")
load("color_factors.robj")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = c(paste0("P_",1:33),
                                            paste0("D_",1:48),
                                            paste0("E_",1:22),
                                            paste0("S_",1:27),
                                            paste0("I_",1:28),
                                            paste0("N_",1)))

table(Idents(KB))
KB <- subset(KB, downsample = 1000)
KB[["RNA"]] <- as(object = KB[["RNA"]], Class = "Assay")

#Get all over-dispersed genes
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200_0424-newData.rda")
pt.od.genes <- sn.od.genes[1:3000]
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200_0424-newData.rda")
dt.od.genes <- sn.od.genes[1:3000]
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Vasculature_k200_0424-newData.rda")
ec.od.genes <- sn.od.genes[1:2000]
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k200_0424-newData.rda")
str.od.genes <- sn.od.genes[1:2000]
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k200_0424-newData.rda")
imm.od.genes <- sn.od.genes[1:2000]
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_SC-NEU_k100_round1_newData-0424.rda")
neu.od.genes <- sn.od.genes[1:10]

save(neu.od.genes,imm.od.genes,str.od.genes,ec.od.genes,pt.od.genes,dt.od.genes, file = "AtlasV2_VariableFeatures.rda")
#od.genes <- unique(c(pt.od.genes,dt.od.genes,imm.od.genes,str.od.genes,ec.od.genes,neu.od.genes))
od.genes <- unique(c(pt.od.genes[1:1000],dt.od.genes[1:1000],imm.od.genes[1:1000],
                     str.od.genes[1:1000], ec.od.genes[1:1000], neu.od.genes[1:10]))
rm(p2)

cl <- Idents(KB)
norm.dat <- KB[["RNA"]]$data

#select.markers
cl.med <- get_cl_medians(norm.dat[od.genes,], cl)

##The preferred order for the leaf nodes.
l.rank <- setNames(1:159,colnames(cl.med))
l.color <- v2.cl.cols[colnames(cl.med)]
dend.result <- build_dend(cl.med,
                          cl.cor=NULL,
                          l.rank, 
                          l.color,
                          nboot = 100)
dend.result$pvclust.result
dend <- dend.result$dend
dend <- dend %>% set("labels_col", v2.cl.cols[labels(dend)])
dend <- dend %>% set("leaves_col", v2.cl.cols[labels(dend)])
dend <- dend %>% pvclust_show_signif_gradient(dend.result$pvclust.result, 
                                              signif_type = "bp", signif_col_fun=colorRampPalette(
                                                c("gray20","gray5","black")))
dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 1) 
plot(dend, main = "Clusters")
order <- labels(dend)

save(dend.result, file = "Plots/Atlas_V2_Cluster_Dendrogram_05082024.rda")
save(order, file = "Plots/Atlas_V2_Cluster_Dendrogram_05082024_label-order.rda")

#Plot dendrogram colored by subclass.l1
load("Plots/Atlas_V2_Cluster_Dendrogram_05082024.rda")
load("color_factors_v2-clusters.robj")
v2.scl1.cols
v2.cl.cols <- setNames(c(rep(as.character(v2.scl1.cols["POD"]),2),
                         rep(as.character(v2.scl1.cols["PEC"]),1),
                         rep(as.character(v2.scl1.cols["PT"]),22),
                         rep(as.character(v2.scl1.cols["DTL"]),6),
                         rep(as.character(v2.scl1.cols["ATL"]),2),
                         rep(as.character(v2.scl1.cols["TAL"]),18),
                         rep(as.character(v2.scl1.cols["DCT"]),5),
                         rep(as.character(v2.scl1.cols["CNT"]),4),
                         rep(as.character(v2.scl1.cols["PC"]),11),
                         rep(as.character(v2.scl1.cols["PapE"]),1),
                         rep(as.character(v2.scl1.cols["IC"]),9),
                         rep(as.character(v2.scl1.cols["EC"]),22),
                         rep(as.character(v2.scl1.cols["FIB"]),19),
                         rep(as.character(v2.scl1.cols["VSM/P"]),7),
                         rep(as.character(v2.scl1.cols["Ad"]),1),
                         rep(as.character(v2.scl1.cols["Lymphoid"]),11),
                         rep(as.character(v2.scl1.cols["Myeloid"]),15),
                         rep(as.character(v2.scl1.cols["Lymphoid"]),1),
                         rep(as.character(v2.scl1.cols["Myeloid"]),1),
                         rep(as.character(v2.scl1.cols["NEU"]),1)),
                       c(paste0("P_",1:33),
                         paste0("D_",1:48),
                         paste0("E_",1:22),
                         paste0("S_",1:27),
                         paste0("I_",1:28),
                         paste0("N_",1))
)


dend <- dend.result$dend
dend <- dend %>% set("labels_col", v2.cl.cols[labels(dend)])
dend <- dend %>% set("leaves_col", v2.cl.cols[labels(dend)])
dend <- dend %>% pvclust_show_signif_gradient(dend.result$pvclust.result, 
                                              signif_type = "bp", signif_col_fun=colorRampPalette(
                                                c("gray20","gray5","black")))
dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 1) 

pdf(file = "Plots/Atlas_V2_Cluster_Dendrogram_05082024_3.pdf", width = 18, height = 4)
par(mar = c(10, 2, 2, 2))
plot(dend, main = "Clusters")
dev.off()



###Plot Cluster level QC
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
load("Plots/Atlas_V2_Cluster_Dendrogram_05082024_label-order.rda")
order

Idents(object = KB) <- "v2.clusters"
Idents(object = KB) <- factor(Idents(object = KB), levels = order)

load("color_factors_v2-clusters.robj")
load("color_factors.robj")

#Violin plots
pdf(file='QC_Plots/AtlasV2_Cluster_Violin_plots.pdf',width=16,height=12)
VlnPlot(object = KB, features = c("nCount_RNA", "nFeature_RNA","percent.er","percent.mt"), 
        ncol = 1, pt.size = -1, cols = v2.cl.cols[levels(Idents(KB))])
dev.off()

#Barplots 1
pdf(file='QC_Plots/AtlasV2_Cluster_Barplots.pdf',width=8,height=12)
layout(matrix(c(1,2,3,4,5), nrow = 5, ncol = 1, byrow = TRUE))
prop1 <- prop.table(table(KB$condition_level2, KB$v2.clusters), margin = 2)[,order]
barplot(prop1, main = "Condition Proportions", names.arg = colnames(prop1), las=2, border = NA, 
        cex.names = 0.5, col = as.character(condl2.cols[rownames(prop1)]))
prop2 <- prop.table(table(KB$patient, KB$v2.clusters), margin = 2)[,order]
barplot(prop2, main = "Patient Proportions", names.arg = colnames(prop2), border = NA,
        cex.names = 0.5, las=2, col = as.character(patient.cols[rownames(prop2)]))
prop3 <- prop.table(table(KB$region_level2, KB$v2.clusters), margin = 2)[,order]
barplot(prop3,main = "Region Proportions", names.arg = colnames(prop3), las=2, border = NA,
        cex.names = 0.5, col = as.character(region.l2.cols[rownames(prop3)]))
prop4 <- prop.table(table(KB$sex, KB$v2.clusters)[c("Male","Female"),], margin = 2)[,order]
barplot(prop4,main = "Sex Proportions", names.arg = colnames(prop3), las=2, border = NA,
        cex.names = 0.5, col = as.character(sex.cols[rownames(prop4)]))

tab1 <- table(KB$v2.clusters)[order]
barplot(tab1, col = "gray", main = "Cluster Size",  border = NA,
        names.arg = order,las=2, cex.names = 0.5)
dev.off()

unique(KB$region_level3)
pdf(file='QC_Plots/AtlasV2_Cluster_Barplots_regl3.pdf',width=8,height=2.4)
prop5 <- prop.table(table(KB$region_level3, KB$v2.clusters)[c("C","CMJ + OM","OM","IM","Papilla"),], margin = 2)[,order]
barplot(prop5,main = "Region Proportions", names.arg = colnames(prop5), las=2, border = NA,
        cex.names = 0.5, col = c("#508f5b","#3D8ACC","#3D8ACC","#7014CC","#B3369F"))
dev.off()

stats <- rbind(prop1,prop3,prop4,prop5,tab1)[,c(paste0("P_",1:33),
                                                paste0("D_",1:48),
                                                paste0("E_",1:22),
                                                paste0("S_",1:27),
                                                paste0("I_",1:28),
                                                paste0("N_",1))]

write.table(stats, file="QC_Plots/AtlasV2_Cluster_region-condition-sex_Stats_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
