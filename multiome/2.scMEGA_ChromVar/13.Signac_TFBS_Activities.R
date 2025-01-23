library(Signac)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(dplyr)
library(JASPAR2022)
library(TFBSTools)

set.seed(1234)
setwd("~/Projects/Human_Kidney/Atlas_V2")



###Cortical Fibroblasts
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_05162024.rda")
KB.STR <- subset(KRAC.nonEpi, idents = paste0("S_", c(9:14,16:19)))

#Link motifs
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
KB.STR <- AddMotifs(
  object = KB.STR,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

KB.STR <- RunChromVAR(
  object = KB.STR,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

KB.STR
#An object of class Seurat 
#624293 features across 17085 samples within 3 assays 
#Active assay: ATAC (586864 features, 0 variable features)
#2 layers present: counts, data
#2 other assays present: RNA, chromvar
#2 dimensional reductions calculated: pca, umap

save(KB.STR, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Stroma_Dual_Seurat_Cortex_Only_05172024.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Stroma_Dual_Seurat_Cortex_Only_05172024.rda")



#TFBS activities Interstitial
KB.STR <- subset(KB.STR, idents = paste0("S_", c(9:14)))
DefaultAssay(KB.STR) <- "chromvar"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = c("S_9","S_10","S_11","S_13","S_12","S_14"))
table(Idents(KB.STR))
rownames(KB.STR@assays$chromvar@data) <- ConvertMotifID(KB.STR,
                                                        id=rownames(KB.STR@assays$chromvar@data),
                                                        assay='ATAC')
tf.markers <- FindAllMarkers(
  object = KB.STR,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

#subset to TFs expressed in at least 10% of cells in the cluster
library(chromfunks)
DefaultAssay(KB.STR) <- "RNA"
clusters <- levels(Idents(KB.STR))
rna.counts <- KB.STR[["RNA"]]$counts
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KB.STR), binarize = T)
min.rna.cl.frac <- 0.1
tf.markers$motif.name <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters


tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
tf.markers.exp <- na.omit(tf.markers.exp)
DefaultAssay(KB.STR) <- 'chromvar'

save(tf.markers.exp, file = "regulatory/TFBS_Activities_Interstitial.rda")

#Plot TFBS Activities
tf.mark <- tf.markers.exp[tf.markers.exp$p_val_adj < 0.01,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 

DotPlot(KB.STR, features = unique(top10$motif.name), dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()



#TFBS activities perivascular
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Stroma_Dual_Seurat_Cortex_Only_05172024.rda")
KB.STR <- subset(KB.STR, idents = paste0("S_", c(16:19)))
DefaultAssay(KB.STR) <- "chromvar"
table(Idents(KB.STR))
rownames(KB.STR@assays$chromvar@data) <- ConvertMotifID(KB.STR,
                                                        id=rownames(KB.STR@assays$chromvar@data),
                                                        assay='ATAC')
tf.markers <- FindAllMarkers(
  object = KB.STR,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

#subset to TFs expressed in at least 10% of cells in the cluster
library(chromfunks)
DefaultAssay(KB.STR) <- "RNA"
clusters <- levels(Idents(KB.STR))
rna.counts <- KB.STR[["RNA"]]$counts
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KB.STR), binarize = T)
min.rna.cl.frac <- 0.1
tf.markers$motif.name <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters


tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
tf.markers.exp <- na.omit(tf.markers.exp)
DefaultAssay(KB.STR) <- 'chromvar'

save(tf.markers.exp, file = "regulatory/TFBS_Activities_Perivascular.rda")

#Plot TFBS Activities
tf.mark <- tf.markers.exp[tf.markers.exp$p_val_adj < 0.05,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 

DotPlot(KB.STR, features = unique(top10$motif.name), dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()



### Proximal Tubules
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_PT-TL_Dual_Seurat_all-peaks_05162024.rda")
KB.PT <- subset(KRAC.PEpi, idents = paste0("P_", c(4,8,10,12:22)))

#Link motifs
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
KB.PT <- AddMotifs(
  object = KB.PT,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

KB.PT <- RunChromVAR(
  object = KB.PT,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

KB.PT
#An object of class Seurat 
#624293 features across 17085 samples within 3 assays 
#Active assay: ATAC (586864 features, 0 variable features)
#2 layers present: counts, data
#2 other assays present: RNA, chromvar
#2 dimensional reductions calculated: pca, umap

save(KB.PT, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_PT_Dual_Seurat_PT-aPT_Only_05202024.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_PT_Dual_Seurat_PT-aPT_Only_05202024.rda")


#TFBS activities - All PT subclass.l3
DefaultAssay(KB.PT) <- "chromvar"
Idents(KB.PT) <- "v2.subclass.l3"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = c("aPT2","aPT1","aPT-S1/S2","PT-S1","PT-S2",
                                                  "frPT-S1/S2",
                                                  "aPT-S3","PT-S3","frPT-S3"))
table(Idents(KB.PT))
rownames(KB.PT@assays$chromvar@data) <- ConvertMotifID(KB.PT,
                                                        id=rownames(KB.PT@assays$chromvar@data),
                                                        assay='ATAC')
tf.markers <- FindAllMarkers(
  object = KB.PT,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

#subset to TFs expressed in at least 10% of cells in the cluster
library(chromfunks)
DefaultAssay(KB.PT) <- "RNA"
clusters <- levels(Idents(KB.PT))
rna.counts <- KB.PT[["RNA"]]$counts
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KB.PT), binarize = T)
min.rna.cl.frac <- 0.1
tf.markers$motif.name <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters


tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
tf.markers.exp <- na.omit(tf.markers.exp)
DefaultAssay(KB.PT) <- 'chromvar'

save(tf.markers.exp, file = "regulatory/TFBS_Activities_PT-aPT.rda")
load("regulatory/TFBS_Activities_PT-aPT.rda")

#Plot TFBS Activities
tf.mark <- tf.markers.exp[tf.markers.exp$p_val_adj < 0.01,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 

DotPlot(KB.PT, features = unique(top10$motif.name), dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()





#TFBS activities - All PT subclass.l2
DefaultAssay(KB.PT) <- "chromvar"
Idents(KB.PT) <- "v2.subclass.l2"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = c("PT-S1","PT-S2","PT-S3","aPT","frPT"))
table(Idents(KB.PT))

tf.markers <- FindAllMarkers(
  object = KB.PT,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

#subset to TFs expressed in at least 10% of cells in the cluster
DefaultAssay(KB.PT) <- "RNA"
clusters <- levels(Idents(KB.PT))
rna.counts <- KB.PT[["RNA"]]$counts
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KB.PT), binarize = T)
min.rna.cl.frac <- 0.1
tf.markers$motif.name <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters


tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
tf.markers.exp <- na.omit(tf.markers.exp)
DefaultAssay(KB.PT) <- 'chromvar'

save(tf.markers.exp, file = "regulatory/TFBS_Activities_PT-aPT_l2.rda")

#Plot TFBS Activities
tf.mark <- tf.markers.exp[tf.markers.exp$p_val_adj < 0.01,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 

DotPlot(KB.PT, features = unique(top10$motif.name), dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()




### Thick Ascending Limb
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_05162024.rda")
KB.DT <- subset(KRAC.D1Epi, idents = paste0("D_", c(5:8,11:13,15:17)))

#Link motifs
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
KB.DT <- AddMotifs(
  object = KB.DT,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

KB.DT <- RunChromVAR(
  object = KB.DT,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

KB.DT
#An object of class Seurat
#624293 features across 113921 samples within 3 assays
#Active assay: ATAC (586864 features, 0 variable features)
#2 layers present: counts, data
#2 other assays present: RNA, chromvar
#2 dimensional reductions calculated: pca, umap

save(KB.DT, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_PT_Dual_Seurat_TAL-aTAL_Only_05212024.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_PT_Dual_Seurat_TAL-aTAL_Only_05212024.rda")


#TFBS activities - All TAL subclass.l3
DefaultAssay(KB.DT) <- "chromvar"
Idents(KB.DT) <- "v2.subclass.l3"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = c("aTAL2","C-TAL-A","C/M-TAL-A","C/M-TAL-B","C-TAL-B","aTAL1","frTAL"))
table(Idents(KB.DT))
rownames(KB.DT@assays$chromvar@data) <- ConvertMotifID(KB.DT,
                                                       id=rownames(KB.DT@assays$chromvar@data),
                                                       assay='ATAC')
tf.markers <- FindAllMarkers(
  object = KB.DT,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

#subset to TFs expressed in at least 10% of cells in the cluster
library(chromfunks)
DefaultAssay(KB.DT) <- "RNA"
clusters <- levels(Idents(KB.DT))
rna.counts <- KB.DT[["RNA"]]$counts
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KB.DT), binarize = T)
min.rna.cl.frac <- 0.1
tf.markers$motif.name <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters


tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
tf.markers.exp <- na.omit(tf.markers.exp)
DefaultAssay(KB.DT) <- 'chromvar'

save(tf.markers.exp, file = "regulatory/TFBS_Activities_TAL-aTAL.rda")
load("regulatory/TFBS_Activities_TAL-aTAL.rda")

#Plot TFBS Activities
tf.mark <- tf.markers.exp[tf.markers.exp$p_val_adj < 0.01,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 

DotPlot(KB.DT, features = unique(top10$motif.name), dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()




#TFBS activities - All TAL subclass.l2
DefaultAssay(KB.DT) <- "chromvar"
Idents(KB.DT) <- "v2.subclass.l2"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = c("aTAL","C-TAL","C/M-TAL","frTAL"))
table(Idents(KB.DT))

tf.markers <- FindAllMarkers(
  object = KB.DT,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

#subset to TFs expressed in at least 10% of cells in the cluster
DefaultAssay(KB.DT) <- "RNA"
clusters <- levels(Idents(KB.DT))
rna.counts <- KB.DT[["RNA"]]$counts
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KB.DT), binarize = T)
min.rna.cl.frac <- 0.1
tf.markers$motif.name <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters


tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
tf.markers.exp <- na.omit(tf.markers.exp)
DefaultAssay(KB.DT) <- 'chromvar'

save(tf.markers.exp, file = "regulatory/TFBS_Activities_TAL-aTAL_l2.rda")

#Plot TFBS Activities
tf.mark <- tf.markers.exp[tf.markers.exp$p_val_adj < 0.01,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 

DotPlot(KB.DT, features = unique(top10$motif.name), dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()






### Immune
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_05162024.rda")
KB.IMM <- subset(KRAC.nonEpi, v2.subclass.l3 %in% c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                                    "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                                    "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))

#Link motifs
pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
KB.IMM <- AddMotifs(
  object = KB.IMM,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

KB.IMM <- RunChromVAR(
  object = KB.IMM,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

KB.IMM
#An object of class Seurat 
#624293 features across 43276 samples within 3 assays 
#Active assay: ATAC (586864 features, 0 variable features)
#2 layers present: counts, data
#2 other assays present: RNA, chromvar
#2 dimensional reductions calculated: pca, umap

save(KB.IMM, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Immune_Dual_Seurat_10012024.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Immune_Dual_Seurat_10012024.rda")

#TFBS activities
DefaultAssay(KB.IMM) <- "chromvar"
rownames(KB.IMM@assays$chromvar@data) <- ConvertMotifID(KB.IMM,
                                                        id=rownames(KB.IMM@assays$chromvar@data),
                                                        assay='ATAC')
tf.markers <- FindAllMarkers(
  object = KB.IMM,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

#subset to TFs expressed in at least 10% of cells in the cluster
library(chromfunks)
DefaultAssay(KB.IMM) <- "RNA"
clusters <- levels(Idents(KB.IMM))
rna.counts <- KB.IMM[["RNA"]]$counts
rna.avg.counts <- GetClusterAvg(rna.counts, Idents(KB.IMM), binarize = T)
min.rna.cl.frac <- 0.1
tf.markers$motif.name <- tf.markers$gene
tf.markers$gene <- gsub("\\s*\\([^\\)]+\\)","",tf.markers$gene) #remove " (var*)"
tf.markers$gene <- toupper(tf.markers$gene)

tf.markers.exp <- lapply(clusters, function(cl) {
  
  cl.expr <- rna.avg.counts[,cl]
  cl.expr <- names(cl.expr[cl.expr > min.rna.cl.frac])
  tf.markers.cl <- tf.markers[tf.markers$cluster == cl,]
  tf.markers.cl <- tf.markers.cl[tf.markers.cl$gene %in% cl.expr,]
  tf.markers.cl
})
names(tf.markers.exp) <- clusters


tf.markers.exp <- do.call("rbind", lapply(tf.markers.exp, as.data.frame)) 
tf.markers.exp <- na.omit(tf.markers.exp)
DefaultAssay(KB.IMM) <- 'chromvar'

save(tf.markers.exp, file = "regulatory/TFBS_Activities_Immune.rda")

#Plot TFBS Activities
tf.mark <- tf.markers.exp[tf.markers.exp$p_val_adj < 0.01,]
tf.mark <- distinct(tf.mark, gene, .keep_all = TRUE)
tf.mark %>% group_by(cluster) %>% top_n(10, avg_diff) -> top10 

DotPlot(KB.IMM, features = unique(top10$motif.name), dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
###




