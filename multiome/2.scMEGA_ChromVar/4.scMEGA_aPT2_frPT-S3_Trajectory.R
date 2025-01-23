library(ArchR)
library(Seurat)
library(Signac)
library(scMEGA)
library(harmony)
library(Nebulosa)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(igraph)
library(ggraph)

set.seed(1234)
setwd("~/Projects/Human_Kidney/Atlas_V2")

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_Lineage5_0424-newData.rda") 
PT.R <- KB.PT
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
cols <- colnames(PT.R@meta.data)[colnames(PT.R@meta.data) %in% colnames(KB@meta.data)]
PT.R@meta.data[,cols] <- KB@meta.data[rownames(PT.R@meta.data),cols]

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_PT-TL_Dual_Seurat_all-peaks_05162024.rda")
PT.AC <- KRAC.PEpi
rm(KB.PT, KRAC.PEpi)


PT.AC <- subset(PT.AC, cells = colnames(PT.AC)[colnames(PT.AC) %in% colnames(PT.R)])
PT.AC@meta.data <- PT.R@meta.data[rownames(PT.AC@meta.data),]
embeddings <- Embeddings(PT.R, reduction = "umap")[rownames(PT.AC@meta.data),]
PT.AC[["umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "umap_", assay = DefaultAssay(PT.AC))
embeddings <- Embeddings(PT.R, reduction = "pca")[rownames(PT.AC@meta.data),]
PT.AC[["pca"]] <- CreateDimReducObject(embeddings = embeddings, key = "pca_", assay = DefaultAssay(PT.AC))

obj <- PT.AC
table(obj$v2.subclass.l3)
cols <- ArchR::paletteDiscrete(obj@meta.data[, "v2.subclass.l3"])
Idents(obj) <- "v2.subclass.l3"
v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#FE8F42","#DC5E93","#A10300","#858567","#783FC1","#FFD300","#9A4D42"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

DimPlot(obj, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.4,
        label.size = 6, repel = TRUE, pt.size = 1) + ggtitle("Subclass l3"
        ) + scale_color_manual(values = c("#1F9698","#B1CC71","#A10300")
        ) + NoLegend()
CellPropPlot(obj,
             group.by = "v2.subclass.l3",
             prop.in = "patient",
             cols = cols)
CellPropPlot(obj,
             group.by = "v2.subclass.l3",
             prop.in = "condition_level3",
             cols = cols)

obj <- RunDiffusionMap(obj, reduction = "pca", dims = 1:20, k = 20)
DimPlot(obj, reduction = "dm", label = TRUE, raster=FALSE,alpha = 0.4,
        label.size = 6, repel = TRUE, pt.size = 1) + ggtitle("Subclass l3"
        ) + scale_color_manual(values = c("#1F9698","#B1CC71","#A10300")
        ) + NoLegend()


###Use Slingshot Trajectory
table(obj$v2.subclass.l3)
TrajectoryPlot(object = obj, 
               reduction = "umap",
               trajectory = "pseudotime.lineage5",
               continuousSet = "blueYellow",
               size = 1,
               addArrow = FALSE)

###TF and gene selection
pfm <- getMatrixSet(
  x = JASPAR2020, 
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm,
  assay = "ATAC"
)

obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "ATAC"
)


###Add trajectory based on ArchR
obj <- AddTrajectory(object = obj, 
                     trajectory = c("aPT2","aPT1","frPT-S3"),
                     group.by = "v2.subclass.l3", 
                     reduction = "umap",
                     use.all = FALSE)

obj <- obj[, !is.na(obj$Trajectory)]
p <- TrajectoryPlot(object = obj, 
                    reduction = "umap",
                    continuousSet = "blueYellow",
                    size = 1,
                    addArrow = FALSE)


pdf(file='regulatory/scMEGA_aPT-frPTS3_l5_ArchR-Trajectory.pdf',width=8,height=8)
print(p)
dev.off()



###TF selection
res <- SelectTFs(object = obj, return.heatmap = TRUE, p.cutoff = 0.01, cor.cutoff = 0.3)
df.cor <- res$tfs
ht <- res$heatmap

pdf(file='regulatory/scMEGA_SelectTFs_aPT-frPTS3_l5_Heatmap.pdf',width=8,height=8)
draw(ht)
dev.off()

write.table(df.cor, file="regulatory/scMEGA_SelectTFs_aPT-frPTS3_l5_Heatmap.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Select Genes
res <- SelectGenes(object = obj,
                   labelTop1 = 0,
                   labelTop2 = 0)

df.p2g <- res$p2g
ht <- res$heatmap

pdf(file='regulatory/scMEGA_SelectGenes_aPT-frPTS3_l5_Heatmap.pdf',width=8,height=8)
draw(ht) 
dev.off()

write.table(df.p2g, file="regulatory/scMEGA_SelectGenes_aPT-frPTS3_l5_Heatmap.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Gene regulatory network inference
tf.gene.cor <- GetTFGeneCorrelation(object = obj, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")


ht <- GRNHeatmap(tf.gene.cor, 
                 tf.timepoint = df.cor$time_point)

pdf(file='regulatory/scMEGA_SelectTFs_aPT-frPTS3_l5_GRNHeatmap.pdf',width=18,height=6)
ht
dev.off()


#associate genes to TFs
motif.matching <- obj@assays$ATAC@motifs@data
colnames(motif.matching) <- obj@assays$ATAC@motifs@motif.names

motif.matching <-
  motif.matching[unique(df.p2g$peak), unique(tf.gene.cor$tf)]


df.grn <- GetGRN(motif.matching = motif.matching, 
                 df.cor = tf.gene.cor, 
                 df.p2g = df.p2g)


write.table(df.grn, file="regulatory/scMEGA_TF2Peak2Gene_aPT-frPTS3_l5_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


#Visualize network
# define colors for nodes representing TFs (i.e., regulators)
df.cor <- df.cor[order(df.cor$time_point), ]
tfs.timepoint <- df.cor$time_point
names(tfs.timepoint) <- df.cor$tfs

# plot the graph, here we can highlight some genes
df.grn2 <- df.grn %>%
  subset(correlation > 0.4) %>%
  select(c(tf, gene, correlation)) %>%
  dplyr::rename(weights = correlation)


write.table(df.grn2, file="regulatory/scMEGA_TF2Peak2Gene_aPT-frPTS3_l5_Table_cor0.4.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

p <- GRNPlot(df.grn2, 
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = FALSE,
             min.importance = 2,
             remove.isolated = FALSE)

options(repr.plot.height = 20, repr.plot.width = 20)

pdf(file='regulatory/scMEGA_SelectTFs_aPT-frPTS3_l5_GRNPlot.pdf',width=10,height=10)
print(p)
dev.off()

save(obj, df.cor, df.grn, df.grn2, df.p2g, tf.gene.cor,
     file = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-frPTS3_l5_Trajectory_GRN_0424-newData.RDA")
#load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-frPTS3_l5_Trajectory_GRN_0424-newData.RDA")

obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay")

###GRN visualization
obj <- AddTargetAssay(object = obj, df.grn = df.grn2)

dfgrn <- df.grn2[df.grn2$weights > 0.6,]
netobj <- graph_from_data_frame(dfgrn,directed = TRUE)
V(netobj)$type <- ifelse(V(netobj)$name %in% dfgrn$tf,"TF/Gene","Gene")

pdf(file='regulatory/scMEGA_SelectTFs_aPT-frPTS3_l5_Central-TFs.pdf',width=10,height=10)
#Network topological measures embedding: PC1 and PC2
p <- TopEmbGRN(df.grn=netobj)
#MEIS2, HOXB7, EHF, POU3F3, PBX3, GRHL2,CEBPD,EMX2,RUNX1,RUNX2,KLF2

#Network topological measures embedding: PC2 and PC3
p <- TopEmbGRN(df.grn=netobj,axis=c(2,3))
#HOXB7,EMX2,NR2F1,CEBPD,EHF

dev.off()

###
p1 <- PseudotimePlot(object = obj, tf.use = "EHF")
p2 <- PseudotimePlot(object = obj, tf.use = "SOX9")

p1 + p2

p3 <- PseudotimePlot(object = obj, tf.use = "EMX2")
p4 <- PseudotimePlot(object = obj, tf.use = "MEIS2")

p3 + p4

p5 <- PseudotimePlot(object = obj, tf.use = "HOXB7")
p6 <- PseudotimePlot(object = obj, tf.use = "CEBPD")

p5 + p6

p7 <- PseudotimePlot(object = obj, tf.use = "HNF4A")
p8 <- PseudotimePlot(object = obj, tf.use = "NR2F1")

p7 + p8

pdf(file='regulatory/scMEGA_SelectTFs_aPT-frPTS3_l5_select-TFs_Pseudotime.pdf',width=20,height=7.5)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
dev.off()






###Network analysis
load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-frPTS3_l5_Trajectory_GRN_0424-newData.RDA")

df.grn2 <- df.grn %>%
  subset(correlation > 0.4) %>%
  select(c(tf, gene, correlation)) %>%
  dplyr::rename(weights = correlation)

netobj <- graph_from_data_frame(df.grn2,directed = TRUE)
V(netobj)$type <- ifelse(V(netobj)$name %in% df.grn2$tf,"TF/Gene","Gene")

outdegree <- igraph::degree(netobj, mode = c("out"))
centrality <- igraph::betweenness(netobj, directed = TRUE)
indegree <- igraph::degree(netobj, mode = c("in"))
eigen.centrality <- eigen_centrality(netobj)$vector

df.cor$outdegree <- outdegree[match(df.cor$tfs, names(outdegree))]
df.cor$indegree <- indegree[match(df.cor$tfs, names(indegree))]
df.cor$between.centrality <- centrality[match(df.cor$tfs, names(centrality))]
df.cor$eigen.centrality <- eigen.centrality[match(df.cor$tfs, names(eigen.centrality))]

write.table(df.cor, file="regulatory/scMEGA_SelectTFs_aPT-frPTS3_l5_Heatmap_Network_Scores_Cor0-4.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Cluster Umaps
Idents(obj) <- "v2.clusters"
v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#FE8F42","#DC5E93","#A10300","#858567","#783FC1","#FFD300","#9A4D42"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

pdf(file='regulatory/scMEGA_aPT-frPTS3_l5_ArchR-Trajectory_clusters_labeled.pdf',width=8,height=6)
DimPlot(obj, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.4, pt.size = 1,
        label.size = 4, repel = TRUE) + ggtitle("Clusters"
        ) + scale_color_manual(values = alpha(v2.fn.cols[levels(Idents(obj))], 0.1), name = "Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) #+ NoLegend()
dev.off()
pdf(file='regulatory/scMEGA_aPT-frPTS3_l5_ArchR-Trajectory_clusters_unlabeled.pdf',width=8,height=6)
DimPlot(obj, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.4, pt.size = 1,
        label.size = 4, repel = TRUE) + ggtitle("Clusters"
        ) + scale_color_manual(values = alpha(v2.fn.cols[levels(Idents(obj))], 0.1), name = "Clusters"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) #+ NoLegend()
dev.off()
