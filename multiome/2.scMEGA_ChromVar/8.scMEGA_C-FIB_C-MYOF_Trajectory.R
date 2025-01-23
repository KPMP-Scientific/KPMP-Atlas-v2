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

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_Lineage4_0424-newData.rda") 
STR.R <- KB.STR
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
cols <- colnames(STR.R@meta.data)[colnames(STR.R@meta.data) %in% colnames(KB@meta.data)]
STR.R@meta.data[,cols] <- KB@meta.data[rownames(STR.R@meta.data),cols]

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_04272024_B.rda")
STR.AC <- KRAC.nonEpi
rm(KB.STR, KRAC.nonEpi)


STR.AC <- subset(STR.AC, cells = colnames(STR.AC)[colnames(STR.AC) %in% colnames(STR.R)])
STR.AC@meta.data <- STR.R@meta.data[rownames(STR.AC@meta.data),]
embeddings <- Embeddings(STR.R, reduction = "umap")[rownames(STR.AC@meta.data),]
STR.AC[["umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "umap_", assay = DefaultAssay(STR.AC))
embeddings <- Embeddings(STR.R, reduction = "pca")[rownames(STR.AC@meta.data),]
STR.AC[["pca"]] <- CreateDimReducObject(embeddings = embeddings, key = "pca_", assay = DefaultAssay(STR.AC))

obj <- STR.AC
table(obj$v2.subclass.l3)
cols <- ArchR::paletteDiscrete(obj@meta.data[, "v2.subclass.l3"])
Idents(obj) <- "v2.subclass.l3"
v2.fn.cols <- setNames(adjustcolor(c("#9A4D42","#02AD24","#0000FF","#009FFF","#00479E","#FF0000"), alpha.f = 0.3),
                       c("S_9","S_10","S_11","S_12","S_13","S_14"))

DimPlot(obj, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.4,
        label.size = 6, repel = TRUE, pt.size = 1) + ggtitle("Subclass l3"
        ) + scale_color_manual(values = c("#0000FF4D","#9A4D424D","#FF0000")
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
        ) + scale_color_manual(values = c("#0000FF4D","#9A4D424D","#FF0000")
        ) + NoLegend()


###Use Slingshot Trajectory
table(obj$v2.subclass.l3)
TrajectoryPlot(object = obj, 
               reduction = "umap",
               trajectory = "pseudotime.lineage4",
               continuousSet = "blueYellow",
               size = 1,
               addArrow = FALSE)

###TF and gene selection
pfm <- getMatrixSet(
  x = JASPAR2020, #Note that JASPAR2022 led to duplicate motif ids that caused errors downstream
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
                     trajectory = c("S_9", "S_11", "S_14"),
                     group.by = "v2.clusters", 
                     reduction = "umap",
                     use.all = FALSE)

obj <- obj[, !is.na(obj$Trajectory)]
p <- TrajectoryPlot(object = obj, 
                    reduction = "umap",
                    continuousSet = "blueYellow",
                    size = 2,
                    addArrow = FALSE)


pdf(file='regulatory/scMEGA_S9-S14_ArchR-Trajectory.pdf',width=8,height=8)
print(p)
dev.off()



###TF selection
res <- SelectTFs(object = obj, return.heatmap = TRUE, p.cutoff = 0.01, cor.cutoff = 0.3)
df.cor <- res$tfs
ht <- res$heatmap

pdf(file='regulatory/scMEGA_SelectTFs_S9-S14_Heatmap.pdf',width=8,height=8)
draw(ht)
dev.off()

write.table(df.cor, file="regulatory/scMEGA_SelectTFs_S9-S14_Heatmap.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Select Genes
res <- SelectGenes(object = obj,
                   labelTop1 = 0,
                   labelTop2 = 0)

df.p2g <- res$p2g
ht <- res$heatmap

pdf(file='regulatory/scMEGA_SelectGenes_S9-S14_Heatmap.pdf',width=8,height=8)
draw(ht) 
dev.off()

write.table(df.p2g, file="regulatory/scMEGA_SelectGenes_S9-S14_Heatmap.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Gene regulatory network inference
tf.gene.cor <- GetTFGeneCorrelation(object = obj, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")


ht <- GRNHeatmap(tf.gene.cor, 
                 tf.timepoint = df.cor$time_point)

pdf(file='regulatory/scMEGA_SelectTFs_S9-S14_GRNHeatmap.pdf',width=18,height=6)
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


write.table(df.grn, file="regulatory/scMEGA_TF2Peak2Gene_S9-S14_Table.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


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
  

write.table(df.grn2, file="regulatory/scMEGA_TF2Peak2Gene_S9-S14_Table_cor0.4.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

p <- GRNPlot(df.grn2, 
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = FALSE,
             min.importance = 2,
             remove.isolated = FALSE)

options(repr.plot.height = 20, repr.plot.width = 20)

pdf(file='regulatory/scMEGA_SelectTFs_S9-S14_GRNPlot.pdf',width=10,height=10)
print(p)
dev.off()

save(obj, df.cor, df.grn, df.grn2, df.p2g, tf.gene.cor,
     file = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S9-S14_Trajectory_GRN_0424-newData.RDA")
#load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S9-S14_Trajectory_GRN_0424-newData.RDA")

obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay")

###GRN visualization
obj <- AddTargetAssay(object = obj, df.grn = df.grn2)

dfgrn <- df.grn2[df.grn2$weights > 0.5,]
netobj <- graph_from_data_frame(dfgrn,directed = TRUE)
V(netobj)$type <- ifelse(V(netobj)$name %in% dfgrn$tf,"TF/Gene","Gene")

pdf(file='regulatory/scMEGA_SelectTFs_S9-S14_Central-TFs.pdf',width=10,height=10)
#Network topological measures embedding: PC1 and PC2
p <- TopEmbGRN(df.grn=netobj)
#SOX4,RUNX1,GLI2,KLF3,LEF1,BACH1,TBX2,FOXP1,STAT1,ZEB1,ETV1,SMAD3,THRB,PBX3,PRRX1,RUNX2,FOXK1,PKNOX2

#Network topological measures embedding: PC2 and PC3
p <- TopEmbGRN(df.grn=netobj,axis=c(2,3))
#SOX4,RUNX1,ZEB1,TCF4,LEF1,ETV1,STAT1,BACH1,TBX3,FOXP1,TBX2,THRB,FOXK1,GLI2,SMAD3

dev.off()

###
p1 <- PseudotimePlot(object = obj, tf.use = "RUNX1")
p2 <- PseudotimePlot(object = obj, tf.use = "SOX4")

p1 + p2

p3 <- PseudotimePlot(object = obj, tf.use = "GLI2")
p4 <- PseudotimePlot(object = obj, tf.use = "ZEB1")

p3 + p4

p5 <- PseudotimePlot(object = obj, tf.use = "TCF4")
p6 <- PseudotimePlot(object = obj, tf.use = "STAT1")

p5 + p6

p7 <- PseudotimePlot(object = obj, tf.use = "LEF1")
p8 <- PseudotimePlot(object = obj, tf.use = "SMAD3")

p7 + p8

pdf(file='regulatory/scMEGA_SelectTFs_S9-S14_select-TFs_Pseudotime.pdf',width=20,height=7.5)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
dev.off()



###Network analysis
load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S9-S14_Trajectory_GRN_0424-newData.RDA")

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

write.table(df.cor, file="regulatory/scMEGA_SelectTFs_S9-S14_Heatmap_Network_Scores_Cor0-4.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
