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

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_Lineage1_0424-newData_TALA.rda")
TAL.R <- KB.TAL
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
cols <- colnames(TAL.R@meta.data)[colnames(TAL.R@meta.data) %in% colnames(KB@meta.data)]
TAL.R@meta.data[,cols] <- KB@meta.data[rownames(TAL.R@meta.data),cols]

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_05162024.rda")
TAL.AC <- KRAC.D1Epi
rm(KB.TAL, KRAC.D1Epi)


TAL.AC <- subset(TAL.AC, cells = colnames(TAL.AC)[colnames(TAL.AC) %in% colnames(TAL.R)])
TAL.AC@meta.data <- TAL.R@meta.data[rownames(TAL.AC@meta.data),]
embeddings <- Embeddings(TAL.R, reduction = "umap")[rownames(TAL.AC@meta.data),]
TAL.AC[["umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "umap_", assay = DefaultAssay(TAL.AC))
embeddings <- Embeddings(TAL.R, reduction = "pca")[rownames(TAL.AC@meta.data),]
TAL.AC[["pca"]] <- CreateDimReducObject(embeddings = embeddings, key = "pca_", assay = DefaultAssay(TAL.AC))

obj <- TAL.AC
table(obj$v2.subclass.l3)
cols <- ArchR::paletteDiscrete(obj@meta.data[, "v2.subclass.l3"])
Idents(obj) <- "v2.subclass.l3"
v2.fn.cols <- setNames(c("#005300","#B1CC71","#886C00","#02AD24","#009FFF","#14F9FF","#1F9698",
                         "#DC5E93","#FFD300","#FE8F42"),
                       c("D_5","D_6","D_7","D_8","D_11","D_12","D_13",
                         "D_15","D_16","D_17"))

DimPlot(obj, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.4,
        label.size = 6, repel = TRUE, pt.size = 1) + ggtitle("Subclass l3"
        ) + scale_color_manual(values = c("#FE8F42","#005300","#FFD300")
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
        ) + scale_color_manual(values = c("#FE8F42","#005300","#FFD300")
        ) + NoLegend()


###Use Slingshot Trajectory
table(obj$v2.subclass.l3)
TrajectoryPlot(object = obj, 
               reduction = "umap",
               trajectory = "pseudotime.lineage1",
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
                     trajectory = c("aTAL1","aTAL2","C/M-TAL-A"),
                     group.by = "v2.subclass.l3", 
                     reduction = "umap",
                     use.all = FALSE)

obj <- obj[, !is.na(obj$Trajectory)]
p <- TrajectoryPlot(object = obj, 
                    reduction = "umap",
                    continuousSet = "blueYellow",
                    size = 1,
                    addArrow = FALSE)


pdf(file='regulatory/scMEGA_aTAL-C-TAL-A_ArchR-Trajectory_TALA.pdf',width=8,height=8)
print(p)
dev.off()




###TF selection
res <- SelectTFs(object = obj, return.heatmap = TRUE, p.cutoff = 0.01, cor.cutoff = 0.3)
df.cor <- res$tfs
ht <- res$heatmap

pdf(file='regulatory/scMEGA_SelectTFs_aTAL-C-TAL-A_Heatmap_TALA.pdf',width=8,height=8)
draw(ht)
dev.off()

write.table(df.cor, file="regulatory/scMEGA_SelectTFs_aTAL-C-TAL-A_Heatmap_TALA.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Select Genes
res <- SelectGenes(object = obj,
                   labelTop1 = 0,
                   labelTop2 = 0)

df.p2g <- res$p2g
ht <- res$heatmap

pdf(file='regulatory/scMEGA_SelectGenes_aTAL-C-TAL-A_Heatmap_TALA.pdf',width=8,height=8)
draw(ht) 
dev.off()

write.table(df.p2g, file="regulatory/scMEGA_SelectGenes_aTAL-C-TAL-A_Heatmap_TALA.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Gene regulatory network inference
tf.gene.cor <- GetTFGeneCorrelation(object = obj, 
                                    tf.use = df.cor$tfs, 
                                    gene.use = unique(df.p2g$gene),
                                    tf.assay = "chromvar", 
                                    gene.assay = "RNA",
                                    trajectory.name = "Trajectory")


ht <- GRNHeatmap(tf.gene.cor, 
                 tf.timepoint = df.cor$time_point)

pdf(file='regulatory/scMEGA_SelectTFs_aTAL-C-TAL-A_GRNHeatmap_TALA.pdf',width=18,height=6)
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


write.table(df.grn, file="regulatory/scMEGA_TF2Peak2Gene_aTAL-C-TAL-A_Table_TALA.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


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


write.table(df.grn2, file="regulatory/scMEGA_TF2Peak2Gene_aTAL-C-TAL-A_Table_cor0.4_TALA.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

p <- GRNPlot(df.grn2, 
             tfs.timepoint = tfs.timepoint,
             show.tf.labels = TRUE,
             seed = 42, 
             plot.importance = FALSE,
             min.importance = 2,
             remove.isolated = FALSE)

options(repr.plot.height = 20, repr.plot.width = 20)

pdf(file='regulatory/scMEGA_SelectTFs_aTAL-C-TAL-A_GRNPlot_TALA.pdf',width=10,height=10)
print(p)
dev.off()

save(obj, df.cor, df.grn, df.grn2, df.p2g, tf.gene.cor,
     file = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aTAL-C-TAL-A_Trajectory_GRN_0424-newData_TALA.RDA")
#load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aTAL-C-TAL-A_Trajectory_GRN_0424-newData_TALA.RDA")

obj[["RNA"]] <- as(obj[["RNA"]], Class = "Assay")

###GRN visualization
obj <- AddTargetAssay(object = obj, df.grn = df.grn2)

dfgrn <- df.grn2[df.grn2$weights > 0.8,]
netobj <- graph_from_data_frame(dfgrn,directed = TRUE)
V(netobj)$type <- ifelse(V(netobj)$name %in% dfgrn$tf,"TF/Gene","Gene")

pdf(file='regulatory/scMEGA_SelectTFs_aTAL-C-TAL-A_Central-TFs_TALA.pdf',width=10,height=10)
#Network topological measures embedding: PC1 and PC2
p <- TopEmbGRN(df.grn=netobj)
#BACH1,FOSL2, IRF1,HNF1B,GLIS3,ETS1,TGIF1,ELK3,PKNOX1,IRF3,RELB,CREM,PPARG,CEBPB,
#HMBOX1,ATF7,CREM,ETS2,RFX7,RFX2,CUX1,ESR2,KLF6,ETV5,TEAD4,STAT3

#Network topological measures embedding: PC2 and PC3
p <- TopEmbGRN(df.grn=netobj,axis=c(2,3))
#STAT1,ESR2,RFX2,CEBPB,SMAD3,HNF1B,ETS1,ELF1,TGIF1,SOX4,ELK3,PPARG,NFAT5,ELK3,ATF7,TEF
dev.off()

###
p1 <- PseudotimePlot(object = obj, tf.use = "CREM")
p2 <- PseudotimePlot(object = obj, tf.use = "STAT1")

p1 + p2

p3 <- PseudotimePlot(object = obj, tf.use = "ELF3")
p4 <- PseudotimePlot(object = obj, tf.use = "PPARG")

p3 + p4

p5 <- PseudotimePlot(object = obj, tf.use = "IRF1")
p6 <- PseudotimePlot(object = obj, tf.use = "SOX9")

p5 + p6

p7 <- PseudotimePlot(object = obj, tf.use = "ESR2")
p8 <- PseudotimePlot(object = obj, tf.use = "ESRRB")

p7 + p8


pdf(file='regulatory/scMEGA_SelectTFs_aTAL-C-TAL-A_select-TFs_Pseudotime_TALA.pdf',width=20,height=7.5)
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
dev.off()




###Network analysis
load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aTAL-C-TAL-A_Trajectory_GRN_0424-newData_TALA.RDA")

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

write.table(df.cor, file="regulatory/scMEGA_SelectTFs_aTAL-C-TAL-A_Heatmap_Network_Scores_Cor0-4_TALA.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
