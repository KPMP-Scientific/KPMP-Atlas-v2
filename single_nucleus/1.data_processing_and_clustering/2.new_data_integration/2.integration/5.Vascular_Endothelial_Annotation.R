library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

nKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_Atlas_V2_04-2024-subset_Object_global-predictions.Rds")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_SCT_Reference.Rda")

nKB <- subset(nKB, predicted.v2.subclass.l1 %in% c("EC")) #subset to stromal cells
nKB <- subset(nKB, predicted.v2.subclass.l1.score > 0.8) #subset to high confidence stromal cells


###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KB,
  query = nKB,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

nKB <- MapQuery(
  anchorset = anchors,
  query = nKB,
  reference = KB,
  refdata = list(
    v2.subclass.l3 = "v2.subclass.l3",
    v2.clusters = "str.v2.clusters"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(nKB, reduction = "ref.umap", group.by = "predicted.v2.subclass.l3", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(nKB, reduction = "ref.umap", group.by = "predicted.v2.clusters", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()


hist(nKB$predicted.v2.subclass.l3.score)
hist(nKB$predicted.v2.clusters.score)


#remove prediction assays
nKB[["prediction.score.v2.subclass.l3"]] <- NULL
nKB[["prediction.score.v2.clusters"]] <- NULL
nKB[["pca"]] <- nKB[["ref.pca"]]
nKB[["umap"]] <- nKB[["ref.umap"]]
nKB[["ref.pca"]] <- NULL
nKB[["ref.umap"]] <- NULL


##Merge with full object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
#update cluster metadata
meta <- KB.EC@meta.data
meta$ec.v2.clusters <- paste0("E_",meta$ec.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$ec.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.EC@meta.data <- meta
KB.EC$v2.clusters <- KB$ec.v2.clusters
KB.EC$group <- "ref"

nKB$group <- "query"
KB <- merge(x = KB.EC, y = nKB, merge.dr = TRUE)
KB <- JoinLayers(KB)
KB <- NormalizeData(KB)
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "group", repel = TRUE) + NoLegend()


##Clustering in Pagoda2 using integrated PCs
#Filtered Count Matrix from Seurat
#Filtered Count Matrix from Seurat
countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#2987 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 2987, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2987]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_ec <- k100infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_ec", repel = TRUE) + NoLegend()

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Vasculature_k100_round1_newData-0424.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Vasculature_k100_round1_newData-0424.rda")

###Identify ambiguous low quality - annotate cluster identities
Idents(KB) <- "pagoda_k100_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v2.clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.v2.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v2.clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v2.clusters)),
                                query = prop.table(table(seurat.obj.sub$predicted.v2.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$predicted.v2.subclass.l3)))
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

max.idents <- getMaxIdents(KB, celltype)
max.idents
write.table(max.idents, file="~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Vasculature_k100_round1_newData-0424_overlaps.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Vasculature_round1_newData-0424_k100-mean-gene_vlnplot.pdf',width=8,height=4)
VlnPlot(KB, features = c("nFeature_RNA"), pt.size = 0)
dev.off()

###Check for tissue specific markers
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
                 "PHACTR1","SLC14A2",             #IMCD
                 "TP63","GPX2",           #M-CBC
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
                 "PTPRC",                                                        #IMM
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


pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Vasculature_round1_newData-0424_k100-marker-gene_dotplot.pdf',width=16,height=8)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
dev.off()

CL35_markers <- FindMarkers(object = KB, ident.1 = c("CL35"), 
                            only.pos = TRUE)


###Ambiguous clusters (checking mean genes, markers, DEGs)
#CL 10 - PT Multiplet
#CL 13 - PC Multiplet
#CL 18 - FIB Multiplet
#CL 21 - TAL Multiplet
#CL 28 - PT Multiplet
#CL 35 - T Multiplet
#CL 41 - TAL Multiplet
#CL 42 - PC Multiplet


###Subset multiplets
to.remove <- c("10","13","18","21","28","35","41","42")
to.remove.query <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_ec %in% to.remove &
                                           KB@meta.data$group %in% "query",])
KB <- subset(KB, cells = to.remove.query, invert = TRUE)
KB
#36588 features across 83782 samples

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")



###Repeat subclustering
#Filtered Count Matrix from Seurat
countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#2896 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 2896, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2896]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_ec <- k100infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k200_infomap_ec <- k200infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k500_infomap_ec <- k500infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.clusters", repel = TRUE) + NoLegend()
FeaturePlot(KB, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Vasculature_k100_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Vasculature_k200_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Vasculature_k500_0424-newData.rda")


save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")




###Annotation of k200 clusters
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")
Idents(KB) <- "pagoda_k200_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))
#remove CL29 of 1 cell
KB <- subset(KB, idents = "29", invert = TRUE)

levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#CL 1: E_2
#CL 2: E_13
#CL 3: E_15
#CL 4: E_6
#CL 5: E_5
#CL 6: E_21
#CL 7: E_15
#CL 8: E_1
#CL 9: E_1
#CL 10: E_14
#CL 11: E_15
#CL 12: E_4
#CL 13: E_16
#CL 14: E_7
#CL 15: E_17
#CL 16: E_11
#CL 17: E_14
#CL 18: E_17
#CL 19: E_22
#CL 20: E_18
#CL 21: E_19
#CL 22: E_19
#CL 23: E_20
#CL 24: E_9
#CL 25: E_12
#CL 26: E_9
#CL 27: E_8
#CL 28: E_10
#CL 29: E_3 (k100 cluster 21)

#missing cluster E_3, check k100
Idents(KB) <- "pagoda_k100_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#Check clusters for patient/region/cluster distribution
table(KB$v2.clusters ,KB$pagoda_k200_infomap_ec)[c("E_3"),] 
table(KB$pagoda_k200_infomap_ec,KB$pagoda_k100_infomap_ec)[c("16"),] #16 splits to 2 clusters at k100 (17,21)
table(KB$pagoda_k100_infomap_ec,KB$patient)[c(21),]
table(KB$pagoda_k100_infomap_ec,KB$region_level2)[c(21),]
table(KB$v2.clusters,KB$pagoda_k100_infomap_ec)[c("E_3"),] #K100 21 = E_3

#Add k100 cluster to k200 clusters (E_3)
cl29 <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_ec %in% "21",])
KB@meta.data$pagoda_k200_infomap_ec <- as.character(KB@meta.data$pagoda_k200_infomap_ec)
KB@meta.data[cl29,]$pagoda_k200_infomap_ec <- "29"
Idents(KB) <- "pagoda_k200_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))


VlnPlot(KB, features = c("nFeature_RNA"), pt.size = 0)
VlnPlot(KB, features = c("percent.mt"), pt.size = 0)
VlnPlot(KB, features = c("percent.er"), pt.size = 0)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = c('PTCHD4',"ZMAT3","B2M", "TMSB4X", "TMSB10", "HLA-B", "IGFBP5"),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = c("EDN1", "GJA4", "CLDN5",                              #Afferent arterioles
                         "KCNN4", "S1PR1", "CXCL12",
                         "KLF4", "SLC6A6"                                     #Efferent arterioles
                         ),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()


#Re-assign cluster numbers for query cells
KB$v2.clusters <- KB$ec.v2.clusters
meta <- KB@meta.data[KB@meta.data$group %in% "query",]
meta$v2.clusters <- meta$pagoda_k200_infomap_ec

cls <- c(1:29)
new.cls <- c(2,13,15,6,5,21,15,1,1,14,15,4,16,7,17,11,14,17,22,18,19,19,20,9,12,9,8,10,3)
meta$v2.clusters <- plyr::mapvalues(meta$v2.clusters, from = cls, to = new.cls)
table(as.character(meta$v2.clusters))

KB@meta.data[rownames(meta),]$v2.clusters <- paste0("E_",meta$v2.clusters)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = paste0("E_",1:22))

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_ec", repel = TRUE) + NoLegend()

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "ec.v2.clusters", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.clusters", repel = TRUE) + NoLegend()

###Update metadata
meta <- KB@meta.data
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB@meta.data <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","set",                           
                        "v1_barcode","source","assay","experiment_long","experiment","patient",                       
                        "specimen","condition_level3","condition_level2","condition_level1","condition",                                      
                        "percent_cortex","percent_medulla","region_level2","region_level1","age","age_binned","sex",                           
                        "race","location","laterality","protocol","tissue_type_full","tissue_type",                   
                        "atlas_version","v1clusters","v1.subclass.l3","v1.subclass.l2","v1.subclass.l1",
                        "v1.state.l1","v1.state.l2","v1.class","v1.structure","group","predicted.v2.subclass.l1.score",
                        "predicted.v2.subclass.l1","predicted.v2.subclass.l3.score","predicted.v2.subclass.l3",
                        "predicted.v2.clusters.score","predicted.v2.clusters","pagoda_k100_infomap_ec",
                        "pagoda_k200_infomap_ec","pagoda_k500_infomap_ec","v2.clusters","v2.subclass.full","v2.subclass.l3",
                        "v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1",
                        "v2.class","v2.structure")]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")
KB
#36588 features across 83781 samples


###Remove unnecessary slots
KB[["RNA"]]$scale.data.1 <- NULL
KB[["RNA"]]$scale.data <- NULL
KB[["ref.pca"]] <- NULL
KB[["ref.umap"]] <- NULL

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered_0424-newData.rda")
KB

