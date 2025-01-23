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
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_SCT_Reference_full.Rda")

nKB.STR <- subset(nKB, predicted.v2.subclass.l1 %in% c("FIB","VSM/P","Ad")) #subset to stromal cells
nKB.STR <- subset(nKB.STR, predicted.v2.subclass.l1.score > 0.8) #subset to high confidence stromal cells


###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KB,
  query = nKB.STR,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

nKB.STR <- MapQuery(
  anchorset = anchors,
  query = nKB.STR,
  reference = KB,
  refdata = list(
    v2.subclass.l3 = "v2.subclass.l3",
    v2.clusters = "str.v2.clusters"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(nKB.STR, reduction = "ref.umap", group.by = "predicted.v2.subclass.l3", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(nKB.STR, reduction = "ref.umap", group.by = "predicted.v2.clusters", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()


hist(nKB.STR$predicted.v2.subclass.l3.score)
hist(nKB.STR$predicted.v2.clusters.score)


#remove prediction assays
nKB.STR[["prediction.score.v2.subclass.l3"]] <- NULL
nKB.STR[["prediction.score.v2.clusters"]] <- NULL
nKB.STR[["pca"]] <- nKB.STR[["ref.pca"]]
nKB.STR[["umap"]] <- nKB.STR[["ref.umap"]]
nKB.STR[["ref.pca"]] <- NULL
nKB.STR[["ref.umap"]] <- NULL


##Merge with full object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
#update cluster metadata
meta <- KB.STR@meta.data
meta$str.v2.clusters <- paste0("S_",meta$str.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$str.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.STR@meta.data <- meta
KB.STR$v2.clusters <- KB$str.v2.clusters
KB.STR$group <- "ref"

nKB.STR$group <- "query"
KB <- merge(x = KB.STR, y = nKB.STR, merge.dr = TRUE)
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
#3792 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3792, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3792]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_str <- k100infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_str", repel = TRUE) + NoLegend()

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k100_round1_newData-0424.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k100_round1_newData-0424.rda")

###Identify ambiguous low quality - annotate cluster identities
Idents(KB) <- "pagoda_k100_infomap_str"
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
write.table(max.idents, file="~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k100_round1_newData-0424_overlaps.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Stroma_round1_newData-0424_k100-mean-gene_vlnplot.pdf',width=8,height=4)
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


pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Stroma_round1_newData-0424_k100-marker-gene_dotplot.pdf',width=16,height=8)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
dev.off()

CL35_markers <- FindMarkers(object = KB, ident.1 = c("CL35"), 
                            only.pos = TRUE)


###Ambiguous clusters (checking mean genes, markers, DEGs)
#CL 17 - PT multiplet
#CL 20 - MAC multiplet
#CL 26 - EC multiplet
#CL 29 - T Multiplet
#CL 30 - epi multiplet
#CL 44 - PC multiplet
#CL 48 - B multiplet
#CL 52 - epi multiplet

###Subset multiplets
to.remove <- c("17","20","26","29","30","44","48","52")
to.remove.query <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_str %in% to.remove &
                                           KB@meta.data$group %in% "query",])
KB <- subset(KB, cells = to.remove.query, invert = TRUE)
KB
#36588 features across 102071 samples

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")



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
#3629 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3629, maxit=1000)

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
sn.od.genes <- rownames(var.info)[1:3629]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_str <- k100infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k200_infomap_str <- k200infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k500_infomap_str <- k500infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k100_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k200_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k500_0424-newData.rda")


save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")




###Annotation of k200 clusters
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
Idents(KB) <- "pagoda_k200_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))
#remove CL32 of 2 cells
KB <- subset(KB, idents = "32", invert = TRUE)

levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#CL 1: S_11
#CL 2: S_10
#CL 3: S_16
#CL 4: S_12
#CL 5: S_18
#CL 6: S_7
#CL 7: S_9
#CL 8: S_14
#CL 9: S_25
#CL 10: S_15
#CL 11: S_24
#CL 12: S_20
#CL 13: S_13
#CL 14: S_9
#CL 15: S_21
#CL 16: S_19
#CL 17: S_22
#CL 18: S_26
#CL 19: S_17
#CL 20: S_23
#CL 21: no query cells - ambiguous mixture
#CL 22: S_5
#CL 23: Imm multiplet
#CL 24: S_1
#CL 25: S_27
#CL 26: Ambiguous (indistinct marker genes)
#CL 27: S_2
#CL 28: S_4
#CL 29: S_1
#CL 30: S_8
#CL 31: S_8
#CL 32: S_3 (k100 cl 38)
#CL 33: S_6 (k100 cl 45)

#missing S_3 and S_6

#Check K100 clusters
Idents(KB) <- "pagoda_k100_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
max.idents <- getMaxIdents(KB, celltype)
max.idents


Idents(KB) <- "pagoda_k200_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
VlnPlot(KB, features = c("nFeature_RNA"), pt.size = 0)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

#Check clusters 10 and 24 for patient/region distribution
table(KB$pagoda_k200_infomap_str,KB$patient)[c(10,24,26),]
table(KB$pagoda_k200_infomap_str,KB$region_level2)[c(10,24),]
#Check for S_3 and S_6
table(KB$v2.clusters ,KB$pagoda_k200_infomap_str)[c("S_3","S_6"),]
table(KB$pagoda_k200_infomap_str ,KB$pagoda_k100_infomap_str)[c(22,24,28),]
table(KB$v2.clusters ,KB$pagoda_k100_infomap_str)[c("S_3"),] #K100 cl38 = S_3
table(KB$pagoda_k200_infomap_str ,KB$pagoda_k100_infomap_str)[c(10),]
table(KB$v2.clusters ,KB$pagoda_k100_infomap_str)[c("S_6"),] #K100 cl45 = S_6
table(KB$pagoda_k100_infomap_str,KB$region_level2)[c(38,45),]

#Assign S_3 and S_6
cl32 <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_str %in% "38",])
KB@meta.data$pagoda_k200_infomap_str <- as.character(KB@meta.data$pagoda_k200_infomap_str)
KB@meta.data[cl32,]$pagoda_k200_infomap_str <- "32"
cl33 <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_str %in% "45",])
KB@meta.data[cl33,]$pagoda_k200_infomap_str <- "33"
Idents(KB) <- "pagoda_k200_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))



### Remove ambiguous clusters
#remove cl 21, 23, 26
Idents(KB) <- "pagoda_k200_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
KB <- subset(KB, idents = c(21, 23, 26), invert = TRUE)

#Re-assign cluster numbers for query cells
KB$v2.clusters <- KB$str.v2.clusters
meta <- KB@meta.data[KB@meta.data$group %in% "query",]
meta$v2.clusters <- meta$pagoda_k200_infomap_str

cls <- c(1:20,22,24,25,27:33)
new.cls <- c(11,10,16,12,18,7,9,14,25,15,24,20,13,9,21,19,22,26,17,23,5,1,27,2,4,1,8,8,3,6)
meta$v2.clusters <- plyr::mapvalues(meta$v2.clusters, from = cls, to = new.cls)
table(as.character(meta$v2.clusters))

KB@meta.data[rownames(meta),]$v2.clusters <- paste0("S_",meta$v2.clusters)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = paste0("S_",1:27))

KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.1, return.model = T)

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_str", repel = TRUE) + NoLegend()

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "str.v2.clusters", repel = TRUE) + NoLegend()
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
                       "predicted.v2.clusters.score","predicted.v2.clusters","pagoda_k100_infomap_str",
                       "pagoda_k200_infomap_str","pagoda_k500_infomap_str","v2.clusters","v2.subclass.full","v2.subclass.l3",
                       "v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1",
                       "v2.class","v2.structure")]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
KB
#36588 features across 101818 samples



###Remove unnecessary slots
KB[["RNA"]]$scale.data.1 <- NULL
KB[["RNA"]]$scale.data <- NULL
KB[["ref.pca"]] <- NULL
KB[["ref.umap"]] <- NULL

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_0424-newData.rda")
KB




###Check stromal markers
#Fibroblasts
fib.markers <- c(
  "DCN","C7",                                                  #Pan FIB markers
  "SYT1","TNC","CA8",                                          #Pan M-FIB markers
  "MCTP2","SNAP25","BDKRB2",                                   #M-FIB
  "IGFBP2","RARRES2",                                          #M-FIB
  "GABRG3","ROBO1","PTCHD4",                                   #M-FIB
  "FREM1","KCNQ3","LOXHD1",                                    #M-FIB
  "AKR1B1","S100A6","SPP1","TIMP1",                            #dM-FIB
  "COL1A2","COL3A1","COL1A1",
  "ADAMTSL1","KCNIP1","ADARB2",
  "ZNF536","SEMA3D","ACTG2",                                   #M-MYOF
  
  "NEGR1","LAMA2","MEG3",                                      #C-FIB
  "LONRF1","SAMHD1","PDGFRA",
  "RIPOR3","ABCA8","GRIN2A","EMID1",
  
  "OSMR", #OSM in inf MAC to OSMR in inf FIB https://doi.org/10.1038/nm.4307.
  "SOD2","UGCG","IL1R1","CCL2",
  "RELB","CXCL10","CCL19", #Buechler (lymph node and spleen); Korsunsky et al cross-tissue T lymphocyte niche inflammatory fibroblast
  "CXCL12",'SELENOP','GGT5',"CCN2","ECRG4",
  "SULF1","NTM","INHBA","FAP","POSTN","SPARC", #Myofib markers
  "BGN","VIM",
  "FLRT2",'IGF1',"FGF14","COL12A1",'ADAMTS3',
  "RSPO3","WNT5B","WNT2B","CHODL",#"OGN","SFRP4", #Smillie et al - epithelial and stem-cell niche supportive 
  "C3","EBF2","SFRP2",
  "CD34","PI16", #Buechler et al - plasticity to give rise to other fib clusters, adventitial, pan-tissue
  'GPC6','PLD5',"MGAT4C","EPHA3","CNTN5",
  
  "MYH11","SYNPO2","ACTA2",'MACROD2','KCNMA1',"PRUNE2"
  #"HMCN2","LAMA3","MCTP1","ZNF804A",
)

pdf(file='stroma/Immune_Cluster_1-19_Stromal_Marker_Dotplot_0424-newData.pdf',width=20,height=6)
DotPlot(KB, features = fib.markers, idents = paste0("S_", c(1:19))) + RotatedAxis()
dev.off()

#Vascular smooth muscle
vsm.markers <- c(
  "PDGFRB","ITGA8",                                       #Pan VSM markers
  "SAMD12","ROBO1","PIEZO2","GATA3",
  "DCC","POSTN","IL1RL1",
  "REN","ROBO2","KCNQ5","SLCO2A1","SPOCK3",
  "NOTCH3","MYH11",
  'RGS6','KCNMA1',"ADRA1A",
  "RGS5",'PLCB4',"ADGRB3","SLC38A11",
  "DNAH11","SLC6A1","TRPC6",
  "HS6ST3","ANO3","MPPED2","SCN3A",
  'PDE1C',"STEAP4","RXFP1",
  'FLNA', 'TAGLN',"ACTA2","VIM","ACTB",
  "CD36","LPL","ADIPOQ","FABP4",
  "PLIN2" #lipofibroblast-like cells
)

pdf(file='stroma/Immune_Cluster_20-27_Stromal_Marker_Dotplot_0424-newData.pdf',width=12,height=5)
DotPlot(KB, features = vsm.markers, idents = paste0("S_", c(20:27))) + RotatedAxis()
dev.off()

table(Idents(KB))

