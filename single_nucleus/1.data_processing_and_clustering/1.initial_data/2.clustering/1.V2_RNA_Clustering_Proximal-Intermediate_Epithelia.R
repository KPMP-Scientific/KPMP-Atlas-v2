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
setwd("~/data/net/home/Human_Kidney/Atlas_V2")

options(future.globals.maxSize = 1e9)
source("misc/utils.R")

options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")
#load("color_factors.robj")

KB.sub <- subset(KB, idents = c(1:100))
KB.sub <- subset(KB.sub, library %in% c("KM150", "KB54"), invert = TRUE)
KB.sub
#36588 features across 1240618 samples within 1 assay



##Subset to vasculature only
KB.PT <- subset(KB.sub, v1.subclass.l1 %in% c("POD","PEC","PT","DTL","ATL"))
KB.PT
#36588 features across 431408 samples


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.PT[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#6186 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3000, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k200infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3000]
VariableFeatures(KB.PT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.PT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.PT))
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                 min.dist = 0.2, return.model = T)
KB.PT@meta.data$pagoda_k200_infomap_pt <- k200infomap[rownames(KB.PT@meta.data)]

DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB.PT, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.PT, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200_round1.rda")

Idents(KB.PT)
KB.PT[["RNA"]] <- as(object = KB.PT[["RNA"]], Class = "Assay")
KB.PT <- ScaleData(KB.PT)

###Identify ambiguous low quality - annotate cluster identities
Idents(KB.PT) <- "pagoda_k200_infomap_pt"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:length(levels(Idents(KB.PT))))
levels(Idents(KB.PT)) <- paste("CL", levels(Idents(KB.PT)), sep = "")
celltype <- Idents(KB.PT)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v1clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$v1.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v1clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v1clusters)),
                                query = prop.table(table(seurat.obj.sub$v1.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$v1.subclass.l3)))
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

max.idents <- getMaxIdents(KB.PT, celltype)
max.idents

###Check Mean Genes
VlnPlot(KB.PT, features = c("nFeature_RNA"), pt.size = 0)

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
                 "PHACTR1","SLC14A2","SOX5","PCDH7",  #IMCD
                 "TP63","GPX2", "FXYD3", "KRT5",           #pape
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
                 "PTPRC",                                                        #ec
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


DotPlot(KB.PT, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.PT, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()


###Ambiguous clusters (checking markers for multiplets)
#CL2 - PT/TAL multiplet
#CL19 - PT/EC multiplet
#CL20 - PT/TAL multiplet
#CL29 - PT/IC multiplet
#CL32 - PT/FIB multiplet
#CL33 - PT/PC multiplet
#CL40 - MYOF/POD/MD/TAL multiplet
#CL42 - PT/IMM multiplet
#CL45 - POD/EC multiplet
#CL53 - ATL/TAL multiplet

                                             
###Subset multiplets
Idents(KB.PT) <- "pagoda_k200_infomap_pt"
to.remove <- c(2,19,20,29,32,33,40,42,45,53)
KB.PT <- subset(KB.PT, idents = to.remove, invert = TRUE)
KB.PT
#36588 features across 403165 samples




###Repeat subclustering
##Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.PT[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#6043 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3000, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately for each k value

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3000]
VariableFeatures(KB.PT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.PT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.PT))
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                 min.dist = 0.1, return.model = T)
KB.PT@meta.data$pagoda_k100_infomap_pt <- k100infomap[rownames(KB.PT@meta.data)]
KB.PT@meta.data$pagoda_k200_infomap_pt <- k200infomap[rownames(KB.PT@meta.data)]
KB.PT@meta.data$pagoda_k500_infomap_pt <- k500infomap[rownames(KB.PT@meta.data)]

DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB.PT, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.PT, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k100.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k500.rda")

save(KB.PT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")


###Annotation of k200 clusters
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
Idents(KB.PT) <- "pagoda_k200_infomap_pt"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:length(levels(Idents(KB.PT))))
levels(Idents(KB.PT)) <- paste("CL", levels(Idents(KB.PT)), sep = "")
celltype <- Idents(KB.PT)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v1clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$v1.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v1clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v1clusters)),
                                query = prop.table(table(seurat.obj.sub$v1.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$v1.subclass.l3)))
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

max.idents <- getMaxIdents(KB.PT, celltype)
max.idents

##Cluster Assessment:
##Check batch, any clusters associated with only 1-2 individuals
##Check marker genes, sufficient markers to be distinct sub-population? merge at K500?
##check for splitting at k100 - splits found primarily for PT clusters, will not split these further
##check for overlapping cell type markers - mark as multiplets if cells show overlapping cell type marker profiles
##check for injury features (high %ERT, high %MT, low genes, injury markers)
##Mark clusters as ambiguous if they show no distinct DEGs, have very low genes detected and very low marker expression

#CL 1: PT-S1 (atlas v1 cl 5 - weak) (no distinct DEGs) (low genes) - merge to cl 16
#CL 2: aPT (atlas v1 cl 11 - weak) 
#CL 3: aPT (atlas v1 cl 9) (1 distinct DEG) - no obvious merging pattern
#CL 4: aPT (atlas v1 cl 11) 
#CL 5: aPT (atlas v1 cl 11)
#CL 6: PT-S1 (atlas v1 cl 5) (no distinct DEGs) - merge to cl 16
#CL 7: PT-S2 (atlas v1 cl 6) (no distinct DEGs)
#CL 8: cycPT (no atlas v1) (no distinct DEGs) - split? weak cyc markers 
#CL 9: dPT (atlas v1 cl 18) (no distinct DEGs) (merges with cl 19 at k500) (high %MT, high %ERT, low genes, injury markers) - merge with cl19
#CL 10: PT-S3 (atlas v1 cl 7)
#CL 11: aPT (atlas v1 cl 11)
#CL 12: dPT (atlas v1 cl 17) (no distinct DEGs) (heavy contribution from pt 31-10035) (very low genes, indistinguishable markers) - ambiguous
#CL 13: PT-S1 (atlas v1 cl 5 - weak) (no distinct DEGs) (merges with cl 1 at k500) (low genes) - merge to cl 16
#CL 14: dPT (atlas v1 cl 17 - weak) (no distinct DEGs) (low genes, indistinguishable markers) - ambiguous
#CL 15: dPT (atlas v1 cl 15 - weak) (no distinct DEGs) (injury markers)
#CL 16: PT-S1 (atlas v1 cl 5)
#CL 17: aPT (atlas v1 cl 8)
#CL 18: POD (atlas v1 cl 1)
#CL 19: dPT/DTL (atlas v1 cl 19) (no distinct DEGs) (merges with cl 9 at k500) (high %MT, high %ERT, low genes)
#CL 20: PEC (atlas v1 cl 3)
#CL 21: PT-S1 (atlas v1 cl 5 - weak) (no distinct DEGs) - merge to cl 16
#CL 22: dPT (atlas v1 cl 17 - weak) (1 distinct DEG) (low genes, indistinguishable markers) - ambiguous
#CL 23: ?dPT? (no atlas v1) - erythrocyte markers - multiplet
#CL 24: PT-S1 (atlas v1 cl 5 - weak) (1 distinct DEG) (merges with cl 16 at k500) - merge to cl 16
#CL 25: dPT-S3 (atlas v1 cl 7 - weak) (high %MT, high %ERT, injury markers)
#CL 26: dPOD (atlas v1 cl 2) (high %ERT, injury markers)
#CL 27: dPT (no atlas v1) (low genes)
#CL 28: cycPT (no atlas v1) - weak cyc markers
#CL 29: aPT (no atlas v1) (unstable at k500)
#CL 30: aPT (atlas v1 cl 11 - weak)
#CL 31: PT-S1 (atlas v1 cl 5 - weak) (1 distinct DEG) - merge to cl 16
#CL 32: DTL2 (atlas v1 cl 20)
#CL 33: DTL1 (atlas v1 cl 21)
#CL 34: ATL (atlas v1 cl 24)
#CL 35: PT-S2 (atlas v1 cl 6) (merges with cl 21 at k500) - merge to cl 7
#CL 36: dDTL1 (atlas v1 cl 21 - weak) (high %ERT, injury markers)
#CL 37: aDTL1 (atlas v1 cl 21 - weak) (1 distinct DEG) (merges with cl 38 at k500) (high %ERT, injury markers)
#CL 38: aPT (atlas v1 cl 10) (high %ERT, injury markers)
#CL 39: aPT (atlas v1 cl 11)
#CL 40: PT-S3 (atlas v1 cl 7 - weak) (merges with cl 10 at k500) - merge to cl 10
#CL 41: dPT-S2 (no atlas v1) (unstable at k500) (1 distinct DEG)(low genes)
#CL 42: dDTL1 (atlas v1 cl 21 - weak) (high %ERT)(1 distinct DEG)
#CL 43: dPT (no atlas v1) (unstable at k500) (no distinct DEGs)(low genes, indistinguishable markers) - ambiguous
#CL 44: aPT (no atlas v1) (unstable at k500)
#CL 45: aPT (no atlas v1) (unstable at k500)
#CL 46: DTL3 (atlas v1 cl 22 - weak) 
#CL 47: dATL (atlas v1 cl 25) (merges with cl 36 at k500) (high %ERT, injury markers)
#CL 48: PT-S2 (atlas v1 cl 6) (unstable at k500) PT/PL multiplet
#CL 49: aPT (atlas v1 cl 11) - only associated with one reference sample - need to remove SDKZ0001 - ambiguous
#CL 50: aPT (atlas v1 cl 11) (unstable at k500) - only associated with one reference sample, need to remove SDKZ0001 - ambiguous



#Check for marker genes K200
pt.markers <- FindAllMarkers(KB.PT, only.pos = TRUE, max.cells.per.ident = 2000,
                             logfc.threshold = 0.25, min.pct = 0.25)
write.table(pt.markers, file="Intermediate_Objects/Kidney_v2_proximal-tubules_k200_Cluster_markers_06282023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(pt.markers, file="Intermediate_Objects/Kidney_v2_proximal-tubules_k200_Cluster_markers_06282023.robj")
#load("Intermediate_Objects/Kidney_v2_proximal-tubules_k200_Cluster_markers_06282023.robj")

cl.mark <- pt.markers[pt.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="Intermediate_Objects/Kidney_v2_proximal-tubules_k200_Cluster_markers_06282023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.PT, features = top5$gene) + RotatedAxis()




###Check clusters for patient distribution
Idents(KB.PT) <- "pagoda_k200_infomap_pt"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:length(levels(Idents(KB.PT))))
load("color_factors.robj")

row.order <- levels(Idents(KB.PT))
prop1 <- prop.table(table(KB.PT$patient, KB.PT$pagoda_k200_infomap_pt), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop1 <- prop.table(table(KB.PT$region_level2, KB.PT$pagoda_k200_infomap_pt), margin = 2)[,row.order]
barplot(prop1,main = "Region Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(region.l2.cols[rownames(prop1)]))
prop1 <- prop.table(table(KB.PT$condition_level1, KB.PT$pagoda_k200_infomap_pt), margin = 2)[,row.order]
barplot(prop1,main = "Condition Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(cond.l1.cols[rownames(prop1)]))

VlnPlot(object = KB.PT, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)

table(Idents(KB.PT))


table(KB.PT$pagoda_k200_infomap_pt,KB.PT$patient)[c(12,39,49,50),]
#patient 31-10035 is showing 19,973 cells in cluster 12!?
table(KB.PT$pagoda_k200_infomap_pt,KB.PT$patient)[,"31-10035"]


table(KB.PT$pagoda_k200_infomap_ec,KB.PT$condition_level1)


###compare cluster resolutions
Idents(KB.PT) <- "pagoda_k200_infomap_pt"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:length(levels(Idents(KB.PT))))

library("scrattch.hicat")

#Pagoda2 K200 and K500
final<-Idents(KB.PT)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.PT$pagoda_k500_infomap_pt
names(prop)<-rownames(KB.PT@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k500_proximal-tubules.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k500", labels=compare.result$cl.id.map$old)
dev.off()


#Pagoda2 K200 and K100
final<-Idents(KB.PT)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.PT$pagoda_k100_infomap_pt
names(prop)<-rownames(KB.PT@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k100_proximal-tubules.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k100", labels=compare.result$cl.id.map$old)
dev.off()


###Check cell type markers
Idents(KB.PT) <- "pagoda_k200_infomap_pt"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:length(levels(Idents(KB.PT))))
DotPlot(KB.PT, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.PT, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

markers <- c(
"NPHS1","NPHS2","PODXL",                              #POD
"CDKN1C","SPOCK2",                                    #dPOD
"CLDN1", "CFH","ALDH1A2",                             #PEC
  

"SLC5A12","SLC22A6",                                  #S1/S2
"PRODH2","SLC5A2","SLC22A8","SLC7A8",                 #S1
"SLC34A1","SLC5A10",                                  #S2                                  
"SLC5A11","SLC7A13",                                  #S3

"ITGB8","CDH6","HAVCR1","VCAM1","PROM1",              #aEpi

"CRYAB","TACSTD2","SLC44A5",                          #TL
"AQP1", "UNC5D",                                      #DTL1
"ADGRL3","ID1",                                       #DTL2
"SMOC2",                                              #DTL3
"AKR1B1","SH3GL3",                                    #DTL3/ATL
"PROX1",                                              #ATL

"EMCN","PECAM1","PDGFRB","ACTA2","C7","DCN","COL1A1", #EC/stromal
"PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
"CD96","NKG7",       
"MS4A2","CPA3", "KIT","S100A9",
"APOE","DEFB1","CST3","GATM","ALDOB",                 #Injury
"IGFBP7","CLU","SPP1",
"UMOD","SLC12A1","SLC12A3","SLC8A1","AQP2",           #Distal tubules
"HSD11B2","SLC26A7",
"MKI67","TOP2A"
)

DotPlot(KB.PT, features = unique(markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()




### Remove ambiguous clusters
#remove cl  - multiplets/ambiguous: 12,14,22,23,43,48,49,50
Idents(KB.PT) <- "pagoda_k200_infomap_pt"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:length(levels(Idents(KB.PT))))
KB.PT <- subset(KB.PT, pagoda_k200_infomap_pt %in% c(12,14,22,23,43,48,49,50), invert = TRUE)
KB.PT <- RunUMAP(object = KB.PT, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                 min.dist = 0.1, return.model = T)
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_pt", repel = TRUE) + NoLegend()

save(KB.PT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")



### Re-order and Merge clusters
#Strategy: merge clusters with no distinct markers (<2) to the largest cluster of that subclass, except if:
#clusters show features of altered states (injury markers, high %ERT/MT, cycling or repair markers)
#and only merge when there is an obvious merging pattern. Merging is also guided by merging found at k500.


Idents(KB.PT) <- "pagoda_k200_infomap_pt"
cls <- c(18,26,20,1,6,13,31,24,21,16,15,9,19,7,35,41,10,40,25,3,39,45,17,30,2,
         11,44,38,29,4,5,8,28,27,32,33,37,36,42,46,34,47)
new.cls <- c(1,2,3,4,4,4,4,4,4,4,5,6,6,7,7,8,9,9,10:33)
Idents(KB.PT) <- plyr::mapvalues(Idents(KB.PT), from = cls, to = new.cls)
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:33)
KB.PT$pt.v2.clusters <- Idents(KB.PT)

DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pt.v2.clusters", repel = TRUE) + NoLegend()

save(KB.PT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
KB.PT
#36588 features across 360023 samples


###Re-Order
Idents(KB.PT) <- "pt.v2.clusters" 

cls <- c(1:5,13,6:10,
         14,12,17,19,11,21,22,
         15,16,18,20,
         23:26,28,27,31,29,30,32,33)
new.cls <- c(1:33)
Idents(KB.PT) <- plyr::mapvalues(Idents(KB.PT), from = cls, to = new.cls)
Idents(KB.PT) <- factor(Idents(KB.PT), levels = 1:33)
KB.PT$pt.v2.clusters <- Idents(KB.PT)

DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pt.v2.clusters", repel = TRUE) + NoLegend()

save(KB.PT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")






###Re-identify cluster markers
pt.markers <- FindAllMarkers(KB.PT, only.pos = TRUE, max.cells.per.ident = 2000,
                             logfc.threshold = 0.25, min.pct = 0.25)
write.table(pt.markers, file="proximal/Kidney_v2_proximal-tubules_Cluster_markers_06302023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(pt.markers, file="proximal/Kidney_v2_proximal-tubules_Cluster_markers_06302023.robj")
#load("Intermediate_Objects/Kidney_v2_proximal-tubules_Cluster_markers_06302023.robj")

cl.mark <- pt.markers[pt.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="proximal/Kidney_v2_proximal-tubules_Cluster_markers_06302023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.PT, features = top5$gene) + RotatedAxis()



###Run NSForest for markers
#Save anndata
expression_matrix <- KB.PT[["RNA"]]$counts
umap <- Embeddings(KB.PT, reduction = "umap")
dfobs <- KB.PT@meta.data
dfobs$pt.v2.clusters <- paste0("CL",dfobs$pt.v2.clusters)
dfvar <- KB.PT@assays$RNA@meta.features
adata <- AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  obsm=list('X_umap'=umap)
)
write_h5ad(adata, "Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.h5ad")



###Re-organize and assign markers
#PT-TL markers
pt.markers <- c(
  "NPHS1","NPHS2","ST6GALNAC3","PODXL",                 #POD
  "PTGDS","BST2","TPPP3","CHI3L1",                      #dPOD
  "ALDH1A2","CFH","FAM155A","CLDN1",                    #PEC
  
  "LRP2","CUBN",                                        #PT
  "SLC5A12","SLC22A6",                                  #S1/S2
  "PRODH2","SLC22A8","PCDH15","RALYL",                  #S1
  "FGA","FGB",                                          #13
  
  "SLC34A1","ANKS1B","SLC5A10","SLC5A11",               #S2                                  
  "SLC22A24","SLC7A13","GPM6A","SLC5A8",                #S3
  
  "ITGB8","CDH6","HAVCR1",                              #aEpi
  "KCTD16","SPON1",                                     #14
  "BTBD11","SYTL3","ZNF250","ITGB3","DLGAP1",           #12
  "PROM1",
  "HDAC9","RAPGEF5","SLC35F3",                          #17
  "SOX4","VIM","IL32","MMP7",                           #19
  "VCAM1",
  "MEG3",                                               #11
  "APBB1IP","LSAMP","KCNIP1","NRXN3",                   #21
  "NRG3","FAM189A1","KITLG","GRM8",                     #22
  
  "CCL2","BIRC3","ICAM1","FSTL3","RELB",                #15
  "DCC","SLC9A9","FAM135B","COL23A1",                   #16
  "PTCHD4",                                             #18
  "ZFPM2","GREB1L","PPFIA2","PDE8B",                    #20
  "TOP2A","MKI67",                                      #cycling
  "DTL","BRIP1","CENPK","EZH2",
  "EGR1","FOS",                                         #dPT
  
  "PAX8-AS1","SLC44A5","CRYAB","TACSTD2",               #TL
  
  "AQP1", "UNC5D","LRRC4C","SLCO1A2","DENND2A",         #DTL2
  
  "SIM2","ADGRL3","JAG1","SMAD9","ID1",                 #DTL1
  "SLC14A2","FAM169A","SMOC2",                          #DTL3
  "ABCA4","BCL6","AKR1B1","SH3GL3",                     #DTL3/ATL
  "CLCNKA","PROX1","CACNA1C",                           #ATL
  "SOD3","HSPA2","CSRP2",
  
  "DEFB1","CST3","APOE","GATM","ALDOB",                 #Injury
  "SPP1","IGFBP7","CLU"
)

DotPlot(KB.PT, features = pt.markers) + RotatedAxis()



###Update Experiment metadata
###Add in experiment metadata
meta <- KB.PT@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race","location","laterality","protocol",
         "tissue_type_full","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB.PT@meta.data <- meta
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "condition_level1", repel = TRUE) 
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 
DimPlot(KB.PT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 

save(KB.PT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")




###Barplots on stats
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
#Generate color factors
load("color_factors.robj")
meta <- KB.PT@meta.data
meta <- meta[!duplicated(meta$pt.v2.clusters),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "pt.v2.clusters", sort_label = F, colorset = "varibow")
pt.cl.cols <- meta$pt.v2.clusters_color; names(pt.cl.cols) <- meta$pt.v2.clusters_label

meta <- KB.PT@meta.data
meta <- meta[,c("assay","condition_level3","condition_level2","condition_level1")]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "assay", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level3", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level2", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level1", sort_label = F, colorset = "varibow")

assay.cols <- meta$assay_color; names(assay.cols) <- meta$assay_label; assay.cols <- factor(na.omit(assay.cols))
assay.cols <- factor(assay.cols[!duplicated(assay.cols)])
condl3.cols <- meta$condition_level3_color; names(condl3.cols) <- meta$condition_level3_label; condl3.cols <- factor(na.omit(condl3.cols))
condl3.cols <- factor(condl3.cols[!duplicated(condl3.cols)])
condl2.cols <- meta$condition_level2_color; names(condl2.cols) <- meta$condition_level2_label; condl2.cols <- factor(na.omit(condl2.cols))
condl2.cols <- factor(condl2.cols[!duplicated(condl2.cols)])
condl1.cols <- meta$condition_level1_color; names(condl1.cols) <- meta$condition_level1_label; condl1.cols <- factor(na.omit(condl1.cols))
condl1.cols <- factor(condl1.cols[!duplicated(condl1.cols)])

save(assay.cols,condl3.cols,condl2.cols,condl1.cols,pt.cl.cols,
     file = "proximal/color_factors_pt.robj")


#Condition (level 3)
Idents(object = KB.PT) <- "condition_level3"
pdf(file='proximal/Kidney_PT-TL_Condition_level3_umap.pdf',width=9,height=6)
DimPlot(KB.PT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 3"
        ) + scale_color_manual(values = alpha(condl3.cols[levels(Idents(KB.PT))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Condition (level 1)
Idents(object = KB.PT) <- "condition_level1"
pdf(file='proximal/Kidney_PT-TL_Condition_level1_umap.pdf',width=9,height=6)
DimPlot(KB.PT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.PT))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='proximal/Kidney_PT-TL_Condition_level1_umap_split.pdf',width=24,height=12)
DimPlot(KB.PT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1, split.by = "condition_level1", ncol = 3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.PT))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()




#Assay
Idents(object = KB.PT) <- "assay"
pdf(file='proximal/Kidney_PT-TL_Assay_umap.pdf',width=9,height=6)
DimPlot(KB.PT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Assay"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(KB.PT))], 0.3), name = "Assay"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Clusters
Idents(object = KB.PT) <- "pt.v2.clusters"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = c(1:33))
pdf(file='proximal/Kidney_PT-TL_Clusters_umap.pdf',width=8,height=6)
DimPlot(KB.PT, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(pt.cl.cols[levels(Idents(KB.PT))], 0.3), name = "Clusters"
        ) + NoLegend()
dev.off()


pdf(file='proximal/Cluster_Contributions_Barplot.pdf',width=6,height=24)
layout(matrix(c(1,1:7), nrow = 8, ncol = 1, byrow = TRUE))
col.order <- levels(Idents(KB.PT))
prop1 <- prop.table(table(KB.PT$patient,KB.PT$pt.v2.clusters), margin = 2)[,col.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop2 <- prop.table(table(KB.PT$assay,KB.PT$pt.v2.clusters), margin = 2)[,col.order]
barplot(prop2,main = "Assay Proportions", cex.names = 0.5, 
        names.arg = colnames(prop2),las=2, col = as.character(assay.cols[rownames(prop2)]))
#red = snRNA-seq
prop3 <- prop.table(table(KB.PT$condition_level3,KB.PT$pt.v2.clusters), margin = 2)[,col.order]
barplot(prop3,main = "Condition l3 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(condl3.cols[rownames(prop3)]))

prop4 <- prop.table(table(KB.PT$condition_level1,KB.PT$pt.v2.clusters), margin = 2)[,col.order]
barplot(prop4,main = "Condition l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop4),las=2, col = as.character(condl1.cols[rownames(prop4)]))
prop5 <- prop.table(table(KB.PT$region_level2,KB.PT$pt.v2.clusters), margin = 2)[,col.order]
barplot(prop5,main = "Region l2 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop5),las=2, col = as.character(region.l2.cols[rownames(prop5)]))
batch.entropy<-table(KB.PT$patient, KB.PT$pt.v2.clusters)
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Experiment Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(KB.PT$pt.v2.clusters)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)
dev.off()






###Comparison with published reference data sets
library(corrplot)
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
KB.PT <- ScaleData(KB.PT, features = rownames(KB.PT))

##Li et al 2022
#Humphreys' data (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139107)
load("~/data/pod/blake_LTS/public/Li2022/Li_Humphreys_2022_mouse_seurat.rda")

hki
#29935 features across 202058 samples

table(Idents(hki))
VariableFeatures(hki)

#Intersect variable genes
common.human <- VariableFeatures(KB.PT)[which(VariableFeatures(KB.PT)[1:2000] %in% rownames(hki))]

ave.KB.PT<-AverageExpression(KB.PT, features = common.human, slot = "scale.data")
ave.hki<-AverageExpression(hki, features = common.human, slot = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.PT$RNA),as.data.frame(ave.hki$RNA)))
ave.cor<-ave.cor[1:33,c("Pod","PT","PT-AcInj","PT-Inj","PT-R","PT-FR","DTL-ATL")]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='proximal/PT-IT_Clusters_Li2022_Mouse.pdf',width=4,height=12)
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()



