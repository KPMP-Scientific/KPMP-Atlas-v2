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
KB.DT <- subset(KB.sub, v1.subclass.l1 %in% c("TAL","DCT","CNT","PC","IC","PapE"))
KB.DT
#36588 features across 543193 samples


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.DT[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#6162 overdispersed genes

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
VariableFeatures(KB.DT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.DT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.DT))
KB.DT <- RunUMAP(object = KB.DT, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                 min.dist = 0.2, return.model = T)
KB.DT@meta.data$pagoda_k200_infomap_dt <- k200infomap[rownames(KB.DT@meta.data)]

DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.DT, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200_round1.rda")
#load("Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200_round1.rda")

Idents(KB.DT)
KB.DT[["RNA"]] <- as(object = KB.DT[["RNA"]], Class = "Assay")
KB.DT <- ScaleData(KB.DT)

###Identify ambiguous low quality - annotate cluster identities
Idents(KB.DT) <- "pagoda_k200_infomap_dt"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:length(levels(Idents(KB.DT))))
levels(Idents(KB.DT)) <- paste("CL", levels(Idents(KB.DT)), sep = "")
celltype <- Idents(KB.DT)

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

max.idents <- getMaxIdents(KB.DT, celltype)
max.idents

###Check Mean Genes
VlnPlot(KB.DT, features = c("nFeature_RNA"), pt.size = 0) + NoLegend()

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


DotPlot(KB.DT, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.DT, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

###Ambiguous clusters (checking markers for multiplets)
#CL2 - EC multiplet
#CL19 - PT multiplet
#CL20 - PT/IC multiplet
#CL25 - PT/PC multiplet
#CL39 - MAC multiplet
#CL42 - FIB multiplet
#CL44 - TL/TAL multiplet
#CL47 - PT/IC multiplet
#CL52 - PT/DCT multiplet
#CL63 - EC/FIB multiplet
#CL67-72 - only one cell each

                                             
###Subset multiplets
Idents(KB.DT) <- "pagoda_k200_infomap_dt"
to.remove <- c(2,19,20,25,39,42,44,47,52,63,67:72)
KB.DT <- subset(KB.DT, idents = to.remove, invert = TRUE)
KB.DT
#36588 features across 516591 samples




###Repeat subclustering
##Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.DT[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#6056 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3000, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 800, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately for each k value

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap
k800infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3000]
VariableFeatures(KB.DT) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.DT[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.DT))
KB.DT <- RunUMAP(object = KB.DT, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                 min.dist = 0.1, return.model = T)
KB.DT@meta.data$pagoda_k100_infomap_dt <- k100infomap[rownames(KB.DT@meta.data)]
KB.DT@meta.data$pagoda_k200_infomap_dt <- k200infomap[rownames(KB.DT@meta.data)]
KB.DT@meta.data$pagoda_k500_infomap_dt <- k500infomap[rownames(KB.DT@meta.data)]
KB.DT@meta.data$pagoda_k800_infomap_dt <- k800infomap[rownames(KB.DT@meta.data)]

DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k800_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.DT, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k100.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k500.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k800.rda")

save(KB.DT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")




###Update Experiment metadata
###Add in experiment metadata
meta <- KB.DT@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race","location","laterality","protocol",
         "tissue_type_full","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB.DT@meta.data <- meta
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "condition_level1", repel = TRUE) 
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 

save(KB.DT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")





###Annotation of k200 clusters
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
Idents(KB.DT) <- "pagoda_k200_infomap_dt"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:length(levels(Idents(KB.DT))))
levels(Idents(KB.DT)) <- paste("CL", levels(Idents(KB.DT)), sep = "")
celltype <- Idents(KB.DT)

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

max.idents <- getMaxIdents(KB.DT, celltype)
max.idents

##Cluster Assessment:
##Check batch, any clusters associated with only 1-2 individuals
##Check marker genes, sufficient markers to be distinct sub-population? merge at K500?
##check for splitting at k100 - splits found primarily for PT clusters, will not split these further
##check for overlapping cell type markers - mark as multiplets if cells show overlapping cell type marker profiles
##check for injury features (high %ERT, high %MT, low genes, injury markers)
##Mark clusters as ambiguous if they show no distinct DEGs, have very low genes detected and very low marker expression


#1 CL 50: M-TAL (atlas v1 cl 31) (merges to cl2 at k100)
#1 CL 58: M-TAL (atlas v1 cl 31) (unstable at k500) - merge to cl 50
#2 CL 61: dM-TAL (atlas v1 cl 31 - weak) (no distinct DEGs) (unstable at k500) (injury markers)
#3 CL 15: dM-TAL (atlas v1 cl 32 - weak) (no distinct DEGs) (high %MT, high %ERT, injury markers)
#4 CL 14: dM-TAL (atlas v1 cl 31 - weak) (1 distinct DEG)

#5 CL 2: C/M-TAL-A (atlas v1 cl 31 - weak) (no distinct DEGs)
#5 CL 47: C/M-TAL-A (atlas v1 cl 31 - weak) (unstable at k500, merges to cl2 at k100) - merge to cl2 based on region
#5 CL 57: C/M-TAL-A (atlas v1 cl 31 - weak) (unstable at k500) - merge to cl2 based on marker genes
#6 CL 34: C-TAL-A (atlas v1 cl 37 - weak) 
#7 CL 32: C/M-TAL-A (atlas v1 cl 31 - weak) (no distinct DEGs) - no clear merging
#8 CL 16: C-TAL-A (atlas v1 cl 36)
#8 CL 48: C-TAL-A (atlas v1 cl 37 - weak) (1 distinct DEG) (unstable at k500) - merge to cl16 based on marker expression

#9 CL 19: M-TAL-B (atlas v1 cl 31 - weak) (no distinct DEGs) (merges with cl11 at k500)
#10 CL 11: C-TAL-B (atlas v1 cl 37)
#11 CL 18: C-TAL-B (atlas v1 cl 37)

#12 CL 17: MD (atlas v1 cl 39)
#13 CL 6: MD (atlas v1 cl 39)

#14 CL 12: dC-TAL (atlas v1 cl 36 - weak) (no distinct DEGs) (low genes) 
#15 CL 3: aTAL1 (C/M) (atlas v1 cl 26) 
#16 CL 24: aTAL1 (C/M) (atlas v1 cl 30 - weak) 
#17 CL 30: aTAL2 (C/M) (atlas v1 cl 28 - weak)

#18 CL 21: cycTAL (no atlas v1) (no distinct DEGs) (split at k100?)

#19 CL 1: DCT1 (atlas v1 cl 40)
#20 CL 23: DCT1 (atlas v1 cl 40) 
#19 CL 29: DCT1 (atlas v1 cl 40) (no distinct DEGs) (low genes) (merges with cl1 at k500) merge to cl 1
#21 CL 8: DCT2 (atlas v1 cl 41) 
#22 CL 36: dDCT (atlas v1 cl 42 - weak) (no distinct DEGs) (high %MT, high %ERT, injury markers)
#23 CL 25: aDCT (atlas v1 cl 40) (no distinct DEGs) (merges with cl1 at k500) adaptive markers
#24 CL 42: aDCT (atlas v1 cl 40) (1 distinct DEG) adaptive markers

#25 CL 7: CNT (atlas v1 cl 44 - weak) 
#25 CL 26: CNT (atlas v1 cl 44 - weak) (no distinct DEGs) (merges with cl7 at k800) - merge to cl7
#26 CL 10: CNT (atlas v1 cl 46 - weak)
#27 CL 35: dCNT (atlas v1 cl 46) (no distinct DEGs) (high %MT, high %ERT, injury markers)
#28 CL 20: aCNT (atlas v1 cl 44) - adaptive markers
#29 CL 31: CNT-PC (atlas v1 cl 45 - weak)
#30 CL 27: dCNT-PC (atlas v1 cl 48 - weak) (no distinct DEGs) (low genes) no clear merging pattern

#31 CL 37: CCD-PC (atlas v1 cl 48) 
#32 CL 9: CCD-PC (atlas v1 cl 48) 
#33 CL 44: OMCD-PC (atlas v1 cl 50 - weak)
#34 CL 33: dOMCD-PC (atlas v1 cl 50) (1 distinct DEG) (high %MT, high %ERT, injury markers)
#35 CL 40: dOMCD-PC (atlas v1 cl 50 - weak) (no distinct DEGs) (no clear merging pattern) (merges with cl 14 at k500, similar markers?)
#36 CL 54: dOMCD-PC (atlas v1 cl 52 - weak) (no distinct DEGs) (merges with cl26 at k500)

#38 CL 41: IMCD (atlas v1 cl 52) 
#39 CL 56: IMCD (atlas v1 cl 52) 
#40 CL 55: dIMCD (atlas v1 cl 53) (high %ERT, injury markers)

#41 CL 46: PapE (atlas v1 cl 61) (high %ERT) (injury markers)

#42 CL 5: CCD-IC-A (atlas v1 cl 55)
#43 CL 13: CCD-IC-A (atlas v1 cl 55) 
#44 CL 39: dCCD-IC-A (atlas v1 cl 55) (injury markers)
#45 CL 4: dCCD-IC-A (atlas v1 cl 57) (high %MT, high %ERT, injury markers)
#46 CL 43: dCCD-IC-A (atlas v1 cl 57) (1 distinct DEG) (high %MT, high %ERT, injury markers) (merges with cl53 at k500)

#47 CL 38: OMCD-IC-A (atlas v1 cl 58)
#48 CL 53: dOMCD-IC-A (atlas v1 cl 57) (high %MT, high %ERT, injury markers)
#49 CL 52: tOMCD-PC-IC (atlas v1 cl 59) (no distinct DEGs)

#51 CL 22: IC-B (atlas v1 cl 60) 


#CL 28: dC-TAL (no atlas v1) (no distinct DEGs) (high %MT, high %ERT) PT multiplet
#CL 45: C-TAL (C/M) (no atlas v1) (1 distinct DEG) - erythrocyte markers - mutliplet
#CL 49: aTAL1 (atlas v1 cl 30) - contributed mostly by SDKZ0001 - ambiguous low quality
#CL 51: CNT (atlas v1 cl 44) (merges with cl10 at k500) (mostly contributed by 2 ref samples) - ambiguous
#37 CL 59: OMCD-PC (atlas v1 cl 49) (unstable at k500) TL multiplet
#50 CL 60: OMCD-IC-A (atlas v1 cl 58) (unstable at k500) TL multiplet


#Check for marker genes K200
dt.markers <- FindAllMarkers(KB.DT, only.pos = TRUE, max.cells.per.ident = 2000,
                             logfc.threshold = 0.25, min.pct = 0.25)
write.table(dt.markers, file="Intermediate_Objects/Kidney_v2_distal-tubules_k200_Cluster_markers_06282023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(dt.markers, file="Intermediate_Objects/Kidney_v2_distal-tubules_k200_Cluster_markers_06282023.robj")
#load("Intermediate_Objects/Kidney_v2_distal-tubules_k200_Cluster_markers_06282023.robj")

cl.mark <- dt.markers[dt.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="Intermediate_Objects/Kidney_v2_distal-tubules_k200_Cluster_markers_06282023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.DT, features = top5$gene) + RotatedAxis()




###Check clusters for patient distribution
Idents(KB.DT) <- "pagoda_k200_infomap_dt"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:length(levels(Idents(KB.DT))))
load("color_factors.robj")

row.order <- c(50,61,58,14,15,2,47,57,32,34,16,48,12,3,24,30,49,19,11,18,17,6,
               21,1,23,29,25,8,36,42,7,26,10,51,20,31,27,35,9,37,44,33,40,54,59,
               56,41,55,5,13,39,4,43,38,53,60,52,22,46,28,45)

prop1 <- prop.table(table(KB.DT$patient, KB.DT$pagoda_k200_infomap_dt), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop1 <- prop.table(table(KB.DT$region_level2, KB.DT$pagoda_k200_infomap_dt), margin = 2)[,row.order]
barplot(prop1,main = "Region Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(region.l2.cols[rownames(prop1)]))
prop1 <- prop.table(table(KB.DT$condition_level1, KB.DT$pagoda_k200_infomap_dt), margin = 2)[,row.order]
barplot(prop1,main = "Condition Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(cond.l1.cols[rownames(prop1)]))

VlnPlot(object = KB.DT, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)

table(Idents(KB.DT))


table(KB.DT$pagoda_k200_infomap_dt,KB.DT$patient)[c(49,51),]




###compare cluster resolutions
Idents(KB.DT) <- "pagoda_k200_infomap_dt"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:length(levels(Idents(KB.DT))))

library("scrattch.hicat")

#Pagoda2 K200 and K500
final<-Idents(KB.DT)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.DT$pagoda_k500_infomap_dt
names(prop)<-rownames(KB.DT@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k500_distal-tubules.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k500", labels=compare.result$cl.id.map$old)
dev.off()

#Pagoda2 K200 and K800
final<-Idents(KB.DT)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.DT$pagoda_k800_infomap_dt
names(prop)<-rownames(KB.DT@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k800_distal-tubules.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k800", labels=compare.result$cl.id.map$old)
dev.off()



#Pagoda2 K200 and K100
final<-Idents(KB.DT)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.DT$pagoda_k100_infomap_dt
names(prop)<-rownames(KB.DT@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k100_distal-tubules.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k100", labels=compare.result$cl.id.map$old)
dev.off()


###Check cell type markers
Idents(KB.DT) <- "pagoda_k200_infomap_dt"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = row.order)
DotPlot(KB.DT, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.DT, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

markers <- c(
  "CASR","SLC12A1","UMOD","EGF",                       #TAL
  "CLDN14", "KCTD16","ANK2",                           #M-TAL
  "ENOX1","TMEM207","CLDN16",                          #C-TAL
  "NOS1","ROBO2",                                      #MD
  
  "SLC12A3","TRPM6",                                   #DCT
  "ADAMTS17",                                          #DCT1
  "SLC8A1","SCN2A","HSD11B2","CALB1",                  #DCT2 / CNT
  "PCDH7",                                             #CNT
  "SCNN1G","SCNN1B",                                   #CNT                   
  "GATA3","AQP2","AQP3",                               #PC
  "SLC4A1", "SLC26A7",                                 #IC
  "SLC4A9", "SLC26A4",                                 #IC-B
  "SLC14A2","FXYD4",                                   #IMCD
  "TP63", "KRT5",                                      #PapE
  
  
  "ITGB8","CDH6","HAVCR1","VCAM1","PROM1",              #aEpi
  
   
  "EMCN","PECAM1","PDGFRB","ACTA2","C7","DCN","COL1A1", #EC/stromal
  "PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
  "CD96","NKG7",      
  "MS4A2","CPA3", "KIT","S100A9",
  "APOE","DEFB1","CST3","GATM","ALDOB",                 #Injury
  "IGFBP7","CLU","SPP1",
  "NPHS1","CLDN1","LRP2","SLC5A12","CRYAB","TACSTD2",          #proximal tubules
  "AKR1B1","SH3GL3",
  "MKI67","TOP2A"
)

DotPlot(KB.DT, features = unique(markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4,
        idents = c(50,61,58,14,15,2,47,57,32,48,16,12,3,24,30,49,19,11,34,18,17,6,
                   21,1,23,29,25,8,36,42,7,26,10,51,20,31,27,35,9,37,44,33,40,54,59,
                   56,41,55,5,13,39,4,43,38,53,60,52,22,46)
) + RotatedAxis()



### Remove ambiguous clusters
#remove cl  - multiplets/ambiguous: 28,45,49,51
Idents(KB.DT) <- "pagoda_k200_infomap_dt"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:length(levels(Idents(KB.DT))))
KB.DT <- subset(KB.DT, pagoda_k200_infomap_dt %in% c(28,45,49,51), invert = TRUE)
KB.DT <- RunUMAP(object = KB.DT, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                 min.dist = 0.1, return.model = T)
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_dt", repel = TRUE) + NoLegend()

save(KB.DT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")



### Re-order and Merge clusters
#Strategy: merge clusters with no distinct markers (<2) to the largest cluster of that subclass, except if:
#clusters show features of altered states (injury markers, high %ERT/MT, cycling or repair markers)
#and only merge when there is an obvious merging pattern. Merging is also guided by merging found at k500/k800.


Idents(KB.DT) <- "pagoda_k200_infomap_dt"
cls <- c(50,58,61,15,14,2,47,57,34,32,16,48,19,11,18,17,6,12,3,24,30,21,1,23,29,
         8,36,25,42,7,26,10,35,20,31,27,37,9,44,33,40,54,59,41,56,55,46,5,13,39,
         4,43,38,53,52,60,22)
new.cls <- c(1,1,2:5,5,5,6:8,8:20,19,21:25,25:51)
Idents(KB.DT) <- plyr::mapvalues(Idents(KB.DT), from = cls, to = new.cls)
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:51)
KB.DT$dt.v2.clusters <- Idents(KB.DT)

DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "dt.v2.clusters", repel = TRUE) + NoLegend()

save(KB.DT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
KB.DT
#36588 features across 506447 samples


### Remove ambiguous clusters and re-order again
#remove cl  - multiplets/ambiguous: 37, 50
Idents(KB.DT) <- "dt.v2.clusters"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:length(levels(Idents(KB.DT))))
table(Idents(KB.DT))
KB.DT <- subset(KB.DT, dt.v2.clusters %in% c(37,50), invert = TRUE)
cls <- c(1:8,12,13,9,10,11,14:30,32,31,33:36,38:41,43,42,44:49,51)
new.cls <- c(1:19,19:48)
Idents(KB.DT) <- plyr::mapvalues(Idents(KB.DT), from = cls, to = new.cls)
Idents(KB.DT) <- factor(Idents(KB.DT), levels = 1:48)
KB.DT$dt.v2.clusters <- Idents(KB.DT)

DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "dt.v2.clusters", repel = TRUE) + NoLegend()

save(KB.DT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")


###Re-identify cluster markers
dt.markers <- FindAllMarkers(KB.DT, only.pos = TRUE, max.cells.per.ident = 2000,
                             logfc.threshold = 0.25, min.pct = 0.25)
write.table(dt.markers, file="distal/Kidney_v2_distal-tubules_Cluster_markers_07102023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(dt.markers, file="distal/Kidney_v2_distal-tubules_Cluster_markers_07102023.robj")
#load("Intermediate_Objects/Kidney_v2_distal-tubules_Cluster_markers_07102023.robj")

cl.mark <- dt.markers[dt.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="distal/Kidney_v2_distal-tubules_Cluster_markers_07102023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.DT, features = top5$gene) + RotatedAxis()



###Run NSForest for markers
#Save anndata
KB.DT.sub <- subset(KB.DT, downsample = 10000)
expression_matrix <- KB.DT.sub[["RNA"]]$counts
umap <- Embeddings(KB.DT.sub, reduction = "umap")
dfobs <- KB.DT.sub@meta.data
dfobs$dt.v2.clusters <- paste0("CL",dfobs$dt.v2.clusters)
dfvar <- KB.DT.sub@assays$RNA@meta.features
adata <- AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  obsm=list('X_umap'=umap)
)
write_h5ad(adata, "Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_downsampled.h5ad")
rm(KB.DT.sub)
gc(reset = TRUE)


###Re-organize and assign markers
#Distal to collecting tubule markers
dt.markers <- c(
  "CASR","SLC12A1","UMOD","EGF",                       #TAL
  
  "ANK2","KCTD16","CLCNKA","CLDN14",                   #M-TAL
  "BMP6","RIPOR3",
  "IGKC",
  
  #"CACNB2",
  "KCNMB2","RGS6",                                     #C-TAL-A
  "HOMER1","ATF3","FOSB","LDLR",
  "ABCC9",
  "ENOX1","CALCR","CNTN1",
  "CABP1","RBM20",
  "PDE3A",
  
  "BBOX1","NOS1","ROBO2",                              #MD
  "TMPRSS4","CAPN8","CD44",
  
  "TENM4","FGF14","PXDNL","GRM8",                      #C-TAL-B
  "CLDN16","TMEM207",#"TMEM52B",                          
  "COL8A1","LINGO2","ELMO1","IFITM10","OLR1",
  
  "SLC12A3","KLHL3","CNNM2","TRPM6",                    #DCT
  
  "ADAMTS17","GRAMD1B",                                 #DCT1
  
  "SLC8A1","UNC5C","HMCN1","TRPV5",                     #DCT2 / CNT
  
  
  "HSD11B2","CALB1",                                    #CNT
  "SNTG1","ANGPT1","ATP1B3","STC1",                     
  "CTSD",                                               #dCNT
  
  "SCNN1G","SCNN1B","SGCD",                             #PC
  "GATA3","AQP2","PAPPA",
  "FGF12","KCNK13","PTCSC3","CNTNAP3B",                 #CCD-PC
  
  "MCTP1","CPAMD8","GABRA2",                            #OMCD-PC
  
  "GREB1L","SLC14A2","HS3ST5",                          #IMCD
  "TENM3","TGFBR3","AKR1C1","FAT3",
  "AOC1","TFPI2","LCN2",
  
  "ADIRF","SNCG","DHRS2","TP63",                        #PapE                              #PapE
  
  "ATP6V0D2", "ATP6V1C2", "CLNK",                       #IC
  "SLC26A7", "SLC4A1",                                  #IC-A                                   
  "HS6ST3","ABR","TMEM255A","PYY",                      #CCD-IC-A
  "NRXN3", "NXPH2", "LEF1", 
  "FAM184B","ADTRP","AQP6","STAP1",                     #OMCD-IC-A
  
  "SLC4A9", "SLC26A4", "INSRR", "TLDC2",                #IC-B
  
                       
  
  "ITGB8","PROM1",                                       #aEpi
  "ARHGAP26","RNF144B","RHEX",
  "CREB5","HIVEP3","COL4A1","CYTOR","DIAPH3","SYT16",
  "ADAMTS1","DNAH5","HAVCR1",#"TFPI","NRP1","ITGB6",
  "RAPGEF5","DLGAP1","BIRC3",
  "EZH2","MKI67","CENPP",#"TOP2A",                       #cycling
  
  
  "CKB","PEBP1","COX8A","UQCRB",                         #Injury
  "CST3","DEFB1","SPP1","IGFBP7","CLU"
)
DotPlot(KB.DT, features = unique(dt.markers)) + RotatedAxis()



###Update Experiment metadata
###Add in experiment metadata
meta <- KB.DT@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race","location","laterality","protocol",
         "tissue_type_full","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB.DT@meta.data <- meta
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "condition_level1", repel = TRUE) 
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 
DimPlot(KB.DT, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 

save(KB.DT, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")




###Barplots on stats
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
#Generate color factors
load("color_factors.robj")
meta <- KB.DT@meta.data
meta <- meta[!duplicated(meta$dt.v2.clusters),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "dt.v2.clusters", sort_label = F, colorset = "varibow")
dt.cl.cols <- meta$dt.v2.clusters_color; names(dt.cl.cols) <- meta$dt.v2.clusters_label

meta <- KB.DT@meta.data
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

save(assay.cols,condl3.cols,condl2.cols,condl1.cols,dt.cl.cols,
     file = "distal/color_factors_dt.robj")


#Condition (level 3)
Idents(object = KB.DT) <- "condition_level3"
pdf(file='distal/Kidney_DT_Condition_level3_umap.pdf',width=9,height=6)
DimPlot(KB.DT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 3"
        ) + scale_color_manual(values = alpha(condl3.cols[levels(Idents(KB.DT))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Condition (level 1)
Idents(object = KB.DT) <- "condition_level1"
pdf(file='distal/Kidney_DT_Condition_level1_umap.pdf',width=9,height=6)
DimPlot(KB.DT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.DT))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='distal/Kidney_DT_Condition_level1_umap_split.pdf',width=24,height=12)
DimPlot(KB.DT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1, split.by = "condition_level1", ncol = 3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.DT))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()




#Assay
Idents(object = KB.DT) <- "assay"
pdf(file='distal/Kidney_DT_Assay_umap.pdf',width=9,height=6)
DimPlot(KB.DT, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Assay"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(KB.DT))], 0.3), name = "Assay"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Clusters
Idents(object = KB.DT) <- "dt.v2.clusters"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = c(1:48))
pdf(file='distal/Kidney_DT_Clusters_umap.pdf',width=8,height=6)
DimPlot(KB.DT, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(dt.cl.cols[levels(Idents(KB.DT))], 0.3), name = "Clusters"
        ) + NoLegend()
dev.off()


pdf(file='distal/Cluster_Contributions_Barplot.pdf',width=6,height=24)
layout(matrix(c(1,1:7), nrow = 8, ncol = 1, byrow = TRUE))
col.order <- levels(Idents(KB.DT))
prop1 <- prop.table(table(KB.DT$patient,KB.DT$dt.v2.clusters), margin = 2)[,col.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop2 <- prop.table(table(KB.DT$assay,KB.DT$dt.v2.clusters), margin = 2)[,col.order]
barplot(prop2,main = "Assay Proportions", cex.names = 0.5, 
        names.arg = colnames(prop2),las=2, col = as.character(assay.cols[rownames(prop2)]))
#red = snRNA-seq
prop3 <- prop.table(table(KB.DT$condition_level3,KB.DT$dt.v2.clusters), margin = 2)[,col.order]
barplot(prop3,main = "Condition l3 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(condl3.cols[rownames(prop3)]))

prop4 <- prop.table(table(KB.DT$condition_level1,KB.DT$dt.v2.clusters), margin = 2)[,col.order]
barplot(prop4,main = "Condition l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop4),las=2, col = as.character(condl1.cols[rownames(prop4)]))
prop5 <- prop.table(table(KB.DT$region_level2,KB.DT$dt.v2.clusters), margin = 2)[,col.order]
barplot(prop5,main = "Region l2 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop5),las=2, col = as.character(region.l2.cols[rownames(prop5)]))
batch.entropy<-table(KB.DT$patient, KB.DT$dt.v2.clusters)
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Experiment Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(KB.DT$dt.v2.clusters)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)
dev.off()



