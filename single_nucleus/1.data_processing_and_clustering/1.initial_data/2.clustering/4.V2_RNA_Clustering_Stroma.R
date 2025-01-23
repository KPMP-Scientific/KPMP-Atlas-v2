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
source("misc/utils.R")

setwd("~/data/net/home/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")
#load("color_factors.robj")

KB.sub <- subset(KB, idents = c(1:100))
KB.sub <- subset(KB.sub, library %in% c("KM150", "KB54"), invert = TRUE)
KB.sub
#36588 features across 1240618 samples within 1 assay



##Subset to stroma only
KB.STR <- subset(KB.sub, v1.class %in% "stroma cells")
KB.STR
#36588 features across 101268 samples


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.STR[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#3641 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3641, maxit=1000)

# Generate K-nearest neighbour graph
#p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
#p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3641]
VariableFeatures(KB.STR) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.STR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.STR))
KB.STR <- RunUMAP(object = KB.STR, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.2, return.model = T)
KB.STR@meta.data$pagoda_k100_infomap_str <- k100infomap[rownames(KB.STR@meta.data)]

DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB.STR, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()

FeaturePlot(KB.STR, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k100_round1.rda")

Idents(KB.STR)
KB.STR[["RNA"]] <- as(object = KB.STR[["RNA"]], Class = "Assay")
KB.STR <- ScaleData(KB.STR)

###Identify ambiguous low quality - annotate cluster identities
Idents(KB.STR) <- "pagoda_k100_infomap_str"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = 1:length(levels(Idents(KB.STR))))
levels(Idents(KB.STR)) <- paste("CL", levels(Idents(KB.STR)), sep = "")
celltype <- Idents(KB.STR)

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

max.idents <- getMaxIdents(KB.STR, celltype)
max.idents

                                             
###Check for mean genes                                             
VlnPlot(KB.STR, features = c("nFeature_RNA"), pt.size = 0)

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


DotPlot(KB.STR, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.STR, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

###Check for DEGs
CL48_markers <- FindMarkers(object = KB.STR, ident.1 = c("CL48"), 
                            only.pos = TRUE)


###Ambiguous clusters (checking mean genes, markers, DEGs)
#CL 9 - TAL multiplet
#CL 13 - EC multiplet
#CL 14 - PT multiplet
#CL 17 - IC multiplet
#CL 21 - PC multiplet
#CL 22 - PT multiplet
#CL 25 - TL multiplet
#CL 49 - MAC multiplet


###Subset multiplets
to.remove <- c("CL9","CL13","CL14","CL17","CL21","CL22","CL25","CL49")
KB.STR <- subset(KB.STR, idents = to.remove, invert = TRUE)
KB.STR
#36588 features across 91693 samples




###Repeat subclustering
##Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.STR[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#1795 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3447, maxit=1000)

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
sn.od.genes <- rownames(var.info)[1:3447]
VariableFeatures(KB.STR) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.STR[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.STR))
KB.STR <- RunUMAP(object = KB.STR, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.1, return.model = T)
KB.STR@meta.data$pagoda_k100_infomap_str <- k100infomap[rownames(KB.STR@meta.data)]
KB.STR@meta.data$pagoda_k200_infomap_str <- k200infomap[rownames(KB.STR@meta.data)]
KB.STR@meta.data$pagoda_k500_infomap_str <- k500infomap[rownames(KB.STR@meta.data)]

DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB.STR, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_str", repel = TRUE) + NoLegend()
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.STR, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k100.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k200.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Stroma_k500.rda")


save(KB.STR, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")



###Add in prior analysis of Stroma subclusters to guide cell type annotations
load("~/data/pod/blake_LTS/Atlas_V2/Stroma/Kidney_BUKMAP-KPMP_05042022_v2_Seurat_FIB-only_b.rda")
meta.str <- KBR.FIB@meta.data
table(meta.str[rownames(KB.STR@meta.data),]$fib.clusters)
KB.STR@meta.data$previous.str.clusters <- meta.str[rownames(KB.STR@meta.data),]$fib.clusters
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "previous.str.clusters", repel = TRUE) + NoLegend()

save(KB.STR, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
rm(KBR.FIB)
gc(reset = TRUE)


###Update Experiment metadata
###Add in experiment metadata
meta <- KB.STR@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race","location","laterality","protocol",
         "tissue_type_full","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB.STR@meta.data <- meta
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "condition_level1", repel = TRUE) 
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 

save(KB.STR, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")



###Annotation of k200 clusters
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
Idents(KB.STR) <- "pagoda_k200_infomap_str"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = 1:length(levels(Idents(KB.STR))))
levels(Idents(KB.STR)) <- paste("CL", levels(Idents(KB.STR)), sep = "")
celltype <- Idents(KB.STR)

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

max.idents <- getMaxIdents(KB.STR, celltype)
max.idents

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$previous.str.clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$v1.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$previous.str.clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$previous.str.clusters)),
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

max.idents <- getMaxIdents(KB.STR, celltype)
max.idents


#CL 1: C-FIB-PATH (previous cl 7)
#CL 2: C-FIB (previous cl 8) (no distinct)
#CL 3: (low genes, no distinct markers) epithelial multiplet
#CL 4: aFIB-RSPO3 (previous cl 14)
#CL 5: C-FIB-OSMRhi (previous cl 10)
#CL 6: C-FIB (previous cl 8) (low genes, no distinct markers) Merge with CL1 or 2?
#CL 7: newC-FIB (previous cl 12)
#CL 8: (no distinct markers) (top markers are TAL) multiplet
#CL 9: C/M-FIB (previous cl 6)
#CL 10: C-MYOF (previous cl 16)
#CL 11: VSMC/P (v1 cl 74)
#CL 12: dFIB (lower genes, no distinct markers)
#CL 13: VSMC/P (v1 cl 74)
#CL 14: MC (v1 cl 71)
#CL 15: C-FIB-Path (previous cl 7)
#CL 16: (high genes) Immune multiplet
#CL 17: REN (previous cls 71/72), REN+
#CL 18: VSMC/P (v1 cl 74)
#CL 19: C/M-MYOF (previous cl 17)
#CL 20: VSMC (v1 cl 73)
#CL 21: dVSMC (v1 cl 75)
#CL 22: aFIB-CD34 (previous cl 15)
#CL 23: VSMC/P (v1 cl 74) (no distinct markers) (top markers are TAL) multiplet
#CL 24: dM-FIB (previous cl 5) (no distinct markers)
#CL 25: nM-FIB (previous cl 1)
#CL 26: Ad (previous cl 13)
#CL 27: M-FIB (previous cl 2)
#CL 28: nM-FIB (previous cl 1)
#CL 29: M-FIB (previous cl 4)
#CL 30: dM-FIB (previous cl 5) (low genes detected) (mostly from one patient) merge with cl 24
#CL 31: M-MYOF (previous cl 18)


##To do:
##Check batch, any clusters associated with only 1-2 individuals
##Check marker genes, sufficient markers to be distinct subpopulation? merge at K500?
##check for splitting at k100
##check region composition

VlnPlot(KB.STR, features = c("nFeature_RNA"), pt.size = 0)
DotPlot(KB.STR, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.STR, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

#Check for marker genes
str.markers <- FindAllMarkers(KB.STR, only.pos = TRUE, max.cells.per.ident = 2000,
                              logfc.threshold = 0.25, min.pct = 0.25)
write.table(str.markers, file="Intermediate_Objects/Kidney_v2_Stroma_k200_Cluster_markers_05082023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(str.markers, file="Intermediate_Objects/Kidney_v2_Stroma_k200_Cluster_markers_05082023.robj")
#load("Intermediate_Objects/Kidney_v2_Stroma_k200_Cluster_markers_05082023.robj")

cl.mark <- str.markers[str.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="Intermediate_Objects/Kidney_v2_Stroma_k200_Cluster_markers_05082023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.STR, features = top5$gene) + RotatedAxis()

#Check clusters 31 and 17 for patient distribution
table(KB.STR$pagoda_k200_infomap_str,KB.STR$patient)[c(3,6,8,30),]
table(KB.STR$pagoda_k200_infomap_str,KB.STR$region_level2)[30,]


###compare cluster resolutions
Idents(KB.STR) <- "pagoda_k200_infomap_str"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = 1:length(levels(Idents(KB.STR))))

library("scrattch.hicat")

#Pagoda2 K200 and K500
final<-Idents(KB.STR)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.STR$pagoda_k500_infomap_str
names(prop)<-rownames(KB.STR@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k500_Stroma.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k500", labels=compare.result$cl.id.map$old)
dev.off()


#Pagoda2 K200 and K100
final<-Idents(KB.STR)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.STR$pagoda_k100_infomap_str
names(prop)<-rownames(KB.STR@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k100_Stroma.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k100", labels=compare.result$cl.id.map$old)
dev.off()



### Remove ambiguous clusters
#remove cl 3, 8, 16, 23 - multiplets
Idents(KB.STR) <- "pagoda_k200_infomap_str"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = 1:length(levels(Idents(KB.STR))))
KB.STR <- subset(KB.STR, idents = c(3,8,16,23), invert = TRUE)


### Re-order clusters
#1 CL 25: nM-FIB (previous cl 1)
#2 CL 28: nM-FIB (previous cl 1)
#3 CL 27: M-FIB (previous cl 2)
#4 CL 29: M-FIB (previous cl 4)
#5 CL 24: dM-FIB (previous cl 5) (no distinct markers)
#6 CL 30: dM-FIB (previous cl 5) (low genes detected) (mostly from one patient) merge with cl 24
#7 CL 9: C/M-FIB (previous cl 6)
#8 CL 31: M-MYOF (previous cl 18)
#9 CL 2: C-FIB (previous cl 8) (no distinct)
#10 CL 6: C-FIB (previous cl 8) (low genes, no distinct markers) Merge with CL1 or 2?
#11 CL 1: C-FIB-PATH (previous cl 7)
#12 CL 5: C-FIB-OSMRhi (previous cl 10)
#13 CL 15: C-FIB-Path (previous cl 7)
#14 CL 10: C-MYOF (previous cl 16)
#15 CL 12: dFIB (lower genes, no distinct markers)
#16 CL 4: aFIB-RSPO3 (previous cl 14)
#17 CL 22: aFIB-CD34 (previous cl 15)
#18 CL 7: newC-FIB (previous cl 12)
#19 CL 19: C/M-MYOF (previous cl 17)
#20 CL 14: MC (v1 cl 71)
#21 CL 17: REN (previous cls 71/72), REN+
#22 CL 20: VSMC (v1 cl 73)
#23 CL 18: VSMC/P (v1 cl 74)
#24 CL 13: VSMC/P (v1 cl 74)
#25 CL 11: VSMC/P (v1 cl 74)
#26 CL 21: dVSMC (v1 cl 75)
#27 CL 26: Ad (previous cl 13)

cls <- c(25,28,27,29,24,30,9,31,2,6,1,5,15,10,12,4,22,7,19,14,17,20,18,13,11,21,26)
new.cls <- c(1:length(cls))
Idents(KB.STR) <- plyr::mapvalues(Idents(KB.STR), from = cls, to = new.cls)
Idents(KB.STR) <- factor(Idents(KB.STR), levels = 1:length(cls))
KB.STR$str.v2.clusters <- Idents(KB.STR)


KB.STR <- RunUMAP(object = KB.STR, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.1, return.model = T)

DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_str", repel = TRUE) + NoLegend()

DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "str.v2.clusters", repel = TRUE) + NoLegend()
DimPlot(KB.STR, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
table(KB.STR$v1.subclass.l3)

FeaturePlot(KB.STR, features = "TOP2A", reduction = "umap")

#Cycling cells not separating out
save(KB.STR, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
KB.STR
#36588 features across 77378 samples



###Re-identify FIB cluster markers
str.markers <- FindAllMarkers(KB.STR, only.pos = TRUE, max.cells.per.ident = 2000,
                              logfc.threshold = 0.25, min.pct = 0.25)
write.table(str.markers, file="stroma/Kidney_v2_Stroma_Cluster_markers_05092023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(str.markers, file="stroma/Kidney_v2_Stroma_Cluster_markers_05092023.robj")
#load("stroma/Kidney_v2_Immune_Cluster_markers_05092023.robj")

cl.mark <- str.markers[str.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="stroma/Kidney_v2_Stroma_Cluster_markers_05092023_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.STR, features = top5$gene) + RotatedAxis()



###Barplots on stats
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.rda")
#Generate color factors
load("color_factors.robj")
load("immune/color_factors_imm.robj")
meta <- KB.STR@meta.data
meta <- meta[!duplicated(meta$str.v2.clusters),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "str.v2.clusters", sort_label = F, colorset = "varibow")
str.cl.cols <- meta$str.v2.clusters_color; names(str.cl.cols) <- meta$str.v2.clusters_label

save(str.cl.cols,
     file = "stroma/color_factors_str.robj")


#Condition (level 3)
Idents(object = KB.STR) <- "condition_level3"
pdf(file='stroma/Kidney_Stroma_Condition_level3_umap.pdf',width=8,height=6)
DimPlot(KB.STR, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 3"
        ) + scale_color_manual(values = alpha(condl3.cols[levels(Idents(KB.STR))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Condition (level 1)
Idents(object = KB.STR) <- "condition_level1"
pdf(file='stroma/Kidney_Stroma_Condition_level1_umap.pdf',width=8,height=6)
DimPlot(KB.STR, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.STR))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='stroma/Kidney_Stroma_Condition_level1_umap_split.pdf',width=24,height=12)
DimPlot(KB.STR, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3, split.by = "condition_level1", ncol = 3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.STR))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Region Level 2
Idents(object = KB.STR) <- "region_level2"
pdf(file='stroma/Kidney_Stroma_Region_level2_umap.pdf',width=7,height=6)
DimPlot(KB.STR, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(region.l2.cols[levels(Idents(KB.STR))], 0.3), name = "Region"
        ) 
dev.off()

#Region Level 1
Idents(object = KB.STR) <- "region_level1"
pdf(file='stroma/Kidney_Stroma_Region_level1_umap.pdf',width=7,height=6)
DimPlot(KB.STR, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(region.l1.cols[levels(Idents(KB.STR))], 0.3), name = "Region"
        ) 
dev.off()


#Assay
Idents(object = KB.STR) <- "assay"
pdf(file='stroma/Kidney_Immune_Assay_umap.pdf',width=8,height=6)
DimPlot(KB.STR, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Assay"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(KB.STR))], 0.3), name = "Assay"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Clusters
Idents(object = KB.STR) <- "str.v2.clusters"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = c(1:28))
pdf(file='stroma/Kidney_Immune_Clusters_umap.pdf',width=7,height=6)
DimPlot(KB.STR, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.2,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(str.cl.cols[levels(Idents(KB.STR))], 0.3), name = "Clusters"
        ) + NoLegend()
dev.off()


pdf(file='stroma/Cluster_Contributions_Barplot.pdf',width=6,height=18)
layout(matrix(c(1,1:6), nrow = 7, ncol = 1, byrow = TRUE))
col.order <- levels(Idents(KB.STR))
prop1 <- prop.table(table(KB.STR$patient,KB.STR$str.v2.clusters), margin = 2)[,col.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop2 <- prop.table(table(KB.STR$assay,KB.STR$str.v2.clusters), margin = 2)[,col.order]
barplot(prop2,main = "Assay Proportions", cex.names = 0.5, 
        names.arg = colnames(prop2),las=2, col = as.character(assay.cols[rownames(prop2)]))
#red = snRNA-seq
prop3 <- prop.table(table(KB.STR$condition_level3,KB.STR$str.v2.clusters), margin = 2)[,col.order]
barplot(prop3,main = "Condition l3 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(condl3.cols[rownames(prop3)]))

prop4 <- prop.table(table(KB.STR$condition_level1,KB.STR$str.v2.clusters), margin = 2)[,col.order]
barplot(prop4,main = "Condition l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop4),las=2, col = as.character(condl1.cols[rownames(prop4)]))

prop5 <- prop.table(table(KB.STR$region_level2,KB.STR$str.v2.clusters), margin = 2)[,col.order]
barplot(prop5,main = "Region L2 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop5),las=2, col = as.character(region.l2.cols[rownames(prop5)]))

batch.entropy<-table(KB.STR$patient, KB.STR$str.v2.clusters)
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Experiment Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(KB.STR$str.v2.clusters)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)
dev.off()







###Correlation with public data
library(corrplot)
###Buechler & Pradhan, ..., Bourgon, Müller, Turley published in Nature, Vol 593, 2021 . 
##Download seurat objects from https://www.fibroxplorer.com/home

##Mouse Steady State
mSS <- readRDS("~/data/pod/blake_LTS/public/Buechler2021/Mouse_SS_Fibro.RDS")
mSS
#21087 features across 120583 samples
Idents(mSS) <- "ClustName"
table(Idents(mSS))
#Ccl19  Cxcl12    Comp Col15a1    Coch    Pi16   Fbln1    Npnt    Hhip    Bmp4 
#15075    8636    7075   44025    6116   23557    4905    7847    1699    1648 

#Intersect variable genes
common.human <- VariableFeatures(KB.STR)[which(VariableFeatures(KB.STR)[1:2000] %in% toupper(VariableFeatures(mSS)[1:2000]))]
common.mouse <- VariableFeatures(mSS)[which(toupper(VariableFeatures(mSS)[1:2000]) %in% (VariableFeatures(KB.STR)[1:2000]))]

KB.STR <- ScaleData(KB.STR, features = rownames(KB.STR))
ave.KB.STR<-AverageExpression(KB.STR, features = common.human, slot = "scale.data")
ave.mSS<-AverageExpression(mSS, features = common.mouse, slot = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.STR$RNA),as.data.frame(ave.mSS$RNA)))
ave.cor<-ave.cor[1:27,28:37]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='stroma/Stromal_Clusters_Buechler2021_SteadyStateMouse.pdf',width=6,height=12)
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()


##Mouse Perturbed State
mPS <- readRDS("~/data/pod/blake_LTS/public/Buechler2021/Mouse_PS_Fibro.RDS")
mPS
#21087 features across 99596 samples
Idents(mPS) <- "ClustName"
table(Idents(mPS))
#Lrrc15     Npnt  Col15a1     Comp     Pi16   Cxcl12    Cxcl5    Ccl19 Adamdec1     Hhip 
#14955     9485    23191     5158    12453     9669     4425    15593     1411     3256

#Intersect variable genes
common.human <- VariableFeatures(KB.STR)[which(VariableFeatures(KB.STR)[1:2000] %in% toupper(VariableFeatures(mPS)[1:2000]))]
common.mouse <- VariableFeatures(mPS)[which(toupper(VariableFeatures(mPS)[1:2000]) %in% (VariableFeatures(KB.STR)[1:2000]))]

ave.KB.STR<-AverageExpression(KB.STR, features = common.human, slot = "scale.data")
ave.mPS<-AverageExpression(mPS, features = common.mouse, slot = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.STR$RNA),as.data.frame(ave.mPS$RNA)))
ave.cor<-ave.cor[1:27,28:37]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='stroma/Stromal_Clusters_Buechler2021_PerturbedStateMouse.pdf',width=6,height=12)
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()



##Human Perturbed State
hPS <- readRDS("~/data/pod/blake_LTS/public/Buechler2021/Human_PS_Fibro.RDS")
hPS
#17410 features across 10355 samples
Idents(hPS) <- "ClustName"
table(Idents(hPS))
#LRRC15 ADAMDEC1   COL3A1    CCL19     PI16     NPNT 
#2130     1625     3247      357     1979     1017 

#Intersect all variable genes
common.genes <- VariableFeatures(KB.STR)[(VariableFeatures(KB.STR) %in% VariableFeatures(hPS))]

ave.KB.STR<-AverageExpression(KB.STR, features = common.genes, slot = "scale.data")
ave.hPS<-AverageExpression(hPS, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.STR$RNA),as.data.frame(ave.hPS$RNA)))
ave.cor<-ave.cor[1:27,28:33]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='stroma/Stromal_Clusters_Buechler2021_PerturbedStateHuman.pdf',width=6,height=8)
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()



###Korsunsky et al., 2022. https://sandbox.zenodo.org/record/772596#.Y0X0TOzMLJx
##Human Perturbed State
counts <- ReadMtx(mtx = "~/data/pod/blake_LTS/public/Korsunsky2022/umis.mtx",
                  cells = "~/data/pod/blake_LTS/public/Korsunsky2022/umis_barcodes.tsv",
                  features = "~/data/pod/blake_LTS/public/Korsunsky2022/umis_genes.tsv",
                  feature.column = 1)
meta <- read.table("~/data/pod/blake_LTS/public/Korsunsky2022/metaData.txt",sep="\t",header=TRUE,row.names=1)
umap <- read.table("~/data/pod/blake_LTS/public/Korsunsky2022/umap_integrated.txt",sep="\t",header=TRUE,row.names=1)

counts <- counts[,rownames(meta)]
fib.atlas <- CreateSeuratObject(counts, meta.data = meta)
fib.atlas
#19952 features across 79968 samples

fib.atlas <- NormalizeData(fib.atlas)
fib.atlas <- FindVariableFeatures(fib.atlas)
fib.atlas <- ScaleData(fib.atlas, features = rownames(fib.atlas))
fib.atlas <- RunPCA(fib.atlas)
fib.atlas <- RunUMAP(fib.atlas, reduction = "pca", dims = 1:30, return.model = T, verbose = F)
Idents(fib.atlas) <- "cell_type_integrated"
DimPlot(fib.atlas, group.by = "cell_type_integrated", reduction = "umap")
colnames(umap) <- c("umap_1", "umap_2")
umap <- umap[rownames(Embeddings(fib.atlas, reduction = "umap")),]
embeddings <- matrix(data = c(as.numeric(umap$dim_1),as.numeric(umap$dim_2)), ncol = 2)
rownames(embeddings) <- rownames(umap)
colnames(embeddings) <- colnames(umap)

fib.atlas[["ref.umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "refumap_", assay = DefaultAssay(fib.atlas))
DimPlot(fib.atlas, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "cell_type_integrated", repel = TRUE) + NoLegend()


table(Idents(fib.atlas))
table(fib.atlas$sample_type)
table(fib.atlas$donor_id)

#Intersect variable genes
common.genes <- VariableFeatures(fib.atlas)[VariableFeatures(fib.atlas) %in% (VariableFeatures(KB.STR))]
common.genes <- VariableFeatures(fib.atlas)[VariableFeatures(fib.atlas) %in% rownames(KB.STR)]
common.genes <- VariableFeatures(fib.atlas)[1:1000][VariableFeatures(fib.atlas)[1:1000] %in% rownames(KB.STR)]
common.genes <- VariableFeatures(fib.atlas)[1:1000][VariableFeatures(fib.atlas)[1:1000] %in% (VariableFeatures(KB.STR))]

ave.KB.STR<-AverageExpression(KB.STR, features = common.genes, slot = "scale.data")
ave.fib.atlas<-AverageExpression(fib.atlas, features = common.genes, slot = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.STR$RNA),as.data.frame(ave.fib.atlas$RNA)))
ave.cor<-ave.cor[1:27,28:41]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='stroma/Stromal_Clusters_Korsunsky2022_Integrated.pdf',width=6,height=12)
corrplot(ave.cor, col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()


###Run NSForest for markers
#Save anndata
expression_matrix <- KB.STR[["RNA"]]$counts
umap <- Embeddings(KB.STR, reduction = "umap")
dfobs <- KB.STR@meta.data
dfobs$str.v2.clusters <- paste0("CL",dfobs$str.v2.clusters)
dfvar <- KB.STR@assays$RNA@meta.features
adata <- AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  obsm=list('X_umap'=umap)
)
write_h5ad(adata, "Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered.h5ad")



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

pdf(file='stroma/Immune_Cluster_1-19_Stromal_Marker_Dotplot.pdf',width=20,height=6)
DotPlot(KB.STR, features = fib.markers, idents = c(1:19)) + RotatedAxis()
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

pdf(file='stroma/Immune_Cluster_20-27_Stromal_Marker_Dotplot.pdf',width=12,height=5)
DotPlot(KB.STR, features = vsm.markers, idents = c(20:27)) + RotatedAxis()
dev.off()


