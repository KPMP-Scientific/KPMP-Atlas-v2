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
KB.sub <- subset(KB, idents = c(1:100))
KB.sub <- subset(KB.sub, library %in% c("KM150", "KB54"), invert = TRUE)
KB.sub
#36588 features across 1240618 samples within 1 assay



##Subset to vasculature only
KB.EC <- subset(KB.sub, v1.subclass.l1 %in% "EC")
KB.EC
#36588 features across 99540 samples


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.EC[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#3293 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3293, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3293]
VariableFeatures(KB.EC) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.EC[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.EC))
KB.EC <- RunUMAP(object = KB.EC, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.2, return.model = T)
KB.EC@meta.data$pagoda_k100_infomap_ec <- k100infomap[rownames(KB.EC@meta.data)]

DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.EC, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "Intermediate_Objects/Kidney_AtlasV2_Pagoda2_vasculature_k100_round1.rda")

Idents(KB.EC)
KB.EC[["RNA"]] <- as(object = KB.EC[["RNA"]], Class = "Assay")
KB.EC <- ScaleData(KB.EC)

###Identify ambiguous low quality - annotate cluster identities
Idents(KB.EC) <- "pagoda_k100_infomap_ec"
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:length(levels(Idents(KB.EC))))
levels(Idents(KB.EC)) <- paste("CL", levels(Idents(KB.EC)), sep = "")
celltype <- Idents(KB.EC)

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

max.idents <- getMaxIdents(KB.EC, celltype)
max.idents

###Check Mean Genes
VlnPlot(KB.EC, features = c("nFeature_RNA"), pt.size = 0)

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


DotPlot(KB.EC, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.EC, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

###Check for DEGs
CL10_markers <- FindMarkers(object = KB.EC, ident.1 = c("CL10"), 
                            only.pos = TRUE)
CL12_markers <- FindMarkers(object = KB.EC, ident.1 = c("CL12"), 
                            only.pos = TRUE)

###Ambiguous clusters (checking mean genes, markers, DEGs)
#CL3 - PT multiplet
#CL7 - high genes, PC multiplet
#CL13 - TAL multiplet
#CL16 - VSM multiplet
#CL20 - low genes
#CL22 - high genes, TL multiplet
#CL25 - PT multiplet
#CL27 - high genes, IC multiplet
#CL29 - PT multiplet
#CL31 - high genes, fibroblast multiplet
#CL34 - immune multiplet
#CL35 - VSM/MYOF multiplet

###Subset multiplets
Idents(KB.EC) <- "pagoda_k100_infomap_ec"
to.remove <- c(3,7,13,16,20,22,25,27,29,31,34,35)
KB.EC <- subset(KB.EC, idents = to.remove, invert = TRUE)
KB.EC
#36588 features across 87602 samples




###Repeat subclustering
##Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.EC[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#2978 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 2978, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')  #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')  #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')  #ran downstream steps separately for each k value
p2$makeKnnGraph(k = 800, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')  #ran downstream steps separately for each k value

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
sn.od.genes <- rownames(var.info)[1:2978]
VariableFeatures(KB.EC) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.EC[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.EC))
KB.EC <- RunUMAP(object = KB.EC, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.1, return.model = T)
KB.EC@meta.data$pagoda_k100_infomap_ec <- k100infomap[rownames(KB.EC@meta.data)]
KB.EC@meta.data$pagoda_k200_infomap_ec <- k200infomap[rownames(KB.EC@meta.data)]
KB.EC@meta.data$pagoda_k500_infomap_ec <- k500infomap[rownames(KB.EC@meta.data)]
KB.EC@meta.data$pagoda_k800_infomap_ec <- k800infomap[rownames(KB.EC@meta.data)]

DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k800_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.EC, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_vasculature_k100.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_vasculature_k200.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_vasculature_k500.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_vasculature_k800.rda")


save(KB.EC, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")



###Update Experiment metadata
###Add in experiment metadata
meta <- KB.EC@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race","location","laterality","protocol",
         "tissue_type_full","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB.EC@meta.data <- meta
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "condition_level1", repel = TRUE) 
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 

save(KB.EC, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")



###Annotation of k200 clusters
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
Idents(KB.EC) <- "pagoda_k200_infomap_ec"
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:length(levels(Idents(KB.EC))))
levels(Idents(KB.EC)) <- paste("CL", levels(Idents(KB.EC)), sep = "")
celltype <- Idents(KB.EC)

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

max.idents <- getMaxIdents(KB.EC, celltype)
max.idents



#CL 1: EC-PTC (no distinct) (stable at k500) (low genes detected) (top markers are for epithelia) multiplet
#CL 2: EC-GC (distinct markers, merges with cl 10 at k500)
#CL 3: EC-PTC (distinct markers, merges with CL21 at k500)
#CL 4: EC-PTC (no distinct markers, stable at k500)
#CL 5: EC-GC (distinct markers, merges with cl10 at k500)
#CL 6: EC-PTC (distinct marker, stable at k500)
#CL 7: EC-DVR (distinct markers, stable at k500)
#CL 8: EC-AVR (distinct markers, stable at k500, merges with cl28)
#CL 9: EC-AEA (distinct markers, stable at k500)
#CL 10: EC-GC (distinct marker, merges with cl 5 & 2 at k500)
#CL 11: EC-LYM (distinct markers, stable at k500)
#CL 12: d?EC-PTC (no distinct markers, stable at k500) (low genes detected) 
#CL 13: EC-AVR (distinct markers, stable at k500)
#CL 14: EC-PTC (distinct markers, stable at k500)
#CL 15: EC-PTC (distinct markers, stable at k500)
#CL 16: EC-PTC (no distinct markers) (top markers are epithelial) multiplet
#CL 17: EC-DVR (distinct markers, stable at k500)
#CL 18: ?EC-PTC? (distinct markers, splits at k500)
#CL 19: ?EC-PTC? (distinct markers, stable at k500)
#CL 20: EC-AVR (distinct markers, stable at k500)
#CL 21: EC-PTC (no distinct markers, merges with cl 3 at k500) merge
#CL 22: cycEC (distinct markers, stable at k500)
#CL 23: EC-PTC (distinct markers, stable at k500)
#CL 24: d?EC-DVR (distinct markers, splits to 16 and 17 at k500) degen markers and high ER? weak epithelial markers
#CL 25: EC-PTC (distinct markers, stable at k500)
#CL 26: EC-PTC (distinct markers, merges with CL 25 at k500)
#CL 27: d?EC-AVR (no distinct markers, unstable at k500) weak epi signatures
#CL 28: EC-AVR (distinct markers, merges with CL 8 at k500)
#CL 29: EC-AVR (distinct markers, unstable at k500)
#CL 30: EC-AVR (distinct markers, unstable at k500)
#CL 31: EC-PTC (distinct markers, unstable at k500)

##To do:
##Check batch, any clusters associated with only 1-2 individuals
##Check marker genes, sufficient markers to be distinct subpopulation? merge at K500?
##check for splitting at k100

#Check for marker genes K200
ec.markers <- FindAllMarkers(KB.EC, only.pos = TRUE, max.cells.per.ident = 2000,
                              logfc.threshold = 0.25, min.pct = 0.25)
write.table(ec.markers, file="Intermediate_Objects/Kidney_v2_vasculature_k200_Cluster_markers_05122023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(ec.markers, file="Intermediate_Objects/Kidney_v2_vasculature_k200_Cluster_markers_05122023.robj")
#load("Intermediate_Objects/Kidney_v2_vasculature_k200_Cluster_markers_05122023.robj")

cl.mark <- ec.markers[ec.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="Intermediate_Objects/Kidney_v2_vasculature_k200_Cluster_markers_05012023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.EC, features = top5$gene) + RotatedAxis()


#Check for marker genes K100
Idents(KB.EC) <- "pagoda_k100_infomap_ec"
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:length(levels(Idents(KB.EC))))
ec.markers <- FindAllMarkers(KB.EC, only.pos = TRUE, max.cells.per.ident = 2000,
                             logfc.threshold = 0.25, min.pct = 0.25)
write.table(ec.markers, file="Intermediate_Objects/Kidney_v2_vasculature_k100_Cluster_markers_05122023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(ec.markers, file="Intermediate_Objects/Kidney_v2_vasculature_k100_Cluster_markers_05122023.robj")
#load("Intermediate_Objects/Kidney_v2_vasculature_k100_Cluster_markers_05122023.robj")

cl.mark <- ec.markers[ec.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="Intermediate_Objects/Kidney_v2_vasculature_k100_Cluster_markers_05012023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.EC, features = top5$gene) + RotatedAxis()


#Check clusters for patient distribution
table(KB.EC$pagoda_k200_infomap_ec,KB.EC$patient)[c(1,12,16,24,27),]
table(KB.EC$pagoda_k200_infomap_ec,KB.EC$condition_level1)


###compare cluster resolutions
Idents(KB.EC) <- "pagoda_k200_infomap_ec"
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:length(levels(Idents(KB.EC))))

library("scrattch.hicat")

#Pagoda2 K200 and K500
final<-Idents(KB.EC)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.EC$pagoda_k500_infomap_ec
names(prop)<-rownames(KB.EC@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k500_vasculature.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k500", labels=compare.result$cl.id.map$old)
dev.off()


#Pagoda2 K200 and K100
final<-Idents(KB.EC)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.EC$pagoda_k100_infomap_ec
names(prop)<-rownames(KB.EC@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k100_vasculature.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k100", labels=compare.result$cl.id.map$old)
dev.off()


#Pagoda2 K200 and K800
final<-Idents(KB.EC)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.EC$pagoda_k800_infomap_ec
names(prop)<-rownames(KB.EC@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k800_vasculature.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k800", labels=compare.result$cl.id.map$old)
dev.off()



Idents(KB.EC) <- "pagoda_k100_infomap_ec"
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:length(levels(Idents(KB.EC))))
DotPlot(KB.EC, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.EC, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()


### Remove ambiguous clusters
#remove cl 1, 16, 24 - multiplets
Idents(KB.EC) <- "pagoda_k200_infomap_ec"
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:length(levels(Idents(KB.EC))))
KB.EC <- subset(KB.EC, idents = c(1,16,24), invert = TRUE)
KB.EC <- RunUMAP(object = KB.EC, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                 min.dist = 0.1, return.model = T)
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_ec", repel = TRUE) + NoLegend()
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_ec", repel = TRUE) + NoLegend()

save(KB.EC, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")

### Re-order and Merge clusters
#1 CL 10: EC-GC (distinct marker, merges with cl 5 & 2 at k500)
#3 CL 2: EC-GC (weak distinct markers, merges with cl 10 at k500)
#4 newCL 32 k100 cluster 25 from k200 cluster 20

#5 CL 9: EC-AEA (distinct markers, stable at k500)
#2 CL 5: EC-PTC (distinct markers, merges with cl10 at k500)

#6 CL 7: EC-DVR (distinct markers, stable at k500)
#7 CL 17: EC-DVR (distinct markers, stable at k500)
#8 CL 14: M-EC-PTC (distinct markers, stable at k500)

#10 CL 30: EC-AVR (distinct markers, unstable at k500)
#11 CL 28: EC-AVR (distinct markers, merges with CL 8 at k500)
#13 CL 8: EC-AVR (distinct markers, stable at k500, merges with cl28)
#12 CL 20: EC-AVR (distinct markers, stable at k500) - split at k100 to 22 (AVR) & 25 (GC)
#14 CL 29: EC-AVR (distinct markers, unstable at k500)
#16 CL 13: EC-AVR (distinct markers, stable at k500)

#17 CL 23: EC-PTC (distinct markers, stable at k500)
#18 CL 3: EC-PTC (distinct markers, merges with CL21 at k500)
#18 CL 21: EC-PTC (no distinct markers, merges with cl 3 at k500) merge
#19 CL 4: EC-PTC (no distinct markers, stable at k500)
#20 CL 6: EC-PTC (distinct marker, stable at k500)
#21 CL 15: EC-PTC (distinct markers, stable at k500)
#22 CL 12: d?EC-PTC (no distinct markers, stable at k500) (low genes detected) 
#23 CL 18: ?EC-PTC? (distinct markers, splits at k500)
#24 CL 25: EC-PTC (distinct markers, stable at k500)
#25 CL 26: EC-PTC (distinct markers, merges with CL 25 at k500)
#26 CL 31: EC-PTC (distinct markers, unstable at k500)

#27 CL 11: EC-LYM (distinct markers, stable at k500)
#28 CL 22: cycEC (distinct markers, stable at k500)

#9 CL 19: ?EC-PTC? (distinct markers, stable at k500) - erythrocyte and epithelial markers, very few DEGs - ambiguous...
#15 CL 27: d?EC-AVR (no distinct markers, unstable at k500) epi signatures


Idents(KB.EC) <- "pagoda_k200_infomap_ec"
newcl32 <- rownames(KB.EC@meta.data[KB.EC@meta.data$pagoda_k100_infomap_ec == 25,])
KB.EC <- SetIdent(KB.EC, cells = newcl32, value = '32')
table(Idents(KB.EC))

cls <- c(10,5,2,32,9,7,17,14,19,30,28,20,8,29,27,13,23,3,21,4,6,15,12,18,25,26,31,11,22)
new.cls <- c(1:17,18,18,19:28)
Idents(KB.EC) <- plyr::mapvalues(Idents(KB.EC), from = cls, to = new.cls)
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:28)
KB.EC$ec.v2.clusters <- Idents(KB.EC)

DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "ec.v2.clusters", repel = TRUE) + NoLegend()

save(KB.EC, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
KB.EC
#36588 features across 71002 samples


##re-adjust clusters (merge/remove)
Idents(KB.EC) <- "ec.v2.clusters"
KB.EC <- subset(KB.EC, idents = c(9, 15), invert = TRUE) #multiplets

cls <- c(1,3,4,5,2,6:8,10,11,13,12,14,16:28)
new.cls <- c(1:9,10,10,11:25) #merge clusters 11 & 13 based on k500, re-order 
Idents(KB.EC) <- plyr::mapvalues(Idents(KB.EC), from = cls, to = new.cls)
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:25)
KB.EC$ec.v2.clusters <- Idents(KB.EC)


KB.EC <- RunUMAP(object = KB.EC, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                 min.dist = 0.1, return.model = T)
DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "ec.v2.clusters", repel = TRUE) + NoLegend()

save(KB.EC, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
KB.EC
#36588 features across 66920 samples

###Re-identify ec cluster markers
ec.markers <- FindAllMarkers(KB.EC, only.pos = TRUE, max.cells.per.ident = 2000,
                              logfc.threshold = 0.25, min.pct = 0.25)
write.table(ec.markers, file="vasculature/Kidney_v2_vasculature_Cluster_markers_05222023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(ec.markers, file="vasculature/Kidney_v2_vasculature_Cluster_markers_05222023.robj")
#load("vasculature/Kidney_v2_vasculature_Cluster_markers_05222023.robj")

cl.mark <- ec.markers[ec.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="vasculature/Kidney_v2_vasculature_Cluster_markers_05222023_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.EC, features = top5$gene) + RotatedAxis()



###Finalize Cluster resolution (merging and annotating)
#Check for distinct DEGs
#Check for injury/state change (injury/repair markers, %genes / %mt / % ribosomal transcripts)
#Check for patient entropy
#merge if no distinct DEGs (<2) and no obvious state change, merge to the largest associated subclass cluster 
#taking in account merging between k200 and k500/k800



###Clusters
### Re-order and Merge clusters
#1 CL 10: EC-GC (distinct marker, merges with cl 5 & 2 at k500)
#2 CL 2: EC-GC (weak distinct markers, merges with cl 10 at k500) - High %ERT
#3 newCL 32 k100 cluster 25 from k200 cluster 20

#5 CL 5: C-EC-PTC (distinct markers, merges with cl10 at k500)
#4 CL 9: EC-AEA (distinct markers, stable at k500)
#6 CL 7: EC-DVR (distinct markers, stable at k500)
#7 CL 17: EC-DVR (distinct markers, stable at k500)
#8 CL 14: M-EC-PTC (distinct markers, stable at k500)

#10 CL 8/28: EC-AVR (distinct markers)
#9 CL 30: EC-AVR (distinct markers, unstable at k500)
#11 CL 20: EC-AVR (distinct markers, stable at k500) - split at k100 to 22 (AVR) & 25 (GC)
#12 CL 29: EC-AVR (distinct markers, unstable at k500)
#13 CL 13: EC-AVR (distinct markers, stable at k500)

#15 CL 3/21: EC-PTC (1 distinct marker) - merges to 13 (EC-DVR) and 16 (EC-PTC) at k800... unclear merging
#16 CL 4: EC-PTC (1 distinct marker, stable at k500)
#17 CL 6: EC-PTC (no distinct markers, stable at k500) - merges with cluster 16 at k800
#18 CL 15: EC-PTC (distinct markers, stable at k500)
#19 CL 12: EC-PTC (no distinct markers, stable at k500) (low genes detected) - merges with cluster 16 at k800 - degenerative?
#20 CL 18: EC-PTC (1 distinct marker, splits at k500) show similar expression signature as clusters 6/7
#14 CL 23: EC-PTC (distinct markers, stable at k500) - High %ERT
#21 CL 25: EC-PTC (distinct markers, stable at k500)
#22 CL 26: EC-PTC (distinct markers, merges with CL 25 at k500)

#24 CL 11: EC-LYM (distinct markers, stable at k500)
#25 CL 22: cycEC (distinct markers, stable at k500)

#23 CL 31: EC-PTC (distinct markers, unstable at k500) only in one ref sample - High %ERT - only found in one reference sample - remove


###Check how clusters merge at k800

#Pagoda2 K200 and K800
final<-Idents(KB.EC)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("ec-cl_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.EC$pagoda_k800_infomap_ec
names(prop)<-rownames(KB.EC@meta.data)
prop <- factor(prop)

pdf("Intermediate_Objects/Pagoda_Cluster_comparision_EC-CL_vs_k800_vasculature.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k800", labels=compare.result$cl.id.map$old)
dev.off()





###Check cell type markers
#vascular markers
ec.markers <- c("PECAM1","PTPRB","FLT1",
                "EMCN",
                "HECW2","EHD3","RXFP1","PLAT", "ITGA8",                                 #EC-GC
                "FILIP1","H19","ESM1","SLC45A4","GJA5","RIPOR2",
                'PTCHD4',"ZMAT3",
                
                        
                "SULF1", "PDE3A","NKAIN2","KIAA1549L","NOS1",             #EC-AEA
                "SERPINE2","FBLN5","CXCL12",
                "VEGFC","SOX17",
                                       
                
                "ADAMTS6","PALMD","AQP1","TM4SF1",
                
                "MCTP1","SLC14A1","ENPP2",                              #EC-DVR
                "LYPD6B",                  
                'KCTD16','CELF2','SLIT3',"EPHA3",
                "ITIH5", "VWF","STXBP6","CP",
                
                "CEACAM1","PLVAP","DNASE1L3", 
                "COL15A1","FABP5",'RGCC',"ALPL",
                "HBB","HBA2","HBA1",
                
                "CACNA2D1","GPM6A","EDIL3","TLL1","ZNF385D","NR2F2",            #EC-AVR
                'CCN2', 'RELN', "BMP6","EFEMP1",
                "RIPOR3","DOK5",
                "VCAN",
                "MX2","RSAD2","IFIT1","ISG15","IFIT3",
                "AQP2", "AKR1B1", "FXYD4",
                "PITPNC1", "GRB10", "SLCO2A1", "RAPGEF4",             #EC-PTC
                'ADAMTSL1',"CPE","TRABD2B",
                'CAVIN2', "RNASE1", "IFITM3","CCL14","CA4",
                'FLRT2', 'KCNIP4', "CMKLR1",
                'NAV3',"LNX1","ASIC2",
                "AKAP12",'SLC2A3', 'PDLIM1', 'ADAMTS9', 'SPRY1',
                'MYO1B', 'ARHGAP18', 'CSGALNACT1',
                "NETO2","PDE4B","PCSK6",
                "ARGLU1",
                'SLC6A6', 'CHRM3',"AFF3",
                'ICAM1', 'CCL2',"TNFAIP3",'SELE',"CXCL2",
                "CXCL10","CXCL11","GBP1","PLA1A","CTSS",
                'HSPA4', 'HSPD1','HSPA6','DNAJB1','HSPA1B','DNAJA4',"BAG3",
                "MMRN1","CD36","TBX1","PKHD1L1","PROX1",
                "TOP2A","MKI67",                                          #cycling
                "B2M", "TMSB4X", "TMSB10", "HLA-B", "IGFBP5"              #degenerative
)


DotPlot(KB.EC, features = ec.markers, idents = c(1:28)) + RotatedAxis()
load("color_factors.robj")

row.order <- levels(Idents(KB.EC))
prop1 <- prop.table(table(KB.EC$patient, KB.EC$ec.v2.clusters), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop1 <- prop.table(table(KB.EC$region_level2, KB.EC$ec.v2.clusters), margin = 2)[,row.order]
barplot(prop1,main = "Region Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(region.l2.cols[rownames(prop1)]))
prop1 <- prop.table(table(KB.EC$condition_level1, KB.EC$ec.v2.clusters), margin = 2)[,row.order]
barplot(prop1,main = "Condition Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(cond.l2.cols))

VlnPlot(object = KB.EC, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)

table(Idents(KB.EC))


###Organize final clusters
Idents(KB.EC) <- "ec.v2.clusters"
KB.EC <- subset(KB.EC, idents = c(23), invert = TRUE) #ambiguous

cls <- c(1,2,3,5,4,6,7,8,10,9,11,12,13,15,16,17,18,19,20,14,21,22,24,25)
new.cls <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15,16,15,17,18,19,20,21,22) #merge clusters 16, 17, 19 based on k800, re-order 
Idents(KB.EC) <- plyr::mapvalues(Idents(KB.EC), from = cls, to = new.cls)
Idents(KB.EC) <- factor(Idents(KB.EC), levels = 1:22)
KB.EC$ec.v2.clusters <- Idents(KB.EC)

DimPlot(KB.EC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "ec.v2.clusters", repel = TRUE) + NoLegend()

save(KB.EC, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")





###Re-identify ec cluster markers
ec.markers <- FindAllMarkers(KB.EC, only.pos = TRUE, max.cells.per.ident = 2000,
                             logfc.threshold = 0.25, min.pct = 0.25)
write.table(ec.markers, file="vasculature/Kidney_v2_vasculature_Cluster_markers_06192023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(ec.markers, file="vasculature/Kidney_v2_vasculature_Cluster_markers_06192023.robj")
#load("vasculature/Kidney_v2_vasculature_Cluster_markers_06192023.robj")

cl.mark <- ec.markers[ec.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="vasculature/Kidney_v2_vasculature_Cluster_markers_06192023_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.EC, features = top5$gene) + RotatedAxis()




###Run NSForest for markers
#Save anndata
expression_matrix <- KB.EC[["RNA"]]$counts
umap <- Embeddings(KB.EC, reduction = "umap")
dfobs <- KB.EC@meta.data
dfobs$ec.v2.clusters <- paste0("CL",dfobs$ec.v2.clusters)
dfvar <- KB.EC@assays$RNA@meta.features
adata <- AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  obsm=list('X_umap'=umap)
)
write_h5ad(adata, "Intermediate_Objects/Kidney_AtlasV2_Seurat_Vasculature_subset_filtered.h5ad")



###Re-organize and assign vascular markers
#vascular markers
ec.markers <- c("PECAM1","PTPRB","FLT1",
                "EMCN","HECW2","ITGA8","EHD3","RXFP1",                     #EC-GC      
                "MGP","SOST",
                'PTCHD4',"ZMAT3",
                "AQP1","FILIP1","H19","ESM1","SLC45A4",
                
                "SULF1", "PDE3A","NKAIN2","NOS1",                          #EC-AEA
                
                "ADAMTS6","MCTP1","PALMD","SLC14A1",                       #EC-DVR  
                "LYPD6B","EPHA3",
                "ITIH5","STXBP6","CP",
                
                "CEACAM1","PLVAP","DNASE1L3",                              #EC-AVR
                "COL15A1","FABP5","ALPL",
                "GPM6A","NR2F2","ZNF385D","RANBP3L",
                "EDIL3","COL3A1",
                "CCN2","RELN","BMP6","EFEMP1",
                "MX2","RSAD2","ISG15","IFIT1",
                
                "SLCO2A1",                                                 #EC-PTC
                "RYR3","ADGRG6","CPE","TRABD2B","ADAMTSL1",
                "CMKLR1","DOK6",
                "NAV3","SYCP2L",
                "USP31","MYO1B","LAMA4","NETO2",
                "SLC6A6","FUT8","ATP13A3","AFF3",
                "IFITM3","HLA-DRA","HLA-DRB1","CAVIN2","CCL14","CA4",

                'ICAM1',"TNFAIP3",'SELE',"CXCL2",'CCL2',"VCAM1",            #Injury
                "CXCL10","GBP1","CXCL11","CTSS",
                "MMRN1","CD36","TBX1","PROX1",
                "TOP2A","MKI67","CENPF"                                     #cycling
                )


DotPlot(KB.EC, features = ec.markers, idents = c(1:28)) + RotatedAxis()





###Barplots on stats
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_vasculature_subset_filtered.rda")
#Generate color factors
load("color_factors.robj")
meta <- KB.EC@meta.data
meta <- meta[!duplicated(meta$ec.v2.clusters),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "ec.v2.clusters", sort_label = F, colorset = "varibow")
ec.cl.cols <- meta$ec.v2.clusters_color; names(ec.cl.cols) <- meta$ec.v2.clusters_label

meta <- KB.EC@meta.data
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

save(assay.cols,condl3.cols,condl2.cols,condl1.cols,ec.cl.cols,
     file = "vasculature/color_factors_ec.robj")


#Condition (level 3)
Idents(object = KB.EC) <- "condition_level3"
pdf(file='vasculature/Kidney_vasculature_Condition_level3_umap.pdf',width=9,height=6)
DimPlot(KB.EC, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 3"
        ) + scale_color_manual(values = alpha(condl3.cols[levels(Idents(KB.EC))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Condition (level 1)
Idents(object = KB.EC) <- "condition_level1"
pdf(file='vasculature/Kidney_vasculature_Condition_level1_umap.pdf',width=9,height=6)
DimPlot(KB.EC, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.EC))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='vasculature/Kidney_vasculature_Condition_level1_umap_split.pdf',width=24,height=12)
DimPlot(KB.EC, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3, split.by = "condition_level1", ncol = 3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.EC))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()




#Assay
Idents(object = KB.EC) <- "assay"
pdf(file='vasculature/Kidney_vasculature_Assay_umap.pdf',width=9,height=6)
DimPlot(KB.EC, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Assay"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(KB.EC))], 0.3), name = "Assay"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Clusters
Idents(object = KB.EC) <- "ec.v2.clusters"
Idents(KB.EC) <- factor(Idents(KB.EC), levels = c(1:22))
pdf(file='vasculature/Kidney_vasculature_Clusters_umap.pdf',width=8,height=6)
DimPlot(KB.EC, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(ec.cl.cols[levels(Idents(KB.EC))], 0.3), name = "Clusters"
        ) + NoLegend()
dev.off()


pdf(file='vasculature/Cluster_Contributions_Barplot.pdf',width=6,height=24)
layout(matrix(c(1,1:7), nrow = 8, ncol = 1, byrow = TRUE))
col.order <- levels(Idents(KB.EC))
prop1 <- prop.table(table(KB.EC$patient,KB.EC$ec.v2.clusters), margin = 2)[,col.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop2 <- prop.table(table(KB.EC$assay,KB.EC$ec.v2.clusters), margin = 2)[,col.order]
barplot(prop2,main = "Assay Proportions", cex.names = 0.5, 
        names.arg = colnames(prop2),las=2, col = as.character(assay.cols[rownames(prop2)]))
#red = snRNA-seq
prop3 <- prop.table(table(KB.EC$condition_level3,KB.EC$ec.v2.clusters), margin = 2)[,col.order]
barplot(prop3,main = "Condition l3 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(condl3.cols[rownames(prop3)]))

prop4 <- prop.table(table(KB.EC$condition_level1,KB.EC$ec.v2.clusters), margin = 2)[,col.order]
barplot(prop4,main = "Condition l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop4),las=2, col = as.character(condl1.cols[rownames(prop4)]))
prop5 <- prop.table(table(KB.EC$region_level2,KB.EC$ec.v2.clusters), margin = 2)[,col.order]
barplot(prop5,main = "Region l2 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop5),las=2, col = as.character(region.l2.cols[rownames(prop5)]))
batch.entropy<-table(KB.EC$patient, KB.EC$ec.v2.clusters)
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Experiment Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(KB.EC$ec.v2.clusters)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)
dev.off()

