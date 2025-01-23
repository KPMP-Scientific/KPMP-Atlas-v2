library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(pagoda2)
library(BPCells)
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

###Load Combined Data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_C.Rds")
table(KB$v2.subclass.l3)
KB <- subset(KB, v2.subclass.l3 %in% c("M-FIB","C-FIB","C-FIB","C-FIB","MYOF",
                                       "pvFIB","VSM/P"))

DefaultAssay(KB) <- "RNA"
KB[["sketch"]] <- NULL
KB[["integrated.rpca"]] <- NULL
KB[["integrated.rpca.full"]] <- NULL
KB[["umap.full"]] <- NULL
KB[["pca"]] <- NULL
KB[["umap"]] <- NULL
KB

##Integrate and Cluster
KB <- JoinLayers(KB)
KB[["RNA"]] <- as(KB[["RNA"]], Class = "Assay")
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- ScaleData(KB)
KB <- RunPCA(KB)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)
KB <- FindNeighbors(KB, reduction = "integrated.rpca", dims = 1:50)
KB <- FindClusters(KB, resolution = 2, cluster.name = "rpca_clusters")

KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:50)

DimPlot(
  KB,
  reduction = "umap",
  group.by = "v2.subclass.l3",
  combine = TRUE, label.size = 2
)
DimPlot(
  KB,
  reduction = "umap",
  group.by = "source",
  combine = TRUE, label.size = 2
)


###cluster using Pagoda2 instead
KB[["RNA"]] <- JoinLayers(KB[["RNA"]])

countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)
#13431

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_str")
DimPlot(KB, group.by = "rpca_k100_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Stroma-Meta.Rds")

#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k100_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(object = KB)) <- paste("CL", levels(Idents(object = KB)), sep = "")
celltype <- Idents(object = KB)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v2.subclass.l1))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.hs.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v2.subclass.l1))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v2.subclass.l1)),
                                query = prop.table(table(seurat.obj.sub$predicted.hs.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$predicted.hs.subclass.l3)))
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
    colnames(Idents.called) <- c("subclassl1", "subclassl1.Freq","subclassl1.Total", "subclassl3",
                                 "subclassl3.Freq","subclassl3.Total","cluster")
    rownames(Idents.called) <- paste(Idents.called$cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(KB, celltype)

###Identify ambiguous low quality
#check for overlapping identities
max.idents

#Check for overlapping cell type markers
Idents(KB) <- "rpca_k100_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

epi.markers <- c(
  "Nphs2","Podxl",#"Robo2",                              #POD
  "Cp","Nkain3",                                         #PEC
  "Lrp2",#"Cubn",                                        #PT
  "Hnf4a",  #"Slc34a1",                                  #S1/S2
  "Slc5a12","Prodh2",#"Slc5a2","Slc7a8",                 #S1
  "Slc22a8","Slc13a3","Slco1a1", #"Cyp2e1",              #Female S2
  #"Slc22a6","Nat8",                                     #Male S2
  "Cyp7b1",#"Slc22a28",                                  #S3
  #"Slc7a13", "Slc5a8","Acsm3",                          #Male S3
  "Slc7a12",#"Acsm1",                                    #Female S3
  #"Itgb8","Itgb6",              #Repair
  "Havcr1","Vcam1","Sox9",#"Cdh6",
  #"Pou5f1","Sox2",
  #"Mki67","Top2a",                                      #cycling
  #"Pax2","Pitx2",                                       #TL
  "Slc4a11","Cdh13",                                    #OM DTL2
  #"Podn","Aqp1", "Cdh6",                                #IM DTL2
  #"Tspan8",                                             #DTL1 (JM nephron)
  "Jag1","Id3","Corin",#"Slc14a2",
  "Sptssb","Wnt5a",#"Fosl2",#"Slc12a2",                 #DTL3
  "Prox1","Clcnka","Cldn10",                            #ATL
  "Slc12a1", "Egf","Umod",                              #MTAL
  "Enox1",                                              #CTAL
  #"Nos1", "Pappa2",                             #MD
  #"Prom1",#"Spp1","Nrp1","Plscr1",                     #Repair
  "Dcdc2a",                                             #Repair
  "Slc12a3",                                            #DCT
  "Trpm6",                                              #DCT1
  #"Pgam2","Pvalb",                                      #DCT2
  "Slc8a1","Calb1","Egfem1",#"S100g",                    #CNT
  "Scnn1g", "Scnn1b",                                   #CDT-PC
  "Fxyd4","Aqp2",#"Aqp3","Ache", #"Gata3","St6gal1",      #PC
  "Ctnnd2","Prickle2",                 #Repair
  #"Mcoln3","Btc",                                       #OMCD
  
  "Aqp4","Slc14a2","Aldh1a3",                           #IMCD
  "Upk1b","Fxyd3", "Krt5",                              #PapE
  "Atp6v1c2","Atp6v0d2",                                #IC
  "Slc4a1","Slc26a7",                                   #ICA
  #"Aqp6","Kit",  
  "Slc4a9", "Slc26a4",                                  #ICB
  "Igfbp5",#"Rhbg",                                      #nonA nonB
  "Sh2d4b"#,"Hmx2"                                      #CNT-CCD-IC-B
)

int.markers <- c("Pecam1", "Ptprb", #"Meis2", "Flt1",                   #EC
                 "Emcn","Kdr","Plat","Hecw2","Ehd3",                    #EC-GC
                 "Plvap",#"Slco3a1",                                    #EC-PTC
                 "Tm4sf1","Vegfc","Sox17","Aqp1",                       #EC-AEA/DVR
                 #"Fbln5",                                              #EC-AEA
                 "Slc14a1",                                             #EC-DVR
                 "Tll1",#"Nr2f2",                                       #EC-AVR
                 #"Aplnr",                                              #Progenitor
                 "Mmrn1", "Prox1","Tbx1",                               #EC-Lym  
                 
                 "Igfbp5","Tnc",#"Kcnk2","Adamtsl1",                    #M-Fib
                 "Cfh","Pdgfra","Lama2","C7",                           #Fib
                 "Osmr", "Cxcl10", "Relb","Ccl2", "Ccl19",              #Inf Fib
                 "Sparc","Col1a2","Col1a1","Col15a1","Sulf1",           #MYOF                     #aStr
                 
                 "Dcn","Meg3",                                          #pvFIB
                 "Flrt2","Igf1","C3","Pi16","Rspo3",   
                 
                 "Piezo2","Gata3","Itga8",                              #MC
                 "Ren1",                                                #Ren 
                 "Rgs6","Notch3",
                 #"Pdgfrb",
                 "Myh11","Rgs5","Acta2",#"Mcam",                        #VSM/P
                 
                 "Ptprc",
                 "Bank1","Cd79a","Ms4a1",                               #B cells
                 "Igha","Jchain","Xbp1",#"Igkc",                        #Plasma Cells
                 "Cd247","Il7r","Camk4",#"Thy1","Cd3e","Cd96","Cd4",    #T cells
                 "Runx3","Nkg7","Ccl5",                                 #NK
                 "Adgre1",                                              #Mon/Mac
                 "Mrc1", "Stab1",                                       #MAC-M2 
                 "C1qc","Cd14","Fcgr4",                                 #Mac - MDC, Fcgr4 = FCGR3A(CD16)
                 "Adgre4","Ace","Ly6c2","Chil3",#"Pglyrp1", #"Tcf7l2",  #MON
                 
                 "Lgals3","Gpnmb","Lipa","Trem2",                      #moFAM
                 "Clec10a","Cd209a",                                   #cDC2
                 "Flt3", "Zbtb46", "Clec9a",                           #cDC1
                 "Itgae","Xcr1",
                 #"Ms4a2", "Cpa3", "Kit",                              #Mast
                 #"S100a9", "S100a8"#,"Ifitm2",                          #N
                 
                 #"Cdh19", "Nrxn1"                                     #SC/Neu
                 #"Cd36","Plin1"                                        #Adipocytes
                 "Mki67","Top2a"    #cycling
)

DotPlot(KB, features = unique(epi.markers)) + RotatedAxis()
DotPlot(KB, features = unique(int.markers)) + RotatedAxis()

###Ambiguous (mixed identity) clusters, plot markers to confirm multiplets
mult <- c(8,12,16,22,24,26)

##Remove multiplets and re-run integration
KB <- subset(KB, rpca_k100_infomap_str %in% mult, invert = TRUE)
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- ScaleData(KB)
KB <- RunPCA(KB)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)
KB <- FindNeighbors(KB, reduction = "integrated.rpca", dims = 1:50)
KB <- FindClusters(KB, resolution = 2, cluster.name = "rpca_clusters")

KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:50)

DimPlot(
  KB,
  reduction = "umap",
  group.by = "v2.subclass.l3",
  combine = TRUE, label.size = 2
)
DimPlot(
  KB,
  reduction = "umap",
  group.by = "source",
  combine = TRUE, label.size = 2
)
DimPlot(
  KB,
  reduction = "umap",
  group.by = "rpca_k100_infomap_str",
  combine = TRUE, label.size = 2
)


###cluster using Pagoda2 instead
KB[["RNA"]] <- JoinLayers(KB[["RNA"]])

countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)
#12244

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 25, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap
k25infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_str")
KB <- AddMetaData(KB, metadata = k50infomap, col.name = "rpca_k50_infomap_str")
KB <- AddMetaData(KB, metadata = k25infomap, col.name = "rpca_k25_infomap_str")
KB <- AddMetaData(KB, metadata = k500infomap, col.name = "rpca_k500_infomap_str")
DimPlot(KB, group.by = "rpca_k100_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k25_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Stroma-Meta_2.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Stroma_K50.Rds")
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset.Rds")


##Remove multiplets and re-run integration
Idents(KB) <- "rpca_k50_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
DotPlot(KB, features = unique(epi.markers)) + RotatedAxis()
DotPlot(KB, features = unique(int.markers)) + RotatedAxis()

mult <- c(18, 6, 30, 23)

KB <- subset(KB, rpca_k50_infomap_str %in% mult, invert = TRUE)
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- ScaleData(KB)
KB <- RunPCA(KB)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)

KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:50)

DimPlot(
  KB,
  reduction = "umap",
  group.by = "v2.subclass.l3",
  combine = TRUE, label.size = 2
)

###cluster using Pagoda2
KB[["RNA"]] <- JoinLayers(KB[["RNA"]])

countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)
#11199

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 25, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap
k25infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_str")
KB <- AddMetaData(KB, metadata = k50infomap, col.name = "rpca_k50_infomap_str")
KB <- AddMetaData(KB, metadata = k25infomap, col.name = "rpca_k25_infomap_str")
KB <- AddMetaData(KB, metadata = k500infomap, col.name = "rpca_k500_infomap_str")
DimPlot(KB, group.by = "rpca_k100_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k25_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Stroma-Meta_2.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Stroma_K50.Rds")
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset.Rds")





###Annotate Clusters
#Identify max predicted subclass per cluster
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset.Rds")

Idents(KB) <- "rpca_k50_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DotPlot(KB, features = unique(epi.markers)) + RotatedAxis()
DotPlot(KB, features = unique(int.markers)) + RotatedAxis()

table(KB$rpca_k50_infomap_str, KB$rpca_k100_infomap_str)
table(KB$rpca_k50_infomap_str, KB$rpca_k500_infomap_str)


#add k25 clusters to k50
cl27 <- names(KB$rpca_k25_infomap_str[KB$rpca_k25_infomap_str %in% "36"])
cl28 <- names(KB$rpca_k25_infomap_str[KB$rpca_k25_infomap_str %in% "32"])
cl29 <- names(KB$rpca_k25_infomap_str[KB$rpca_k25_infomap_str %in% "14"])
KB@meta.data$rpca_k50_infomap_str <- as.character(KB@meta.data$rpca_k50_infomap_str)
KB@meta.data[cl27,]$rpca_k50_infomap_str <- "27"
KB@meta.data[cl28,]$rpca_k50_infomap_str <- "28"
KB@meta.data[cl29,]$rpca_k50_infomap_str <- "29"


Idents(KB) <- "rpca_k50_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(object = KB)) <- paste("CL", levels(Idents(object = KB)), sep = "")
celltype <- Idents(object = KB)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v2.subclass.l1))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.hs.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v2.subclass.l1))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v2.subclass.l1)),
                                query = prop.table(table(seurat.obj.sub$predicted.hs.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$predicted.hs.subclass.l3)))
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
    colnames(Idents.called) <- c("subclassl1", "subclassl1.Freq","subclassl1.Total", "subclassl3",
                                 "subclassl3.Freq","subclassl3.Total","cluster")
    rownames(Idents.called) <- paste(Idents.called$cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(KB, celltype)


###Annotate based on human predicted identities, published annotations, marker expression
#check for overlapping identities
max.idents

#Check for overlapping cell type markers
Idents(KB) <- "rpca_k50_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

table(KB$rpca_k50_infomap_str, KB$rpca_k100_infomap_str)
table(KB$rpca_k50_infomap_str, KB$rpca_k500_infomap_str)

###Annotations
reorder <- c(
  22, #IM-FIB
  12, #OM-FIB
  26, #OM-FIB
  10, #C/M-FIB? (also possible pvFIB?)
  3, #C-FIB; no distinct markers (merges with 2 at k500)
  2, #C-FIB; 
  1, #C-FIB; no distinct markers
  17, #C-FIB; no distinct markers (merges with 2 at k500)
  8, #C-FIB-PATH; no distinct DEGs (merges with 1 at k500)
  16, #C-FIB; no distinct DEGs, unstable at k500
  15, #C-FIB-OSMRhi
  20, #C-FIB-OSMRhi?
  23, #C-MYOF
  5, #C-MYOF
  25, #cycMYOF
  
  29, #pvFIB-Rspo3+
  13, #pvFIB-Pi16+
  21, #pvFIB; 1 distinct DEG
  14, #pvFIB; no distinct DEGs (merges with cl 19 at k100)
  19, #pvFIB; 1 distinct DEG
  6, #pvFIB; 2 distinct DEGs
  28, #pvMYOF
  
  
  7, #MC
  11, #REN
  4, #VSMC
  24, #VSMC/P
  18, #VSMC/P
  9, #Ad
  27 #SC/Neu 
  
)

DimPlot(KB, group.by = "rpca_k100_infomap_str", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_str", reduction = "umap", label = TRUE) + NoLegend()

Idents(KB) <- "rpca_k50_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = reorder)

##Check markers, injury features and meta and add details to reordered list above
str.markers <- c("Cfh","Pdgfra","Lama2",                           #Fib
                 "Igfbp5","Tnc",                                   #M-Fib
                 "Vcan","Car8","Kcnk2",                            #IM-FIB
                 "Ranbp3l","Cck",                                  #IM/OM-FIB
                 "Adamtsl1","Frem1","Negr1","Cxcl12",              #OM & C/M-FIB
                 "Npy","Mt2","Npas3",                              #OM-FIB
                 
                 "C7",                                             #Pan cortical FIB
                 "Smoc2",
                 "Ccn2",                                 
                 "Nrxn3","Aoah","Clca3a1",
                 "Itga5","Pakap","Adamts6","Spp1","Vim","Timp1","Relb",  #Activated FIB
                 #"Tnc","Ccl2",
                 
                 "Osmr",
                 "Selenop","Cxcl14","Prelp",                      #C-FIB-OSMRlo
                 "Iigp1","Oasl2","Parp14","Stat1","Cxcl10","Ccl2", #C-FIB-OSMRhi
                 "Prr16","Dok6","Tfap2a","Fndc1","Kif26b","Inhba",
                 "Sparc","Col1a2","Col1a1","Col15a1","Sulf1",           #MYOF                     #aStr
                 "Unc13c","Ryr2","Ryr3","Wdr17","Gpm6a",  
                 
                 
                 "Top2a","Mki67",    #cycling
                 
                 "Dcn","Flrt2","Meg3",                             #Pan pvFIB
                 "Fap","Igfbp6","Rspo3","Gli2","Rspo1",            #pvFIB-RSPO3+
                 "Postn","Pdzrn4","Igf1",
                 "C3","Pi16","Ebf2","Sfrp2","Cd34",                #pvFIB-PI16+
                 "Clec3b","Scara5","Sfrp4",
                 "Tmeff2","Itgbl1","Epha3","Plxna4","Hpgd",        #pvFIB
                 "Adgrb3",#"Piezo2",
                 "Myh11","Acta2",'Kcnma1',"Pcdh7",                #pvMYOF
                 "Myocd","Synpo2",                     #pvMYOF
                 "Dpp6","Pgm5","Map1b","Myocd","Col23a1",
                 
                 "Piezo2","Gata3","Kcnq5","Coro2a",#"Mbnl3","Slco2a1", #MC
                 
                 "Ren1","Rgs6","Ptp4a3","Akr1b7","Ephb1",#"Ccdc3",       #REN
                 
                 "Notch3",
                 "Myh11","Rgs5","Acta2","Mcam",                        #VSM/P
                 "Col19a1","Grip2","Ldb3","Prkg2",#"Entpd1",
                 "Rergl","Tenm2",#"Dgkg","Slc38a11","Map3k7cl",
                 "Parp8","Sgip1","Ntm","Grip1",
                 "Cd36","Plin1","Lpl","Adipoq",#"Fabp4",                     #Ad
                 "Nrxn1", "Grik2", "Cdh19"          #SC/NEU
                 
                 )

DotPlot(KB, features = unique(str_to_title(str.markers))) + RotatedAxis()

row.order <- reorder
prop1 <- prop.table(table(KB$sex, KB$rpca_k50_infomap_str), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k50_infomap_str), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder)


#Find Distinct marker genes
markers <- FindAllMarkers(KB, only.pos = TRUE, max.cells.per.ident = 2000,
                          logfc.threshold = 0.25, min.pct = 0.25)
write.table(markers, file="mouse_IRI/Integrated_Mouse_Stroma_Cluster_Markers_11062024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_Stroma_Cluster_Markers_11062024.robj")
#load("mouse_IRI/Integrated_Mouse_Stroma_Cluster_Markers_11062024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_Stroma_Cluster_Markers_11062024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


###Merge clusters and re-assess DEGs
neworder <- c(
  1, #22, #IM-FIB
  2, #12, #OM-FIB
  3, #26, #OM-FIB
  4, #10, #C/M-FIB?; 1 distinct DEG
  5, #3, #C-FIB; no distinct markers (merges with 2 at k500, but mixed at k100)
  6, #2, #C-FIB; 
  7, #1, #C-FIB; no distinct markers
  6, #17, #C-FIB; no distinct markers (merges with 2 at k500)
  7, #8, #C-FIB-PATH; no distinct DEGs (merges with 1 at k500)
  8, #16, #C-FIB; no distinct DEGs, unstable at k500
  9, #15, #C-FIB-Inf 
  10, #20, #C-FIB-Inf
  11, #23, #C-MYOF
  12, #5, #C-MYOF
  13, #25, #cycFIB
  
  14, #29, #pvFIB-Rspo3+
  15, #13, #pvFIB-Pi16+
  16, #21, #pvFIB; 1 distinct DEG
  17, #14, #pvFIB; no distinct DEGs (merges with cl 19 at k100)
  17, #19, #pvFIB; 1 distinct DEG
  18, #6, #pvFIB; 2 distinct DEGs
  19, #28, #pvMYOF
  
  
  20, #7, #MC
  21, #11, #REN
  22, #4, #VSMC
  23, #24, #VSMC/P
  24, #18, #VSMC/P
  25, #9, #Ad
  26  #27 #SC/Neu 
  
)

Idents(KB) <- "rpca_k50_infomap_str"
Idents(KB) <- factor(Idents(KB), levels = reorder)
current.cluster.ids <- reorder
new.cluster.ids <- neworder
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
KB$v2.clusters <- Idents(KB)


DotPlot(KB, features = unique(str_to_title(str.markers))) + RotatedAxis()

row.order <- reorder
prop1 <- prop.table(table(KB$condition_level1, KB$v2.clusters), margin = 2)
barplot(prop1,main = "Condition Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$Group, KB$v2.clusters), margin = 2)
barplot(prop1,main = "Group", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder)

#Find Distinct marker genes
markers <- FindAllMarkers(KB, only.pos = TRUE, max.cells.per.ident = 2000,
                          logfc.threshold = 0.25, min.pct = 0.25)
write.table(markers, file="mouse_IRI/Integrated_Mouse_Stroma_V2-Cluster_Markers_11052024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_Stroma_V2-Cluster_Markers_11052024.robj")
#load("mouse_IRI/Integrated_Mouse_Stroma_Cluster_Markers_11052024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.12,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_Stroma_V2-Cluster_Markers_11052024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


table(KB$v2.clusters,KB$condition_level1)
table(KB$v2.clusters,KB$source)
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_B.Rds")


###Update metadata
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_B.Rds")

KB@meta.data[KB@meta.data$source == "Kirita2020",]$library <- KB@meta.data[KB@meta.data$source == "Kirita2020",]$orig.ident

#Add in experiment metadata
meta <- KB@meta.data
exp.meta <- read.delim("mouse_IRI/Mouse_Experiment_Metadata_11122024.txt")
emc <- c("library","source","assay","experiment","patient","source_ID",       
         "specimen","injury_condition","condition_level3","condition_level2","condition_level1","condition",      
         "age_months","age_group","sex","genotype","strain","protocol","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB@meta.data <- meta
order <- c("library","nCount_RNA", "nFeature_RNA","percent.er","percent.mt","doubletdetection",
           "source","assay","experiment","patient",       
           "specimen","injury_condition","condition_level3","condition_level2","condition_level1","condition",      
           "age_months","age_group","sex","genotype","strain","protocol","tissue_type",
           "predicted.hs.subclass.l1.score", "predicted.hs.subclass.l1", "predicted.hs.subclass.l2.score",
           "predicted.hs.subclass.l2", "predicted.hs.subclass.l3.score", "predicted.hs.subclass.l3",
           "Gerhardt2023.Celltype", "celltype","rpca_k100_infomap_str", "rpca_k50_infomap_str", "rpca_k500_infomap_str",
           "rpca_k25_infomap_str","v2.clusters")
KB@meta.data <- KB@meta.data[,order]

DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")


###Re-order clusters
table(KB$v2.clusters,KB$condition_level3)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = 1:26)
current.cluster.ids <- c(1:4,6,7,8,10,5,9,11:26)
new.cluster.ids <- c(paste0("S_", c(1:25)),"N_1")
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
Idents(KB) <- factor(Idents(KB), levels = c(paste0("S_", c(1:25)),"N_1"))
KB$v2.clusters <- Idents(KB)

DotPlot(KB, features = unique(str_to_title(str.markers))) + RotatedAxis()

DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
