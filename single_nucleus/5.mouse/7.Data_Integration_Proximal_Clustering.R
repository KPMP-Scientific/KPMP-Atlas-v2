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
KB <- subset(KB, v2.subclass.l3 %in% c("POD","PEC","PT-S1","PT-S1/S2","dPT-S1/S2","PT-S2","dPT-S2",
                                       "PT-S3","dPT-S3","aPT-S1/S2","aPT","aPT-S3/DTL","cycPT",
                                       "DTL2","DTL1","aDTL","DTL3","ATL"))

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
#158120

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k200infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k200infomap, col.name = "rpca_k200_infomap_pt")
DimPlot(KB, group.by = "rpca_k200_infomap_pt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "v2.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Pagoda2-PT-TL-Meta.Rds")

#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k200_infomap_pt"
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
Idents(KB) <- "rpca_k200_infomap_pt"
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
mult <- c(16,31,37,38,52)

cycDT <- WhichCells(KB, idents = 52)
saveRDS(cycDT, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_CycDT_Cells.Rds")

##Remove multiplets and re-run integration
KB <- subset(KB, rpca_k200_infomap_pt %in% mult, invert = TRUE)
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
#155257

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_pt")
KB <- AddMetaData(KB, metadata = k200infomap, col.name = "rpca_k200_infomap_pt")
KB <- AddMetaData(KB, metadata = k500infomap, col.name = "rpca_k500_infomap_pt")
DimPlot(KB, group.by = "rpca_k100_infomap_pt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k200_infomap_pt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_pt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "v2.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-PT-TL-Meta_2.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-PT-TL_K500.Rds")
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset.Rds")


###Annotate Clusters
#Identify max predicted subclass per cluster
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset.Rds")

Idents(KB) <- "rpca_k500_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DotPlot(KB, features = unique(epi.markers)) + RotatedAxis()
DotPlot(KB, features = unique(int.markers)) + RotatedAxis()

levels(Idents(object = KB)) <- paste("CL", levels(Idents(object = KB)), sep = "")
celltype <- Idents(object = KB)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v2.subclass.l3))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.hs.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v2.subclass.l3))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v2.subclass.l3)),
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
    colnames(Idents.called) <- c("subclassl3", "subclassl3.Freq","subclassl3.Total", "hs.subclassl3",
                                 "hs.subclassl3.Freq","hs.subclassl3.Total","cluster")
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
Idents(KB) <- "rpca_k500_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

table(KB$rpca_k200_infomap_pt, KB$rpca_k500_infomap_pt)[1:20,]
table(KB$rpca_k200_infomap_pt, KB$rpca_k500_infomap_pt)[21:40,]
table(KB$rpca_k200_infomap_pt, KB$rpca_k500_infomap_pt)[41:68,]

###Annotations
reorder <- c(
  14, #POD
  
  16, #PEC
  
  4, #PT-S1
  15, #PT-S1; no distinct markers (merge to 4)
  31, #PT-S1; no distinct markers (merge to 4)
  6, #PT-S1; no distinct markers (merge to 4)
  7, #PT-S1; no distinct markers (merge to 4)
  17, #PT-S1; 1 distinct marker (merge to 4)
  
  9, #PT-S1/S2
  30, #PT-S1/S2; no distinct markers (merge to 9)
  23, #PT-S1/S2; no distinct markers (merge to 9)
  
  
  1, #PT-S2
  41, #PT-S2; no distinct markers (merge to 1) 
  5, #PT-S2; 2 distinct markers (merge to 1)
  39, #PT-S2; 2 distinct markers (merge to 1)
  36, #PT-S2; no distinct markers (merge to 1) 
  21, #PT-S2; 2 distinct markers (merge to 1)
  40, #PT-S2; Mostly Female; 2 distinct markers (merge to 1)
  37, #PT-S2; 1 distinct marker (merge to 1)
  3, #dPT-S2; High %MT, High %ERT; Injury Markers
  
  8, #PT-S3
  38, #PT-S3; 2 distinct markers (merge to 8)
  33, #PT-S3
  45, #PT-S3
  32, #PT-S3
  43, #PT-S3; no distinct markers (merge to 8)
  28, #dPT-S3; High %MT, High %ERT; Injury Markers
  35, #dPT-S3; Injury Markers
  
  26, #aPT
  12, #aPT
  42, #aPT
  
  27, #dPT
  25, #dPT-Havcr1+
  22, #dPT-Havcr1+
  29, #dPT-Havcr1+
  34, #dPT-Havcr1+
  
  20, #cycPT
  44, #cycPT
  
  2, #DTL1
  18, #DTL1
  13, #DTL2
  11, #DTL2
  24, #aDTL
  19, #DTL3
  10  #ATL
  
)

DimPlot(KB, group.by = "rpca_k200_infomap_pt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_pt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k100_infomap_pt", reduction = "umap", label = TRUE) + NoLegend()

Idents(KB) <- "rpca_k500_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = reorder)
table(Idents(KB))

##Check markers, injury features and meta and add details to reordered list above
pt.markers <- c("Nphs2","Podxl","Robo2",                              #POD
                "Cp","Tnc","Apba1","Nkain3",#"Kcnmb2",                #PEC
                "Lrp2","Cubn","Hnf4a",                                #PT
                "Slc34a1",                                            #S1/S2
                "Slc5a12","Prodh2","Slc5a2","Slc7a8",                 #S1
                "Slc22a8","Slc13a3","Cyp2e1",                         #S2
                "Slc22a6","Nat8","Slco1a1",                           #S2
                "Slc7a13", "Slc5a8","Acsm3","Cyp7b1",                 #S3
                "Cdh6","Havcr1","Sox4","Sox9",                        #aPT2?
                "Dock10","Kcnip4","Kcnh8","Vcam1","Sorcs1",           #aPT
                "Ccl2","Dlgap1",                                      #aPT
                "Nlrc5","Parp14","Igtp","Iigp1","Irgm2","Gda",        #aPT-S1/S2
                "Scin","Epha7","Prom1",
                
                "Havcr1",
                "Top2a","Mki67",                                      #cycling
                
                "Pax2","Tbc1d4",                                      #TL
                "Slc14a2","Htr4","Corin","Tshr",                      #DTL1 (JM nephron)
                "Jag1","Id3","Tspan8",                                #DTL1 (JM nephron)
                "Cdh13","Slc4a11","Fst","Stk32a","Pitx2",             #OM DTL2
                "Aqp1", "Cdh6","Tiam1","Susd4","Fam78b","Kirrel3",    #IM DTL2
                
                "Lhfp","Ankrd1","Runx1",                              #aDTL
                
                                                              
                "Frmpd4","Ephb2",                                     #DTL3
                "Rbm20","Akr1b3","Sptssb",                            #DTL3/ATL
                "Sgcz","Atp10b","Clcnka", "Cldn10", "Prox1",          #ATL
                "Nrip3","Bsnd",
                
                "Igfbp7","Spp1","Itgb6","Cdh6","Havcr1",
                "Tpt1","Aldob","Gpx3","Wfdc2",        #injury
                "Cst3","Apoe","Clu",
                "Fabp1","Lcn2"
                )

DotPlot(KB, features = unique(str_to_title(pt.markers))) + RotatedAxis()

row.order <- reorder
prop1 <- prop.table(table(KB$sex, KB$rpca_k200_infomap_pt), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k200_infomap_pt), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)
prop1 <- prop.table(table(KB$source, KB$rpca_k200_infomap_pt), margin = 2)[,row.order]
barplot(prop1,main = "Source Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder)


#Find Distinct marker genes
markers <- FindAllMarkers(KB, only.pos = TRUE, max.cells.per.ident = 2000,
                          logfc.threshold = 0.25, min.pct = 0.25)
write.table(markers, file="mouse_IRI/Integrated_Mouse_PT-TL_Cluster_Markers_11082024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_PT-TL_Cluster_Markers_11082024.robj")
#load("mouse_IRI/Integrated_Mouse_PT-TL_Cluster_Markers_11082024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_PT-TL_Cluster_Markers_11082024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


###Merge clusters and re-assess DEGs
neworder <- c(
  1, #14, #POD
  
  2, #16, #PEC
  
  3, #4, #PT-S1
  3, #15, #PT-S1; no distinct markers (merge to 4)
  3, #31, #PT-S1; no distinct markers (merge to 4)
  3, #6, #PT-S1; no distinct markers (merge to 4)
  3, #7, #PT-S1; no distinct markers (merge to 4)
  3, #17, #PT-S1; 1 distinct marker (merge to 4)
  
  4, #9, #PT-S1/S2
  4, #30, #PT-S1/S2; no distinct markers (merge to 9)
  4, #23, #PT-S1/S2; no distinct markers (merge to 9)

  5, #1, #PT-S2
  5, #41, #PT-S2; no distinct markers (merge to 1) 
  5, #5, #PT-S2; 2 distinct markers (merge to 1)
  5, #39, #PT-S2; 2 distinct markers (merge to 1)
  5, #36, #PT-S2; no distinct markers (merge to 1) 
  5, #21, #PT-S2; 2 distinct markers (merge to 1)
  5, #40, #PT-S2; Mostly Female; 2 distinct markers (merge to 1)
  5, #37, #PT-S2; 1 distinct marker (merge to 1)
  6, #3, #dPT-S2; High %MT, High %ERT; Injury Markers
  
  7, #8, #PT-S3
  7, #38, #PT-S3; 2 distinct markers (merge to 8)
  8, #33, #PT-S3
  9, #45, #PT-S3
  10, #32, #PT-S3
  7, #43, #PT-S3; no distinct markers (merge to 8)
  11, #28, #dPT-S3; High %MT, High %ERT; Injury Markers
  12, #35, #dPT-S3; Injury Markers
  
  13, #26, #aPT
  14, #12, #aPT
  15, #42, #aPT
  
  16, #27, #dPT
  17, #25, #dPT-Havcr1+
  18, #22, #dPT-Havcr1+
  19, #29, #dPT-Havcr1+
  20, #34, #dPT-Havcr1+
  
  21, #20, #cycPT
  22, #44, #cycPT
  
  23, #2, #DTL1
  24, #18, #DTL1
  25, #13, #DTL2
  26, #11, #DTL2
  27, #24, #aDTL
  28, #19, #DTL3
  29  #10  #ATL
)

Idents(KB) <- "rpca_k500_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = reorder)
current.cluster.ids <- reorder
new.cluster.ids <- neworder
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
KB$v2.clusters <- Idents(KB)


DotPlot(KB, features = unique(str_to_title(pt.markers))) + RotatedAxis()

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
write.table(markers, file="mouse_IRI/Integrated_Mouse_PT-TL_V2-Cluster_Markers_11082024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_PT-TL_V2-Cluster_Markers_11082024.robj")
#load("mouse_IRI/Integrated_Mouse_PT-TL_Cluster_Markers_11082024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.12,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_PT-TL_V2-Cluster_Markers_11082024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


table(KB$v2.clusters,KB$condition_level1)
table(KB$v2.clusters,KB$source)
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_B.Rds")


###Update metadata
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_B.Rds")

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
           "Gerhardt2023.Celltype", "celltype","rpca_k100_infomap_pt", "rpca_k200_infomap_pt", "rpca_k500_infomap_pt",
           "v2.clusters")
KB@meta.data <- KB@meta.data[,order]

DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "rpca_k100_infomap_pt", reduction = "umap", label = TRUE) #+ NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")


### Re-order clusters
table(KB$v2.clusters,KB$condition_level3)

#pull out aPT subset rpca_k200_infomap_pt cluster 63 and 62... 
Idents(KB) <- "rpca_k200_infomap_pt"
CL62 <- WhichCells(KB, idents = 62)
CL63 <- WhichCells(KB, idents = 63)

Idents(KB) <- "v2.clusters"
KB <- SetIdent(KB, value = 30, cells = CL62)
KB <- SetIdent(KB, value = 31, cells = CL63)


Idents(KB) <- factor(Idents(KB), levels = 1:31)
current.cluster.ids <- c(1:12,30,14,13,31,15:29)
new.cluster.ids <- paste0("P_", c(1:31))
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
Idents(KB) <- factor(Idents(KB), levels = paste0("P_", c(1:31)))
KB$v2.clusters <- Idents(KB)

DotPlot(KB, features = unique(str_to_title(pt.markers))) + RotatedAxis()

DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")
