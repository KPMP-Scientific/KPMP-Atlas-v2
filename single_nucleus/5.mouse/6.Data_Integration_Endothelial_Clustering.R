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
KB <- subset(KB, v2.subclass.l3 %in% c("EC","cycEC"))

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
#24177

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
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_ec")
DimPlot(KB, group.by = "rpca_k100_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Pagoda2-Endothelial-Meta.Rds")

#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k100_infomap_ec"
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
Idents(KB) <- "rpca_k100_infomap_ec"
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
mult <- c(12,16,25,26,29,30)

##Remove multiplets and re-run integration
KB <- subset(KB, rpca_k100_infomap_ec %in% mult, invert = TRUE)
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
#22191

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
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_ec")
KB <- AddMetaData(KB, metadata = k50infomap, col.name = "rpca_k50_infomap_ec")
KB <- AddMetaData(KB, metadata = k25infomap, col.name = "rpca_k25_infomap_ec")
KB <- AddMetaData(KB, metadata = k500infomap, col.name = "rpca_k500_infomap_ec")
DimPlot(KB, group.by = "rpca_k100_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k25_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Endothelial-Meta_2.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Endothelial_K500.Rds")
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset.Rds")


##Remove multiplets and re-run integration
Idents(KB) <- "rpca_k50_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
DotPlot(KB, features = unique(epi.markers)) + RotatedAxis()
DotPlot(KB, features = unique(int.markers)) + RotatedAxis()

mult <- c(20,22)

KB <- subset(KB, rpca_k50_infomap_ec %in% mult, invert = TRUE)
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
#20610

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
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_ec")
KB <- AddMetaData(KB, metadata = k50infomap, col.name = "rpca_k50_infomap_ec")
KB <- AddMetaData(KB, metadata = k25infomap, col.name = "rpca_k25_infomap_ec")
KB <- AddMetaData(KB, metadata = k500infomap, col.name = "rpca_k500_infomap_ec")
DimPlot(KB, group.by = "rpca_k100_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k25_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Endothelial-Meta_2.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Endothelial_K500.Rds")
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset.Rds")





###Annotate Clusters
#Identify max predicted subclass per cluster
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset.Rds")

Idents(KB) <- "rpca_k100_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DotPlot(KB, features = unique(epi.markers)) + RotatedAxis()
DotPlot(KB, features = unique(int.markers)) + RotatedAxis()

table(KB$rpca_k100_infomap_ec, KB$rpca_k500_infomap_ec)

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
Idents(KB) <- "rpca_k50_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

table(KB$rpca_k100_infomap_ec, KB$rpca_k500_infomap_ec)

###Annotations
reorder <- c(
  17, #EC-GC
  22, #EC-AA
  3, #EC-DVR
  8, #EC-EA
  6, #M-EC-PTC
  1, #EC-PTC
  16, #EC-PTC
  13, #EC-PTC
  4, #EC-PTC
  9, #EC-PTC
  21, #EC-PTC
  23, #EC-PTC
  10, #EC-PTC
  11, #EC-V
  7, #EC-PCV
  5, ##EC-AVR
  12, #EC-AVR
  2, #EC-AVR
  15, #EC-AVR
  14, #EC-AVR
  24, #EC-AVR (papilla)
  18, #cycEC
  20, #EC-LYM
  19  #EC-PTC
  
)

DimPlot(KB, group.by = "rpca_k100_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k25_infomap_ec", reduction = "umap", label = TRUE) + NoLegend()

Idents(KB) <- "rpca_k100_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = reorder)

##Check markers, injury features and meta and add details to reordered list above


###vascular markers - Dumas et al., 2021. Nature Reviews | Nephrology
ec.markers <- c("PLAT",                                               #gRECs
                "TSPAN7", "GATA5", "TBX3",
                "KDR", "EHD3", "SMAD6",                               #Capillaries
                "VEGFA", 
                "GJA5", "CXCR4",                                      #AA
                "FBLN5","SOX17", "PI16",                              #Arterioles
                "EDN1", "GJA4", "CLDN5",                              #Distant afferent arterioles
                "SLC6A6"                                              #Distal efferent arterioles
)

DotPlot(KB, features = unique(str_to_title(ec.markers))) + RotatedAxis()

#cortical EC
ec.markers <- c("IGFBP3","NPR3",                                      #cRECs
                "SOX17", "SEMA3G", "FBLN5",                           #Arteries
                "ELN", "LTBP4","EDN1", "GJA5",                        #Large arteries
                "SLC8A1", "CLDN5", "JAG1", 
                "EDN1", "GJA4", "CLDN5",                              #Afferent arterioles
                "KCNN4", "S1PR1", "CXCL12",
                "KLF4", "SLC6A6",                                     #Efferent arterioles
                "NR2F2", "CD9", "GAS6", "PLVAP",                      #Veins
                "KDR", "NR2F2", "TNXB", "JUP",                        #Postcapillary venules
                "KDR", "FLT1", "NRP1", "PLVAP", "INSR",               #Peritubular capillaries
                "APOE", "PLPP3", "THRSP",                            #Type I cappillaries (or Type II - low Apoe)
                "ISG15", "IFIT1", "IFIT3", "IFI16",                   #Interferon activated
                "IFIT2", "IRF7", 
                "ESM1", "COL4A1", "COL4A2", "APLNR", "APLN",          #Angiogenic
                "TP53I11", "PLK2", "FSCN1",
                "LYVE1", "PROX1", "FLT4", "PDPN"                      #Lymphatics
                
                
)

DotPlot(KB, features = unique(str_to_title(ec.markers))) + RotatedAxis()



#medullary EC
ec.markers <- c("IGFBP7","CD36",                                         #mRECs
                
                "SOX17", "GJA4", "FBLN5",                                #Arteries
                
                "SLC14A1", "AQP1", "SCIN", "CXCL12", "HPGD", "EDN1",     #Descending vasa recta
                "ADIPOR2", "CLDN5", "SLC5A3", 
                
                "S100A4", "S100A6", "AKR1B1",                            #Descending vasa recta papilla
                
                "NR2F2", "PLVAP",                                        #Veins
                
                "TEK", "GAS6",                                           #Ascending vasa recta
                
                "CRYAB", "FXYD2", "CD9", "S100A6", "AKR1B1",             #Ascending vasa recta papilla
                "GAPDH", "LDHA",
                
                "ISG15", "IFIT1", "IFIT3", "IFI16",                      #Ascending vasa recta Interferon activated
                "IFIT2", "IRF7", 
                
                "KDR", "NR2F2", "TNXB", "JUP",                           #Postcapillary venules
                
                "KDR", "FLT1", "NRP1", "PLVAP",                          #Capillaries
                
                "PLPP3", "CD36",                                         #Medullary capillaries
                
                "ISG15", "IFIT1", "IFIT3", "IFI16",                      #Interferon activated
                "IFIT2", "IRF7", 
                
                "ESM1", "COL4A1", "COL4A2", "APLNR", "APLN",             #Angiogenic
                "TP53I11", "PLK2", "FSCN1"
                
                
)

DotPlot(KB, features = unique(str_to_title(ec.markers))) + RotatedAxis()


ec.markers <- c("Pecam1","Ptprb","Flt1",                                #Broad EC
                "Emcn","Hecw2","Plat","Kdr",                            #EC-GC
                "Sema5a","Ehd3","Tbx3",#"Tspan7","Epha6",
                
                "Pi16","Smad6",                                         #EC-GC/EC-AA
                "Efna5","Ptprj","Eln","Pde3a","Ltbp4","Sulf1",          #EC-AA
                #"Cldn5","Jag1","Slc8a1","Gja5","Edn1","Edn1","Gja4",                                  #Distant afferent arterioles
                  
                "Fbln5","Sema3g","Sox17","Adipor2",                     #Arteries
                
                "Lef1","Slc14a1","Prom1","Aqp1","Scin",                 #EC-DVR  
                #"Hpgd","Slc5a3",
                
                "Slc6a6","Atp13a3","Zbtb7c","Fut8","Klf4",              #EC-EA
                
                "Kdr", "Nrp1","Plpp3","Plvap","Insr",                   #EC-PTC
                "Plpp3", "Cd36","Col15a1","Fabp5",                      #M-EC-PTC
                
                "Npr3","Esm1",                                          #EC-PTC
                
                "Kit","Nav3","Adam19","Plk2","Usp31","Aplnr",           #angEC-PTC
                "Apln","Fscn1",
                "Inhbb","Osmr","Gda","Galnt15",                         #angEC-PTC
                "Pakap","Akap12","Spry4",                               #infEC-PTC
                "Gbp4","Ifit1","Ifit2",'Icam1',                         #iaEC-PTC (Interferon activated)
                #"Isg15","Ifit3","Irf7","Cxcl10",
                
                "Nr2f2","Cd9","Gas6","Plvap",                           #Veins
                "Pcdh7","Slco2a1","Adamtsl1","Adgrg6","Nrxn3",          #EC-V
                #"Bnc2","Vwf",
                
                "Kdr", "Nr2f2", "Tnxb","Nhs","Epb41l4a",                #EC-PCV
                
                
                "Tek", "GAS6","Tll1",                                     #EC-AVR
                "Slc24a3","C1qtnf7","Slco1a4","Unc5c",
                "Nlgn1","Igf1","Thsd4","Flrt2","Lepr",
                "Ppfia2","Tox","Cyria","Brip1","Ppip5k2",
                "Tln2","Sec22b","Col4a3",
                
                "Mrps6","Akr1b3","Cryab","Pax2","Ldha",                   #EC-AVR papilla
                
                
                "TOP2A","MKI67","CENPF",                                  #cycling
                
                "Ccl21a","Mmrn1","Prox1","Flt4","Tbx1"#,"Lyve1","Pdpn"              #EC-LYM
                
)


DotPlot(KB, features = unique(str_to_title(ec.markers))) + RotatedAxis()

row.order <- reorder
prop1 <- prop.table(table(KB$sex, KB$rpca_k100_infomap_ec), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k100_infomap_ec), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder)


#Find Distinct marker genes
markers <- FindAllMarkers(KB, only.pos = TRUE, max.cells.per.ident = 2000,
                          logfc.threshold = 0.25, min.pct = 0.25)
write.table(markers, file="mouse_IRI/Integrated_Mouse_Endothelial_Cluster_Markers_11072024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_Endothelial_Cluster_Markers_11072024.robj")
#load("mouse_IRI/Integrated_Mouse_Endothelial_Cluster_Markers_11072024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_Endothelial_Cluster_Markers_11072024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


###Merge clusters and re-assess DEGs
neworder <- c(
  1, #17, #EC-GC
  2, #22, #EC-AA
  3, #3, #EC-DVR
  4, #8, #EC-EA
  5, #6, #M-EC-PTC; No distinct DEGs
  6, #1, #EC-PTC
  6, #16, #EC-PTC; 1 distinct DEG (merges with 1 at k500)
  7, #13, #EC-PTC; no distinct DEGs (merges with 7 at k500)
  8, #4, #EC-PTC; no distinct DEGs (merges with 19 at k500)
  9, #9, #angEC-PTC-Kit+
  10, #21, #angEC-PTC-Spry4+
  11, #23, #EC-PTC-Inhbb+
  12, #10, #iaEC-PTC
  13, #11, #EC-V
  14, #7, #EC-PCV
  15, #5, #EC-AVR, no distinct DEGs (no clear merging at k500)
  16, #12, #EC-AVR
  16, #2, #EC-AVR; 1 distinct DEG (merges with 12 at k500)
  17, #15, #EC-AVR
  18, #14, #EC-AVR; 1 distinct DEG (distinct at k500) (papilla)
  19, #24, #EC-AVR (papilla) possible epi multiplet?
  20, #18, #cycEC
  21, #20,  #EC-LYM
  22  #19 #EC-PTC/Ery multiplet
)

Idents(KB) <- "rpca_k100_infomap_ec"
Idents(KB) <- factor(Idents(KB), levels = reorder)
current.cluster.ids <- reorder
new.cluster.ids <- neworder
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
KB$v2.clusters <- Idents(KB)


DotPlot(KB, features = unique(str_to_title(ec.markers))) + RotatedAxis()

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
write.table(markers, file="mouse_IRI/Integrated_Mouse_Endothelial_V2-Cluster_Markers_11072024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_Endothelial_V2-Cluster_Markers_11072024.robj")
#load("mouse_IRI/Integrated_Mouse_Endothelial_Cluster_Markers_11072024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.12,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_Endothelial_V2-Cluster_Markers_11072024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


table(KB$v2.clusters,KB$condition_level1)
table(KB$v2.clusters,KB$source)
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_B.Rds")


###Update metadata
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_B.Rds")

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
           "Gerhardt2023.Celltype", "celltype","rpca_k100_infomap_ec", "rpca_k50_infomap_ec", "rpca_k500_infomap_ec",
           "rpca_k25_infomap_ec","v2.clusters")
KB@meta.data <- KB@meta.data[,order]

#remove multiplets
KB <- subset(KB, v2.clusters %in% c(22), invert = TRUE)
table(KB$v2.clusters)
KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:50)
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_C.Rds")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_C.Rds")



###Re-order clusters
table(KB$v2.clusters,KB$condition_level3)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = 1:21)
current.cluster.ids <- c(1:9,11,10,12:14,16:19,15,20,21)
new.cluster.ids <- paste0("E_", c(1:21))
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
Idents(KB) <- factor(Idents(KB), levels = paste0("E_", c(1:21)))
KB$v2.clusters <- Idents(KB)

DotPlot(KB, features = unique(str_to_title(ec.markers))) + RotatedAxis()

DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_C.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_C.Rds")
