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
cycDT <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_CycDT_Cells.Rds")

KB$v2.subclass.l3 <- KB$v2.subclass.l3 %>% as.character()
KB$v2.subclass.l3[cycDT] <- "cycDT"

table(KB$v2.subclass.l3)
KB <- subset(KB, v2.subclass.l3 %in% c("M-TAL","dM-TAL","C/M-TAL","C-TAL","dC/M-TAL","dC-TAL",
                                       "aTAL","DCT1","dDCT","aDCT","DCT2","CNT","CNT-PC","PC",
                                       "IMCD","PapE","IC-A","IC-B","cycDT"))

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
#121607

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
KB <- AddMetaData(KB, metadata = k200infomap, col.name = "rpca_k200_infomap_dt")
DimPlot(KB, group.by = "rpca_k200_infomap_dt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "v2.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Pagoda2-TAL-CD-Meta.Rds")

#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k200_infomap_dt"
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
Idents(KB) <- "rpca_k200_infomap_dt"
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
mult <- c(17,24,45,46)

##Remove multiplets and re-run integration
KB <- subset(KB, rpca_k200_infomap_dt %in% mult, invert = TRUE)
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
#119763

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
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_dt")
KB <- AddMetaData(KB, metadata = k200infomap, col.name = "rpca_k200_infomap_dt")
KB <- AddMetaData(KB, metadata = k500infomap, col.name = "rpca_k500_infomap_dt")
DimPlot(KB, group.by = "rpca_k100_infomap_dt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k200_infomap_dt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_dt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "v2.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Pagoda2-TAL-CD-Meta_2.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Pagoda2-TAL-CD_K500.Rds")
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset.Rds")


###Annotate Clusters
#Identify max predicted subclass per cluster
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset.Rds")

Idents(KB) <- "rpca_k200_infomap_dt"
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
Idents(KB) <- "rpca_k500_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

table(KB$rpca_k200_infomap_dt, KB$rpca_k500_infomap_dt)
table(KB$rpca_k200_infomap_dt, KB$rpca_k500_infomap_dt)[20:50,]
table(KB$rpca_k200_infomap_dt, KB$rpca_k500_infomap_dt)[41:68,]

###Annotations
reorder <- c(
  2, #M-TAL; 1 distinct DEG
  22, #M-TAL; no distinct DEGs (merge to cluster 2)
  10, #M-TAL; 2 distinct DEGs (merge to cluster 2)
  27, #M-TAL; no distinct DEGs (merge to cluster 2)
  39, #M-TAL; no distinct DEGs (merge to cluster 2)
  28, #M-TAL
  3, #dM-TAL; High %MT & High % ERT; Injury markers
  15, #dM-TAL; High %MT & High % ERT; Injury markers
  58, #dM-TAL; very low genes
  
  12, #C/M-TAL-A; 2 distinct DEGs
  33, #C-TAL-A; 
  40, #MD
  32, #C/M-TAL-B; 1 distinct DEG
  13, #C-TAL-B
  44, #dC-TAL; High %MT & High % ERT; Injury markers
  52, #dC-TAL; High %MT; Injury markers
  43, #aTAL1
  45, #aTAL
  18, #aTAL
  51, #aTAL; Injury Markers; High %MT & High % ERT
  38, #ifgTAL; interferon gamma response
  57, #cycDT
  
  25, #DCT1
  5, #DCT1; no distinct DEGs (merge to cluster 25)
  11, #DCT1
  21, #DCT1; no distinct DEGs (merge to cluster 25)
  41, #dDCT1; High %MT & High % ERT; Injury markers
  50, #dDCT1; High %MT & High % ERT; Injury markers
  17, #DCT2
  19, #DCT2
  47, #dDCT2; Injury markers
  31, #aDCT
  46, #aDCT
  
  14, #CNT
  6, #CNT-PC
  36, #dCNT-PC; High %MT & High % ERT; Injury markers
  49, #dCNT-PC; Injury markers
  
  
  4, #PC
  7, #PC
  8, #PC; no distinct DEGs (merge to cluster 7 based on k500)
  23, #dPC; High %MT & High % ERT; Injury markers
  35, #dPC; Injury markers
  16, #IMCD
  42, #IMCD
  56, #IMCD
  
  30, #PapE
  26, #PapE, TL multiplet?
  37, #PapE
  48, #PapE
  
  
  1, #CCD-IC-A
  24, #OMCD-IC-A
  20, #OMCD-IC-A
  34, #OMCD-IC-A; no distinct DEGs (merge to cluster 20 based on k500)
  
  9, #IC-B
  29, #IC-B

  53,54,55 #PT multiplets
)

DimPlot(KB, group.by = "rpca_k200_infomap_dt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_dt", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k100_infomap_dt", reduction = "umap", label = TRUE) + NoLegend()

Idents(KB) <- "rpca_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = reorder)
table(Idents(KB))

##Check markers, injury features and meta and add details to reordered list above
dt.markers <- c("Slc12a1", "Egf","Umod",                              #TAL
                "Enpp2","Ano4","Kcnt1","Prox1","Ank2","Clcnka",       #M-TAL
                
                "Slco3a1",                                            #C-TAL-A
                "Tenm4","Slc5a1","Enox1",                             #C-TAL
                "Thsd4",                                              #C-TAL-A
                "Robo2","Nos1","Ndst4","Auts2","Pappa2",              #MD
                
                
                "Tmem207","Tinag","Cldn16","Tmem52b","Kcnj10",        #C-TAL-B
                
                "Itgb6",                                              #aTAL
                "Nrg1","Samd4","Creb5","Hmga2","Runx1","Igf2bp2",     #aTAL1
                "Itga3","Havcr1","Cd44",#"Map1b",
                "Dysf","Dcdc2a","Asic2","Rapgef4",                    #aTAL2
                "Adamts1",                                            #aTAL2
                "Itgb8","Prom1",                                      #frTAL
                "Ifi44","Iigp1","Ifit1","Gbp3",                       #ifnTAL
                "Top2a","Cdca2","Mki67",                              #cycling
                
                "Slc12a3","Cnnm2","Klhl3","Trpm6",                    #DCT
                "Fgf13",                                              #DCT1
                "Trpv5","Pgam2","Pvalb",                              #DCT2
                "Cdh3","Aff2","Camk2b","Scg5",                        #aDCT
                "Fam129a","Ltc4s","Rab11fip5",
                
                "Slc8a1","Calb1","S100g","Egfem1",                    #CNT
                "Adamts16","Depdc7",
                
                "Scnn1g", "Scnn1b",                                   #CNT-PC
                "Gata3","St6gal1","Aqp2","Fxyd4","Aqp3","Ache",       #PC
                "Mcoln3","Btc","Fanca",                               
                "Il1rapl1","Spock3",
                
                "Aqp4","Slc14a2","Aldh1a3",                           #IMCD
                "Tenm3",
                "Fxyd3", "Krt5","Upk1b",                              #PapE
                "Atp6v1c2","Atp6v0d2","Clnk",                         #IC
                "Slc26a7","Slc4a1",                                   #ICA
                "Hs6st3",                                             #CCD-IC-A
                "Aqp6","Kit",                                         #OMCD-IC-A 
                "Slc4a9", "Slc26a4",                                  #ICB
                "Hmx2","Sh2d4b",                                      
                "Igfbp7","Tpt1","Aldob","Gpx3","Wfdc2","Spp1",        #injury
                "Cst3","Apoe","Fabp1","Lcn2","Clu"
                )

DotPlot(KB, features = unique(str_to_title(dt.markers))) + RotatedAxis()

row.order <- reorder
prop1 <- prop.table(table(KB$sex, KB$rpca_k200_infomap_dt), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k200_infomap_dt), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)
prop1 <- prop.table(table(KB$source, KB$rpca_k200_infomap_dt), margin = 2)[,row.order]
barplot(prop1,main = "Source Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder)


#Find Distinct marker genes
markers <- FindAllMarkers(KB, only.pos = TRUE, max.cells.per.ident = 2000,
                          logfc.threshold = 0.25, min.pct = 0.25)
write.table(markers, file="mouse_IRI/Integrated_Mouse_TAL-CD_Cluster_Markers_11112024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_TAL-CD_Cluster_Markers_11112024.robj")
#load("mouse_IRI/Integrated_Mouse_TAL-CD_Cluster_Markers_11112024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_TAL-CD_Cluster_Markers_11112024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


###Merge clusters and re-assess DEGs
neworder <- c(
  1, #2, #M-TAL; 1 distinct DEG
  1, #22, #M-TAL; no distinct DEGs (merge to cluster 2)
  1, #10, #M-TAL; 2 distinct DEGs (merge to cluster 2)
  1, #27, #M-TAL; no distinct DEGs (merge to cluster 2)
  1, #39, #M-TAL; no distinct DEGs (merge to cluster 2)
  2, #28, #M-TAL
  3, #3, #dM-TAL; High %MT & High % ERT; Injury markers
  4, #15, #dM-TAL; High %MT & High % ERT; Injury markers
  5, #58, #dM-TAL; very low genes
  
  6, #12, #C/M-TAL-A; 2 distinct DEGs
  7, #33, #C-TAL-A; 
  8, #40, #MD
  9, #32, #C/M-TAL-B; 1 distinct DEG
  10, #13, #C-TAL-B
  11, #44, #dC-TAL; High %MT & High % ERT; Injury markers
  12, #52, #dC-TAL; High %MT; Injury markers
  13, #43, #aTAL1
  14, #45, #aTAL
  15, #18, #aTAL
  16, #51, #dTAL; Injury Markers; High %MT & High % ERT
  17, #38, #ifnTAL; interferon gamma response
  18, #57, #cycTAL
  
  19, #25, #DCT1
  19, #5, #DCT1; no distinct DEGs (merge to cluster 25)
  20, #11, #DCT1
  19, #21, #DCT1; no distinct DEGs (merge to cluster 25)
  21, #41, #dDCT1; High %MT & High % ERT; Injury markers
  22, #50, #dDCT1; High %MT & High % ERT; Injury markers
  23, #17, #DCT2
  24, #19, #DCT2
  25, #47, #dDCT2; Injury markers
  26, #31, #aDCT
  27, #46, #aDCT
  
  28, #14, #CNT
  29, #6, #CNT-PC
  30, #36, #dCNT-PC; High %MT & High % ERT; Injury markers
  31, #49, #dCNT-PC; Injury markers
  
  
  32, #4, #PC (Aqp4+, OMCD-PC?)
  33, #7, #PC
  33, #8, #PC; no distinct DEGs (merge to cluster 7 based on k500)
  34, #23, #dPC; High %MT & High % ERT; Injury markers
  35, #35, #dPC; Injury markers
  36, #16, #IMCD
  37, #42, #IMCD
  38, #56, #IMCD
  
  39, #30, #PapE
  40, #26, #PapE, TL multiplet?
  41, #37, #PapE
  42, #48, #PapE
  
  
  43, #1, #OMCD-IC-A
  44, #24, #OMCD-IC-A
  45, #20, #CCD-IC-A
  45, #34, #CCD-IC-A; no distinct DEGs (merge to cluster 20 based on k500)
  
  46, #9, #IC-B
  47, #29, #IC-B
  
  48, #53 #PT multiplets
  48, #54 #PT multiplets
  48 #55 #PT multiplets
)

Idents(KB) <- "rpca_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = reorder)
current.cluster.ids <- reorder
new.cluster.ids <- neworder
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
KB$v2.clusters <- Idents(KB)


DotPlot(KB, features = unique(str_to_title(dt.markers))) + RotatedAxis()

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
write.table(markers, file="mouse_IRI/Integrated_Mouse_TAL-CD_V2-Cluster_Markers_11112024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_TAL-CD_V2-Cluster_Markers_11112024.robj")
#load("mouse_IRI/Integrated_Mouse_TAL-CD_Cluster_Markers_11112024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.12,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_TAL-CD_V2-Cluster_Markers_11112024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


table(KB$v2.clusters,KB$condition_level1)
table(KB$v2.clusters,KB$source)
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_B.Rds")




###Update metadata
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_B.Rds")

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
           "Gerhardt2023.Celltype", "celltype","rpca_k100_infomap_dt", "rpca_k200_infomap_dt", "rpca_k500_infomap_dt",
           "v2.clusters")
KB@meta.data <- KB@meta.data[,order]

###Remove Multiplets
KB <- subset(KB, v2.clusters %in% c(48), invert = TRUE)
table(KB$v2.clusters)
KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:50)
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")




### Re-order clusters
table(KB$v2.clusters,KB$condition_level3)
table(KB$v2.clusters,KB$Gerhardt2023.Celltype)

Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = 1:47)
current.cluster.ids <- c(1:12,16,13:15,17:31,33,32,34:47)
new.cluster.ids <- paste0("D_", c(1:47))
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
Idents(KB) <- factor(Idents(KB), levels = paste0("D_", c(1:47)))
KB$v2.clusters <- Idents(KB)

DotPlot(KB, features = unique(str_to_title(dt.markers))) + RotatedAxis()

DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")
