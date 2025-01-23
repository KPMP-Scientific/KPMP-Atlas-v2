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
KB <- subset(KB, v2.subclass.l3 %in% c("B","PL","T/NK","MAC","FAM","DC","cycMAC"))

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
#13183

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
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_imm")
DimPlot(KB, group.by = "rpca_k100_infomap_imm", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-IMM-Meta.Rds")

#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k100_infomap_imm"
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
Idents(KB) <- "rpca_k100_infomap_imm"
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
mult <- c(17,19,21,23,25,26,28)

##Remove multiplets and re-run integration
KB <- subset(KB, rpca_k100_infomap_imm %in% mult, invert = TRUE)
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
#11873

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
k50infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap_imm")
KB <- AddMetaData(KB, metadata = k50infomap, col.name = "rpca_k50_infomap_imm")
KB <- AddMetaData(KB, metadata = k500infomap, col.name = "rpca_k500_infomap_imm")
DimPlot(KB, group.by = "rpca_k100_infomap_imm", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_imm", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k500_infomap_imm", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Immune-Meta_2.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Immune_K50.Rds")
saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset.Rds")


###Annotate Clusters
#Identify max predicted subclass per cluster
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset.Rds")

Idents(KB) <- "rpca_k100_infomap_imm"
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
Idents(KB) <- "rpca_k50_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DotPlot(KB, features = unique(epi.markers)) + RotatedAxis()
DotPlot(KB, features = unique(int.markers)) + RotatedAxis()

table(KB$rpca_k50_infomap_imm, KB$rpca_k100_infomap_imm)
table(KB$rpca_k50_infomap_imm, KB$rpca_k500_infomap_imm)

###Annotations
reorder <- c(
  20, #B
  28, #B - Age-associated B cells?
  12, #PL
  17, #Naïve Th
  16, #T-REG
  25, #Tcm/Naïve cytotoxic T cells
  3, #Tem/Trm CD8+ Cytotoxic T Cell
  6, #Tem/Temra (Terminally Differentiated) CD8+ Cytotoxic T Cells
  13, #NK
  23, #resMAC
  11, #MAC; 1 distinct marker (merge to 2/5 at K500)
  2, #MAC; no distinct markers (merge to 5 at K500)
  5, #MAC; no distinct markers (merge to 2 at K500)
  
  31, #MAC; no distinct immune markers (merge to 10/14 at K500)
  10, #MAC; no distinct markers (merge to 14 at K500)
  14, #MAC; no distinct markers (merge to 10 at K500)
  19, #MAC; no distinct markers (merge to 10/14 at K500)
  
  1, #moMAC
  9, #moMAC
  30, #infMAC
  29, #infMAC
  32, #infMAC
  
  27, #MAC; no distinct markers (merges only partially with 21 at k500)
  21, #MAC
  7, #moFAM

  18, #MON
  
  
  8, #cDC2
  22, #mDC
  4, #cDC1
  26, #pDC
  24, #cycMAC
  15 #cycDC
  
  
)

DimPlot(KB, group.by = "rpca_k100_infomap_imm", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "rpca_k50_infomap_imm", reduction = "umap", label = TRUE) + NoLegend()

Idents(KB) <- "rpca_k50_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = reorder)


##Check markers, injury features and meta and add details to reordered list above
imm.markers <- c("Ptprc",
                 "Bank1","Cd79a","Ms4a1","Cd37",                    #B cells
                 "Fcrl5","Ccr6","Sox5",                             #Age-associated B cells?
                 "Igha","Jchain","Xbp1",#"Igkc",                    #Plasma Cells
                 
                 "Il7r","Bcl11b","Cd247","Camk4",#"Cd96",           #T cells
                 "Cd3d","Cd3e", "Ccr7",                                    #T cells
                 "Ikzf2","Il2ra","Ctla4","Cd4","Foxp3",             #T-REG
                 "Ms4a4b","Themis",
                 "Lef1","Tcf7","Sell",#"Ccr7",                      #Tcm/Naïve cytotoxic T cells
                 "Cd8a",                                            #CD8+
                 "Gzmk",#"Ccl5",                                    #Tem/Trm CD8+ Cytotoxic T Cell
                 "Runx3","Ccl5","Nkg7","Samd3",                     #CD8+ & NK 
                 "Klrk1","Txk",                                     #Tem/Temra (Terminally Differentiated) CD8+ Cytotoxic T Cells
                 "Klrc2","Klrd1","Klre1","Ncr1","Gzma",             #NK
                 
                 "Ms4a7",
                 "Adgre1","Csf1r",                                   #General MAC
                 "Mrc1","Stab1","F13a1","Cd163",#"Lyve1",            #resMAC
                 "H2-Aa","H2-Eb1",#"H2-Ab1",                         #resMAC-MHChi
                 "C1qb",# "C1qc", "C1qa",                             
                 
                 "Cx3cr1","Cd14",                                    #MAC
                 "Adam22","Negr1",                                   #resMAC-MHChi
                 "Fnbp1l","Col14a1","Nadk","Il1b",#"Ccl3",           #resMAC-MHChi
                 
                 "Apbb2","Lhfpl2","Pmepa1","Cpd","Rapgef5",          #MAC-Cpd+
                 
                 "Ifi213","Ifi204","Pde7b","Mx1","Iigp1",            #moMAC-Ifn+
                 "Cxcl9","Cxcl10","Dnase1l3","Cd274",#"Pstpip2",     #moMAC-Cxcl10+
                 "C3",
                 
                 "Gpnmb","Lgals3","Lipa", "Trem2","Clec4d",#"Capg",  #moFAM
                 #"C3","Cx3cr1",                                     #moMAC-C3+
                 #"TCF7L2","COTL1","Fcgr4",                          #ncMON, Fcgr4 = FCGR3A(CD16)
                 
                 "Ccr2","Fn1","Pltp","Fgfr1","Ccl9",                 #moMAC-Fn1+
                 #"Cd200r1","Tnip3","Fcrls",
                 "Itgam","Ccl2","Adam8","S100a6","Arg1",#"Gda",      #moMAC-Arg1+
                 
                 "Plcb1","Gsr","Vcan","Ace","Chil3",                 #MON
                 #"Adgre4","Pglyrp1",
                 
                 "Dpp4","Flt3", "Zbtb46",                            #DC
                 "Kmo","Clec10a","Cd209a",                           #cDC2
                 "Tbc1d4","Slco5a1","Ccr7",                          #mDC
                 "Wdfy4","Clec9a","Itgae","Xcr1",                    #cDC1
                 
                 "Bcl11a","Siglech","Ly6c2",                         #pDC
                 #"S100a9","Fcgr4","S100a8","Ifitm2",                #N
                 
                 "Top2a","Mki67"                                     #cycling                
)

DotPlot(KB, features = unique(str_to_title(imm.markers))) + RotatedAxis()

row.order <- reorder
prop1 <- prop.table(table(KB$sex, KB$rpca_k50_infomap_imm), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k50_infomap_imm), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder)


table(KB$rpca_k50_infomap_imm, KB$rpca_k100_infomap_imm)

#Find Distinct marker genes
markers <- FindAllMarkers(KB, only.pos = TRUE, max.cells.per.ident = 2000,
                                  logfc.threshold = 0.25, min.pct = 0.25)
write.table(markers, file="mouse_IRI/Integrated_Mouse_Immune_Cluster_Markers_11042024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_Immune_Cluster_Markers_11042024.robj")
#load("mouse_IRI/Integrated_Mouse_Immune_Cluster_Markers_11042024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_Immune_Cluster_Markers_11042024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


###Merge clusters and re-assess DEGs
neworder <- c(
  1, #B
  2, #B
  3, #PL
  4, #Naïve Th
  5, #T-REG
  6, #Tcm/Naïve cytotoxic T cells
  7, #Tem/Trm CD8+ Cytotoxic T Cell
  8, #Tem/Temra (Terminally Differentiated) CD8+ Cytotoxic T Cells
  9, #NK
  10, #resMAC
  11, #MAC; 1 distinct marker (merge to 2/5 at K500)
  11, #MAC; no distinct markers (merge to 5 at K500)
  11, #MAC; no distinct markers (merge to 2 at K500)
  
  12, #MAC; no distinct immune markers (merge to 10 at K500)
  12, #MAC; no distinct markers (merge to 14 at K500)
  12, #MAC; no distinct markers (merge to 10 at K500)
  12, #MAC; no distinct markers
  
  13, #moMAC
  14, #moMAC
  15, #infMAC
  16, #infMAC
  17, #infMAC
  
  18, #MAC; no distinct markers (merges with 21 at k500)
  19, #MAC
  20, #moFAM
  
  21, #MON
  
  
  22, #cDC2
  23, #mDC
  24, #cDC1
  25, #pDC
  26, #cycMAC
  27 #cycDC
  
)

Idents(KB) <- "rpca_k50_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = reorder)
current.cluster.ids <- reorder
new.cluster.ids <- neworder
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
KB$v2.clusters <- Idents(KB)


DotPlot(KB, features = unique(str_to_title(imm.markers))) + RotatedAxis()

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
write.table(markers, file="mouse_IRI/Integrated_Mouse_Immune_V2-Cluster_Markers_11052024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(markers, file="mouse_IRI/Integrated_Mouse_Immune_V2-Cluster_Markers_11052024.robj")
#load("mouse_IRI/Integrated_Mouse_Immune_Cluster_Markers_11052024.robj")

cl.mark <- markers[markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.12,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(10, avg_log2FC) -> top10
cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="mouse_IRI/Integrated_Mouse_Immune_V2-Cluster_Markers_11052024_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB, features = top5$gene) + RotatedAxis()


table(KB$v2.clusters,KB$condition_level1)
table(KB$v2.clusters,KB$source)
table(KB$v2.clusters,KB$source)
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_B.Rds")


###Update metadata
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_B.Rds")

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
           "Gerhardt2023.Celltype", "celltype","rpca_k100_infomap_imm", "rpca_k50_infomap_imm", "rpca_k500_infomap_imm",
           "v2.clusters")
KB@meta.data <- KB@meta.data[,order]

DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_C.Rds")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_C.Rds")


### Re-order clusters
table(KB$v2.clusters,KB$condition_level3)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = 1:27)
current.cluster.ids <- c(1:11,13,19,18,15,16,20,12,14,17,21:27)
new.cluster.ids <- paste0("I_", c(1:27))
Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
Idents(KB) <- factor(Idents(KB), levels = paste0("I_", c(1:27)))
KB$v2.clusters <- Idents(KB)

DotPlot(KB, features = unique(str_to_title(imm.markers))) + RotatedAxis()

DimPlot(KB, group.by = "condition_level3", reduction = "umap", label = TRUE) #+ NoLegend()
DimPlot(KB, group.by = "v2.clusters", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_C.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_C.Rds")
