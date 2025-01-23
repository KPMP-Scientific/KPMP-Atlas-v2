library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)
library(tidyverse)
library(patchwork)
library(pagoda2)
library(BPCells)
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

###Load Combined Data
##Mouse Injury data (multiome)
mmKidAt <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object_filtered_B.Rds")

###Kirita data
load("~/hsKidAt/blake_LTS/public/Kirita2020/Humphreys_2019_mouse_seurat.rda")
counts <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/public/Kirita2020/GSE139107/Kirita_2020_mouse_Counts_BP")
counts <- counts[,colnames(hki)]
meta <- hki@meta.data
hki <- CreateSeuratObject(counts = counts, project = "Kirita2020", min.cells = 3, min.features = 200, 
                          meta.data = meta)
Idents(hki) <- "celltype"
table(Idents(hki))

###Merge objects
table(mmKidAt$source)
hki$source <- "Kirita2020"

KB <- merge(x = mmKidAt, y = c(hki))
KB <- JoinLayers(KB)
KB[["RNA"]] <- split(KB[["RNA"]], f = KB$source)
KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- SketchData(
  object = KB,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
DefaultAssay(KB) <- "sketch"

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
  group.by = c("source", "v2.subclass.l1"),
  combine = TRUE, label.size = 2
)

DimPlot(
  KB,
  reduction = "umap",
  group.by = c("rpca_clusters"),
  combine = FALSE, label.size = 2
)
DimPlot(
  KB,
  reduction = "umap",
  group.by = c("v2.subclass.l1"),
  label = TRUE, label.size = 2
) + NoLegend()
DimPlot(
  KB,
  reduction = "umap",
  group.by = c("celltype"),
  label = TRUE, label.size = 2
) + NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024.Rds")


##Join Layers and Annotate
KB[["sketch"]] <- JoinLayers(KB[["sketch"]])
KB[["sketch"]] <- as(object = KB[["sketch"]], Class = "Assay")
KB <- NormalizeData(KB)


###cluster using Pagoda2 instead
countMatrix <- KB[["sketch"]]$counts
dim(countMatrix)
#200000

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
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap")
DimPlot(KB, group.by = "rpca_k100_infomap", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Round1-Meta.Rds")

#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k100_infomap"
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
max.idents[1:50,]
max.idents[51:100,]
max.idents[101:150,]
max.idents[151:200,]
max.idents[201:250,]
max.idents[251:300,]

#CL19 - EC/Lymphoid
#CL34 - Lymphoid/FIB
#CL51 - TAL/Lymphoid
#CL60 - PT/FIB
#CL85 - EC/DCT
#CL87 - TAL/DTL
#CL89 - PT/Myeloid
#CL91 - PT/Lymphoid
#CL92 - PT/DTL
#CL96 - EC/Myeloid
#CL97 - Myeloid/FIB
#CL98 - TAL/Myeloid
#CL102 - DTL/EC
#CL104 - Lymphoid/PT
#CL109 - FIB/TAL
#CL110 - TAL/DCT
#CL112 - TAL/PC
#CL113 - PT/PC
#CL118 - PT/DTL
#CL119 - PT/TAL
#CL120 - PT/PC
#CL121 - PT/TAL
#CL128 - PT/EC
#CL129 - PT/DTL
#CL130 - PT/FIB
#CL131 - PT/FIB
#CL134 - PC/EC
#CL136 - TAL/DTL
#CL138 - PC/PapE
#CL141 - PT/TAL
#CL142 - PT/FIB
#CL143 - DCT/PT

poss.mult <- c(19,34,51,60,85,87,89,91,92,96,97,98,102,104,109,110,112,113,118,119,120,121,128,129,
               130,131,134,136,138,141,142,143)

#Check for overlapping cell type markers
Idents(KB) <- "rpca_k100_infomap"
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

DotPlot(KB, features = unique(epi.markers), idents = poss.mult) + RotatedAxis()
DotPlot(KB, features = unique(int.markers), idents = poss.mult) + RotatedAxis()



###Ambiguous (mixed identity) clusters, plot markers to confirm multiplets
#CL19 - EC/Lymphoid
#CL34 - Lymphoid/FIB
#CL51 - TAL/Lymphoid
#CL85 - EC/DCT
#CL87 - TAL/DTL
#CL89 - PT/Myeloid
#CL91 - PT/Lymphoid
#CL96 - EC/Myeloid
#CL97 - Myeloid/FIB
#CL98 - TAL/Myeloid
#CL102 - DTL/EC
#CL104 - Lymphoid/PT
#CL109 - FIB/TAL
#CL110 - TAL/DCT
#CL112 - TAL/PC
#CL113 - PT/PC
#CL119 - Ambiguous
#CL120 - PT/PC
#CL128 - PT/EC
#CL129 - PT/DTL
#CL130 - PT/FIB
#CL131 - PT/FIB
#CL134 - PC/EC
#CL141 - PT/TAL
#CL142 - PT/FIB
#CL143 - DCT/PT

mult <- c(19,34,51,85,87,89,91,96,97,98,102,104,109,110,112,113,119,120,128,129,
          130,131,134,141,142,143)

##Project to full data set
meta <- KB@meta.data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024.Rds")
KB <- AddMetaData(KB, metadata = meta[,"rpca_k100_infomap"], col.name = "rpca_k100_infomap")
KB <- ProjectIntegration(object = KB, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


KB <- ProjectData(object = KB, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                  full.reduction = "integrated.rpca.full", dims = 1:50, refdata = list(rpca_k100_infomap = "rpca_k100_infomap"))

KB <- RunUMAP(KB, reduction = "integrated.rpca.full", dims = 1:50, reduction.name = "umap.full",
              reduction.key = "UMAP_full_")

Idents(KB) <- "rpca_k100_infomap"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DimPlot(KB, reduction = "umap.full", group.by = "group", alpha = 0.1)
DimPlot(KB, reduction = "umap.full", group.by = "rpca_k100_infomap", label = TRUE,
        alpha = 0.1) + NoLegend()


#Remove multiplets
table(Idents(KB))
KB <- subset(KB, idents = mult, invert = TRUE)
KB
#An object of class Seurat 
#70956 features across 330539 samples within 2 assays 
#Active assay: sketch (35478 features, 2000 variable features)
#9 layers present: counts.Gerhardt2023, counts.WashU, counts.Altos, counts.Kirita2020, data.Gerhardt2023, data.WashU, data.Altos, data.Kirita2020, scale.data
#1 other assay present: RNA
#5 dimensional reductions calculated: pca, integrated.rpca, umap, integrated.rpca.full, umap.full

table(KB$source)

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered.Rds")




###Re-run integration (on existing sketch assay)
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
  group.by = c("source", "v2.subclass.l1"),
  combine = TRUE, label.size = 2
)

DimPlot(
  KB,
  reduction = "umap",
  group.by = c("rpca_clusters"),
  combine = FALSE, label.size = 2
)
DimPlot(
  KB,
  reduction = "umap",
  group.by = c("v2.subclass.l1"),
  label = TRUE, label.size = 2
) + NoLegend()
DimPlot(
  KB,
  reduction = "umap",
  group.by = c("celltype"),
  label = TRUE, label.size = 2
) + NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_B.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_B.Rds")


##Join Layers and Annotate
KB[["sketch"]] <- JoinLayers(KB[["sketch"]])
KB[["sketch"]] <- as(object = KB[["sketch"]], Class = "Assay")
KB <- NormalizeData(KB)


###cluster using Pagoda2 instead
countMatrix <- KB[["sketch"]]$counts
dim(countMatrix)
#200000

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
KB <- AddMetaData(KB, metadata = k200infomap, col.name = "rpca_k200_infomap")
DimPlot(KB, group.by = "rpca_k200_infomap", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.hs.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Round2-Meta.Rds")
saveRDS(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Pagoda2-Round2-k200-Object.Rds")



#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k200_infomap"
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
max.idents[1:50,]
max.idents[51:100,]
max.idents[101:150,]
max.idents[151:200,]

#Check for overlapping cell type markers
Idents(KB) <- "rpca_k200_infomap"
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

DotPlot(KB, features = unique(epi.markers), idents = 1:25) + RotatedAxis()
DotPlot(KB, features = unique(int.markers), idents = 1:25) + RotatedAxis()
DotPlot(KB, features = unique(epi.markers), idents = 26:50) + RotatedAxis()
DotPlot(KB, features = unique(int.markers), idents = 26:50) + RotatedAxis()
DotPlot(KB, features = unique(epi.markers), idents = 51:75) + RotatedAxis()
DotPlot(KB, features = unique(int.markers), idents = 51:75) + RotatedAxis()
DotPlot(KB, features = unique(epi.markers), idents = 76:99) + RotatedAxis()
DotPlot(KB, features = unique(int.markers), idents = 76:99) + RotatedAxis()

table(KB$rpca_k200_infomap, KB$Gerhardt2023.Celltype)[1:25,]
table(KB$rpca_k200_infomap, KB$Gerhardt2023.Celltype)[26:50,]
table(KB$rpca_k200_infomap, KB$Gerhardt2023.Celltype)[51:75,]
table(KB$rpca_k200_infomap, KB$Gerhardt2023.Celltype)[76:99,]

###Annotations
reorder <- c(
  52, #POD
  
  49, #PEC
  
  6, #PT-S1
  16, #PT-S1
  64, #PT-S1
  42, #PT-S1
  33, #PT-S1
  53, #PT-S1
  73, #PT-S1
  54, #PT-S1
  12, #PT-S1
  
  13, #PT-S1/S2
  71, #dPT-S1/S2?, reduced marker genes; low genes
  67, #dPT-S1/S2, injury markers; high %MT; high %ERT
  95, #dPT-S1/S2, injury markers
  29, #dPT-S1/S2?, reduced marker genes; low genes
  
  62, #PT-S2
  2, #PT-S2
  9, #PT-S2
  19, #PT-S2
  18, #PT-S2
  59, #PT-S2
  72, #PT-S2
  75, #dPT-S2, injury markers; high %MT; high %ERT
  
  76, #PT-S3
  69, #PT-S3
  93, #PT-S3
  15, #PT-S3
  70, #PT-S3
  87, #PT-S3
  88, #PT-S3
  91, #dPT-S3, injury markers; high %MT; high %ERT
  84, #PT-S3, 
  90, #PT-S3, female
  89, #dPT-S3, female, injury markers; high %MT; high %ERT
  
  31, #aPT-S1/S2
  23, #aPT
  57, #aPT
  85, #aPT-S3/DTL
  82, #aPT?, HAVCR1/SOX9
  81, #aPT?, HAVCR1, low markers
  86, #aPT?, HAVCR1, low markers
  
  97, #dPT; possible pt multiplet; high genes
  99, #dPT; high %MT; high %ERT
  
  79, #cycPT
  
  41, #DTL2 (OM)
  50, #DTL2 (IM)
  47, #DTL1
  60, #DTL1
  74, #aDTL
  65, #DTL3
  51, #ATL
  
  
  28, #M-TAL
  40, #M-TAL
  98, #dM-TAL; injury markers
  80, #dM-TAL; injury markers
  
  14, #C/M-TAL
  17, #C-TAL
  34, #C-TAL + MD
  61, #dC/M-TAL; low genes
  63, #dC-TAL; high %MT; injury markers
  37, #aTAL
  27, #aTAL
  
  4, #DCT1
  78, #dDCT; high %MT
  44, #aDCT
  22, #DCT2
  
  21, #CNT
  7, #CNT-PC
  
  5, #PC (CCD-PC?)
  3, #PC (OMCD-PC?)
  24, #IMCD
  30, #PapE
  1, #IC-A (OMCD-IC-A?)
  8, #IC-A (OMCD-IC-A?)
  20, #IC-B
  32, #IC-B (CNT-CCD-IC)
  
  11, #EC
  35, #EC
  56, #EC
  83, #EC
  66, #cycEC
  96, #EC - epi multiplet?
  94, #EC - epi multiplet?
  
  58, #M-FIB
  10, #C-FIB
  38, #C-FIB
  43, #C-FIB
  46, #MYOF
  48, #pvFIB
  
  36, #VSM/P
  
  68, #B
  55, #PL
  25, #T/NK
  92, #MAC (multiplet?)
  26, #MAC
  45, #FAM
  39, #DC
  77  #cycMAC
  
)


Idents(KB) <- "rpca_k200_infomap"
Idents(KB) <- factor(Idents(KB), levels = reorder)


##Check Proximal markers, injury features and meta and add details to reordered list above
pt.markers <- c("Nphs2","Podxl","Robo2",                              #POD
                "Cp","Nkain3",                                        #PEC
                "Lrp2","Cubn",                                        #PT
                "Slc34a1","Hnf4a",                                    #S1/S2
                "Slc5a12","Prodh2","Slc5a2","Slc7a8",                 #S1
                "Slc22a8","Slc13a3","Cyp2e1","Slco1a1",               #Female S2
                "Slc22a6","Nat8",                                     #Male S2
                "Cyp7b1",                                             #S3
                "Slc7a13", "Slc5a8","Acsm3",                          #Male S3
                "Slc7a12","Acsm1",                                    #Female S3
                "Igfbp7","Spp1","Itgb6","Cdh6","Havcr1",
                
                "Cst3","Apoe","Clu","Vim",                           #aPT2
                "Sox4","Vcam1","Sox9","Ccl2",                                #aPT2
                
                "Gda",                                                        #aPT-S1/S2
                "Dlgap1","Prom1",                                             #frPTS1/S2
                #aPT-S3
                "Dtna","Kitl",                       #frPT-S3
                "Pax2","Pitx2",                                       #TL
                "Slc4a11","Cdh13",                                    #OM DTL2
                "Podn","Aqp1", "Cdh6",                                #IM DTL2
                "Tspan8", "Id3","Jag1","Corin","Slc14a2",             #DTL1 (JM nephron)
                "Wnt5a","Fosl2","Slc12a2",                            #DTL3
                "Clcnka", "Cldn10", "Prox1", "Sptssb",                #ATL
                "Tpt1","Aldob","Gpx3","Wfdc2",        #injury
                "Fabp1","Lcn2",
                "Mki67","Top2a"                                       #cycling
)

DotPlot(KB, features = unique(pt.markers), idents = reorder[1:53]) + RotatedAxis()

row.order <- reorder[1:53]
prop1 <- prop.table(table(KB$sex, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder[1:53])


##Check Distal markers, injury features and meta and add details to reordered list above
dt.markers <- c("Slc12a1", "Egf","Umod",                              #TAL
                "Ank2","Clcnka",                                      #M-TAL
                "Slco3a1",                                            #C-TAL-A
                "Tmem207","Cldn16",                                   #C-TAL-B
                "Enox1",                                              #C-TAL
                "Nos1","Robo2", "Pappa2",                             #MD
                "Creb5","Havcr1","Itgb6","Nrp1",                      #aTAL1
                "Adamts1",                                            #aTAL2
                "Itgb8","Prom1",                                      #frTAL
                
                "Slc12a3",                                            #DCT
                "Trpm6",                                              #DCT1
                "Pgam2","Pvalb",                                      #DCT2
                "Prickle2", "Cdh3",                                   #aDCT
                "Slc8a1","Calb1","S100g","Egfem1",                    #CNT
                "Scnn1g", "Scnn1b",                                   #CNT-PC
                "Gata3","St6gal1","Aqp2","Fxyd4","Aqp3","Ache",       #PC
                "Mcoln3","Btc",                                       #OMCD
                "Aqp4","Slc14a2","Aldh1a3",                           #IMCD
                "Fxyd3", "Krt5","Upk1b",                              #PapE
                "Atp6v1c2","Atp6v0d2",                                #IC
                "Slc26a7","Slc4a1",                                   #ICA
                "Aqp6","Kit",  
                "Slc4a9", "Slc26a4",                                  #ICB
                "Igfbp5","Rhbg",                                      #nonA nonB
                "Hmx2","Sh2d4b",                                      #CNT-CCD-IC
                "Igfbp7","Tpt1","Aldob","Gpx3","Wfdc2","Spp1",        #injury
                "Cst3","Apoe","Fabp1","Lcn2","Clu","Vcam1",
                "Top2a","Mki67"                                        #cycling
)

DotPlot(KB, features = unique(dt.markers), idents = reorder[53:77]) + RotatedAxis()

row.order <- reorder[53:77]
prop1 <- prop.table(table(KB$sex, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder[53:77])





##Check EC markers, injury features and meta and add details to reordered list above
ec.markers <- c("Pecam1", "Ptprb", "Meis2", "Flt1",                   #EC
                "Emcn","Kdr","Hecw2","Ehd3","Abcc4","Plat",           #EC-GC
                "Esm1","Filip1",                                      #Capilary
                "Aqp1","Sox17","Tm4sf1",                              #EC-AEA/DVR
                "Slc14a1",                                            #EC-DVR
                "Fbln5",                                              #EC-AEA
                "Plvap","Slco3a1",                                    #EC-PTC
                "Tll1","Nr2f2", "Gpm6a",                              #EC-AVR
                "Adamtsl1", "Cmklr1", "Dok6","Nav3",                  #EC-PTC
                "Aplnr",                                              #Progenitor?
                
                "Lrp2","Cubn", "Umod", "Slc12a1","Pax2",              #Tubules
                "Slc12a3","Slc8a1","Aqp2","Fxyd4","Atp6v1c2",         #Tubules
                "Fabp1","C7","Acta2",                                 #Others
                "Apoe", "Pdzk1ip1", "Aldob",
                "Mrc1","Cd163","Cd96","Bank1","Cd247","Ptprc",
                "Timp2","Igfbp7","Clu","Cst3",                        #Injury
                "Top2a","Mki67"                                       #cycling
)


DotPlot(KB, features = unique(ec.markers), idents = reorder[78:84]) + RotatedAxis()

row.order <- reorder[78:84]
prop1 <- prop.table(table(KB$sex, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder[78:84])



##Check FIB markers, injury features and meta and add details to reordered list above
fib.markers <- c("Igfbp5","Tnc",#"Kcnk2","Adamtsl1",                    #M-Fib
                 "Cfh","Pdgfra","Lama2","C7",                           #Fib
                 "Ccn2",                   #C-FIB (interstitial fib)
                 "Selenop","Cxcl12",                  #C-FIB-OSMRlo
                 "Osmr", "Cxcl10", "Relb","Ccl2", "Sod2","Ugcg","Il1r1",#Inf Fib
                 "Sparc","Col1a2","Col1a1","Col15a1","Sulf1",           #MYOF                     #aStr
                 "Inhba","Fap","Postn",
                 
                 "Dcn","Meg3",                                          #pvFIB
                 "Flrt2","Igf1","C3","Pi16","Rspo3",   
                 "Col12a1","Fgf14",                                     #Pan pvFIB
                 "Pdzrn4",'Igf1','Adamts3',"Rspo3","Wnt5b",                #pvFIB-RSPO3+
                 "C3","Ebf2","Sfrp2","Cd34","Pi16",                        #pvFIB-PI16+
                 "Itgbl1","Pld5","Cntn5",                                  #pvFIB & pvMYOF
                 "Mgat4c","Epha3",                                         #pvFIB
                 "Itga8",                                                  #pvMYOF & VSMC
                 "Adgrl3","Myh11","Acta2",'Kcnma1',"Pcdh7",                #pvMYOF
                 "Prune2","Myocd", #"SYNPO2",'MACROD2',                     #pvMYOF
                 
                 "Piezo2","Gata3","Itga8",                              #MC
                 "Ren1",                                                #Ren 
                 "Rgs6","Notch3",
                 "Pdgfrb",
                 "Myh11","Rgs5","Acta2"#"Mcam",                        #VSM/P
)


DotPlot(KB, features = unique(fib.markers), idents = reorder[85:91]) + RotatedAxis()

row.order <- reorder[85:91]
prop1 <- prop.table(table(KB$sex, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder[85:91])

##Check IMM markers, injury features and meta and add details to reordered list above
imm.markers <- c("Ptprc",
                 "Bank1","Cd79a","Ms4a1",                              #B cells
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
                 "Ms4a2", "Cpa3", "Kit",                              #Mast
                 "S100a9", "S100a8",#,"Ifitm2",                          #N
                 
                 "Cdh19", "Nrxn1",                                     #SC/Neu
                 "Cd36","Plin1" ,                                       #Adipocytes
                 
                 "Lrp2","Cubn", "Umod", "Slc12a1","Pax2",              #Tubules
                 "Slc12a3","Slc8a1","Aqp2","Fxyd4","Atp6v1c2",         #Tubules
                 "Meis2", "Flt1","Fabp1",                              #Others
                 "Pdzk1ip1", "Aldob",
                 
                 "Timp2","Igfbp7","Clu","Cst3",                        #Injury
                 "Top2a","Mki67"                                       #cycling
)


DotPlot(KB, features = unique(imm.markers), idents = reorder[92:99]) + RotatedAxis()

row.order <- reorder[92:99]
prop1 <- prop.table(table(KB$sex, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Sex Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = c("red","green"))
prop1 <- prop.table(table(KB$patient, KB$rpca_k200_infomap), margin = 2)[,row.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2)

VlnPlot(object = KB, features = c("percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1, idents = reorder[92:99])


##Project to full data set
meta <- KB@meta.data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_B.Rds")
KB <- AddMetaData(KB, metadata = meta[,"rpca_k200_infomap"], col.name = "rpca_k200_infomap")
KB <- ProjectIntegration(object = KB, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


KB <- ProjectData(object = KB, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                  full.reduction = "integrated.rpca.full", dims = 1:50, refdata = list(rpca_k200_infomap = "rpca_k200_infomap"))

KB <- RunUMAP(KB, reduction = "integrated.rpca.full", dims = 1:50, reduction.name = "umap.full",
              reduction.key = "UMAP_full_")

Idents(KB) <- "rpca_k200_infomap"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DimPlot(KB, reduction = "umap.full", group.by = "group", alpha = 0.1)
DimPlot(KB, reduction = "umap.full", group.by = "rpca_k200_infomap", label = TRUE,
        alpha = 0.1) + NoLegend()


#Add preliminary annotations
table(Idents(KB))

Idents(KB) <- "rpca_k200_infomap"
Idents(KB) <- factor(Idents(KB), levels = reorder)
current.cluster.ids <- reorder
new.cluster.ids <- c("POD","PEC","PT-S1","PT-S1","PT-S1","PT-S1","PT-S1","PT-S1",
                     "PT-S1","PT-S1","PT-S1","PT-S1/S2","dPT-S1/S2","dPT-S1/S2",
                     "dPT-S1/S2","dPT-S1/S2","PT-S2","PT-S2","PT-S2","PT-S2","PT-S2",
                     "PT-S2","PT-S2","dPT-S2","PT-S3","PT-S3","PT-S3","PT-S3","PT-S3",
                     "PT-S3","PT-S3","dPT-S3","PT-S3","PT-S3","dPT-S3","aPT-S1/S2",
                     "aPT","aPT","aPT-S3/DTL","aPT","aPT","aPT","dPT","dPT","cycPT",
                     "DTL2","DTL2","DTL1","DTL1","aDTL","DTL3","ATL","M-TAL","M-TAL",
                     "dM-TAL","dM-TAL","C/M-TAL","C-TAL","C-TAL","dC/M-TAL","dC-TAL",
                     "aTAL","aTAL","DCT1","dDCT","aDCT","DCT2","CNT","CNT-PC","PC",
                     "PC","IMCD","PapE","IC-A","IC-A","IC-B","IC-B","EC","EC","EC",
                     "EC","cycEC","EC","EC","M-FIB","C-FIB","C-FIB","C-FIB","MYOF",
                     "pvFIB","VSM/P","B","PL","T/NK","MAC","MAC","FAM","DC","cycMAC"
                     
)

Idents(KB) <- plyr::mapvalues(Idents(KB), from = current.cluster.ids, to = new.cluster.ids)
table(Idents(KB))
KB$v2.subclass.l3 <- Idents(KB)


DimPlot(KB, reduction = "umap.full", group.by = "v2.subclass.l3", label = TRUE,
        alpha = 0.1, repel = TRUE) + NoLegend()



saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_C.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_C.Rds")


