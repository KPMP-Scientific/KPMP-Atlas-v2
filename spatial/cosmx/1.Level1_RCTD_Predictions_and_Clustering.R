library(Seurat)
library(future)
library(ggplot2)
library(pagoda2)
library(spacexr)
library(pagoda2)

plan("multisession", workers = 10)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(future.globals.maxSize = 8000 * 1024^2)
load("color_factors_v2-clusters.robj")

ts <- c("R5341.S1", "R5446.S2")
nano.list <- lapply(ts, function(x) {
  print(paste("Running for Tissue Section:", x))
  obj <- readRDS(paste0("~/hsKidAt/spatial_cosmx/SMI-0070_SanjayJain_WU/5_Raw_data/",x,"_seurat_object_unannotated.Rds"))
  return(obj)
})

names(nano.list) <- ts


###Spatial Deconvolution using RCTD
#Reference data
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")

ref$v2.subclass.l1 <- gsub("/","",ref$v2.subclass.l1) 
ref$v2.subclass.l1 <- gsub("DTL","TL",ref$v2.subclass.l1) 
ref$v2.subclass.l1 <- gsub("ATL","TL",ref$v2.subclass.l1) 

Idents(object = ref) <- "v2.subclass.l1"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT","TL","TAL","DCT","CNT","PC","PapE","IC","EC","FIB",
  "VSMP","Ad","Lymphoid","Myeloid","NEU"))
ref <- subset(ref, idents = c("PapE","Ad","NEU"), invert = TRUE)
ref.ds <- subset(ref, downsample = 10000)
ref.ds[["RNA"]] <- as(object = ref.ds[["RNA"]], Class = "Assay")

# extract information to pass to the RCTD Reference function
counts <- ref.ds[["RNA"]]$counts
cluster <- as.factor(ref.ds$v2.subclass.l1)
names(cluster) <- colnames(ref.ds)
nUMI <- ref.ds$nCount_RNA
names(nUMI) <- colnames(ref.ds)
reference <- Reference(counts, cluster, nUMI)

nano.list <- lapply(ts, function(x) {
  print(paste("Running for Tissue Section:", x))
  
  #RCTD method with full mode
  obj <- nano.list[[x]]
  counts <- obj[["Nanostring"]]$counts
  coords <- GetTissueCoordinates(obj, which = "centroids")
  rownames(coords) <- coords$cell
  coords <- coords[,c("x", "y")]
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  RCTD <- create.RCTD(query, reference, max_cores = 8)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  
  #add l1 prediction weights
  rctd.cls <- RCTD@results$weights
  rctd.cls <- 100 * sweep(rctd.cls, 1, rowSums(rctd.cls), "/")
  rctd.cls <- as.data.frame(rctd.cls)
  maxWeight.l1 <- as.numeric(apply(rctd.cls, 1, max))
  names(maxWeight.l1) <- rownames(rctd.cls)
  maxWeight.l1.ct <- colnames(rctd.cls)[apply(rctd.cls, 1, which.max)]
  names(maxWeight.l1.ct) <- rownames(rctd.cls)
  
  rctd.cls <- as.data.frame(rctd.cls)[colnames(obj),]
  rctd.cls[is.na(rctd.cls)] <- 0
  rownames(rctd.cls) <- colnames(obj)
  
  obj[["l1.predictions"]] <- CreateAssayObject(counts = t(as.matrix(rctd.cls)))
  obj$maxWeight.l1 <- as.numeric(maxWeight.l1[rownames(obj@meta.data)])
  obj$maxWeight.l1.ct <- as.character(maxWeight.l1.ct[rownames(obj@meta.data)])
  
  return(obj) 
})

names(nano.list) <- ts

nano.obj <- merge(x = nano.list[[1]],
                  y = c(nano.list[[2]]))  
nano.obj


table(nano.obj$maxWeight.l1.ct)
hist(nano.obj$maxWeight.l1)


###Visualize level 1 annotations
#Update color factors
v2.scl1.cols.mod <- v2.scl1.cols
names(v2.scl1.cols.mod) <- gsub("/","",names(v2.scl1.cols.mod))
names(v2.scl1.cols.mod) <- gsub("DTL","TL",names(v2.scl1.cols.mod))

ImageDimPlot(nano.obj, fov = "R5341.S1", group.by = "maxWeight.l1.ct",
             axes = TRUE, cols = v2.scl1.cols.mod, 
             size = 0.05, border.size = 0.001)
ImageDimPlot(nano.obj, fov = "R5446.S2", group.by = "maxWeight.l1.ct",
             axes = TRUE, cols = v2.scl1.cols.mod, 
             size = 0.05, border.size = 0.001)


###Update metadata
obj <- readRDS("~/hsKidAt/spatial_cosmx/SMI-0070_SanjayJain_WU/5_Raw_data/seurat_object.Rds")
meta <- obj@meta.data

nano.obj@meta.data
table(unlist(lapply(rownames(nano.obj@meta.data),function(x) unlist(strsplit(x,"_"))[3])))
nano.obj@meta.data$Run_Tissue_name <- unlist(lapply(rownames(nano.obj@meta.data),function(x) unlist(strsplit(x,"_"))[3]))
current.ids <- c(1:2)
new.ids <- ts
nano.obj@meta.data$Run_Tissue_name <- plyr::mapvalues(nano.obj@meta.data$Run_Tissue_name, 
                                                      from = current.ids, to = new.ids)
table(nano.obj@meta.data$Run_Tissue_name)
nano.obj@meta.data$barcode <- paste0(nano.obj@meta.data$Run_Tissue_name,"_",
                                     unlist(lapply(rownames(nano.obj@meta.data),function(x) unlist(strsplit(x,"_"))[1])),
                                     "_", 
                                     unlist(lapply(rownames(nano.obj@meta.data),function(x) unlist(strsplit(x,"_"))[2])))

meta <- meta[nano.obj@meta.data$barcode,]
rownames(meta) <- rownames(nano.obj@meta.data)
meta <- meta[,4:25]
nano.obj <- AddMetaData(nano.obj, metadata = meta)

saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024.RDS")




###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#84 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference overdispersed genes
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
var.genes <- VariableFeatures(ref)
sn.od.genes <- var.genes[var.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap <- k100infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxWeight.l1.ct") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024.RDS")



###Annotate to subclass.l1
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024.RDS")

##Generate table of cluster cell types
Idents(nano.obj) <- "pagoda_k100_infomap"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxWeight.l1.ct))), 2)
    max.len = 2
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(query = prop.table(table(seurat.obj.sub$maxWeight.l1.ct))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxWeight.l1.ct)))
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(nano.obj, celltype)
max.idents

##Check marker genes
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


Idents(nano.obj) <- "pagoda_k100_infomap"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))

DotPlot(nano.obj, features = epi.markers[epi.markers %in% rownames(nano.obj)], dot.scale = 8) + RotatedAxis()
DotPlot(nano.obj, features = int.markers[int.markers %in% rownames(nano.obj)], dot.scale = 8) + RotatedAxis()


###Ambiguous clusters
#3 - no distinct markers - ambiguous
#8 - Erythroid?
#31 - cycling cells

CL3_markers <- FindMarkers(object = nano.obj, ident.1 = 3, 
                            max.cells.per.ident = 500, only.pos = TRUE)
CL8_markers <- FindMarkers(object = nano.obj, ident.1 = 8, 
                            max.cells.per.ident = 500, only.pos = TRUE)
CL31_markers <- FindMarkers(object = nano.obj, ident.1 = 31, 
                            max.cells.per.ident = 500, only.pos = TRUE)

#Checked ambigous clusters showing poor prediction assignments and poor correlations
#and also showing very low DEGs were labeled as ambiguous and removed


###Re-run clustering without ambiguous
nano.obj <- subset(nano.obj, idents = 3, invert = TRUE)

###Clustering 
countMatrix <- nano.obj[["Nanostring"]]$counts
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[3])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
library(pagoda2)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#90 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the overdispersed genes intersected with reference overdispersed genes
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
var.genes <- VariableFeatures(ref)
sn.od.genes <- var.genes[var.genes %in% rownames(var.info)[1:2000]]
p2$calculatePcaReduction(nPcs = 50, odgenes = sn.od.genes, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
cell.embeddings <- p2$reductions$PCA
nano.obj[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(nano.obj))

nano.obj <- RunUMAP(object = nano.obj, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                    min.dist = 0.2, return.model = T)

nano.obj@meta.data$pagoda_k100_infomap <- k100infomap[rownames(nano.obj@meta.data)]

DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "pagoda_k100_infomap") + NoLegend()
DimPlot(nano.obj, label = T, label.size = 3, reduction = "umap", raster=FALSE,
        group.by = "maxWeight.l1.ct") + NoLegend()


saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")



###Annotate to subclass.l1
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")

##Generate table of cluster cell types
Idents(nano.obj) <- "pagoda_k100_infomap"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))
levels(Idents(nano.obj)) <- paste("CL", levels(Idents(nano.obj)), sep = "")
celltype <- Idents(nano.obj)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    query.top2 <- tail(names(sort(table(seurat.obj.sub$maxWeight.l1.ct))), 2)
    max.len = 2
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(query = prop.table(table(seurat.obj.sub$maxWeight.l1.ct))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$maxWeight.l1.ct)))
    if (is.null(Idents.called$query.Var1) == TRUE) {
      Idents.called$query.Var1 <- rep(NA, 2)
      Idents.called$query.Freq <- rep(NA, 2)
      Idents.called$query.Total <- rep(NA, 2)
      Idents.called <- Idents.called[,c("query.Var1",
                                        "query.Freq","query.Total")]
    } else {
      print(paste("query cell types called for:", ct))
    }
    
    Idents.called$cluster <- rep(ct, length(rownames(Idents.called)))
    colnames(Idents.called) <- c("query.subclass",
                                 "query.Freq","query.Total","int.cluster")
    rownames(Idents.called) <- paste(Idents.called$int.cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(nano.obj, celltype)
max.idents

##Check marker genes
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


Idents(nano.obj) <- "pagoda_k100_infomap"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = 1:length(levels(Idents(nano.obj))))

DotPlot(nano.obj, features = epi.markers[epi.markers %in% rownames(nano.obj)], dot.scale = 8) + RotatedAxis()
DotPlot(nano.obj, features = int.markers[int.markers %in% rownames(nano.obj)], dot.scale = 8) + RotatedAxis()


###Cluster annotations
#16 - POD
#2 - PT
#6 - PT (aPT)
#34 - PEC
#35 - TL
#8 - TAL (aTAL)
#30 - TAL
#33 - TAL
#12 - TAL (aTAL)
#10 - DCT/CNT
#22 - IC
#20 - IC
#19 - PC
#1 - EC
#24 - EC
#32 - EC
#14 - EC 
#36 - EC
#28 - VSMP (MC?)
#3 - VSMP
#25 - VSMP
#4 - FIB
#17 - FIB
#29 - FIB
#15 - FIB
#5 - Myeloid
#13 - Myeloid
#21 - Myeloid (MAST cells)
#9 - Lymphoid
#27 - Lymphoid
#23 - Lymphoid
#11 - Lymphoid

#7 - cycling cells
#18 - Erythrocytes

#31 - Ambiguous
#26 - Ambiguous (no marker genes)
#37 - Ambiguous (1 cell only)


CL7_markers <- FindMarkers(object = nano.obj, ident.1 = 7, 
                           max.cells.per.ident = 500, only.pos = TRUE)
CL18_markers <- FindMarkers(object = nano.obj, ident.1 = 18, 
                           max.cells.per.ident = 500, only.pos = TRUE)
CL31_markers <- FindMarkers(object = nano.obj, ident.1 = 31, 
                            max.cells.per.ident = 500, only.pos = TRUE)
CL21_markers <- FindMarkers(object = nano.obj, ident.1 = 21, 
                            max.cells.per.ident = 500, only.pos = TRUE)
CL26_markers <- FindMarkers(object = nano.obj, ident.1 = 26, 
                            max.cells.per.ident = 500, only.pos = TRUE)

#Checked ambigous clusters showing poor prediction assignments and poor correlations
#and also showing very low DEGs were labeled as ambiguous and removed




###Generate filtered nano.obj object
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")
Idents(nano.obj) <- "pagoda_k100_infomap"
cls <- c(16,2,6,34,35,8,30,33,12,10,22,20,19,1,24,32,14,36,28,3,25,4,17,29,15,5,
         13,21,9,27,23,11,7,18,31,26,37)
new.cls <- c("POD","PT","PT","PEC","TL","TAL","TAL","TAL","TAL","DCT/CNT","IC","IC","PC",
             "EC","EC","EC","EC","EC","VSMP","VSMP","VSMP","FIB","FIB","FIB","FIB","Myeloid",
             "Myeloid","Myeloid","Lymphoid","Lymphoid","Lymphoid","Lymphoid",
             "cycling","ERY","NA","NA","NA")
   
Idents(nano.obj) <- plyr::mapvalues(Idents(nano.obj), from = cls, to = new.cls)
Idents(nano.obj) <- factor(Idents(nano.obj), levels = c("POD","PEC","PT","TL","TAL","DCT/CNT",
                                                        "PC","IC","EC","FIB","VSMP",
                                                        "Lymphoid","Myeloid","ERY","cycling","NA"))

table(Idents(nano.obj))
nano.obj$v2.subclass.l1 <- Idents(nano.obj)

v2.scl1.cols.mod <- v2.scl1.cols
names(v2.scl1.cols.mod) <- gsub("/","",names(v2.scl1.cols.mod))
names(v2.scl1.cols.mod) <- gsub("DTL","TL",names(v2.scl1.cols.mod))
names(v2.scl1.cols.mod) <- gsub("DCT","DCT/CNT",names(v2.scl1.cols.mod))

DimPlot(nano.obj, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,raster=FALSE,
        group.by = "v2.subclass.l1", repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = alpha(v2.scl1.cols.mod[levels(Idents(object = nano.obj))], 0.1), name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))

saveRDS(nano.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")
#nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Filtered.RDS")

