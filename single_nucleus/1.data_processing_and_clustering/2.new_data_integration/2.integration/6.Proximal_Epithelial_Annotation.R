library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

nKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_Atlas_V2_04-2024-subset_Object_global-predictions.Rds")
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Proximal-tubules_subset_filtered_SCT_Reference.Rda")

nKB.PT <- subset(nKB, predicted.v2.subclass.l1 %in% c("POD","PEC","PT","DTL","ATL")) #subset to proximal-intermediate tubule cells
nKB.PT <- subset(nKB.PT, predicted.v2.subclass.l1.score > 0.8) #subset to high confidence cells


###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KB,
  query = nKB.PT,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

nKB.PT <- MapQuery(
  anchorset = anchors,
  query = nKB.PT,
  reference = KB,
  refdata = list(
    v2.subclass.l3 = "v2.subclass.l3",
    v2.clusters = "pt.v2.clusters"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(nKB.PT, reduction = "ref.umap", group.by = "predicted.v2.subclass.l3", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(nKB.PT, reduction = "ref.umap", group.by = "predicted.v2.clusters", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()


hist(nKB.PT$predicted.v2.subclass.l3.score)
hist(nKB.PT$predicted.v2.clusters.score)


#remove prediction assays
nKB.PT[["prediction.score.v2.subclass.l3"]] <- NULL
nKB.PT[["prediction.score.v2.clusters"]] <- NULL
nKB.PT[["pca"]] <- nKB.PT[["ref.pca"]]
nKB.PT[["umap"]] <- nKB.PT[["ref.umap"]]
nKB.PT[["ref.pca"]] <- NULL
nKB.PT[["ref.umap"]] <- NULL


##Merge with full object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered.rda")
#update cluster metadata
meta <- KB.PT@meta.data
meta$pt.v2.clusters <- paste0("P_",meta$pt.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$pt.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.PT@meta.data <- meta
KB.PT$v2.clusters <- KB.PT$pt.v2.clusters
KB.PT$group <- "ref"

nKB.PT$group <- "query"
KB <- merge(x = KB.PT, y = nKB.PT, merge.dr = TRUE)
KB <- JoinLayers(KB)
KB <- NormalizeData(KB)
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "group", repel = TRUE) + NoLegend()


##Clustering in Pagoda2 using integrated PCs
#Filtered Count Matrix from Seurat
#Filtered Count Matrix from Seurat
countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#6914 overdispersed genes

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
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k200_infomap_pt <- k200infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_pt", repel = TRUE) + NoLegend()

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200_round1_newData-0424.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200_round1_newData-0424.rda")

###Identify ambiguous low quality - annotate cluster identities
Idents(KB) <- "pagoda_k200_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v2.clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.v2.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v2.clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v2.clusters)),
                                query = prop.table(table(seurat.obj.sub$predicted.v2.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$predicted.v2.subclass.l3)))
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

max.idents <- getMaxIdents(KB, celltype)
max.idents
write.table(max.idents, file="~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200_round1_newData-0424_overlaps.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/PT_round1_newData-0424_k200-mean-gene_vlnplot.pdf',width=8,height=4)
VlnPlot(KB, features = c("nFeature_RNA"), pt.size = 0)
dev.off()

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


pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/PT_round1_newData-0424_k100-marker-gene_dotplot.pdf',width=16,height=8)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
dev.off()

CL35_markers <- FindMarkers(object = KB, ident.1 = c("CL35"), 
                            only.pos = TRUE)


###Ambiguous clusters (checking mean genes, markers, DEGs)
#CL12 - FIB multiplet
#CL19 - IMM multiplet
#CL25 - EC multiplet
#CL32 - PT/PC multiplet
#CL43 - FIB multiplet
#CL47 - PT/TAL multiplet
#CL50 - PT/IC multiplet


###Subset multiplets
to.remove <- c("12","19","25","32","43","47","50")
to.remove.query <- rownames(KB@meta.data[KB@meta.data$pagoda_k200_infomap_pt %in% to.remove &
                                           KB@meta.data$group %in% "query",])
KB <- subset(KB, cells = to.remove.query, invert = TRUE)
KB
#36588 features across 470913 samples

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")



###Repeat subclustering
#Filtered Count Matrix from Seurat
countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#6815 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 3000, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing
#p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
#k500infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:3000]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_pt <- k100infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k200_infomap_pt <- k200infomap[rownames(KB@meta.data)]
#KB@meta.data$pagoda_k500_infomap_pt <- k500infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_pt", repel = TRUE) + NoLegend()
#DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
#        group.by = "pagoda_k500_infomap_pt", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k100_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k200_0424-newData.rda")
#save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_proximal-tubules_k500_0424-newData.rda")


save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")




###Annotation of k200 clusters
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")
Idents(KB) <- "pagoda_k200_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))
#remove cl43,46:52 (5 or less cells)
KB <- subset(KB, idents = c(43,46:52), invert = TRUE)

levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#CL 1: P_4
#CL 2: P_16
#CL 3: P_17
#CL 4: P_18
#CL 5: P_8
#CL 6: P_23
#CL 7: P_7
#CL 8: P_10
#CL 9: P_14
#CL 10: P_12
#CL 11: P_4
#CL 12: P_1
#CL 13: P_4
#CL 14: P_20
#CL 15: P_3
#CL 16: P_4
#CL 17: P_4
#CL 18: P_4
#CL 19: P_2
#CL 20: P_4
#CL 21: P_5
#CL 22: PT Distal multiplet
#CL 23: P_19
#CL 24: P_24
#CL 25: P_28
#CL 26: P_21
#CL 27: P_25
#CL 28: P_26
#CL 29: P_32
#CL 30: P_30
#CL 31: P_33
#CL 32: P_11
#CL 33: P_10
#CL 34: P_15
#CL 35: Ambiguous
#CL 36: P_31
#CL 37: P_16 (assigned based on marker genes and correlations)
#CL 38: Multiplet
#CL 39: P_13
#CL 40: P_29
#CL 41: P_27
#CL 42: P_20 (assigned based on marker genes and correlations)
#CL 44: P_22
#CL 45: P_29

#missing clusters 6,9, check k100
Idents(KB) <- "pagoda_k100_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#Check clusters for patient/region/cluster distribution
table(KB$pagoda_k200_infomap_pt,KB$patient)[c(21),]
table(KB$pagoda_k200_infomap_pt,KB$region_level2)[c(21),]
table(KB$v2.clusters ,KB$pagoda_k200_infomap_pt)[c("P_6","P_9"),] 
table(KB$v2.clusters ,KB$pagoda_k100_infomap_pt)[c("P_6","P_9"),] 
table(KB$pagoda_k200_infomap_pt,KB$pagoda_k100_infomap_pt)[c(1,37),]

#CL P_6 and P_9 don't seem to be stable with addition of new samples

#Add query/ref correlations to confirm some aPT assignements
marker.genes = VariableFeatures(KB)
KB <- ScaleData(KB, features = marker.genes, assay = "RNA")
Idents(KB) <- "group"
ref <- subset(KB, idents = "ref")
Idents(ref) <- "v2.clusters"
query <- subset(KB, idents = "query")
Idents(query) <- "pagoda_k200_infomap_pt"
levels(Idents(object = query)) <- paste("CL", levels(Idents(object = query)), sep = "")
ave.query <- AverageExpression(query, features = marker.genes, assays = "RNA",
                               layer = "scale.data")
ave.ref <- AverageExpression(ref, features = marker.genes, assays = "RNA",
                             layer = "scale.data")
library("corrplot")
ave.cor <- cor(as.data.frame(ave.ref),as.data.frame(ave.query))
ave.cor
colnames(ave.cor) <- gsub("RNA.","",colnames(ave.cor))

#Top 2 correlations for each cluster
tail(sort(ave.cor[,"CL26"]), 2)
#RNA.P.20  RNA.P.21 
#0.4223030 0.5827525 
tail(sort(ave.cor[,"CL37"]), 2)
#RNA.P.16   RNA.P.6 
#0.4008490 0.4555076
tail(sort(ave.cor[,"CL42"]), 2)
#RNA.P.20  RNA.P.12 
#0.6240661 0.7172918


#Check marker genes
VlnPlot(KB, features = c("nFeature_RNA"), pt.size = 0)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()



### Remove ambiguous clusters
#remove cl 22, 35, 38
Idents(KB) <- "pagoda_k200_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = c(1:42,44:45))
KB <- subset(KB, idents = c(22, 35, 38), invert = TRUE)

#Re-assign cluster numbers for query cells
KB$v2.clusters <- KB$pt.v2.clusters
meta <- KB@meta.data[KB@meta.data$group %in% "query",]
meta$v2.clusters <- meta$pagoda_k200_infomap_pt

cls <- c(1:21,23:34,36:37,39:42,44:45)
new.cls <- c(4,16,17,18,8,23,7,10,14,12,4,1,4,20,3,4,4,4,2,4,5,19,24,28,21,25,26,32,30,33,11,10,15,31,16,13,29,27,20,22,29)
meta$v2.clusters <- plyr::mapvalues(meta$v2.clusters, from = cls, to = new.cls)
table(as.character(meta$v2.clusters))

KB@meta.data[rownames(meta),]$v2.clusters <- paste0("P_",meta$v2.clusters)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = paste0("P_",1:33))

KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_pt", repel = TRUE) + NoLegend()

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pt.v2.clusters", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.clusters", repel = TRUE) + NoLegend()

###Update metadata
meta <- KB@meta.data
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB@meta.data <- meta[,c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","set",                           
                        "v1_barcode","source","assay","experiment_long","experiment","patient",                       
                        "specimen","condition_level3","condition_level2","condition_level1","condition",                                      
                        "percent_cortex","percent_medulla","region_level2","region_level1","age","age_binned","sex",                           
                        "race","location","laterality","protocol","tissue_type_full","tissue_type",                   
                        "atlas_version","v1clusters","v1.subclass.l3","v1.subclass.l2","v1.subclass.l1",
                        "v1.state.l1","v1.state.l2","v1.class","v1.structure","group","predicted.v2.subclass.l1.score",
                        "predicted.v2.subclass.l1","predicted.v2.subclass.l3.score","predicted.v2.subclass.l3",
                        "predicted.v2.clusters.score","predicted.v2.clusters","pagoda_k100_infomap_pt",
                        "pagoda_k200_infomap_pt","v2.clusters","v2.subclass.full","v2.subclass.l3",
                        "v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1",
                        "v2.class","v2.structure")]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()

###Remove unnecessary slots
KB[["RNA"]]$scale.data.1 <- NULL
KB[["RNA"]]$scale.data <- NULL
KB[["ref.pca"]] <- NULL
KB[["ref.umap"]] <- NULL

KB
#36588 features across 463477 samples

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_0424-newData.rda")



###Check  markers
#PT-TL markers
pt.markers <- c(
  "NPHS1","NPHS2","ST6GALNAC3","PODXL",                 #POD
  "PTGDS","BST2","TPPP3","CHI3L1",                      #dPOD
  "ALDH1A2","CFH","FAM155A","CLDN1",                    #PEC
  
  "LRP2","CUBN",                                        #PT
  "SLC5A12","SLC22A6",                                  #S1/S2
  "PRODH2","SLC22A8","PCDH15","RALYL",                  #S1
  "FGA","FGB",                                          #13
  
  "SLC34A1","ANKS1B","SLC5A10","SLC5A11",               #S2                                  
  "SLC22A24","SLC7A13","GPM6A","SLC5A8",                #S3
  
  "ITGB8","CDH6","HAVCR1",                              #aEpi
  "KCTD16","SPON1",                                     #14
  "BTBD11","SYTL3","ZNF250","ITGB3","DLGAP1",           #12
  "PROM1",
  "HDAC9","RAPGEF5","SLC35F3",                          #17
  "SOX4","VIM","IL32","MMP7",                           #19
  "VCAM1",
  "MEG3",                                               #11
  "APBB1IP","LSAMP","KCNIP1","NRXN3",                   #21
  "NRG3","FAM189A1","KITLG","GRM8",                     #22
  
  "CCL2","BIRC3","ICAM1","FSTL3","RELB",                #15
  "DCC","SLC9A9","FAM135B","COL23A1",                   #16
  "PTCHD4",                                             #18
  "ZFPM2","GREB1L","PPFIA2","PDE8B",                    #20
  "TOP2A","MKI67",                                      #cycling
  "DTL","BRIP1","CENPK","EZH2",
  "EGR1","FOS",                                         #dPT
  
  "PAX8-AS1","SLC44A5","CRYAB","TACSTD2",               #TL
  
  "AQP1", "UNC5D","LRRC4C","SLCO1A2","DENND2A",         #DTL2
  
  "SIM2","ADGRL3","JAG1","SMAD9","ID1",                 #DTL1
  "SLC14A2","FAM169A","SMOC2",                          #DTL3
  "ABCA4","BCL6","AKR1B1","SH3GL3",                     #DTL3/ATL
  "CLCNKA","PROX1","CACNA1C",                           #ATL
  "SOD3","HSPA2","CSRP2",
  
  "DEFB1","CST3","APOE","GATM","ALDOB",                 #Injury
  "SPP1","IGFBP7","CLU"
)

pdf(file='proximal/Proximal_Cluster_Marker_Dotplot_0424-newData.pdf',width=14,height=7)
DotPlot(KB, features = pt.markers) + RotatedAxis()
dev.off()
