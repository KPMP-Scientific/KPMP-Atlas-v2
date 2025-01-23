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
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_SCT_Reference_full.Rda")

nKB.IMM <- subset(nKB, predicted.v2.subclass.l1 %in% c("Lymphoid","Myeloid")) #subset to stromal cells
nKB.IMM <- subset(nKB.IMM, predicted.v2.subclass.l1.score > 0.8) #subset to high confidence stromal cells


###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KB,
  query = nKB.IMM,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

nKB.IMM <- MapQuery(
  anchorset = anchors,
  query = nKB.IMM,
  reference = KB,
  refdata = list(
    v2.subclass.l3 = "v2.subclass.l3",
    v2.clusters = "imm.v2.clusters"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(nKB.IMM, reduction = "ref.umap", group.by = "predicted.v2.subclass.l3", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(nKB.IMM, reduction = "ref.umap", group.by = "predicted.v2.clusters", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()


hist(nKB.IMM$predicted.v2.subclass.l3.score)
hist(nKB.IMM$predicted.v2.clusters.score)


#remove prediction assays
nKB.IMM[["prediction.score.v2.subclass.l3"]] <- NULL
nKB.IMM[["prediction.score.v2.clusters"]] <- NULL
nKB.IMM[["pca"]] <- nKB.IMM[["ref.pca"]]
nKB.IMM[["umap"]] <- nKB.IMM[["ref.umap"]]
nKB.IMM[["ref.pca"]] <- NULL
nKB.IMM[["ref.umap"]] <- NULL


##Merge with full object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
#update cluster metadata
meta <- KB.IMM@meta.data
meta$imm.v2.clusters <- paste0("I_",meta$imm.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$imm.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.IMM@meta.data <- meta
KB.IMM$v2.clusters <- KB.IMM$imm.v2.clusters
KB.IMM$group <- "ref"

nKB.IMM$group <- "query"
KB <- merge(x = KB.IMM, y = nKB.IMM, merge.dr = TRUE)
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
#2520 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 2520, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2520]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_imm <- k100infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_imm", repel = TRUE) + NoLegend()

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k100_round1_newData-0424.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k100_round1_newData-0424.rda")

###Identify ambiguous low quality - annotate cluster identities
Idents(KB) <- "pagoda_k100_infomap_imm"
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
write.table(max.idents, file="~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k100_round1_newData-0424_overlaps.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Immune_round1_newData-0424_k100-mean-gene_vlnplot.pdf',width=8,height=4)
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


pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Immune_round1_newData-0424_k100-marker-gene_dotplot.pdf',width=16,height=8)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
dev.off()

CL35_markers <- FindMarkers(object = KB, ident.1 = c("CL35"), 
                            only.pos = TRUE)


###Ambiguous clusters (checking mean genes, markers, DEGs)
#CL5 - PT Multiplet
#CL33 - FIB Multiplet
#CL44 - FIB Multiplet
#CL51 - EC Multiplet
#CL53 - TAL Multiplet
#CL54 - Epi Multiplet
#CL55 - FIB Multiplet
#CL59 - PC Multiplet
#CL62 - Epi Multiplet
#CL63 - Epi Multiplet
#CL65 - Epi Multiplet


###Subset multiplets
to.remove <- c("5","33","44","51","53","54","55","59","62","63","65")
to.remove.query <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_imm %in% to.remove &
                                           KB@meta.data$group %in% "query",])
KB <- subset(KB, cells = to.remove.query, invert = TRUE)
KB
#36588 features across 74361 samples

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")



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
#2158 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 2158, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #done sucessively with downstream P2 processing

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2158]
VariableFeatures(KB) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB))
KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)
KB@meta.data$pagoda_k100_infomap_imm <- k100infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k200_infomap_imm <- k200infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k500_infomap_imm <- k500infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k100_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k200_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k500_0424-newData.rda")


save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")




###Annotation of k200 clusters
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")
Idents(KB) <- "pagoda_k200_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))
#remove CL37 of 2 cells
KB <- subset(KB, idents = "37", invert = TRUE)

levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#CL 1: I_14
#CL 2: I_8
#CL 3: I_3
#CL 4: I_7
#CL 5: I_1
#CL 6: I_6
#CL 7: I_25
#CL 8: I_12
#CL 9: I_14
#CL 10: I_14
#CL 11: I_3
#CL 12: I_20
#CL 13: I_4
#CL 14: I_19
#CL 15: I_11
#CL 16: I_14
#CL 17: I_13
#CL 18: I_15
#CL 19: I_23
#CL 20: I_2
#CL 21: I_22
#CL 22: I_24
#CL 23: I_18
#CL 24: I_11
#CL 25: I_26
#CL 26: I_4
#CL 27: I_17
#CL 28: I_14
#CL 29: I_20
#CL 30: I_21
#CL 31: I_28 (needs splitting? cycT and cycMAC) sep out k100 cl 38 as new cluster 37 = I_27
#CL 32: I_8
#CL 33: I_5
#CL 34: I_14
#CL 35: I_8
#CL 36: PC Multiplet
#CL 37: I_27 (k100 cluster 38)
#CL 38: I_10 (K100 cluster 45)
#CL 39: I_16 (k100 cluster 39)

#missing clusters 10,16, check k100
Idents(KB) <- "pagoda_k100_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#Check clusters for patient/region/cluster distribution
table(KB$pagoda_k200_infomap_imm,KB$patient)[c(31,35),]
table(KB$pagoda_k200_infomap_imm,KB$region_level2)[c(31,35),]
table(KB$pagoda_k200_infomap_imm,KB$pagoda_k100_infomap_imm)[c(31,35),]
table(KB$pagoda_k100_infomap_imm,KB$v2.clusters)[c(38,46),] #38 = cycT, 46 = cycMAC
table(KB$v2.clusters ,KB$pagoda_k200_infomap_imm)[c("I_10"),] 
table(KB$pagoda_k200_infomap_imm,KB$pagoda_k100_infomap_imm)[c("24"),] #24 splits to 2 clusters at k100 (29,45)
table(KB$v2.clusters,KB$pagoda_k100_infomap_imm)[c("I_10","I_11"),] #K100 29 = I_11, 45 = I_10
table(KB$v2.clusters ,KB$pagoda_k200_infomap_imm)[c("I_16"),] 
table(KB$pagoda_k200_infomap_imm,KB$pagoda_k100_infomap_imm)[c("16"),] #24 splits to several clusters at k100 (inc 39)
table(KB$v2.clusters,KB$pagoda_k100_infomap_imm)[c("I_16"),] #K100 39 = I_16


#Add k100 cluster to k200 clusters (cycT cycMAC split)
cl37 <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_imm %in% "38",])
KB@meta.data[cl37,]$pagoda_k200_infomap_imm <- "37"
Idents(KB) <- "pagoda_k200_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))
cl38 <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_imm %in% "45",])
KB@meta.data$pagoda_k200_infomap_imm <- as.character(KB@meta.data$pagoda_k200_infomap_imm)
KB@meta.data[cl38,]$pagoda_k200_infomap_imm <- "38"
cl39 <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap_imm %in% "39",])
KB@meta.data[cl39,]$pagoda_k200_infomap_imm <- "39"
Idents(KB) <- "pagoda_k200_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))


VlnPlot(KB, features = c("nFeature_RNA"), pt.size = 0)
DotPlot(KB, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()



### Remove ambiguous clusters
#remove cl 36
Idents(KB) <- "pagoda_k200_infomap_imm"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
KB <- subset(KB, idents = c(36), invert = TRUE)

#Re-assign cluster numbers for query cells
KB$v2.clusters <- KB$imm.v2.clusters
meta <- KB@meta.data[KB@meta.data$group %in% "query",]
meta$v2.clusters <- meta$pagoda_k200_infomap_imm

cls <- c(1:35,37:39)
new.cls <- c(14,8,3,7,1,6,25,12,14,14,3,20,4,19,11,14,13,15,23,2,22,24,18,11,26,4,17,14,20,21,28,8,5,14,8,27,10,16)
meta$v2.clusters <- plyr::mapvalues(meta$v2.clusters, from = cls, to = new.cls)
table(as.character(meta$v2.clusters))

KB@meta.data[rownames(meta),]$v2.clusters <- paste0("I_",meta$v2.clusters)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = paste0("I_",1:28))

KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_imm", repel = TRUE) + NoLegend()

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "imm.v2.clusters", repel = TRUE) + NoLegend()
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
                        "predicted.v2.clusters.score","predicted.v2.clusters","pagoda_k100_infomap_imm",
                        "pagoda_k200_infomap_imm","pagoda_k500_infomap_imm","v2.clusters","v2.subclass.full","v2.subclass.l3",
                        "v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1",
                        "v2.class","v2.structure")]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")
KB
#36588 features across 72454 samples



###Check Immune markers
#lymphoid
lym.markers <- c("BANK1","MS4A1","CD79A",                            #B Cells
                 "IGKC","JCHAIN","MZB1",                             #PL Cells
                 "IL7R","LEF1","TRAC","CD4","CD3D","CCR7","SELL",    #T Cells
                 "SLC4A10","KLRB1","CCR6",                           #MAIT
                 "TOX2","KIT","RORC",                                #ILC3
                 "RTKN2","IL2RA","CTLA4","FOXP3",                    #T-REG
                 "GZMK","CD8A","CCL5","CCL4","GZMA",                 #T-CYT
                 "NKG7","GNLY","GZMB","CX3CR1",                      #tdT-CYT/NK
                 "KLRF1","PDGFD","NCAM1",                            #NK
                 
                 "TOP2A","MKI67"                                     #cycling
)


pdf(file='immune/Immune_Cluster_Lymphoid_Marker_Dotplot_0424-newData.pdf',width=10,height=4.5)
DotPlot(KB, features = lym.markers, idents = paste0("I_",c(1:11,27))) + RotatedAxis()
dev.off()


#Myeloid markers
myel.markers <- c("HBB","HBA2","HBA1",                                #Erythrocyte
                  "IL18R1","CPA3","TPSB2","TPSAB1",                   #MAST
                  "MRC1","CD163","MSR1","CSF1R", "LYVE1","CD14",      #MAC-LYVE1
                  "HLA-DPA1","C1QA","FGL2","CD68","HMOX1","TREM2",    #MAC-MHC
                  "ITGAX","HIF1A","HBEGF","LUCAT1",                   #MAC-INF
                  "TIMP1","OSM", "IL1B",                              #MAC-INF
                  "CCL2","CCL3","CXCL10","CXCL9",                     #MAC-INF
                  "GPNMB","SPP1","PLA2G7","APOC1","CAPG",             #SPP1+ (SAM/LAM)
                  "C3","KCNQ3","ADGRB3","VASH1","CX3CR1",             #MDC
                  "CLEC10A","FCER1A","CD1C",                          #cDC2
                  "TCF7L2","FCGR3A","FCN1",                           #ncMON
                  "LYZ", "CD36",                                      #MON
                  "CCR7","EBI3","CCL19",                              #mDC
                  "CADM1","CLEC9A","BATF3",                           #cDC1
                  "IL3RA","LILRA4","PLD4",                            #pDC
                  "S100A9","FCGR3B","S100A8","IFITM2",                #N
                  "TOP2A","MKI67"
)

pdf(file='immune/Immune_Cluster_Myeloid_Marker_Dotplot_0424-newData.pdf',width=14,height=5.5)
DotPlot(KB, features = myel.markers, idents = paste0("I_",c(12:26,28))) + RotatedAxis()
dev.off()

table(Idents(KB))


###Remove unnecessary slots
KB[["RNA"]]$scale.data.1 <- NULL
KB[["RNA"]]$scale.data <- NULL
KB[["ref.pca"]] <- NULL
KB[["ref.umap"]] <- NULL

KB

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_0424-newData.rda")


