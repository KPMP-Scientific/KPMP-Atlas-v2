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
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Distal-tubules_subset_filtered_SCT_Reference.Rda")

nKB.DT <- subset(nKB, predicted.v2.subclass.l1 %in% c("TAL","DCT","CNT","PC","PapE","IC")) #subset to distal to collecting tubule cells
nKB.DT <- subset(nKB.DT, predicted.v2.subclass.l1.score > 0.8) #subset to high confidence cells


###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KB,
  query = nKB.DT,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

nKB.DT <- MapQuery(
  anchorset = anchors,
  query = nKB.DT,
  reference = KB,
  refdata = list(
    v2.subclass.l3 = "v2.subclass.l3",
    v2.clusters = "dt.v2.clusters"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(nKB.DT, reduction = "ref.umap", group.by = "predicted.v2.subclass.l3", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(nKB.DT, reduction = "ref.umap", group.by = "predicted.v2.clusters", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()


hist(nKB.DT$predicted.v2.subclass.l3.score)
hist(nKB.DT$predicted.v2.clusters.score)


#remove prediction assays
nKB.DT[["prediction.score.v2.subclass.l3"]] <- NULL
nKB.DT[["prediction.score.v2.clusters"]] <- NULL
nKB.DT[["pca"]] <- nKB.DT[["ref.pca"]]
nKB.DT[["umap"]] <- nKB.DT[["ref.umap"]]
nKB.DT[["ref.pca"]] <- NULL
nKB.DT[["ref.umap"]] <- NULL


##Merge with full object
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered.rda")
#update cluster metadata
meta <- KB.DT@meta.data
meta$dt.v2.clusters <- paste0("D_",meta$dt.v2.clusters)
cl.meta <- read.delim("Kidney_AtlasV2_Cluster_Metadata_11272023.txt")
emc <- c("v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$dt.v2.clusters, cl.meta$v2.clusters)]
}
colnames(meta)
KB.DT@meta.data <- meta
KB.DT$v2.clusters <- KB.DT$dt.v2.clusters
KB.DT$group <- "ref"

nKB.DT$group <- "query"
KB <- merge(x = KB.DT, y = nKB.DT, merge.dr = TRUE)
KB <- JoinLayers(KB)

kidney.mat <- KB[["RNA"]]$counts
kidney.mat <- convert_matrix_type(kidney.mat, type = "uint32_t")

write_matrix_dir(
  mat = kidney.mat,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/distal-tubules_kidney_count_subset",
  overwrite = TRUE
)

kidney.mat <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/distal-tubules_kidney_count_subset")
KB[["RNA"]]$counts <- kidney.mat
KB <- NormalizeData(KB)
KB[["RNA"]]$scale.data.1 <- NULL
KB[["RNA"]]$scale.data <- NULL
KB[["ref.pca"]] <- NULL
KB[["ref.umap"]] <- NULL

rm(KB.DT, nKB, nKB.DT)
gc(reset = TRUE)

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
#6283 overdispersed genes

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
KB@meta.data$pagoda_k200_infomap_dt <- k200infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_dt", repel = TRUE) + NoLegend()

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200_round1_newData-0424.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200_round1_newData-0424.rda")
save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")

###Identify ambiguous low quality - annotate cluster identities
Idents(KB) <- "pagoda_k200_infomap_dt"
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
write.table(max.idents, file="~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200_round1_newData-0424_overlaps.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/DT_round1_newData-0424_k200-mean-gene_vlnplot.pdf',width=8,height=4)
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

table(Idents(KB))
KB.sub <- subset(KB, downsample = 5000)
KB.sub[["RNA"]]$counts <- as(object = KB.sub[["RNA"]]$counts, Class = "dgCMatrix")
KB.sub <- NormalizeData(KB.sub)

pdf(file='~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/DT_round1_newData-0424_k100-marker-gene_dotplot.pdf',width=16,height=8)
DotPlot(KB.sub, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.sub, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
dev.off()

CL35_markers <- FindMarkers(object = KB, ident.1 = c("CL35"), 
                            only.pos = TRUE)


###Ambiguous clusters (checking mean genes, markers, DEGs)
#CL39 - PT multiplet
#CL45 - IMM multiplet
#CL50 - IMM multiplet
#CL52 - FIB multiplet
#CL56 - TAL/PC multiplet
#CL58 - EC multiplet
#CL60 - PT/TAL multiplet
#CL61 - TAL/IMCD multiplet
#CL64 - TL/PC multiplet


###Subset multiplets
to.remove <- c("39","45","50","52","56","58","60","61","64")
to.remove.query <- rownames(KB@meta.data[KB@meta.data$pagoda_k200_infomap_dt %in% to.remove &
                                           KB@meta.data$group %in% "query",])
KB <- subset(KB, cells = to.remove.query, invert = TRUE)
KB
#36588 features across 666981 samples

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")



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
#6271 overdispersed genes

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
KB@meta.data$pagoda_k100_infomap_dt <- k100infomap[rownames(KB@meta.data)]
KB@meta.data$pagoda_k200_infomap_dt <- k200infomap[rownames(KB@meta.data)]
#KB@meta.data$pagoda_k500_infomap_dt <- k500infomap[rownames(KB@meta.data)]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_dt", repel = TRUE) + NoLegend()
#DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
#        group.by = "pagoda_k500_infomap_dt", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k100_0424-newData.rda")
save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k200_0424-newData.rda")
#save(p2, sn.od.genes, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_distal-tubules_k500_0424-newData.rda")


save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
#load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-distal_subset_filtered_0424-newData.rda")




###Annotation of k200 clusters
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
Idents(KB) <- "pagoda_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))
#remove cl55:56 (only one cell each)
KB <- subset(KB, idents = c(55:56), invert = TRUE)

levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents

#CL 1: D_19
#CL 2: D_5
#CL 3: D_15
#CL 4: D_43
#CL 5: D_41
#CL 6: D_10
#CL 7: D_24
#CL 8: D_20
#CL 9: D_30
#CL 10: D_25
#CL 11: D_12
#CL 12: D_14
#CL 13: D_40
#CL 14: D_5
#CL 15: D_9
#CL 16: D_8
#CL 17: D_13
#CL 18: D_3
#CL 19: D_27
#CL 20: D_18
#CL 21: D_48
#CL 22: D_11
#CL 23: D_19
#CL 24: D_16
#CL 25: D_22
#CL 26: D_19
#CL 27: D_24
#CL 28: D_29
#CL 29: D_19
#CL 30: D_31
#CL 31: D_7
#CL 32: D_26
#CL 33: D_33
#CL 34: D_28
#CL 35: D_30
#CL 36: D_21
#CL 37: D_1
#CL 38: D_17
#CL 39: D_45
#CL 40: D_23
#CL 41: D_38
#CL 42: D_5
#CL 43: D_44
#CL 44: D_32
#CL 45: D_9
#CL 46: D_37
#CL 47: D_39
#CL 48: D_6
#CL 49: D_47
#CL 50: D_35
#CL 51: D_46
#CL 52: D_5
#CL 53: TAL/PC multiplet
#CL 54: TAL/PC multiplet

#Check marker genes
table(Idents(KB))
KB.sub <- subset(KB, downsample = 5000)
KB.sub[["RNA"]]$counts <- as(object = KB.sub[["RNA"]]$counts, Class = "dgCMatrix")
KB.sub <- NormalizeData(KB.sub)

VlnPlot(KB, features = c("nFeature_RNA"), pt.size = 0)
DotPlot(KB.sub, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.sub, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()


#missing clusters 2,4,34,36,42, check k100
Idents(KB) <- "pagoda_k100_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)
table(celltype)

max.idents <- getMaxIdents(KB, celltype)
max.idents[max.idents$ref.cluster %in% paste0("D_",c(2,4,34,36,42)),]
#none of the missing clusters are resolved at k100, therefore unstable with addition of new samples

#Check clusters for patient/region/cluster distribution
table(KB$pagoda_k200_infomap_dt,KB$patient)[c(45,53,54),]




### Remove ambiguous clusters
#remove cl 53,54
Idents(KB) <- "pagoda_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = c(1:54))
KB <- subset(KB, idents = c(53:54), invert = TRUE)

#Re-assign cluster numbers for query cells
KB$v2.clusters <- KB$dt.v2.clusters
meta <- KB@meta.data[KB@meta.data$group %in% "query",]
meta$v2.clusters <- meta$pagoda_k200_infomap_dt

cls <- c(1:52)
new.cls <- c(19,5,15,43,41,10,24,20,30,25,12,14,40,5,9,8,13,3,27,18,48,11,19,16,22,19,24,29,19,31,7,26,33,28,30,21,1,17,45,23,38,5,44,32,9,37,39,6,47,35,46,5)
meta$v2.clusters <- plyr::mapvalues(meta$v2.clusters, from = cls, to = new.cls)
table(as.character(meta$v2.clusters))

KB@meta.data[rownames(meta),]$v2.clusters <- paste0("D_",meta$v2.clusters)
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = paste0("D_",1:48))

KB <- RunUMAP(object = KB, reduction = "pca", dims = 1:50, n.neighbors = 20L,
              min.dist = 0.1, return.model = T)

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_dt", repel = TRUE) + NoLegend()

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "dt.v2.clusters", repel = TRUE) + NoLegend()
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
                        "predicted.v2.clusters.score","predicted.v2.clusters","pagoda_k100_infomap_dt",
                        "pagoda_k200_infomap_dt","v2.clusters","v2.subclass.full","v2.subclass.l3",
                        "v2.subclass.l2","v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1",
                        "v2.class","v2.structure")]

DimPlot(KB, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l3", repel = TRUE) + NoLegend()

KB
#36588 features across 666636 samples

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_0424-newData.rda")



###Check  markers
KB.sub <- subset(KB, downsample = 1000)
KB.sub[["RNA"]]$counts <- as(object = KB.sub[["RNA"]]$counts, Class = "dgCMatrix")
dt.markers <- c(
  "CASR","SLC12A1","UMOD","EGF",                       #TAL
  
  "ANK2","KCTD16","CLCNKA","CLDN14",                   #M-TAL
  "BMP6","RIPOR3",
  "IGKC",
  
  #"CACNB2",
  "KCNMB2","RGS6",                                     #C-TAL-A
  "HOMER1","ATF3","FOSB","LDLR",
  "ABCC9",
  "ENOX1","CALCR","CNTN1",
  "CABP1","RBM20",
  "PDE3A",
  
  "BBOX1","NOS1","ROBO2",                              #MD
  "TMPRSS4","CAPN8","CD44",
  
  "TENM4","FGF14","PXDNL","GRM8",                      #C-TAL-B
  "CLDN16","TMEM207",#"TMEM52B",                          
  "COL8A1","LINGO2","ELMO1","IFITM10","OLR1",
  
  "SLC12A3","KLHL3","CNNM2","TRPM6",                    #DCT
  
  "ADAMTS17","GRAMD1B",                                 #DCT1
  
  "SLC8A1","UNC5C","HMCN1","TRPV5",                     #DCT2 / CNT
  
  
  "HSD11B2","CALB1",                                    #CNT
  "SNTG1","ANGPT1","ATP1B3","STC1",                     
  "CTSD",                                               #dCNT
  
  "SCNN1G","SCNN1B","SGCD",                             #PC
  "GATA3","AQP2","PAPPA",
  "FGF12","KCNK13","PTCSC3","CNTNAP3B",                 #CCD-PC
  
  "MCTP1","CPAMD8","GABRA2",                            #OMCD-PC
  
  "GREB1L","SLC14A2","HS3ST5",                          #IMCD
  "TENM3","TGFBR3","AKR1C1","FAT3",
  "AOC1","TFPI2","LCN2",
  
  "ADIRF","SNCG","DHRS2","TP63",                        #PapE                              #PapE
  
  "ATP6V0D2", "ATP6V1C2", "CLNK",                       #IC
  "SLC26A7", "SLC4A1",                                  #IC-A                                   
  "HS6ST3","ABR","TMEM255A","PYY",                      #CCD-IC-A
  "NRXN3", "NXPH2", "LEF1", 
  "FAM184B","ADTRP","AQP6","STAP1",                     #OMCD-IC-A
  
  "SLC4A9", "SLC26A4", "INSRR", "TLDC2",                #IC-B
  
  
  
  "ITGB8","PROM1",                                       #aEpi
  "ARHGAP26","RNF144B","RHEX",
  "CREB5","HIVEP3","COL4A1","CYTOR","DIAPH3","SYT16",
  "ADAMTS1","DNAH5","HAVCR1",#"TFPI","NRP1","ITGB6",
  "RAPGEF5","DLGAP1","BIRC3",
  "EZH2","MKI67","CENPP",#"TOP2A",                       #cycling
  
  
  "CKB","PEBP1","COX8A","UQCRB",                         #Injury
  "CST3","DEFB1","SPP1","IGFBP7","CLU"
)
DotPlot(KB.sub, features = unique(dt.markers)) + RotatedAxis()



