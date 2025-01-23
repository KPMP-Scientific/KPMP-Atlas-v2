library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
library(anndata)

options(future.globals.maxSize = 1e9)
source("misc/utils.R")

setwd("~/data/net/home/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")
#load("color_factors.robj")

KB.sub <- subset(KB, idents = c(1:100))
KB.sub <- subset(KB.sub, library %in% c("KM150", "KB54"), invert = TRUE)
KB.sub
#36588 features across 1240618 samples within 1 assay



##Subset to immune only
KB.IMM <- subset(KB.sub, v1.subclass.l1 %in% "IMM")
KB.IMM
#36588 features across 61923 samples


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.IMM[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#2180 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 2180, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2180]
VariableFeatures(KB.IMM) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.IMM[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.IMM))
KB.IMM <- FindNeighbors(KB.IMM, dims = 1:50)
KB.IMM <- RunUMAP(object = KB.IMM, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.2, return.model = T)
KB.IMM@meta.data$pagoda_k100_infomap_imm <- k100infomap[rownames(KB.IMM@meta.data)]

DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.IMM, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k100_round1.rda")

Idents(KB.IMM)
KB.IMM[["RNA"]] <- as(object = KB.IMM[["RNA"]], Class = "Assay")
KB.IMM <- ScaleData(KB.IMM)

###Identify ambiguous low quality - annotate cluster identities
Idents(KB.IMM) <- "pagoda_k100_infomap_imm"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = 1:length(levels(Idents(KB.IMM))))
levels(Idents(KB.IMM)) <- paste("CL", levels(Idents(KB.IMM)), sep = "")
celltype <- Idents(KB.IMM)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v1clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$v1.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v1clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v1clusters)),
                                query = prop.table(table(seurat.obj.sub$v1.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$v1.subclass.l3)))
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

max.idents <- getMaxIdents(KB.IMM, celltype)
max.idents

###Check Mean genes
VlnPlot(KB.IMM, features = c("nFeature_RNA"))

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


DotPlot(KB.IMM, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.IMM, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

###Check for DEGs
CL10_markers <- FindMarkers(object = KB.IMM, ident.1 = c("CL10"), 
                            only.pos = TRUE)
CL12_markers <- FindMarkers(object = KB.IMM, ident.1 = c("CL12"), 
                            only.pos = TRUE)

###Ambiguous clusters (checking mean genes, markers, DEGs)
#3 - MACM2/T - high mean genes - TAL markers
#5 - MACM2/T - high mean genes - PT markers
#10 - MACM2/PL - low mean genes - aEpi markers
#17 - MACM2/T - high mean genes - TL markers
#30 - MACM2/B - high mean genes - Stromal markers
#31 - ncMON/T - high mean genes - POD markers
#33 - MACM2/T - high mean genes - PC markers
#39 - B/ncMON - high mean genes - EC markers
#41 - MACM2/T - high mean genes - DCT markers
#45 - MACM2/T - high mean genes - IC markers


                                             
###Subset multiplets
to.remove <- c("CL3","CL5","CL10","CL17","CL30","CL31","CL33","CL39","CL41","CL45")
KB.IMM <- subset(KB.IMM, idents = to.remove, invert = TRUE)
KB.IMM
#36588 features across 48699 samples




###Repeat subclustering
##Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.IMM[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#1795 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 1795, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #ran downstream steps separately

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
k500infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:1795]
VariableFeatures(KB.IMM) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.IMM[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.IMM))
KB.IMM <- RunUMAP(object = KB.IMM, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.1, return.model = T)
KB.IMM@meta.data$pagoda_k100_infomap_imm <- k100infomap[rownames(KB.IMM@meta.data)]
KB.IMM@meta.data$pagoda_k200_infomap_imm <- k200infomap[rownames(KB.IMM@meta.data)]
KB.IMM@meta.data$pagoda_k500_infomap_imm <- k500infomap[rownames(KB.IMM@meta.data)]

DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k100_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k500_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.IMM, features = "nFeature_RNA", reduction = "umap")

save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k100.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k200.rda")
save(p2, sn.od.genes, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Pagoda2_Immune_k500.rda")


save(KB.IMM, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")


###Add in prior analysis of immune subclusters to guide cell type annotations
load("~/data/pod/blake_LTS/Atlas_V2/Immune/Kidney_BUKMAP-KPMP_05042022_v2_Seurat_IMM-only_b.rda")
meta.imm <- KBR.IMM@meta.data
table(meta.imm[rownames(KB.IMM@meta.data),]$imm.clusters)
KB.IMM@meta.data$previous.imm.clusters <- meta.imm[rownames(KB.IMM@meta.data),]$imm.clusters
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "previous.imm.clusters", repel = TRUE) + NoLegend()

save(KB.IMM, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
rm(KBR.IMM)
gc(reset = TRUE)

###Update Experiment metadata
###Add in experiment metadata
meta <- KB.IMM@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race","location","laterality","protocol",
         "tissue_type_full","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB.IMM@meta.data <- meta
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "condition_level1", repel = TRUE) 
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = FALSE, alpha = 0.1,
        group.by = "region_level2", repel = TRUE) 

save(KB.IMM, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")

#Save anndata
expression_matrix <- KB.IMM[["RNA"]]$counts
umap <- Embeddings(KB.IMM, reduction = "umap")
dfobs <- KB.IMM@meta.data
dfvar <- KB.IMM@assays$RNA@meta.features
adata <- AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  obsm=list('X_umap'=umap)
)
write_h5ad(adata, "Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.h5ad")

###Perform CellTypist in python then reload the adata file
adata <- read_h5ad('/home/blake/data/net/home/Human_Kidney/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_celltypist.h5ad', backed = NULL)
adata$obs
KB.IMM@meta.data <- adata$obs[rownames(KB.IMM@meta.data),]

Idents(KB.IMM) <- "pagoda_k200_infomap_imm"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = 1:length(levels(Idents(KB.IMM))))
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_imm", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "predicted_labels", repel = TRUE) + NoLegend()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "majority_voting", repel = TRUE) + NoLegend()

save(KB.IMM, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")


###Annotation of k200 clusters
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
Idents(KB.IMM) <- "pagoda_k200_infomap_imm"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = 1:length(levels(Idents(KB.IMM))))
levels(Idents(KB.IMM)) <- paste("CL", levels(Idents(KB.IMM)), sep = "")
celltype <- Idents(KB.IMM)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$v1clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$v1.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$v1clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$v1clusters)),
                                query = prop.table(table(seurat.obj.sub$v1.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$v1.subclass.l3)))
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

max.idents <- getMaxIdents(KB.IMM, celltype)
max.idents

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$previous.imm.clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$v1.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$previous.imm.clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$previous.imm.clusters)),
                                query = prop.table(table(seurat.obj.sub$v1.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$v1.subclass.l3)))
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

max.idents <- getMaxIdents(KB.IMM, celltype)
max.idents


#CL 1: MAC-M2 (no distinct, cl5) (previous cl 11 - MAC-LYVE1) (merges with cl5 at k500) (weak epithelial signature)
#CL 2: T (previous cl 5 - T-CYT) (Celltypist: TEM/TRM cytotoxic T Cells)
#CL 3: T  (no distinct) (previous cl 4 - T) (merges with cl 11 at k500) (Celltypist: Tcm/Naive helper T Cells)
#CL 4: T (previous cl 6 - T-REG) (Celltypist: Regulatory T Cells)
#CL 5: MAC-M2 (previous cl 11 - MAC-LYVE1) 
#CL 6: B (previous cl 1 - B) (Celltypist: Memory/Naive B Cells) (split?)
#CL 7: ILC3 (previous cl 4 - T) (Celltypist: ILC3)
#CL 8: pDC (previous cl 20 - pDC) (Celltypist: pDC)
#CL 9: ERY (previous cl 9 - ERY) 
#CL 10: LAM (previous cl 14 - LAM)
#CL 11: T (previous cl 4 - T) (Celltypist: TCM/Naive helper T cells)
#CL 12: cDC2 (previous cl 16, cDC2) (Celltypist: DC2)
#CL 13: T (1 distinct, cl11) (previous cl 7 - T-CYT) (Celltypist: TCM/Naive helper T cells)
#CL 14: MDC (previous cl 15 - MDC) 
#CL 15: NKT (previous cl 8 - NKT) (Celltypist: CD16+/CD16- NK cells)
#CL 16: MAST (previous cl 10 - MAST) (Celltypist: Mast cells)
#CL 17: amb (previous cl 24 - ambiguous) 
#CL 18: MAC-MHC (2 distinct) (previous cluster 12 - MAC-MHC) (merges with cl 23 at k500)
#CL 19: mDC (previous cluster 18, mDC) (Celltypist: Migratory DCs)
#CL 20: PL (previous cluster 2, PL) (Celltypist: Plasma cells)
#CL 21: cDC1 (previous cl 19, cDC1) (Celltypist: DC1)
#CL 22: N (previous cluster 21, N) 
#CL 23: MAC-M2 (T cell multiplet?) (no distinct) (some epithelial signatures)
#CL 24: MON (previous cl 17 - ncMON) (Celltypist: Classical Monocytes)
#CL 25: MAC-INF? (weak previous cl 13 - MAC-INF)
#CL 26: ncMON (previous cl 17 - ncMON) (Celltypist: Non-classical monocytes)
#CL 27: cycMNP and cycT (needs splitting) 
#CL 28: MAIT (previous cl 4 - T)
#CL 29: NKT (previous cl 8 - NKT) (Celltypist: Tem/Temra cytotoxic T cells)
#CL 30: MAC-INF (previous cl 13 - MAC-INF)
#CL 31: amb (previous cl 23 - ambiguous)
#CL 32: T (previous cl 7 - T-CYT) (merges more with cl 2 at K500)


##Cluster Assessment:
##Check batch, any clusters associated with only 1-2 individuals
##Check marker genes, sufficient markers to be distinct sub-population? merge at K500?
##check for splitting at k100
##check for overlapping cell type markers - mark as multiplets if cells show overlapping cell type marker profiles

#Check for marker genes
imm.markers <- FindAllMarkers(KB.IMM, only.pos = TRUE, max.cells.per.ident = 2000,
                              logfc.threshold = 0.25, min.pct = 0.25)
write.table(imm.markers, file="Intermediate_Objects/Kidney_v2_Immune_k200_Cluster_markers_05012023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(imm.markers, file="Intermediate_Objects/Kidney_v2_Immune_k200_Cluster_markers_05012023.robj")
#load("Intermediate_Objects/Kidney_v2_Immune_k200_Cluster_markers_05012023.robj")

cl.mark <- imm.markers[imm.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="Intermediate_Objects/Kidney_v2_Immune_k200_Cluster_markers_05012023_top5_distinct.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.IMM, features = top5$gene) + RotatedAxis()

#Check clusters 31 and 17 for patient distribution
table(KB.IMM$pagoda_k200_infomap_imm,KB.IMM$patient)[c(17,31),]
#both predominantly from patient 3691, therefore tagged as ambiguous low quality

table(KB.IMM$pagoda_k200_infomap_imm,KB.IMM$condition_level1)


###compare cluster resolutions
Idents(KB.IMM) <- "pagoda_k200_infomap_imm"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = 1:length(levels(Idents(KB.IMM))))

library("scrattch.hicat")

#Pagoda2 K200 and K500
final<-Idents(KB.IMM)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.IMM$pagoda_k500_infomap_imm
names(prop)<-rownames(KB.IMM@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k500_Immune.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k500", labels=compare.result$cl.id.map$old)
dev.off()


#Pagoda2 K200 and celltypist results
final<-Idents(KB.IMM)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.IMM$majority_voting
names(prop)<-rownames(KB.IMM@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_CellTypist_Immune.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="CellTypist", labels=compare.result$cl.id.map$old)
dev.off()


#Pagoda2 K200 and K100
final<-Idents(KB.IMM)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k200_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.IMM$pagoda_k100_infomap_imm
names(prop)<-rownames(KB.IMM@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k200_vs_k100_Immune.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="k100", labels=compare.result$cl.id.map$old)
dev.off()


#Pagoda2 K100 and celltypist results
Idents(KB.IMM) <- "pagoda_k100_infomap_imm"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = 1:length(levels(Idents(KB.IMM))))
final<-Idents(KB.IMM)
final.df <- data.frame(cluster_id = unique(final),
                       cluster_label = paste0("k100_",unique(final)),
                       cluster_color = rainbow(length(unique(final))))
rownames(final.df) <- final.df$cluster_id
final.df.order<-final.df[order(final.df$cluster_id),]

prop<-KB.IMM$majority_voting
names(prop)<-rownames(KB.IMM@meta.data)
prop <- factor(prop)


pdf("Intermediate_Objects/Pagoda_Cluster_comparision_k100_vs_CellTypist_Immune.pdf",width = 8, height = 8)
compare.result <- compare_annotate(prop, final, final.df.order, reorder = TRUE)
compare.result$g + scale_x_discrete(name ="CellTypist", labels=compare.result$cl.id.map$old)
dev.off()



### Remove ambiguous clusters
#remove cl 17 and cl 31 - both with distinct DEGs yet are specific to a single ref sample
#remove cl 23 - multiplet

Idents(KB.IMM) <- "pagoda_k200_infomap_imm"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = 1:length(levels(Idents(KB.IMM))))
KB.IMM <- subset(KB.IMM, idents = c(17, 23, 31), invert = TRUE)


### Re-order and Merge clusters
#Strategy: merge clusters with no distinct markers (<2) guided by merging found at k500.

#CL 6: B (previous cl 1 - B) (Celltypist: Memory/Naive B Cells) (split?)
#CL 20: PL (previous cluster 2, PL) (Celltypist: Plasma cells)
#CL 3: T  (no distinct) (previous cl 4 - T) (merges with cl 11 at k500) (Celltypist: Tcm/Naive helper T Cells)
#CL 11: T (previous cl 4 - T) (Celltypist: Tcm/Naive helper T cells)
#CL 13: T (1 distinct, cl11) (previous cl 7 - T-CYT) (Celltypist: TCM/Naive helper T cells)
#CL 28: MAIT (previous cl 4 - T) (Celltypist: MAIT)
#CL 7: ILC3 (previous cl 4 - T) (Celltypist: ILC3)
#CL 4: T (previous cl 6 - T-REG) (Celltypist: Regulatory T Cells)
#CL 2: T (previous cl 5 - T-CYT) (Celltypist: TEM/TRM cytotoxic T Cells)
#CL 32: T (previous cl 7 - T-CYT) (merges more with cl 2 at K500)
#CL 29: NKT (previous cl 8 - NKT) (Celltypist: Tem/Temra cytotoxic T cells)
#CL 15: NKT (previous cl 8 - NKT) (Celltypist: CD16+/CD16- NK cells)
#CL 9: ERY (previous cl 9 - ERY) 
#CL 16: MAST (previous cl 10 - MAST) (Celltypist: Mast cells)
#CL 1: MAC-M2 (no distinct, cl5) (previous cl 11 - MAC-LYVE1) (merges with cl5 at k500) (weak epithelial signature)
#CL 5: MAC-M2 (previous cl 11 - MAC-LYVE1) 
#CL 18: MAC-MHC (2 distinct) (previous cluster 12 - MAC-MHC) (merges with cl 23 at k500)
#CL 30: MAC-INF (previous cl 13 - MAC-INF)
#CL 25: MAC-INF? (weak previous cl 13 - MAC-INF)
#CL 10: LAM (previous cl 14 - LAM)
#CL 14: MDC (previous cl 15 - MDC) 
#CL 12: cDC2 (previous cl 16, cDC2) (Celltypist: DC2)
#CL 26: ncMON (previous cl 17 - ncMON) (Celltypist: Non-classical monocytes)
#CL 24: MON (previous cl 17 - ncMON) (Celltypist: Classical Monocytes)
#CL 19: mDC (previous cluster 18, mDC) (Celltypist: Migratory DCs)
#CL 21: cDC1 (previous cl 19, cDC1) (Celltypist: DC1)
#CL 8: pDC (previous cl 20 - pDC) (Celltypist: pDC)
#CL 22: N (previous cluster 21, N) 
#CL 27: cycT (after splitting) 
#CL 27: cycMNP (after splitting) 

cls <- c(6,20,3,11,13,28,7,4,2,32,29,15,9,16,1,5,18,30,25,10,14,12,26,24,19,21,8,22,27)
#merge cl 1 with cl 5
#merge cl 3 with cl 11
new.cls <- c(1:3,3,4:14,14,15:27)
Idents(KB.IMM) <- plyr::mapvalues(Idents(KB.IMM), from = cls, to = new.cls)
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = 1:27)
KB.IMM$imm.v2.clusters <- Idents(KB.IMM)


KB.IMM <- RunUMAP(object = KB.IMM, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                  min.dist = 0.1, return.model = T)

DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "imm.v2.clusters", repel = TRUE) + NoLegend()



###Subset cycling cells
#CL 27 cycMNP and cycT
#Pagoda2
library(pagoda2)
require(parallel)
KB.cyc <- subset(KB.IMM, idents = c(27))
countMatrix <- KB.cyc[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#5 overdispersed genes
p2$calculatePcaReduction(nPcs = 5, n.odgenes = 10, maxit=1000)
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k200infomap <- p2$clusters$PCA$infomap
KB.cyc[["pagoda_k200_infomap_cyc"]] <- k200infomap[rownames(KB.cyc@meta.data)]
KB.cyc <- RunUMAP(object = KB.cyc, reduction = "pca", dims = 1:10, n.neighbors = 30L,
                   min.dist = 0.3)
DimPlot(KB.cyc, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "pagoda_k200_infomap_cyc", repel = TRUE) + NoLegend()
DimPlot(KB.cyc, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
#cl 1 = cyc Lymphocytes
#cl 2 = cyc Myloid cells or cyc MNPs

Idents(KB.cyc) <- "pagoda_k200_infomap_cyc"
cycT <- WhichCells(object = KB.cyc, idents = c(1))
cycM <- WhichCells(object = KB.cyc, idents = c(2))

DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycT, sizes.highlight = 0.1) + NoLegend()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycM, sizes.highlight = 0.1) + NoLegend()
KB.IMM <- SetIdent(KB.IMM, cells = cycT, value = 27)
KB.IMM <- SetIdent(KB.IMM, cells = cycM, value = 28)
KB.IMM$imm.v2.clusters <- Idents(KB.IMM)

Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = c(1:28))
table(Idents(KB.IMM))

DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "imm.v2.clusters", repel = TRUE) + NoLegend()


save(KB.IMM, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
KB.IMM
#36588 features across 45547 samples



###Re-identify IMM cluster markers
imm.markers <- FindAllMarkers(KB.IMM, only.pos = TRUE, max.cells.per.ident = 2000,
                              logfc.threshold = 0.25, min.pct = 0.25)
write.table(imm.markers, file="immune/Kidney_v2_Immune_Cluster_markers_05032023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
save(imm.markers, file="immune/Kidney_v2_Immune_Cluster_markers_05032023.robj")
#load("immune/Kidney_v2_Immune_Cluster_markers_05032023.robj")

cl.mark <- imm.markers[imm.markers$p_val_adj < 0.05,]
cl.mark <- cl.mark[cl.mark$pct.2 < 0.1,]
cl.mark <- cl.mark[cl.mark$pct.1 > 0.3,]
cl.mark <- cl.mark[cl.mark$avg_log2FC > 1,]
cl.mark <- distinct(cl.mark, gene, .keep_all = TRUE)

cl.mark %>% group_by(cluster) %>% top_n(5, avg_log2FC) -> top5

write.table(top5, file="immune/Kidney_v2_Immune_Cluster_markers_05032023_top5_distinct_markers.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

DotPlot(KB.IMM, features = top5$gene) + RotatedAxis()


FeaturePlot(KB.IMM, features = "HLA-DRA", reduction = "umap")
FeaturePlot(KB.IMM, features = "LYVE1", reduction = "umap")
FeaturePlot(KB.IMM, features = "FCGR3A", reduction = "umap") #CD16
FeaturePlot(KB.IMM, features = "CD14", reduction = "umap") 
FeaturePlot(KB.IMM, features = "SPP1", reduction = "umap")
FeaturePlot(KB.IMM, features = "CD9", reduction = "umap")
FeaturePlot(KB.IMM, features = "CD207", reduction = "umap")
FeaturePlot(KB.IMM, features = "MKI67", reduction = "umap")
FeaturePlot(KB.IMM, features = "CD8A", reduction = "umap")


DotPlot(KB.IMM, features = c("CD4","CD8A","CD14","FCGR3A","HLA-DRA","LYVE1")) + RotatedAxis()
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "imm.v2.clusters", repel = TRUE) + NoLegend()




###Comparison with Eraslan et. al., 2022
###Download immune object from Eraslan et. al., 2022
url <- "https://storage.googleapis.com/gtex_analysis_v9/snrna_seq_data/GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.h5ad"
curl::curl_download(url, basename(url))
#GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.h5ad
adata <- read_h5ad('GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.h5ad', backed = NULL)
meta <- adata$obs
features <- adata$var
counts <- t(adata$layers[["counts"]])
counts <- as.matrix(counts)
IMM <- CreateSeuratObject(counts, meta.data = meta)
IMM <- NormalizeData(IMM, normalization.method = "LogNormalize", scale.factor = 10000)
IMM <- ScaleData(IMM)

pca.embeddings <- adata$obsm$X_pca
rownames(pca.embeddings) <- rownames(meta)
IMM[["pca"]] <- CreateDimReducObject(embeddings = pca.embeddings, key = "pca_", assay = DefaultAssay(IMM))

umap.embeddings <- adata$obsm$X_umap
rownames(umap.embeddings) <- rownames(meta)
IMM[["umap"]] <- CreateDimReducObject(embeddings = umap.embeddings, key = "umap_", assay = DefaultAssay(IMM))


colnames(IMM@meta.data) <- c("orig.ident","n_genes","tissue","prep","individual","nGenes",                     
                             "nUMIs","PercentMito","PercentRibo","Age_bin","Sex","Sample.ID",                  
                            "participant_id","Sample.ID.short","Sample.Ischemic.Time.mins",
                             "Tissue.Site.Detail","scrublet","scrublet_score",             
                             "barcode","batch","annotation","broad","granular","leiden",                     
                             "Tissue","LAM_prediction","nCount_RNA","nFeature_RNA")

DimPlot(IMM, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "annotation", repel = TRUE) + NoLegend()
DimPlot(IMM, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "Tissue.Site.Detail", repel = TRUE) + NoLegend()

save(IMM, file = "Eraslan_Science_2022_Immune_subset_Seurat.rda")
load("Eraslan_Science_2022_Immune_subset_Seurat.rda")


###Correlation of Eraslan imm clusters with ours
KB.IMM <- ScaleData(KB.IMM, features = VariableFeatures(KB.IMM))

Idents(object = IMM) <- "annotation"
order <- c("MΦ LYVE1 hi","Mo/MΦ FCGR3A hi","MΦ HLAII hi","Inflammatory MΦ","LAM-like",
           "DC2","Mo/MΦ FCGR3A lo","CD14+ monocyte","CD16+ monocyte","Mature DC","DC1" ,
           "Proliferating MΦ","Langerhans","Lung MΦ" )
Idents(IMM) <- factor(Idents(IMM), levels = order)

select.markers <- intersect(rownames(IMM), VariableFeatures(KB.IMM))
ave.ref<-AverageExpression(IMM, features = select.markers, assays = "RNA",
                           slot = "scale.data")
ave.query<-AverageExpression(KB.IMM, features = select.markers, assays = "RNA",
                             slot = "scale.data")

library("corrplot")
ave.cor<-cor(as.data.frame(ave.ref$RNA),as.data.frame(ave.query$RNA))
write.table(ave.cor, file="immune/Eraslan2022_Cluster_corr.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='immune/Eraslan2022_Cluster_corr_Plot.pdf',width=8,height=8)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor[1:12,14:28], col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()




###Check cell type markers

#lymphocyte markers
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


pdf(file='immune/Immune_Cluster_Lymphoid_Marker_Dotplot.pdf',width=10,height=4.5)
DotPlot(KB.IMM, features = lym.markers, idents = c(1:11,27)) + RotatedAxis()
dev.off()

#Myeloid markers
myel.markers <- c("HBB","HBA2","HBA1",                                #Erythrocyte
                  "IL18R1","CPA3","TPSB2","TPSAB1",                   #MAST
                  "MRC1","CD163","MSR1","CSF1R", "LYVE1","CD14",      #MAC-LYVE1
                  "HLA-DPA1","C1QA","FGL2","CD68","HMOX1","TREM2",    #MAC-MHC
                  "ITGAX","HIF1A","HBEGF","LUCAT1",                   #MAC-INF
                  "TIMP1","OSM", "IL1B",                              #MAC-INF
                  "CCL2","CCL3","CXCL10","CXCL9",                     #MAC-INF
                  "GPNMB","PLA2G7","APOC1","HTRA4",                   #LAM
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

pdf(file='immune/Immune_Cluster_Myeloid_Marker_Dotplot.pdf',width=14,height=5.5)
DotPlot(KB.IMM, features = myel.markers, idents = c(12:26,28)) + RotatedAxis()
dev.off()


#"HLA-DQA1","C1QB","TMEM176B","TMEM176A","CD81","APOE", #Tissue resident MAC - Zimmerman et al. 2019 https://doi.org/10.1681/ASN.2018090931; Fu et al 2019 https://doi.org/10.1681/ASN.2018090896
#"TYROBP","C1QC",

#LAM markers: "C1QA","CD68","TREM2","LGALS3","CD9","APOE",#"SPP1","LIPA","FABP5","PLA2G7","CTSB","CHIT1","CHI3L1", #LAM markers Eraslan et al
DotPlot(KB.IMM, features = c("C1QA","CD68","TREM2","LGALS3","CD9","APOE","SPP1","LIPA","FABP5","PLA2G7","CTSB","CHIT1","CHI3L1"), idents = c(12:26,28)) + RotatedAxis()


###Re-generate celltypist plot
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
#Save anndata
expression_matrix <- KB.IMM[["RNA"]]$counts
umap <- Embeddings(KB.IMM, reduction = "umap")
dfobs <- KB.IMM@meta.data
dfvar <- KB.IMM@assays$RNA@meta.features
adata <- AnnData(
  X=t(expression_matrix),
  obs=dfobs,
  var=dfvar,
  obsm=list('X_umap'=umap)
)
write_h5ad(adata, "Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.h5ad")




###Barplots on stats
load("~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered.rda")
#Generate color factors
load("color_factors.robj")
meta <- KB.IMM@meta.data
meta <- meta[!duplicated(meta$imm.v2.clusters),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "imm.v2.clusters", sort_label = F, colorset = "varibow")
imm.cl.cols <- meta$imm.v2.clusters_color; names(imm.cl.cols) <- meta$imm.v2.clusters_label

meta <- KB.IMM@meta.data
meta <- meta[,c("assay","condition_level3","condition_level2","condition_level1","majority_voting")]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "assay", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level3", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level2", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "condition_level1", sort_label = F, colorset = "varibow")
meta <- annotate_cat(meta, col = "majority_voting", sort_label = F, colorset = "varibow")

assay.cols <- meta$assay_color; names(assay.cols) <- meta$assay_label; assay.cols <- factor(na.omit(assay.cols))
assay.cols <- factor(assay.cols[!duplicated(assay.cols)])
condl3.cols <- meta$condition_level3_color; names(condl3.cols) <- meta$condition_level3_label; condl3.cols <- factor(na.omit(condl3.cols))
condl3.cols <- factor(condl3.cols[!duplicated(condl3.cols)])
condl2.cols <- meta$condition_level2_color; names(condl2.cols) <- meta$condition_level2_label; condl2.cols <- factor(na.omit(condl2.cols))
condl2.cols <- factor(condl2.cols[!duplicated(condl2.cols)])
condl1.cols <- meta$condition_level1_color; names(condl1.cols) <- meta$condition_level1_label; condl1.cols <- factor(na.omit(condl1.cols))
condl1.cols <- factor(condl1.cols[!duplicated(condl1.cols)])
majority_voting.cols <- meta$majority_voting_color; names(majority_voting.cols) <- meta$majority_voting_label; majority_voting.cols <- factor(na.omit(majority_voting.cols))
majority_voting.cols <- factor(majority_voting.cols[!duplicated(majority_voting.cols)])

save(assay.cols,condl3.cols,condl2.cols,condl1.cols,majority_voting.cols,imm.cl.cols,
     file = "immune/color_factors_imm.robj")


#Condition (level 3)
Idents(object = KB.IMM) <- "condition_level3"
pdf(file='immune/Kidney_Immune_Condition_level3_umap.pdf',width=8,height=6)
DimPlot(KB.IMM, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 3"
        ) + scale_color_manual(values = alpha(condl3.cols[levels(Idents(KB.IMM))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()

#Condition (level 1)
Idents(object = KB.IMM) <- "condition_level1"
pdf(file='immune/Kidney_Immune_Condition_level1_umap.pdf',width=8,height=6)
DimPlot(KB.IMM, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.IMM))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
pdf(file='immune/Kidney_Immune_Condition_level1_umap_split.pdf',width=24,height=12)
DimPlot(KB.IMM, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3, split.by = "condition_level1", ncol = 3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Condition Level 1"
        ) + scale_color_manual(values = alpha(condl1.cols[levels(Idents(KB.IMM))], 0.3), name = "Condition"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Celltypist labels
Idents(object = KB.IMM) <- "majority_voting"
pdf(file='immune/Kidney_Immune_Celltypist_Labels_umap.pdf',width=7,height=6)
DimPlot(KB.IMM, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(majority_voting.cols[levels(Idents(KB.IMM))], 0.3), name = "Condition"
        ) + NoLegend()
dev.off()


#Assay
Idents(object = KB.IMM) <- "assay"
pdf(file='immune/Kidney_Immune_Assay_umap.pdf',width=8,height=6)
DimPlot(KB.IMM, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("Assay"
        ) + scale_color_manual(values = alpha(assay.cols[levels(Idents(KB.IMM))], 0.3), name = "Assay"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()


#Clusters
Idents(object = KB.IMM) <- "imm.v2.clusters"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = c(1:28))
pdf(file='immune/Kidney_Immune_Clusters_umap.pdf',width=7,height=6)
DimPlot(KB.IMM, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE, pt.size = 0.5) + ggtitle("V2 Clusters"
        ) + scale_color_manual(values = alpha(imm.cl.cols[levels(Idents(KB.IMM))], 0.3), name = "Clusters"
        ) + NoLegend()
dev.off()


pdf(file='immune/Cluster_Contributions_Barplot.pdf',width=6,height=18)
layout(matrix(c(1,1:6), nrow = 7, ncol = 1, byrow = TRUE))
col.order <- levels(Idents(KB.IMM))
prop1 <- prop.table(table(KB.IMM$patient,KB.IMM$imm.v2.clusters), margin = 2)[,col.order]
barplot(prop1,main = "Patient Proportions", cex.names = 0.5, 
        names.arg = colnames(prop1),las=2, col = as.character(new.patient.cols[rownames(prop1)]))
prop2 <- prop.table(table(KB.IMM$assay,KB.IMM$imm.v2.clusters), margin = 2)[,col.order]
barplot(prop2,main = "Assay Proportions", cex.names = 0.5, 
        names.arg = colnames(prop2),las=2, col = as.character(assay.cols[rownames(prop2)]))
#red = snRNA-seq
prop3 <- prop.table(table(KB.IMM$condition_level3,KB.IMM$imm.v2.clusters), margin = 2)[,col.order]
barplot(prop3,main = "Condition l3 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop3),las=2, col = as.character(condl3.cols[rownames(prop3)]))

prop4 <- prop.table(table(KB.IMM$condition_level1,KB.IMM$imm.v2.clusters), margin = 2)[,col.order]
barplot(prop4,main = "Condition l1 Proportions", cex.names = 0.5, 
        names.arg = colnames(prop4),las=2, col = as.character(condl1.cols[rownames(prop4)]))
batch.entropy<-table(KB.IMM$patient, KB.IMM$imm.v2.clusters)
batch.entropy<-apply(batch.entropy,2,entropy)
barplot(batch.entropy[col.order], col = "gray", main = "Experiment Entropy",  
        names.arg = col.order,las=2,
        cex.names = 0.5, ylim=c(0,1))
tab1 <- table(KB.IMM$imm.v2.clusters)[col.order]
barplot(tab1, col = "gray", main = "Cluster Size", 
        names.arg = col.order,las=2, cex.names = 0.5)
dev.off()
