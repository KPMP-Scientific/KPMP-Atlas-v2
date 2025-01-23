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
setwd("~/data/net/home/Human_Kidney/Atlas_V2")

options(future.globals.maxSize = 1e9)
source("misc/utils.R")

options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")

KB.sub <- subset(KB, idents = c(1:100))
KB.sub <- subset(KB.sub, library %in% c("KM150", "KB54"), invert = TRUE)
KB.sub
#36588 features across 1240618 samples within 1 assay

table(KB.sub$v1.subclass.l3)
##Subset to Neural only
KB.N <- subset(KB.sub, v1.subclass.l3 %in% c("SC/NEU"))
KB.N
#36588 features across 3286 samples


###Pagoda2 Clustering
#Filtered Count Matrix from Seurat
countMatrix <- KB.N[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#28 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 500, maxit=1000)

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

k200infomap <- p2$clusters$PCA$infomap

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:500]
VariableFeatures(KB.N) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.N[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.N))
KB.N <- RunUMAP(object = KB.N, reduction = "pca", dims = 1:50, n.neighbors = 20L,
                 min.dist = 0.2, return.model = T)
KB.N@meta.data$pagoda_k200_infomap_n <- k200infomap[rownames(KB.N@meta.data)]

DimPlot(KB.N, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_n", repel = TRUE) + NoLegend()
DimPlot(KB.N, reduction = "ref.umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "pagoda_k200_infomap_n", repel = TRUE) + NoLegend()
DimPlot(KB.N, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
FeaturePlot(KB.N, features = "nFeature_RNA", reduction = "umap")


Idents(KB.N)
KB.N[["RNA"]] <- as(object = KB.N[["RNA"]], Class = "Assay")
KB.N <- ScaleData(KB.N)

###Identify ambiguous low quality - annotate cluster identities
Idents(KB.N) <- "pagoda_k200_infomap_n"

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
                 "S100A9","S100A8","FCGR3B",                     #NC
                 "CDH19", "NRXN1", "GINS3" #SC/NEU
)


DotPlot(KB.N, features = unique(epi.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()
DotPlot(KB.N, features = unique(int.markers),cols = c("lightgrey", "#225ea8"),dot.scale = 4) + RotatedAxis()

#only clusters 1 and 2 shows SC/NEU marker genes
#cluster 2 shows TAL markers

SC.NEU.cells <- WhichCells(KB.N, idents = 1)
saveRDS(SC.NEU.cells, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_SC-NEU_Cell_subset.rda")

non.SC.NEU.cells <- WhichCells(KB.N, idents = c(2,3))
saveRDS(non.SC.NEU.cells, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_nonSC-NEU_Cell_subset.rda")
                                             
neu.markers <- FindAllMarkers(KB.N, only.pos = TRUE, logfc.threshold = 0.25)
neu.markers[neu.markers$cluster == 3,]$gene
table(Idents(KB.N))

table(KB.N$pagoda_k200_infomap_n,KB.N$patient)

saveRDS(KB.N, file = "~/data/pod/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_SC-NEU_Obj.rds")
