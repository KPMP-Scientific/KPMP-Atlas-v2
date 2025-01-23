library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
library("scrattch.hicat")

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

###Prepare Reference Data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
meta <- KB@meta.data
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_Atlas_V2_07-2023_Object_Downsampled_SCT_Reference_B.Rda")
meta <- meta[rownames(KB@meta.data),]

KB@meta.data$v2.subclass.l1 <- meta$v2.subclass.l1
KB@meta.data$v2.subclass.l2 <- meta$v2.subclass.l2
KB@meta.data$v2.subclass.l3 <- meta$v2.subclass.l3

###Prepare Query Data
mmKidAt <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object.Rds")

#Convert to homologous genes (need to do this for both Human and Mouse)
###Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/Projects/Mouse_Kidney/Atlas_V1/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% rownames(KB),]

rownames(hs.ms.table) <- hs.ms.table$Mouse.gene.name

counts <- mmKidAt[["RNA"]]$counts
counts <- counts[rownames(counts) %in% rownames(hs.ms.table),]
rownames(counts) <- as.character(hs.ms.table[rownames(counts),]$Gene.name)

meta <- mmKidAt@meta.data
mmKidAt <- CreateSeuratObject(counts = counts, project = "Mouse Kidney Atlas V2", min.cells = 3, min.features = 200, 
                          meta.data = meta)


#Free up unused memory
gc(reset = TRUE)


###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KB,
  query = mmKidAt,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  dims = 1:50
)


mmKidAt <- MapQuery(
  anchorset = anchors,
  query = mmKidAt,
  reference = KB,
  refdata = list(
    hs.subclass.l1 = "v2.subclass.l1",
    hs.subclass.l2 = "v2.subclass.l2",
    hs.subclass.l3 = "v2.subclass.l3"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

p1 = DimPlot(mmKidAt, reduction = "ref.umap", group.by = "predicted.hs.subclass.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(mmKidAt, reduction = "ref.umap", group.by = "predicted.hs.subclass.l3", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2

meta <- mmKidAt@meta.data
save(meta, file = "mouse_IRI/Mouse_Kidney_Multiome_Atlas_082024_metadata.rda")
#load("mouse_IRI/Mouse_Kidney_Multiome_Atlas_082024_metadata.rda")

rm(KB)
#Free up unused memory
gc(reset = TRUE)

hist(mmKidAt$predicted.hs.subclass.l1.score)
hist(mmKidAt$predicted.hs.subclass.l2.score)
hist(mmKidAt$predicted.hs.subclass.l3.score)


###add predictions to mouse object
mmKidAt <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object.Rds")
mmKidAt@meta.data <- meta[rownames(mmKidAt@meta.data),]
table(mmKidAt$predicted.hs.subclass.l1)


###Clustering in Pagoda2 - Round 1 for QC
#Filtered Count Matrix from Seurat
countMatrix <- mmKidAt[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#5028 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 5208, maxit=1000)

# Generate K-nearest neighbour graph
#p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

save(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Multiome_Atlas_082024_Pagoda_K200_round1.rda") 
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Multiome_Atlas_082024_Pagoda_K200_round1.rda")

#Add pagoda2 clusters 
#k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap

#mmKidAt[["pagoda_k100_infomap"]] <- k100infomap[rownames(mmKidAt@meta.data)]
mmKidAt[["pagoda_k200_infomap"]] <- k200infomap[rownames(mmKidAt@meta.data)]

#var info and umap from K100 object
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:5208]
VariableFeatures(mmKidAt) <- sn.od.genes
cell.embeddings <- p2$reductions$PCA
mmKidAt[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(mmKidAt))
mmKidAt <- RunUMAP(object = mmKidAt, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                   min.dist = 0.3)
#DimPlot(mmKidAt, group.by = "pagoda_k100_infomap", label = TRUE) + NoLegend()
DimPlot(mmKidAt, group.by = "pagoda_k200_infomap", label = TRUE) + NoLegend()

DimPlot(mmKidAt, group.by = "orig.ident", label = FALSE) 
DimPlot(mmKidAt, group.by = "predicted.hs.subclass.l2", label = TRUE, label.size = 3 ,
        repel = TRUE, raster=FALSE) + NoLegend()


###Identify ambiguous low quality - annotate cluster identities
mmKidAt[["RNA"]] <- as(object = mmKidAt[["RNA"]], Class = "Assay")
mmKidAt <- NormalizeData(mmKidAt)

#Identify max predicted subclass per cluster
Idents(mmKidAt) <- "pagoda_k200_infomap"
Idents(mmKidAt) <- factor(Idents(mmKidAt), levels = 1:length(levels(Idents(mmKidAt))))
levels(Idents(object = mmKidAt)) <- paste("CL", levels(Idents(object = mmKidAt)), sep = "")
celltype <- Idents(object = mmKidAt)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.hs.subclass.l1))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.hs.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$predicted.hs.subclass.l1))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$predicted.hs.subclass.l1)),
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
    colnames(Idents.called) <- c("pred.hs.subclassl1", "pred.hs.subclassl1.Freq","pred.hs.subclassl1.Total", "pred.hs.subclassl2",
                                 "pred.hs.subclassl2.Freq","pred.hs.subclassl2.Total","cluster")
    rownames(Idents.called) <- paste(Idents.called$cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(mmKidAt, celltype)
max.idents[1:50,]
max.idents[51:100,]
max.idents[101:150,]
max.idents[151:200,]


###Check for Ambiguous (mixed identity) clusters, plot markers to confirm multiplets
#CL22 - EC/FIB multiplet
#CL66 - cycling cluster
#CL79 - FIB/PT
#CL83 - PT/EC
#CL84 - TAL/EC

c(22,66,79,83,84)

##Check for multiplets based on overlapping cell type markers
Idents(mmKidAt) <- "pagoda_k200_infomap"
Idents(mmKidAt) <- factor(Idents(mmKidAt), levels = 1:length(levels(Idents(mmKidAt))))

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

DotPlot(mmKidAt, features = unique(epi.markers), idents = c(22,54,62,66,79,83,84,92,94)) + RotatedAxis()
DotPlot(mmKidAt, features = unique(int.markers), idents = c(22,54,62,66,79,83,84,92,94)) + RotatedAxis()

###Remove multiplets and re-cluster
mmKidAt <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object.Rds")
mmKidAt@meta.data <- meta[rownames(mmKidAt@meta.data),]
mmKidAt[["pagoda_k200_infomap"]] <- k200infomap[rownames(mmKidAt@meta.data)]
Idents(mmKidAt) <- "pagoda_k200_infomap"
mmKidAt <- subset(mmKidAt, idents = c(22,79,83,84), invert = TRUE)
saveRDS(mmKidAt, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object_filtered.Rds")


###Clustering in Pagoda2 - Round 2 for Clusters
mmKidAt <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object_filtered.Rds")
#Filtered Count Matrix from Seurat
countMatrix <- mmKidAt[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)

# Adjust the variance
p2$adjustVariance(plot = T, gam.k = 10)
#5027 overdispersed genes

# Calculate a PCA reduction with the number of PCs specified by nPCs
# and using only the n.odgenes overdispersed genes 
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 5027, maxit=1000)

# Generate K-nearest neighbour graph
#p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #Downstream clustering steps run speparately for each k value
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #Downstream clustering steps run speparately for each k value
#p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine') #Downstream clustering steps run speparately for each k value

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

save(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Multiome_Atlas_082024_Pagoda_K200_round2.rda") 
#save(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Multiome_Atlas_082024_Pagoda_K100_round2.rda") 
#save(p2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Multiome_Atlas_082024_Pagoda_K500_round2.rda") 
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Multiome_Atlas_082024_Pagoda_K200_round2.rda")

#var info and umap from K100 object
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:5027]
VariableFeatures(mmKidAt) <- sn.od.genes
cell.embeddings <- p2$reductions$PCA
mmKidAt[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(mmKidAt))
mmKidAt <- RunUMAP(object = mmKidAt, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                   min.dist = 0.3)

#Add pagoda2 clusters 
#k100infomap <- p2$clusters$PCA$infomap
k200infomap <- p2$clusters$PCA$infomap
#k500infomap <- p2$clusters$PCA$infomap

#mmKidAt[["pagoda_k100_infomap"]] <- k100infomap[rownames(mmKidAt@meta.data)]
mmKidAt[["pagoda_k200_infomap"]] <- k200infomap[rownames(mmKidAt@meta.data)]
#mmKidAt[["pagoda_k500_infomap"]] <- k500infomap[rownames(mmKidAt@meta.data)]

#DimPlot(mmKidAt, group.by = "pagoda_k100_infomap", label = TRUE) + NoLegend()
DimPlot(mmKidAt, group.by = "pagoda_k200_infomap", label = TRUE) + NoLegend()
#DimPlot(mmKidAt, group.by = "pagoda_k500_infomap", label = TRUE) + NoLegend()


#update metadata
###Add in experiment metadata
meta <- mmKidAt@meta.data
exp.meta <- read.delim("mouse_IRI/Mouse_Experiment_Metadata_08092024.txt")
emc <- c("library","source","assay","experiment","patient","source_ID",       
         "specimen","injury_condition","condition_level3","condition_level2","condition_level1","condition",      
         "age_months","age_group","sex","genotype","strain","protocol","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$orig.ident, exp.meta$library)]
}
mmKidAt@meta.data <- meta
order <- c("library","nCount_RNA", "nFeature_RNA","percent.er","percent.mt","doubletdetection",
           "source","assay","experiment","patient","source_ID",       
           "specimen","injury_condition","condition_level3","condition_level2","condition_level1","condition",      
           "age_months","age_group","sex","genotype","strain","protocol","tissue_type",
           "predicted.hs.subclass.l1.score", "predicted.hs.subclass.l1", "predicted.hs.subclass.l2.score",
           "predicted.hs.subclass.l2", "predicted.hs.subclass.l3.score", "predicted.hs.subclass.l3",
           "pagoda_k200_infomap")
mmKidAt@meta.data <- mmKidAt@meta.data[,order]



###Add in Gerhardt2023 cell type annotations
gobj <- readRDS("~/hsKidAt/blake_LTS/public/Gerhardt2023/Ki67_SObj_rna_peak_motif_geneactivity_assay.Rds")

Idents(gobj) <- "SampleID"
samples <- c("AKI_GFP+_6months_1","Ctrl_4weeks_2","AKI_GFP+_4weeks_3","AKI_GFP-_4weeks","AKI_GFP+_6months_2","AKI_GFP+_6months_4",
             "AKI_GFP+_6months_3","AKI_GFP+_4weeks_2","Ctrl_4weeks_1","AKI_GFP-_6months","Ctrl_6months_1","AKI_GFP+_4weeks_1")
libs <- c("AKIGFPpos6mo1","Cont4wk2","AKIGFPpos4wk3","AKIGFPneg4wk","AKIGFPpos6mo2","AKIGFPpos6mo4",
          "AKIGFPpos6mo3","AKIGFPpos4wk2","Cont4wk1","AKIGFPneg6mo","Cont6mo","AKIGFPpos4wk1")
Idents(gobj) <- plyr::mapvalues(Idents(gobj), from = samples, to = libs)
gobj$library <- Idents(gobj)

gobj.meta <- gobj@meta.data
rownames(gobj.meta) <- paste(gobj.meta$library, rownames(gobj.meta), sep = "_")
rownames(gobj.meta) <- sub("-.*", "-1", rownames(gobj.meta))
gobj.meta <- gobj.meta[rownames(gobj.meta) %in% rownames(mmKidAt@meta.data),]

meta <- mmKidAt@meta.data
meta$Gerhardt2023.Celltype <- "NA"
meta[rownames(gobj.meta),]$Gerhardt2023.Celltype <- as.character(gobj.meta$Celltype)
mmKidAt@meta.data$Gerhardt2023.Celltype <- meta$Gerhardt2023.Celltype
table(mmKidAt@meta.data$Gerhardt2023.Celltype)
DimPlot(mmKidAt, group.by = "Gerhardt2023.Celltype", cells = rownames(gobj.meta),
        label = TRUE, raster=FALSE) + NoLegend()

DimPlot(mmKidAt, group.by = "library", label = FALSE) 
DimPlot(mmKidAt, group.by = "predicted.hs.subclass.l3", label = TRUE, label.size = 3 ,
        repel = TRUE, raster=FALSE) + NoLegend()
DimPlot(mmKidAt, group.by = "condition_level3", label = TRUE, label.size = 3 ,
        repel = TRUE, raster=FALSE) 

meta <- mmKidAt@meta.data
save(meta, file = "Mouse_Kidney_Multiome_Atlas_08092024_metadata.rda")

saveRDS(mmKidAt, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object_filtered.Rds")






###Annotate to Subclass l1
mmKidAt <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object_filtered.Rds")
#Identify max predicted subclass per cluster
Idents(mmKidAt) <- "pagoda_k200_infomap"
Idents(mmKidAt) <- factor(Idents(mmKidAt), levels = 1:length(levels(Idents(mmKidAt))))
levels(Idents(object = mmKidAt)) <- paste("CL", levels(Idents(object = mmKidAt)), sep = "")
celltype <- Idents(object = mmKidAt)

getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.hs.subclass.l1))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$Gerhardt2023.Celltype))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$predicted.hs.subclass.l1))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$predicted.hs.subclass.l1)),
                                query = prop.table(table(seurat.obj.sub$Gerhardt2023.Celltype))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$Gerhardt2023.Celltype)))
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
    colnames(Idents.called) <- c("pred.hs.subclassl1", "pred.hs.subclassl1.Freq","pred.hs.subclassl1.Total", "Gerhardt.celltype",
                                 "Gerhardt.celltype.Freq","Gerhardt.celltype.Total","cluster")
    rownames(Idents.called) <- paste(Idents.called$cluster, rownames(Idents.called), sep = ".")
    Idents.called
  }))
  return(max.idents)
}

max.idents <- getMaxIdents(mmKidAt, celltype)
max.idents[1:50,]
max.idents[51:100,]
max.idents[101:150,]
max.idents[151:200,]


#check marker gene expression and reassign above as needed
Idents(mmKidAt) <- "pagoda_k200_infomap"
Idents(mmKidAt) <- factor(Idents(mmKidAt), levels = 1:length(levels(Idents(mmKidAt))))
mmKidAt.sub <- subset(mmKidAt, downsample = 1000)
mmKidAt.sub[["RNA"]] <- as(object = mmKidAt.sub[["RNA"]], Class = "Assay")
mmKidAt.sub <- NormalizeData(mmKidAt.sub)
DotPlot(mmKidAt.sub, features = unique(epi.markers), idents = c(1:20)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(int.markers), idents = c(1:20)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(epi.markers), idents = c(21:40)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(int.markers), idents = c(21:40)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(epi.markers), idents = c(41:60)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(int.markers), idents = c(41:60)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(epi.markers), idents = c(61:88)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(int.markers), idents = c(61:88)) + RotatedAxis()


#CL1 - IC
#CL2 - PT
#CL3 - PT
#CL4 - TAL
#CL5 - PT
#CL6 - PC
#CL7 - DCT
#CL8 - CNT
#CL9 - DTL
#CL10 - PT
#CL11 - PC
#CL12 - PT
#CL13 - PT
#CL14 - IC
#CL15 - PT
#CL16 - PT
#CL17 - PT
#CL18 - FIB
#CL19 - EC
#CL20 - TAL
#CL21 - PT
#CL22 - PT
#CL23 - EC
#CL24 - TAL
#CL25 - ATL
#CL26 - EC
#CL27 - Lymphoid
#CL28 - Myeloid
#CL29 - CNT
#CL30 - TAL
#CL31 - PC
#CL32 - DCT
#CL33 - DTL
#CL34 - PT
#CL35 - PT
#CL36 - TAL
#CL37 - IC
#CL38 - PT
#CL39 - PT
#CL40 - PapE
#CL41 - DTL
#CL42 - IC
#CL43 - POD
#CL44 - PEC
#CL45 - VSM/P
#CL46 - DCT
#CL47 - FIB
#CL48 - PT
#CL49 - Myeloid
#CL50 - PC
#CL51 - PT
#CL52 - Myeloid
#CL53 - FIB
#CL54 - PT
#CL55 - DCT
#CL56 - Myeloid
#CL57 - DTL
#CL58 - PT
#CL59 - TAL
#CL60 - Myeloid
#CL61 - PT
#CL62 - cycling
#CL63 - PT
#CL64 - TAL
#CL65 - Lymphoid
#CL66 - cycling
#CL67 - PT
#CL68 - PT
#CL69 - Lymphoid
#CL70 - PT
#CL71 - PT
#CL72 - PT
#CL73 - Myeloid
#CL74 - PT
#CL75 - PT
#CL76 - PT
#CL77 - PT
#CL78 - PT
#CL79 - TAL
#CL80 - PT
#CL81 - PT
#CL82 - PT
#CL83 - PT
#CL84 - PT
#CL85 - EC
#CL86 - NA (ec multiplet)
#CL87 - NA (ec multiplet)
#CL88 - TAL

#Add subclasses
current.ident <- c(1:88)
new.ident <- c("IC","PT","PT","TAL","PT","PC","DCT","CNT","DTL","PT","PC","PT",
  "PT","IC","PT","PT","PT","FIB","EC","TAL","PT","PT","EC","TAL","ATL","EC",
  "Lymphoid","Myeloid","CNT","TAL","PC","DCT","DTL","PT","PT","TAL","IC","PT",
  "PT","PapE","DTL","IC","POD","PEC","VSM/P","DCT","FIB","PT","Myeloid","PC","PT",
  "Myeloid","FIB","PT","DCT","Myeloid","DTL","PT","TAL","Myeloid","PT","cycling",
  "PT","TAL","Lymphoid","cycling","PT","PT","Lymphoid","PT","PT","PT","Myeloid",
  "PT","PT","PT","PT","PT","TAL","PT","PT","PT","PT","PT","EC","NA","NA","TAL"
)
Idents(mmKidAt) <- "pagoda_k200_infomap"
Idents(mmKidAt) <- plyr::mapvalues(Idents(mmKidAt), from = current.ident, to = new.ident)
mmKidAt$v2.subclass.l1 <- Idents(mmKidAt)
table(Idents(mmKidAt))

#Remove doublets
mmKidAt <- subset(mmKidAt, idents = "NA", invert = TRUE)

DimPlot(mmKidAt, group.by = "v2.subclass.l1", label = TRUE) 


###Subset and re-annotate cycling cells
mmKidAt.sub <- subset(mmKidAt, idents = "cycling")
countMatrix <- mmKidAt.sub[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#92 overdispersed genes
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 200, maxit=1000)
p2$makeKnnGraph(k = 50, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k50infomap <- p2$clusters$PCA$infomap
mmKidAt.sub[["pagoda_k50_infomap.cyc"]] <- k50infomap[rownames(mmKidAt.sub@meta.data)]
cell.embeddings <- p2$reductions$PCA
mmKidAt.sub[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(mmKidAt.sub))
mmKidAt.sub <- RunUMAP(object = mmKidAt.sub, reduction = "pca", dims = 1:50, n.neighbors = 30L,
                   min.dist = 0.3)
DimPlot(mmKidAt.sub, group.by = "pagoda_k50_infomap.cyc", label = TRUE) + NoLegend()

Idents(mmKidAt.sub) <- "pagoda_k50_infomap.cyc"
DotPlot(mmKidAt.sub, features = unique(epi.markers)) + RotatedAxis()
DotPlot(mmKidAt.sub, features = unique(int.markers)) + RotatedAxis()

#CL1 - PT
#CL2 - Myeloid
#CL3 - TAL
#CL4 - EC

PT <- names(Idents(mmKidAt.sub)[Idents(mmKidAt.sub) %in% c(1)])
EC <- names(Idents(mmKidAt.sub)[Idents(mmKidAt.sub) %in% c(4)])
Myeloid <- names(Idents(mmKidAt.sub)[Idents(mmKidAt.sub) %in% c(2)])
TAL <- names(Idents(mmKidAt.sub)[Idents(mmKidAt.sub) %in% c(3)])

mmKidAt <- SetIdent(object = mmKidAt, cells = PT, value = "PT")
mmKidAt <- SetIdent(object = mmKidAt, cells = EC, value = "EC")
mmKidAt <- SetIdent(object = mmKidAt, cells = Myeloid, value = "Myeloid")
mmKidAt <- SetIdent(object = mmKidAt, cells = TAL, value = "TAL")

table(Idents(mmKidAt))

mmKidAt$v2.subclass.l1 <- Idents(mmKidAt)
saveRDS(mmKidAt, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object_filtered_B.Rds")

