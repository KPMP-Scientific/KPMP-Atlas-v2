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

setwd("~/data/net/home/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")

###Prepare Reference Data
load("~/data/pod/blake_LTS/Atlas_V1/Kidney_KPMP-Biopsy_Ref_10X-R_092020_Seurat_G_SCT.rda")
colnames(KBR@reductions$pca@feature.loadings) <- paste0("pca_", 1:50)
sn.od.genes <- VariableFeatures(KBR)

###Prepare Query Data (reconfigure seurat object of V1, V2-R and V2-M layers?)
KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")
KB[["RNA"]]$data <- NormalizeData(KB[["RNA"]]$counts)

###Find anchor genes
anchors <- FindTransferAnchors(
  reference = KBR,
  query = KB,
  normalization.method = "SCT",
  k.filter = NA,
  reference.reduction = "pca",
  reference.neighbors = "pca.annoy.neighbors",
  dims = 1:50
)

KB <- MapQuery(
  anchorset = anchors,
  query = KB,
  reference = KBR,
  refdata = list(
    class = "class",
    subclass.l1 = "subclass.l1",
    subclass.l2 = "subclass.l2",
    subclass.l3 = "subclass.l3"
  ),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(KB, reduction = "ref.umap", group.by = "predicted.subclass.l1", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(KB, reduction = "ref.umap", group.by = "predicted.subclass.l3", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()

#remove prediction assays
KB[["prediction.score.class"]] <- NULL
KB[["prediction.score.subclass.l1"]] <- NULL
KB[["prediction.score.subclass.l2"]] <- NULL
KB[["prediction.score.subclass.l3"]] <- NULL

saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_04-2023_Object.Rds",
  destdir = "~/tmp_scratch/kidney_object"
)


###Add in v1 labels to v1 data
KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")
meta <- KB@meta.data
metav1 <- readRDS("Atlas_V1_Metadata_Tables.RDS")

#Check overlap of atlas v1 metadata
length(rownames(meta)[rownames(meta) %in% rownames(metav1)])
metav1.sub <- metav1[!rownames(metav1) %in% rownames(meta), ]
table(metav1.sub$library)
table(metav1$library)[names(table(metav1.sub$library))]
#Needing barcode fixes: dKC1 dKC2 dKC3 KB1 KB2 KB3 KB4 KB5 KB6 KC33 KC34 KC35 KC36

metav1$barcode <- rownames(metav1)
metav1$corrected.barcode <- rownames(metav1)
rownames(meta[meta$orig.ident == "dKC1",])
metav1[metav1$library == "dKC1",]$corrected.barcode <- paste0(metav1[metav1$library == "dKC1",]$barcode,"-1")
rownames(meta[meta$orig.ident == "dKC2",])
metav1[metav1$library == "dKC2",]$corrected.barcode <- paste0(metav1[metav1$library == "dKC2",]$barcode,"-1")
rownames(meta[meta$orig.ident == "dKC3",])
metav1[metav1$library == "dKC3",]$corrected.barcode <- paste0(metav1[metav1$library == "dKC3",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KB1",])
metav1[metav1$library == "KB1",]$corrected.barcode <- paste0(metav1[metav1$library == "KB1",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KB2",])
metav1[metav1$library == "KB2",]$corrected.barcode <- paste0(metav1[metav1$library == "KB2",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KB3",])
metav1[metav1$library == "KB3",]$corrected.barcode <- paste0(metav1[metav1$library == "KB3",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KB4",])
metav1[metav1$library == "KB4",]$corrected.barcode <- paste0(metav1[metav1$library == "KB4",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KB5",])
metav1[metav1$library == "KB5",]$corrected.barcode <- paste0(metav1[metav1$library == "KB5",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KB6",])
metav1[metav1$library == "KB6",]$corrected.barcode <- paste0(metav1[metav1$library == "KB6",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KC33",])
metav1[metav1$library == "KC33",]$corrected.barcode <- paste0(metav1[metav1$library == "KC33",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KC34",])
metav1[metav1$library == "KC34",]$corrected.barcode <- paste0(metav1[metav1$library == "KC34",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KC35",])
metav1[metav1$library == "KC35",]$corrected.barcode <- paste0(metav1[metav1$library == "KC35",]$barcode,"-1")
rownames(meta[meta$orig.ident == "KC36",])
metav1[metav1$library == "KC36",]$corrected.barcode <- paste0(metav1[metav1$library == "KC36",]$barcode,"-1")
rownames(metav1) <- metav1$corrected.barcode

###Update meta tables
meta$clusters <- metav1$clusters[match(rownames(meta), rownames(metav1))]
meta$subclass.l3 <- metav1$subclass.l3[match(rownames(meta), rownames(metav1))]
meta$subclass.l2 <- metav1$subclass.l2[match(rownames(meta), rownames(metav1))]
meta$subclass.l1 <- metav1$subclass.l1[match(rownames(meta), rownames(metav1))]
meta$class <- metav1$class[match(rownames(meta), rownames(metav1))]
meta$v1_barcode <- metav1$barcode[match(rownames(meta), rownames(metav1))]

KB@meta.data <- meta

DimPlot(KB, reduction = "ref.umap", group.by = "clusters", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(KB, reduction = "ref.umap", group.by = "subclass.l1", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()
DimPlot(KB, reduction = "ref.umap", group.by = "predicted.subclass.l1", alpha = 0.1, label = TRUE, #split.by = "assay", ncol = 3, 
        label.size = 3) + NoLegend()

table(KB$predicted.subclass.l1)

#Free up unused memory
rm(KBR)
gc(reset = TRUE)

saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_04-2023_Object.Rds",
  destdir = "~/tmp_scratch/kidney_object"
)



###Clustering in Pagoda2 using integrated PCs
KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")
VariableFeatures(KB)

#Filtered Count Matrix from Seurat
countMatrix <- KB[["RNA"]]$counts[rownames(KB[["RNA"]]$counts) %in% sn.od.genes,]
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)
#28043

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Adjust the variance
#p2$adjustVariance(plot = T, gam.k = 10)
#2393 overdispersed genes

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["ref.pca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 8, distance = 'cosine')

# Identify clusters using the infomap.community method
# on the basis of the reduction called 'PCA' (generated above)
# Save the resulting clustering as 'infomap'
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

save(sn.od.genes, p2, file = "Kidney_V1-V2_Integration_p2_Object_04172023.rda") #k100
save(sn.od.genes, p2, file = "Kidney_V1-V2_Integration_p2_Object_04172023_K500.rda") #k500
#load("Kidney_V1-V2_Integration_p2_Object_04172023.rda")


#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
KB[["pagoda_k100_infomap"]] <- k100infomap[rownames(KB@meta.data)]
k500infomap <- p2$clusters$PCA$infomap
#KB[["pagoda_k500_infomap"]] <- k500infomap[rownames(KB@meta.data)]

KB <- RunUMAP(object = KB, reduction = "ref.pca", dims = 1:50, n.neighbors = 30L,
                    min.dist = 0.2)

DimPlot(KB, group.by = "pagoda_k100_infomap", reduction = "umap", label = TRUE) + NoLegend()
#DimPlot(KB, group.by = "pagoda_k500_infomap", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "clusters", reduction = "umap", label = TRUE) + NoLegend()
table(KB$pagoda_k100_infomap)

saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_04-2023_Object.Rds",
  destdir = "~/tmp_scratch/kidney_object"
)



###Add in experiment metadata
meta <- KB@meta.data
exp.meta <- read.delim("Kidney_AtlasV2_Experiment_Metadata.txt")
colnames(meta)[colnames(meta) == "orig.ident"] <- "library"
emc <- c("source","assay","experiment_long","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "percent_cortex","percent_medulla","region_level2","region_level1",
         "age","age_binned","sex","race","location","laterality","protocol",
         "tissue_type_full","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
KB@meta.data <- meta
saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_04-2023_Object.Rds",
  destdir = "~/tmp_scratch/kidney_object"
)





###Annotate
###Annotate based on sn subclass l3 (including predicted subclass and top corr values)
#Correlate query clusters to reference cell types using variable genes
KB <- readRDS("~/tmp_scratch/kidney_object/Kidney_Atlas_V2_04-2023_Object.Rds")
VariableFeatures(KB)
marker.genes = sn.od.genes
KB <- ScaleData(KB, features = marker.genes, assay = "RNA")
Idents(KB) <- "atlas_version"
ref <- subset(KB, idents = "V1")
Idents(ref) <- "subclass.l3"
query <- subset(KB, idents = "V2")
Idents(query) <- "pagoda_k100_infomap"
levels(Idents(object = query)) <- paste("CL", levels(Idents(object = query)), sep = "")
ave.query <- AverageExpression(query, features = marker.genes, assays = "RNA",
                               slot = "scale.data")
ave.ref <- AverageExpression(ref, features = marker.genes, assays = "RNA",
                             slot = "scale.data")
library("corrplot")
ave.cor <- cor(as.data.frame(ave.ref),as.data.frame(ave.query))
ave.cor
colnames(ave.cor) <- gsub("RNA.","",colnames(ave.cor))

#Generate table
Idents(KB) <- "pagoda_k100_infomap"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(KB)) <- paste("CL", levels(Idents(KB)), sep = "")
celltype <- Idents(KB)

tail(sort(ave.cor[,levels(celltype)[1]]), 2)
data.frame(ref.cor.name = tail(names(sort(ave.cor[,levels(celltype)[1]])), 2),
           ref.cor.value = tail(sort(ave.cor[,levels(celltype)[1]]), 2))



getMaxIdents <- function(seurat.obj, celltype) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    ref.top2 <- tail(names(sort(table(seurat.obj.sub$clusters))), 2)
    query.top2 <- tail(names(sort(table(seurat.obj.sub$predicted.subclass.l3))), 2)
    max.len = 2
    ref.top2  = c(ref.top2 , rep(NA, max.len - length(ref.top2 )))
    query.top2 = c(query.top2 , rep(NA, max.len - length(query.top2 )))
    Idents.called <- data.frame(ref = prop.table(table(seurat.obj.sub$clusters))[ref.top2],
                                ref.Total = sum(table(seurat.obj.sub$clusters)),
                                query = prop.table(table(seurat.obj.sub$predicted.subclass.l3))[query.top2],
                                query.Total = sum(table(seurat.obj.sub$predicted.subclass.l3)))
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

getMaxCor <- function(ave.cor, celltype) {
  max.cor <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    
    if (ct %in% colnames(ave.cor) == FALSE) {
      max.cor <- data.frame(ref.cor.name = c(NA,NA),
                            ref.cor.value = c(NA,NA))
    } else {
      max.cor <- data.frame(ref.cor.name = tail(names(sort(ave.cor[,ct])), 2),
                            ref.cor.value = tail(sort(ave.cor[,ct]), 2))
    } 
    
  }))
  return(max.cor)
}

max.cor <- getMaxCor(ave.cor, celltype)
max.idents$cor.subclass <- max.cor$ref.cor.name
max.idents$cor.value <- max.cor$ref.cor.value

#Write out table
write.table(max.idents, file="Intermediate_Objects/Kidney_Integrated_cluster_annotation_overlaps_04182023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)





### Check max cluster correlations for merging
Idents(ref) <- "clusters"
Idents(ref) <- factor(Idents(ref), levels = 1:119)
ave.ref <- AverageExpression(ref, features = marker.genes, assays = "RNA",
                             slot = "scale.data")
ave.cor <- cor(as.data.frame(ave.ref),as.data.frame(ave.query))
ave.cor
colnames(ave.cor) <- gsub("RNA.","",colnames(ave.cor))
rownames(ave.cor) <- gsub("RNA.","",rownames(ave.cor))

#Generate table
max.cor <- getMaxCor(ave.cor, celltype)

#Write out table
write.table(max.cor, file="Intermediate_Objects/Kidney_Integrated_cluster_annotation_overlaps_04182023_cl-cor.txt", sep = "\t", row.names=TRUE, col.names=TRUE)





###Calculate scores
#degen.score
load("misc/Combined_Degen-State_Conserved_Conditionl1_Markers_Injury-score-set.rda")
count.data = LayerData(KB, layer = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(degen.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KB$degen.score <- aaa[rownames(KB@meta.data)]

#aEpi.score
load("misc/aEpi_state_Conserved_Conditionl1_Markers_score-set.rda")
count.data = GetAssayData(KB, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(adap.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KB$aEpi.score <- aaa[rownames(KB@meta.data)]

#aStr.score
load("misc/aStr_State_Conserved_Conditionl1_Markers_score-set.rda")
count.data = GetAssayData(KB, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(alt.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KB$aStr.score <- aaa[rownames(KB@meta.data)]

#cyc.score
load("misc/Cycling_State_Markers_Combined_score-set.rda")
count.data = GetAssayData(KB, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(cyc.state.genes) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KB$cyc.score <- aaa[rownames(KB@meta.data)]

#matrisome.score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Division %in% c("Core matrisome"))]
count.data = GetAssayData(KB, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KB$matrisome.score <- aaa[rownames(KB@meta.data)]

#collagen.score
matrisome_set = read.table("https://raw.githubusercontent.com/mahmoudibrahim/KidneyMap/master/assets/public/ecm_genes_human.txt", sep = "\t", header = TRUE)
matrisome = matrisome_set$Gene.Symbol[which(matrisome_set$Category %in% c("Collagens"))]
count.data = GetAssayData(KB, slot = "data")
count.data = as.matrix(count.data[rownames(count.data) %in% toupper(matrisome) ,])
combined_matrix = rbind(count.data,average=apply(count.data,2,mean))
cor_matrix = cor(t(combined_matrix))
cor_vector = cor_matrix[,dim(cor_matrix)[1]]
count.data_restricted = count.data[rownames(count.data) %in% names(cor_vector[cor_vector >= 0.1]),]
aaa = apply(count.data_restricted,2,mean)
KB$collagen.score <- aaa[rownames(KB@meta.data)]


###Add assigned v1 cluster IDs/subclasses for each integrated cluster
assigned <- read.delim("Intermediate_Objects/Kidney_Integrated_cluster_annotation_lookup_table_04182023.txt")
KB@meta.data$v1clusters <- assigned$v1.cl[match(KB@meta.data$pagoda_k100_infomap, assigned$int.cl)]

Idents(KB) <- "atlas_version"
query <- subset(KB, idents = "V2")
Idents(query) <- "pagoda_k100_infomap"
Idents(query) <- factor(Idents(query), levels = 1:length(levels(Idents(query))))

pdf(file='Intermediate_Objects/Kidney_v1_integration_scores_04192023.pdf',width=50,height=24)
VlnPlot(object = query, features = c("degen.score", "aEpi.score","aStr.score","cyc.score",
                                     "percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)
dev.off()


###Check marker genes
sum(table(KB$v1clusters))
##POD, PEC
Idents(KB) <- "pagoda_k100_infomap"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

KB.sub <- subset(KB, subset = v1clusters %in% c(1:3))
DotPlot(KB.sub, features = c("NPHS1","NPHS2","PODXL",                                #POD
                         "CDKN1C","SPOCK2",                                          #dPOD
                         "CLDN1", "CFH","ALDH1A2",                             #PEC
                               
                         "LRP2",                                               #PT
                         "SLC5A12","SLC22A6","SLC34A1",                        #S1/S2
                         "SLC5A11","SLC7A13",                                  #S3
                         "CRYAB","TACSTD2","SLC44A5",                          #TL
                         "FABP1","EMCN","PECAM1","C7","ACTA2","DCN",           #EC/stromal
                         "PTPRC","MRC1","BANK1","JCHAIN","IL7R",               #Immune
                         "MS4A2","S100A9",
                         "APOE","DEFB1","CST3","GATM","ALDOB",              #Injury
                         "UMOD","SLC12A1","SLC12A3","SLC8A1","AQP2",            #Distal tubules
                          "HSD11B2","SLC26A7"
), split.by = "atlas_version",
cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()


cells <- WhichCells(KB, idents = 202)
DimPlot(KB, cells.highlight = cells, reduction = "umap", label = FALSE) + NoLegend()



##PT
KB.sub <- subset(KB, subset = v1clusters %in% c(4:19))
length(levels(Idents(KB.sub)))
#91
set1 <- levels(Idents(KB.sub))[1:30]
set2 <- levels(Idents(KB.sub))[31:60]
set3 <- levels(Idents(KB.sub))[61:91]
KB.sub <- subset(KB, idents = set1)
pt.markers <- c("SLC5A12","SLC22A6",                                  #S1/S2
                "PRODH2","SLC5A2","SLC22A8","SLC7A8",                 #S1
                "SLC34A1","SLC5A10",                                  #S2                                  
                "SLC5A11","SLC7A13",                                  #S3
                
                "ITGB8","CDH6","HAVCR1","VCAM1",                      #AS
                
                "PODXL","CLDN1", "CFH", 
                "CRYAB","TACSTD2","SLC44A5", "SH3GL3",                #TL
                "EMCN","PECAM1","PDGFRB","ACTA2","C7","DCN","COL1A1", #EC/stromal
                "PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
                "CD96","NKG7",       
                "MS4A2","CPA3", "KIT","S100A9",
                "APOE","DEFB1","CST3","GATM","ALDOB",                 #Injury
                "UMOD","SLC12A1","SLC12A3","SLC8A1","AQP2",            #Distal tubules
                "HSD11B2","SLC26A7",
                "MKI67","TOP2A")

DotPlot(KB.sub, features = pt.markers,split.by = "atlas_version",
cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

KB.sub <- subset(KB, idents = set2)
DotPlot(KB.sub, features = pt.markers,split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

KB.sub <- subset(KB, idents = set3)
DotPlot(KB.sub, features = pt.markers,split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()


VlnPlot(object = KB.sub, features = c("degen.score", "aEpi.score",
                                     "percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)

cells <- WhichCells(KB, idents = 257)
DimPlot(KB, cells.highlight = cells, reduction = "umap", label = FALSE) + NoLegend()
DimPlot(KB, group.by = "predicted.subclass.l1", reduction = "umap", label = TRUE) + NoLegend()



#TL
KB.sub <- subset(KB, subset = v1clusters %in% c(20:25))
length(levels(Idents(KB.sub)))
#27

loh.markers <- c(
"CRYAB","TACSTD2","SLC44A5",                          #TL
"AQP1", "UNC5D",                                     #DTL1
"ADGRL3","ID1",                                      #DTL2
"SMOC2",                                             #DTL3
"AKR1B1","SH3GL3",                                   #DTL3/ATL
"PROX1",                                             #ATL

"ITGB6","CDH6","HAVCR1","VCAM1","PROM1",              #AS

"SLC5A12","SLC22A6",                                  #PT
"PRODH2","SLC34A1","SLC7A13",
"PODXL","CLDN1", "CFH", 

"EMCN","PECAM1","PDGFRB","ACTA2","C7","DCN","COL1A1", #EC/stromal
"PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
"CD96","NKG7",       
"MS4A2","CPA3", "KIT","S100A9",
"APOE","DEFB1","CST3","SERPINA1","S100A4",         #Injury
"UMOD","SLC12A1","SLC12A3","SLC8A1","AQP2",            #Distal tubules
"HSD11B2","SLC26A7",
"MKI67","TOP2A"
)

DotPlot(KB.sub, features = loh.markers,split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()


VlnPlot(object = KB.sub, features = c("degen.score", "aEpi.score",
                                      "percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)



##TAL
KB.sub <- subset(KB, subset = v1clusters %in% c(26:39))
length(levels(Idents(KB.sub)))
#83
set1 <- levels(Idents(KB.sub))[1:27]
set2 <- levels(Idents(KB.sub))[28:54]
set3 <- levels(Idents(KB.sub))[55:83]

loh.markers <- c(
                 "CASR","SLC12A1","UMOD","EGF",                   #TAL
                 "CLDN14", "KCTD16","ANK2",                       #M-TAL
                 "ENOX1","TMEM207","CLDN16",                      #C-TAL
                 "NOS1","ROBO2",                                  #MD
                 
                 "CRYAB","TACSTD2","SLC44A5",                          #TL
                 "UNC5D","ADGRL3","AKR1B1","SH3GL3",
                 
                 "ITGB6","CDH6","HAVCR1","VCAM1","PROM1",              #AS
                 
                 "SLC5A12","SLC22A6",                                  #PT
                 "PRODH2","SLC34A1","SLC7A13",
                 "PODXL","CLDN1", "CFH", 
                 "EMCN","PECAM1","PDGFRB","ACTA2","C7","DCN","COL1A1", #EC/stromal
                 "PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
                 "CD96","NKG7",       
                 "MS4A2","CPA3", "S100A9",
                 "APOE","DEFB1","CST3","SPP1","WFDC2",                #Injury
                 "GATA3","PCDH7","AQP2","SLC14A2",                           #Distal tubules
                 "HSD11B2","SLC26A7",
                 "MKI67","TOP2A"
                 )

KB.sub <- subset(KB, idents = set1)
DotPlot(KB.sub, features = loh.markers,split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

KB.sub <- subset(KB, idents = set2)
DotPlot(KB.sub, features = loh.markers,split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

KB.sub <- subset(KB, idents = set3)
DotPlot(KB.sub, features = loh.markers,split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()


VlnPlot(object = KB.sub, features = c("degen.score", "aEpi.score",
                                      "percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)


cells <- WhichCells(KB, idents = 340)
DimPlot(KB, cells.highlight = cells, reduction = "umap", label = FALSE) + NoLegend()
DimPlot(KB, group.by = "predicted.subclass.l1", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

table(KB.sub$pagoda_k100_infomap,KB.sub$region_level1)
table(query$pagoda_k100_infomap_coembed,query$region.l2)
table(KB.sub$pagoda_k100_infomap,KB.sub$patient)[c(420),]
table(KB.sub$pagoda_k100_infomap,KB.sub$region_level1)[436,]

table(query$pagoda_k100_infomap_coembed,query$condition.l2)


##Distal Tubules
KB.sub <- subset(KB, subset = v1clusters %in% c(40:61))
length(levels(Idents(KB.sub)))
#104
set1 <- levels(Idents(KB.sub))[1:34]
set2 <- levels(Idents(KB.sub))[35:68]
set3 <- levels(Idents(KB.sub))[69:104]



dt.markers <- c(
  "SLC12A3","TRPM6",                                   #DCT
  "ADAMTS17",                                          #DCT1
  "SLC8A1","SCN2A","HSD11B2","CALB1",                  #DCT2 / CNT
  "PCDH7",                                             #CNT
  "SCNN1G","SCNN1B",                                   #CNT                   
  "GATA3","AQP2","AQP3",                               #PC
  "SLC4A1", "SLC26A7",                                 #IC
  "SLC4A9", "SLC26A4",                                 #IC-B
  "SLC14A2","FXYD4",                                   #IMCD
  "TP63", "KRT5",                                      #PapE
  
  
  "PODXL","CLDN1", "CFH",                          #POD/PEC
  "SLC5A12","SLC22A6",                                  #PT
  "PRODH2","SLC34A1","SLC7A13",
  "CRYAB","TACSTD2","SLC44A5",                          #TL
  "AKR1B1","SH3GL3",
  "CASR","SLC12A1","UMOD",                         #TAL
  
  "EMCN","PECAM1","PDGFRB","ACTA2","C7","DCN","COL1A1", #EC/stromal
  "PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
  "CD96","CD247", "THEMIS","NKG7",       
  "MS4A2","CPA3", "S100A9",
  "APOE","DEFB1","CST3","WFDC2",                #Injury
                            
  "MKI67","TOP2A"
)

KB.sub <- subset(KB, idents = set1)
DotPlot(KB.sub, features = dt.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

KB.sub <- subset(KB, idents = set2)
DotPlot(KB.sub, features = dt.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

KB.sub <- subset(KB, idents = set3)
DotPlot(KB.sub, features = dt.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()



VlnPlot(object = KB.sub, features = c("degen.score", "aEpi.score",
                                      "percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)

VlnPlot(object = KB.sub, idents = c(268), features = c("SLC8A1","SLC14A2","FXYD4"), 
        ncol = 1, pt.size = -1)


cells <- WhichCells(KB, idents = 340)
DimPlot(KB, cells.highlight = cells, reduction = "umap", label = FALSE) + NoLegend()
DimPlot(KB, group.by = "predicted.subclass.l1", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "predicted.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

table(KB.sub$pagoda_k100_infomap,KB.sub$region_level1)
table(query$pagoda_k100_infomap_coembed,query$region.l2)
table(KB.sub$pagoda_k100_infomap,KB.sub$patient)[c(420),]
table(KB.sub$pagoda_k100_infomap,KB.sub$region_level1)[116,]

table(query$pagoda_k100_infomap_coembed,query$condition.l2)




#EC-FIB
KB.sub <- subset(KB, subset = v1clusters %in% c(62:87))
length(levels(Idents(KB.sub)))
#49
set1 <- levels(Idents(KB.sub))[1:24]
set2 <- levels(Idents(KB.sub))[25:49]


st.markers <- c(
  "PECAM1", "CD34",                                     #EC
  "EMCN","HECW2","KDR",                                 #EC-GC
  "BTNL9", "PALMD",                                     #EC-AEA/DVR
  "CEACAM1", "DNASE1L3", "PLVAP",                       #EC-PTC
  "NR2F2",                                              #EC-AVR
  "MMRN1",                                              #EC-Lym
  
  "PDGFRB", "NTRK3", "MYH11",                           #VSMC/P
  "PIP5K1B", "POSTN", "REN",                            #MC
  
  "COL1A1","C7","FBLN5","DCN",                          #FIB 
  "ACTA2","FAP",                                        #MYOF
  "FLRT2", "IGF1",
  
  "PODXL","CLDN1", "CFH",                          #POD/PEC
  "SLC5A12","SLC22A6",                                  #PT
  "PRODH2","SLC34A1","SLC7A13",
  "CRYAB","TACSTD2","SLC44A5",                          #TL
  "AKR1B1","SH3GL3",
  "CASR","SLC12A1","UMOD",                         #TAL
  
  "SLC12A3","TRPM6",                                   #DCT
  "SLC8A1","SCN2A","CALB1",                            #DCT2 / CNT
  "PCDH7",                                             #CNT
  "GATA3","AQP2","AQP3",                               #PC
  "SLC26A7","SLC26A4",                                 #IC
  "SLC14A2","FXYD4",                                   #IMCD
  "TP63",                                              #PapE
  
  
  "PTPRC","MRC1","CD163","MS4A1","JCHAIN","IL7R",       #Immune
  "CD96","CD247", "THEMIS","NKG7",       
  "MS4A2","CPA3", "S100A9",
  "APOE","DEFB1","CST3","WFDC2",                #Injury
  
  "MKI67","TOP2A"
)

KB.sub <- subset(KB, idents = set1)
DotPlot(KB.sub, features = st.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()

KB.sub <- subset(KB, idents = set2)
DotPlot(KB.sub, features = st.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()


VlnPlot(object = KB.sub, features = c("degen.score", "aEpi.score",
                                      "percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)

table(KB.sub$pagoda_k100_infomap,KB.sub$experiment)[c(377),]
table(KB.sub$pagoda_k100_infomap,KB.sub$patient)[c(377),]
table(KB.sub$pagoda_k100_infomap,KB.sub$region_level1)[377,]



##Immune
KB.sub <- subset(KB, subset = v1clusters %in% c(88:100))
length(levels(Idents(KB.sub)))
#21

imm.markers <- c(
  "PTPRC",                                              #Immune
  "BANK1","MS4A1",                                      #B
  "IGKC", "JCHAIN",                                     #PL
  "CD96","CD247","IL7R",                                #T
  "NKG7","GZMA",                                        #NKT
  "MS4A2","CPA3","KIT",                                 #MAST
  "MRC1","CD163","MSR1",                                #MAC
  "ITGAX",                                              #cDC
  "CUX2",                                               #pDC
  "CTSS","FCGR3A",                                      #ncMON
  "S100A9", "FCGR3B",                                   #N
  
  "PODXL","CLDN1", "CFH",                          #POD/PEC
  "SLC5A12","SLC22A6",                                  #PT
  "PRODH2","SLC34A1","SLC7A13",
  "CRYAB","TACSTD2","SLC44A5",                          #TL
  "AKR1B1","SH3GL3",
  "CASR","SLC12A1","UMOD",                         #TAL
  
  "SLC12A3","TRPM6",                                   #DCT
  "SLC8A1","SCN2A","HSD11B2","CALB1",                  #DCT2 / CNT
  "PCDH7",                                             #CNT
  "GATA3","AQP2","AQP3",                               #PC
  "SLC26A7","SLC26A4",                                 #IC
  "SLC14A2","FXYD4",                                   #IMCD
  "TP63",                                              #PapE
  
  "ITGB6","CDH6","HAVCR1","VCAM1","PROM1",             #aEpi
  
  "PECAM1","EMCN","PLVAP","NR2F2","MMRN1",             #EC
  "PDGFRB", "MYH11","POSTN",                           #VSMC/P
  "COL1A1","C7","DCN",                                 #FIB 
  "ACTA2","FAP",                                        #MYOF
  "CDH19",                                             #SC/NEU
  "APOE","DEFB1","CST3","WFDC2",                #Injury
  
  "MKI67","TOP2A"
)

DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()


VlnPlot(object = KB.sub, features = c("degen.score", "aEpi.score",
                                      "percent.er","percent.mt","nFeature_RNA"), 
        ncol = 1, pt.size = -1)

table(KB.sub$pagoda_k100_infomap,KB.sub$experiment)[c(381),]
table(KB.sub$pagoda_k100_infomap,KB.sub$patient)[c(433),]
table(KB.sub$pagoda_k100_infomap,KB.sub$region_level1)[433,]

DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()



##misc
KB.sub <- subset(KB, subset = v1clusters %in% c(101:119))



##Double check ambiguous assigned clusters (assigned 101)
#Add assigned v1 cluster IDs/subclasses for each integrated cluster
assigned <- read.delim("Intermediate_Objects/Kidney_Integrated_cluster_annotation_lookup_table_04212023.txt")
KB@meta.data$v1clusters <- assigned$v1.cl[match(KB@meta.data$pagoda_k100_infomap, assigned$int.cl)]

KB.sub <- subset(KB, subset = v1clusters %in% c(101))
length(levels(Idents(KB.sub)))
#168
set1 <- levels(Idents(KB.sub))[1:25]
set2 <- levels(Idents(KB.sub))[26:50]
set3 <- levels(Idents(KB.sub))[51:75]
set4 <- levels(Idents(KB.sub))[76:100]
set5 <- levels(Idents(KB.sub))[101:125]
set6 <- levels(Idents(KB.sub))[126:150]
set7 <- levels(Idents(KB.sub))[151:168]


KB.sub <- subset(KB, idents = set1)
DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
KB.sub <- subset(KB, idents = set2)
DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
KB.sub <- subset(KB, idents = set3)
DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
KB.sub <- subset(KB, idents = set4)
DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
KB.sub <- subset(KB, idents = set5)
DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
KB.sub <- subset(KB, idents = set6)
DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()
KB.sub <- subset(KB, idents = set7)
DotPlot(KB.sub, features = imm.markers, split.by = "atlas_version",
        cols = c("red", "dodgerblue4"), dot.scale = 8) + RotatedAxis()




###Update metadata
###Add assigned v1 cluster IDs/subclasses for each integrated cluster
assigned <- read.delim("Intermediate_Objects/Kidney_Integrated_cluster_annotation_lookup_table_04212023.txt")
KB@meta.data$v1clusters <- assigned$v1.cl[match(KB@meta.data$pagoda_k100_infomap, assigned$int.cl)]

#Update to actual v1 cluster annotations for v1 data
v1.cells <- rownames(KB@meta.data[KB@meta.data$clusters %in% c(1:119),])
length(v1.cells)
#193880
v2.cells <- rownames(KB@meta.data)[!rownames(KB@meta.data) %in% v1.cells]
length(v2.cells)
#1298752

KB@meta.data[v1.cells,]$v1clusters <- KB@meta.data[v1.cells,]$clusters


metav1 <- readRDS("Atlas_V1_Metadata_Tables.RDS")

###Update meta tables
KB@meta.data$v1.subclass.l3 <- metav1$subclass.l3[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.subclass.l2 <- metav1$subclass.l2[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.subclass.l1 <- metav1$subclass.l1[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.state.l1 <- metav1$state.l1[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.state.l2 <- metav1$state.l2[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.class <- metav1$class[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.structure <- metav1$structure[match(KB@meta.data$v1clusters, metav1$clusters)]

colnames(KB@meta.data)
order <- c("library","nCount_RNA","nFeature_RNA","percent.er","percent.mt","set","v1_barcode",
           "source","assay","experiment_long","experiment","patient","specimen",                   
           "condition_level3","condition_level2","condition_level1","condition","percent_cortex","percent_medulla",
           "region_level2","region_level1","age","age_binned","sex","race","location","laterality","protocol",
           "tissue_type_full","tissue_type","atlas_version","predicted.class.score",      
           "predicted.class","predicted.subclass.l1.score","predicted.subclass.l1","predicted.subclass.l2.score",
           "predicted.subclass.l2","predicted.subclass.l3.score","predicted.subclass.l3",
           "pagoda_k100_infomap","v1clusters","v1.subclass.l3","v1.subclass.l2",
           "v1.subclass.l1","v1.state.l1","v1.state.l2","v1.class","v1.structure")

KB@meta.data <- KB@meta.data[,order]




###Subset REN for query
ren <- WhichCells(object = KB, expression = REN > 3)
Idents(KB) <- "v1.class"
stroma <- WhichCells(object = KB, idents = "stroma cells")
ren <- ren[ren %in% stroma]
ren <- ren[ren %in% v2.cells]

DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = ren) + NoLegend()
Idents(object = KB) <- "v1clusters"
KB <- SetIdent(KB, cells = ren, value = 72)
KB$v1clusters <- Idents(object = KB)




###Subset AEA from DVR
#Pagoda2
library(pagoda2)
require(parallel)
Idents(KB) <- "v1clusters"
KBR.aea <- subset(KB, idents = c(63:64))

countMatrix <- KBR.aea[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#2311 overdispersed genes
p2$calculatePcaReduction(nPcs = 10, n.odgenes = 1000, maxit=1000)
p2$makeKnnGraph(k = 500, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k500infomap <- p2$clusters$PCA$infomap
KBR.aea[["pagoda_k500_infomap_AEA"]] <- k500infomap[rownames(KBR.aea@meta.data)]
cell.embeddings <- p2$reductions$PCA
KBR.aea[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KBR.aea))

KBR.aea <- RunUMAP(object = KBR.aea, reduction = "pca", dims = 1:10, n.neighbors = 30L,
                   min.dist = 0.3)
DimPlot(KBR.aea, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "pagoda_k500_infomap_AEA", repel = TRUE) + NoLegend()
DimPlot(KBR.aea, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "region_level2", repel = TRUE) + NoLegend()
DimPlot(KBR.aea, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
Idents(KBR.aea) <- "pagoda_k500_infomap_AEA"

#Annotate based on overlapping ref annotations
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
marker.genes = rownames(var.info)[1:1000]
ref.KBR.aea <- subset(KBR.aea, cells = v1.cells[v1.cells %in% colnames(KBR.aea)])
Idents(ref.KBR.aea) <- "v1clusters"
ref.KBR.aea <- ScaleData(ref.KBR.aea, features = marker.genes)
query.KBR.aea <- subset(KBR.aea, cells = v2.cells[v2.cells %in% colnames(KBR.aea)])
levels(Idents(object = query.KBR.aea)) <- paste("CL", levels(Idents(object = query.KBR.aea)), sep = "")
query.KBR.aea <- ScaleData(query.KBR.aea, features = marker.genes)

ave.query <- AverageExpression(query.KBR.aea, features = marker.genes, assays = "RNA",
                               slot = "scale.data")
ave.ref <- AverageExpression(ref.KBR.aea, features = marker.genes, assays = "RNA",
                             slot = "scale.data")
library("corrplot")
ave.cor <- cor(as.data.frame(ave.ref),as.data.frame(ave.query))
ave.cor

table(ref.KBR.aea$subclass.l3,ref.KBR.aea$pagoda_k100_infomap_AEA)

#Clusters that overlap and correlate with AEA = 1, 2

AEA <- WhichCells(object = KBR.aea, idents = c(1, 2))
AEA <- AEA[AEA %in% v2.cells]

DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = AEA) + NoLegend()
KB <- SetIdent(KB, cells = AEA, value = 63)
KB$v1clusters <- Idents(object = KB)

Idents(KB) <- factor(Idents(KB), levels = c(1:119))
table(Idents(KB))






###Subset cycling cells
#Pagoda2
KBR.cyc <- subset(KB, idents = c(12,43,47,69,79,94))

countMatrix <- KBR.cyc[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)

## Batch annotations
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#1253 overdispersed genes
p2$calculatePcaReduction(nPcs = 30, n.odgenes = 1253, maxit=1000)
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')
k100infomap <- p2$clusters$PCA$infomap
KBR.cyc[["pagoda_k100_infomap_cyc"]] <- k100infomap[rownames(KBR.cyc@meta.data)]
cell.embeddings <- p2$reductions$PCA
KBR.cyc[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KBR.cyc))

KBR.cyc <- RunUMAP(object = KBR.cyc, reduction = "pca", dims = 1:30, n.neighbors = 30L,
                   min.dist = 0.3)
DimPlot(KBR.cyc, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "pagoda_k100_infomap_cyc", repel = TRUE) + NoLegend()
#DimPlot(KBR.cyc, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "region_level2", repel = TRUE) + NoLegend()
DimPlot(KBR.cyc, reduction = "umap", pt.size = 0.5, label = TRUE, group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()


#Annotate based on overlapping ref annotations
Idents(KBR.cyc) <- "pagoda_k100_infomap_cyc"
table(KBR.cyc@meta.data[rownames(KBR.cyc@meta.data) %in% v1.cells,]$v1.subclass.l3,
      KBR.cyc@meta.data[rownames(KBR.cyc@meta.data) %in% v1.cells,]$pagoda_k100_infomap_cyc)

#cycPT = 1,3,4,5,10
#cycDCT = 7
#cycCNT = 6
#cycEC = 8
#cycMNP = 2
#cycMYOF = 9

cycPT <- WhichCells(object = KBR.cyc, idents = c(1,3,4,5,10))
cycPT <- cycPT[cycPT %in% v2.cells]
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycPT) + NoLegend()

cycDCT <- WhichCells(object = KBR.cyc, idents = c(7))
cycDCT <- cycDCT[cycDCT %in% v2.cells]
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycDCT) + NoLegend()

cycCNT <- WhichCells(object = KBR.cyc, idents = c(6))
cycCNT <- cycCNT[cycCNT %in% v2.cells]
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycCNT) + NoLegend()

cycEC <- WhichCells(object = KBR.cyc, idents = c(8))
cycEC <- cycEC[cycEC %in% v2.cells]
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycEC) + NoLegend()

cycMNP <- WhichCells(object = KBR.cyc, idents = c(2))
cycMNP <- cycMNP[cycMNP %in% v2.cells]
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycMNP) + NoLegend()

cycMYOF <- WhichCells(object = KBR.cyc, idents = c(9))
cycMYOF <- cycMYOF[cycMYOF %in% v2.cells]
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = cycMYOF) + NoLegend()


KB <- SetIdent(KB, cells = cycPT, value = 12)
KB <- SetIdent(KB, cells = cycDCT, value = 43)
KB <- SetIdent(KB, cells = cycCNT, value = 47)
KB <- SetIdent(KB, cells = cycEC, value = 69)
KB <- SetIdent(KB, cells = cycMYOF, value = 79)
KB <- SetIdent(KB, cells = cycMNP, value = 94)
KB$v1clusters <- Idents(object = KB)

Idents(KB) <- factor(Idents(KB), levels = c(1:119))
table(Idents(KB))



###Update meta tables
KB@meta.data$v1.subclass.l3 <- metav1$subclass.l3[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.subclass.l2 <- metav1$subclass.l2[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.subclass.l1 <- metav1$subclass.l1[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.state.l1 <- metav1$state.l1[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.state.l2 <- metav1$state.l2[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.class <- metav1$class[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.structure <- metav1$structure[match(KB@meta.data$v1clusters, metav1$clusters)]


DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE, group.by = "v1clusters", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE, group.by = "v1.subclass.l3", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE, group.by = "v1.subclass.l2", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE, group.by = "v1.subclass.l1", repel = TRUE) + NoLegend()
DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE, group.by = "v1.state.l1", repel = TRUE) + NoLegend()

#Make all multiplet clusters a single cluster 101
KB@meta.data[KB@meta.data$v1clusters %in% c(101:119), ]$v1clusters <- 101
KB@meta.data$v1clusters <- droplevels(KB@meta.data$v1clusters)
table(KB@meta.data$v1clusters)

Idents(KB) <- "v1clusters"
Idents(KB) <- factor(Idents(KB), levels = c(1:101))
table(Idents(KB))



##Fix FIB annotation - int cluster 90 - Fib cluster 81
Idents(KB)
v1.cells <- rownames(KB@meta.data)[rownames(KB@meta.data) %in% rownames(metav1)]
length(v1.cells)
FIB81 <- rownames(KB@meta.data[KB@meta.data$pagoda_k100_infomap == "90",])
FIB81 <- FIB81[!FIB81 %in% v1.cells]

DimPlot(KB, reduction = "ref.umap", pt.size = 0.5, label = TRUE,
        repel = TRUE, cells.highlight = FIB81) + NoLegend()
KB <- SetIdent(KB, cells = FIB81, value = 81)
KB$v1clusters <- Idents(object = KB)

Idents(KB) <- factor(Idents(KB), levels = c(1:101))
table(Idents(KB))

#Update meta tables
KB@meta.data$v1.subclass.l3 <- metav1$subclass.l3[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.subclass.l2 <- metav1$subclass.l2[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.subclass.l1 <- metav1$subclass.l1[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.state.l1 <- metav1$state.l1[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.state.l2 <- metav1$state.l2[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.class <- metav1$class[match(KB@meta.data$v1clusters, metav1$clusters)]
KB@meta.data$v1.structure <- metav1$structure[match(KB@meta.data$v1clusters, metav1$clusters)]

saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_04-2023_Object.Rds",
  destdir = "~/tmp_scratch/kidney_object"
)

