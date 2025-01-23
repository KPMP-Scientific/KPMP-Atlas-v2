library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(BPCells)
options(future.globals.maxSize = 1e10)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.STR <- subset(KB, v2.subclass.l3 %in% c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo",
                                           "C-FIB-OSMRhi","C-MYOF","pvFIB-RSPO3+",
                                           "pvFIB-PI16+","pvFIB","pvMYOF"))
Idents(KB.STR) <- "v2.subclass.l3"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo",
                                                    "C-FIB-OSMRhi","C-MYOF","pvFIB-RSPO3+",
                                                    "pvFIB-PI16+","pvFIB","pvMYOF"))
KB.STR[["RNA"]] <- as(KB.STR[["RNA"]], Class = "Assay")
KB.STR <- FindVariableFeatures(KB.STR, nfeatures = 2000)

###Correlation with public data
library(corrplot)
###Buechler & Pradhan, ..., Bourgon, Müller, Turley published in Nature, Vol 593, 2021 . 
##Download seurat objects from https://www.fibroxplorer.com/home

##Human Perturbed State
hPS <- readRDS("~/hsKidAt/blake_LTS/public/Buechler2021/Human_PS_Fibro.RDS")
hPS
#17410 features across 10355 samples
Idents(hPS) <- "ClustName"
table(Idents(hPS))
#LRRC15 ADAMDEC1   COL3A1    CCL19     PI16     NPNT 
#2130     1625     3247      357     1979     1017 

#Intersect all variable genes
common.genes <- rownames(KB.STR)[(rownames(KB.STR) %in% VariableFeatures(hPS)[1:500])]

KB.STR <- ScaleData(KB.STR, features = common.genes)
ave.KB.STR<-AverageExpression(KB.STR, features = common.genes, layer = "scale.data")
hPS <- ScaleData(hPS, features = common.genes)
ave.hPS<-AverageExpression(hPS, features = common.genes, layer = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.STR$RNA),as.data.frame(ave.hPS$RNA)))
ave.cor<-ave.cor[1:9,10:15]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='stroma/Stromal_Clusters_Buechler2021_PerturbedStateHuman.pdf',width=6,height=8)
corrplot(ave.cor[,c("ADAMDEC1","CCL19","LRRC15","COL3A1","PI16","NPNT")], col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()



###Korsunsky et al., 2022. https://sandbox.zenodo.org/record/772596#.Y0X0TOzMLJx
##Human Perturbed State
counts <- ReadMtx(mtx = "~/hsKidAt/blake_LTS/public/Korsunsky2022/umis.mtx",
                  cells = "~/hsKidAt/blake_LTS/public/Korsunsky2022/umis_barcodes.tsv",
                  features = "~/hsKidAt/blake_LTS/public/Korsunsky2022/umis_genes.tsv",
                  feature.column = 1)
meta <- read.table("~/hsKidAt/blake_LTS/public/Korsunsky2022/metaData.txt",sep="\t",header=TRUE,row.names=1)
umap <- read.table("~/hsKidAt/blake_LTS/public/Korsunsky2022/umap_integrated.txt",sep="\t",header=TRUE,row.names=1)

counts <- counts[,rownames(meta)]
fib.atlas <- CreateSeuratObject(counts, meta.data = meta)
fib.atlas
#19952 features across 79968 samples

fib.atlas <- NormalizeData(fib.atlas)
fib.atlas <- FindVariableFeatures(fib.atlas, nfeatures = 2000)
fib.atlas <- ScaleData(fib.atlas)
fib.atlas <- RunPCA(fib.atlas)
fib.atlas <- RunUMAP(fib.atlas, reduction = "pca", dims = 1:30, return.model = T, verbose = F)
Idents(fib.atlas) <- "cell_type_integrated"
DimPlot(fib.atlas, group.by = "cell_type_integrated", reduction = "umap")
colnames(umap) <- c("umap_1", "umap_2")
umap <- umap[rownames(Embeddings(fib.atlas, reduction = "umap")),]
embeddings <- matrix(data = c(as.numeric(umap$umap_1),as.numeric(umap$umap_2)), ncol = 2)
rownames(embeddings) <- rownames(umap)
colnames(embeddings) <- colnames(umap)

fib.atlas[["umap"]] <- CreateDimReducObject(embeddings = embeddings, key = "umap_", assay = DefaultAssay(fib.atlas))
DimPlot(fib.atlas, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        repel = TRUE) + NoLegend()


table(Idents(fib.atlas))
table(fib.atlas$sample_type)
table(fib.atlas$donor_id)

#Intersect variable genes
common.genes <- rownames(KB.STR)[rownames(KB.STR) %in% VariableFeatures(fib.atlas)[1:500]]

KB.STR <- ScaleData(KB.STR, features = common.genes)
ave.KB.STR<-AverageExpression(KB.STR, features = common.genes, layer = "scale.data")
fib.atlas <- ScaleData(fib.atlas, features = common.genes)
ave.fib.atlas<-AverageExpression(fib.atlas, features = common.genes, layer = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.STR$RNA),as.data.frame(ave.fib.atlas$RNA)))
ave.cor<-ave.cor[1:9,10:23]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='stroma/Stromal_Clusters_Korsunsky2022_Integrated.pdf',width=7,height=8)
corrplot(ave.cor[,c("C6","C10","PTGS2+SEMA4A+ C8","C12","C1","C7","C3","C2","CXCL10+CCL19+ C11",
                    "SPARC+COL3A1+ C4","CD34+MFAP5+ C9","FBLN1+ C5","MYH11+ C13","C0")], col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()

saveRDS(fib.atlas, file = "~/hsKidAt/blake_LTS/public/Korsunsky2022/Korsunsky2022_Seurat.RDS")
#fib.atlas <- readRDS("~/hsKidAt/blake_LTS/public/Korsunsky2022/Korsunsky2022_Seurat.RDS")





#Myeloid Cells
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.mye <- subset(KB, v2.subclass.l3 %in% c("resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                           "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N",
                                           "cycMAC"))
Idents(KB.mye) <- "v2.subclass.l3"
Idents(KB.mye) <- factor(Idents(KB.mye), levels = c("resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                                    "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N",
                                                    "cycMAC"))
KB.mye <- subset(KB.mye, downsample = 1000)
KB.mye[["RNA"]] <- as(KB.mye[["RNA"]], Class = "Assay")
KB.mye <- NormalizeData(KB.mye)
KB.mye <- FindVariableFeatures(KB.mye, nfeatures = 2000)

###Comparison with Eraslan et. al., 2022
#Download immune object from Eraslan et. al., 2022
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

save(IMM, file = "~/hsKidAt/blake_LTS/Atlas_V2/public/Eraslan2022/Eraslan_Science_2022_Immune_subset_Seurat.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/public/Eraslan2022/Eraslan_Science_2022_Immune_subset_Seurat.rda")


###Correlation of Eraslan imm clusters
Idents(object = IMM) <- "annotation"
order <- c("MΦ LYVE1 hi","Mo/MΦ FCGR3A hi","MΦ HLAII hi","Inflammatory MΦ","LAM-like",
           "DC2","Mo/MΦ FCGR3A lo","CD14+ monocyte","CD16+ monocyte","Mature DC","DC1" ,
           "Proliferating MΦ","Langerhans","Lung MΦ" )
Idents(IMM) <- factor(Idents(IMM), levels = order)

#Intersect variable genes
common.genes <- intersect(rownames(IMM), VariableFeatures(KB.mye))

KB.mye <- ScaleData(KB.mye, features = common.genes)
ave.KB.mye<-AverageExpression(KB.mye, features = common.genes, layer = "scale.data")
IMM <- ScaleData(IMM, features = common.genes)
ave.IMM<-AverageExpression(IMM, features = common.genes, layer = "scale.data")

ave.cor<-cor(cbind(as.data.frame(ave.KB.mye$RNA),as.data.frame(ave.IMM$RNA)))
ave.cor<-ave.cor[1:14,15:28]

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='immune/Myeloid_Clusters_Eraslan2022_CorrPlot.pdf',width=7,height=7)
corrplot(ave.cor[,c("MΦ LYVE1 hi","Mo/MΦ FCGR3A hi","MΦ HLAII hi",
                    "Mo/MΦ FCGR3A lo","CD14+ monocyte",
                    "Inflammatory MΦ","LAM-like",
                    "CD16+ monocyte","DC2","DC1" ,"Mature DC",
                    "Proliferating MΦ","Langerhans","Lung MΦ" )], col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()




###Comparison with Conway et. al., 2020
IMM <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/public/Conway2020/myeloid_sc_object.rds")
IMM <- UpdateSeuratObject(IMM)
IMM <- NormalizeData(IMM, normalization.method = "LogNormalize", scale.factor = 10000)
IMM <- ScaleData(IMM)
IMM <- RunUMAP(object = IMM, reduction = "pca", dims = 1:20, n.neighbors = 20L,
               min.dist = 0.2, return.model = T)
DimPlot(IMM, reduction = "umap", pt.size = 0.5, label = TRUE, 
        group.by = "final_classification", repel = TRUE) + NoLegend()
DimPlot(IMM, reduction = "umap", pt.size = 0.5, label = TRUE, 
        group.by = "Phenotype", repel = TRUE) 

table(Idents(IMM))

Idents(IMM) <- factor(Idents(IMM), levels = c("Mrc1+ macrophage","Resident macrophage","Ccr2+ macrophage",
                                              "Arg1+ monocyte","Mmp12+ macrophage","IFN induced macrophage",
                                              "cDC2","Patrolling monocyte","Ly6c2+ monocyte","Ccr7+ monocyte",
                                              "cDC1","Proliferating cell"))

save(IMM, file = "~/hsKidAt/blake_LTS/Atlas_V2/public/Conway2020/Conway_JASN_2020_Myeloid_subset_Seurat.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/public/Conway2020/Conway_JASN_2020_Myeloid_subset_Seurat.rda")


##using homologous genes
###Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% rownames(KB.mye),]
rownames(hs.ms.table) <- hs.ms.table$Mouse.gene.name
hs.ms.table <- hs.ms.table[rownames(hs.ms.table) %in% rownames(IMM),]

#Intersect variable genes
IMM <- FindVariableFeatures(IMM, nfeatures = 1000)
IMM <- ScaleData(IMM, features = VariableFeatures(IMM))
ave.IMM<-AverageExpression(IMM, features = VariableFeatures(IMM), assays = "RNA",
                           slot = "scale.data")
KB.mye <- ScaleData(KB.mye, features = VariableFeatures(KB.mye))
ave.KB.mye<-AverageExpression(KB.mye, features = VariableFeatures(KB.mye), assays = "RNA",
                             slot = "scale.data")

ave.IMM$RNA <- ave.IMM$RNA[rownames(ave.IMM$RNA) %in% rownames(hs.ms.table),]
rownames(ave.IMM$RNA) <- as.character(hs.ms.table[rownames(ave.IMM$RNA),]$Gene.name)
select.markers <- intersect(rownames(ave.IMM$RNA), rownames(ave.KB.mye$RNA))

ave.cor<-cor(as.data.frame(ave.KB.mye$RNA[select.markers,]),as.data.frame(ave.IMM$RNA[select.markers,]))

col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
pdf(file='immune/Myeloid_Clusters_Conway2020_CorrPlot.pdf',width=7,height=7)
corrplot(ave.cor[,c("Mrc1+ macrophage","Resident macrophage","Ccr2+ macrophage",
                    "Ly6c2+ monocyte","Arg1+ monocyte","Mmp12+ macrophage","IFN induced macrophage",
                    "Patrolling monocyte","cDC2","cDC1","Ccr7+ monocyte",
                    "Proliferating cell")], col = col2(500),  method="color",order = "original", cl.lim=c(-1,1))
dev.off()
