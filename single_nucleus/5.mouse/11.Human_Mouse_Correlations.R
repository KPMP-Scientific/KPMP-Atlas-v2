library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(pagoda2)
library(BPCells)
library(corrplot)
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
mouse.meta <- KB@meta.data


##PT-TL
#load mouse object
mKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_PT-TL-Subset_C.Rds")
mKB@meta.data <- mouse.meta[rownames(mKB@meta.data),]
mKB[["RNA"]] <- as(object = mKB[["RNA"]], Class = "Assay")
Idents(mKB) <- "v2.subclass.l3"
order <- c("POD","PEC","PT-S1","PT-S1/2","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT",
           "aPT1-frPT","aPT-S1/S2","aPT-S3","dPT","cycPT","DTL1","DTL2","aDTL",
           "DTL3","ATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B",
           "C-TAL-B","dC-TAL","aTAL1","aTAL2","frTAL","infTAL","cycTAL","DCT1",
           "dDCT1","DCT2","dDCT2","aDCT","CNT","CNT-PC","dCNT-PC","CCD-PC","OMCD-PC",
           "dPC","IMCD","PapE","CCD-IC-A","OMCD-IC-A","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-EA","M-EC-PTC","C-EC-PTC","dEC-PTC","angEC-PTC","infEC-PTC","iaEC-PTC",
           "EC-V","EC-PCV","EC-AVR","P-EC-AVR","dEC-AVR","cycEC","EC-LYM","IM-FIB",
           "OM-FIB","dM-FIB","C/M-FIB","C-FIB","C-FIB-PATH","eaC-FIB","C-FIB-Cxcl14+",
           "C-FIB-Cxcl10+","C-MYOF","cycFIB","pvFIB-Rspo3+","pvFIB-Pi16+","pvFIB",
           "pvMYOF","MC","REN","VSMC","VSMC/P","Ad","B","PL","T","T-REG","Tcm/Naive",
           "CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","resMAC","resMAC-H2hi","resMAC-Cpd+",
           "MAC","moMAC-IA","moMAC-Cxcl10+","moFAM","moMAC","moMAC-Ccr2+","moMAC-Agr1+",
           "MON","cDC2","mDC","cDC1","pDC","cycMAC","cycDC","SC/NEU")
Idents(mKB) <- factor(Idents(mKB), levels = order)

mKB <- subset(mKB, downsample = 500)

#load human object
hKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(hKB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

hKB <- subset(hKB, idents = c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
                              "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
                              "dDTL","ATL","dATL"))

hKB <- subset(hKB, downsample = 500)
hKB[["RNA"]] <- as(hKB[["RNA"]], Class = "Assay")

#Convert to homologous genes
#Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/Projects/Mouse_Kidney/Atlas_V1/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% rownames(hKB),]

rownames(hs.ms.table) <- hs.ms.table$Gene.name

counts <- hKB[["RNA"]]$counts
counts <- counts[rownames(counts) %in% rownames(hs.ms.table),]
rownames(counts) <- as.character(hs.ms.table[rownames(counts),]$Mouse.gene.name)

meta <- hKB@meta.data
hKB <- CreateSeuratObject(counts = counts, project = "Human Kidney Atlas V2", min.cells = 3, min.features = 200, 
                          meta.data = meta)
Idents(hKB) <- "v2.subclass.l3"
order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","dPT","cycPT","DTL1","DTL2","aDTL2","DTL3",
           "dDTL","ATL","dATL")
Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

load("AtlasV2_VariableFeatures.rda")
od.genes <- unique(c(pt.od.genes[1:500],dt.od.genes[1:500],imm.od.genes[1:500],
                     str.od.genes[1:500],ec.od.genes[1:500]))
sn.od.genes <- str_to_title(pt.od.genes)[str_to_title(pt.od.genes) %in% VariableFeatures(mKB)]

#hKB.IMM[["RNA"]] <- as(object = hKB.IMM[["RNA"]], Class = "Assay")
hKB <- NormalizeData(hKB)
hKB <- ScaleData(hKB, features = sn.od.genes, assay = "RNA")

mKB <- NormalizeData(mKB)
mKB <- ScaleData(mKB, features = sn.od.genes, assay = "RNA")

ave.hKB <- AverageExpression(hKB, features = sn.od.genes, assays = "RNA", slot = "scale.data")
ave.mKB <- AverageExpression(mKB, features = sn.od.genes, assay = "RNA", slot = "scale.data")
to.use <- rownames(ave.hKB$RNA)[rownames(ave.hKB$RNA) %in% rownames(ave.mKB$RNA)]
ave.cor<-cor(as.data.frame(ave.mKB$RNA[to.use,]),as.data.frame(ave.hKB$RNA[to.use,]))
colnames(ave.cor) <- gsub("RNA.", "", colnames(ave.cor))
rownames(ave.cor) <- gsub("RNA.", "", rownames(ave.cor))
write.table(ave.cor, file="mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_PT-TL_corr.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_PT-TL_corr_Plot.pdf',width=14,height=14)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))
dev.off()





##TAL-CD
#load mouse object
mKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_11-2024_Filtered_TAL-CD-Subset_C.Rds")
mKB@meta.data <- mouse.meta[rownames(mKB@meta.data),]
mKB[["RNA"]] <- as(object = mKB[["RNA"]], Class = "Assay")
Idents(mKB) <- "v2.subclass.l3"
order <- c("POD","PEC","PT-S1","PT-S1/2","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT",
           "aPT1-frPT","aPT-S1/S2","aPT-S3","dPT","cycPT","DTL1","DTL2","aDTL",
           "DTL3","ATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B",
           "C-TAL-B","dC-TAL","aTAL1","aTAL2","frTAL","infTAL","cycTAL","DCT1",
           "dDCT1","DCT2","dDCT2","aDCT","CNT","CNT-PC","dCNT-PC","CCD-PC","OMCD-PC",
           "dPC","IMCD","PapE","CCD-IC-A","OMCD-IC-A","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-EA","M-EC-PTC","C-EC-PTC","dEC-PTC","angEC-PTC","infEC-PTC","iaEC-PTC",
           "EC-V","EC-PCV","EC-AVR","P-EC-AVR","dEC-AVR","cycEC","EC-LYM","IM-FIB",
           "OM-FIB","dM-FIB","C/M-FIB","C-FIB","C-FIB-PATH","eaC-FIB","C-FIB-Cxcl14+",
           "C-FIB-Cxcl10+","C-MYOF","cycFIB","pvFIB-Rspo3+","pvFIB-Pi16+","pvFIB",
           "pvMYOF","MC","REN","VSMC","VSMC/P","Ad","B","PL","T","T-REG","Tcm/Naive",
           "CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","resMAC","resMAC-H2hi","resMAC-Cpd+",
           "MAC","moMAC-IA","moMAC-Cxcl10+","moFAM","moMAC","moMAC-Ccr2+","moMAC-Agr1+",
           "MON","cDC2","mDC","cDC1","pDC","cycMAC","cycDC","SC/NEU")
Idents(mKB) <- factor(Idents(mKB), levels = order)

mKB <- subset(mKB, downsample = 500)

#load human object
hKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(hKB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

hKB <- subset(hKB, idents = c("M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
                              "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
                              "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
                              "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B"))

hKB <- subset(hKB, downsample = 500)
hKB[["RNA"]] <- as(hKB[["RNA"]], Class = "Assay")

#Convert to homologous genes
#Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/Projects/Mouse_Kidney/Atlas_V1/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% rownames(hKB),]

rownames(hs.ms.table) <- hs.ms.table$Gene.name

counts <- hKB[["RNA"]]$counts
counts <- counts[rownames(counts) %in% rownames(hs.ms.table),]
rownames(counts) <- as.character(hs.ms.table[rownames(counts),]$Mouse.gene.name)

meta <- hKB@meta.data
hKB <- CreateSeuratObject(counts = counts, project = "Human Kidney Atlas V2", min.cells = 3, min.features = 200, 
                          meta.data = meta)
Idents(hKB) <- "v2.subclass.l3"
order <- c("M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","dDCT","DCT2","aDCT","frDCT","CNT","CNT-PC","dCNT","aCNT",
           "dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B")
Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

load("AtlasV2_VariableFeatures.rda")
od.genes <- unique(c(pt.od.genes[1:500],dt.od.genes[1:500],imm.od.genes[1:500],
                     str.od.genes[1:500],ec.od.genes[1:500]))
sn.od.genes <- str_to_title(dt.od.genes)[str_to_title(dt.od.genes) %in% VariableFeatures(mKB)]

hKB <- NormalizeData(hKB)
hKB <- ScaleData(hKB, features = sn.od.genes, assay = "RNA")

mKB <- NormalizeData(mKB)
mKB <- ScaleData(mKB, features = sn.od.genes, assay = "RNA")

ave.hKB <- AverageExpression(hKB, features = sn.od.genes, assays = "RNA", slot = "scale.data")
ave.mKB <- AverageExpression(mKB, features = sn.od.genes, assay = "RNA", slot = "scale.data")
to.use <- rownames(ave.hKB$RNA)[rownames(ave.hKB$RNA) %in% rownames(ave.mKB$RNA)]
ave.cor<-cor(as.data.frame(ave.mKB$RNA[to.use,]),as.data.frame(ave.hKB$RNA[to.use,]))
colnames(ave.cor) <- gsub("RNA.", "", colnames(ave.cor))
rownames(ave.cor) <- gsub("RNA.", "", rownames(ave.cor))
write.table(ave.cor, file="mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_TAL-CD_corr.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_TAL-CD_corr_Plot.pdf',width=14,height=14)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))
dev.off()




###EC
#load mouse object
mKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Endothelial-Subset_C.Rds")
mKB@meta.data <- mouse.meta[rownames(mKB@meta.data),]
mKB[["RNA"]] <- as(object = mKB[["RNA"]], Class = "Assay")
Idents(mKB) <- "v2.subclass.l3"
order <- c("POD","PEC","PT-S1","PT-S1/2","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT",
           "aPT1-frPT","aPT-S1/S2","aPT-S3","dPT","cycPT","DTL1","DTL2","aDTL",
           "DTL3","ATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B",
           "C-TAL-B","dC-TAL","aTAL1","aTAL2","frTAL","infTAL","cycTAL","DCT1",
           "dDCT1","DCT2","dDCT2","aDCT","CNT","CNT-PC","dCNT-PC","CCD-PC","OMCD-PC",
           "dPC","IMCD","PapE","CCD-IC-A","OMCD-IC-A","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-EA","M-EC-PTC","C-EC-PTC","dEC-PTC","angEC-PTC","infEC-PTC","iaEC-PTC",
           "EC-V","EC-PCV","EC-AVR","P-EC-AVR","dEC-AVR","cycEC","EC-LYM","IM-FIB",
           "OM-FIB","dM-FIB","C/M-FIB","C-FIB","C-FIB-PATH","eaC-FIB","C-FIB-Cxcl14+",
           "C-FIB-Cxcl10+","C-MYOF","cycFIB","pvFIB-Rspo3+","pvFIB-Pi16+","pvFIB",
           "pvMYOF","MC","REN","VSMC","VSMC/P","Ad","B","PL","T","T-REG","Tcm/Naive",
           "CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","resMAC","resMAC-H2hi","resMAC-Cpd+",
           "MAC","moMAC-IA","moMAC-Cxcl10+","moFAM","moMAC","moMAC-Ccr2+","moMAC-Agr1+",
           "MON","cDC2","mDC","cDC1","pDC","cycMAC","cycDC","SC/NEU")
Idents(mKB) <- factor(Idents(mKB), levels = order)

mKB <- subset(mKB, downsample = 500)

#load human object
hKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(hKB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

hKB <- subset(hKB, idents = c("EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
                              "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
                              "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC"))

hKB <- subset(hKB, downsample = 500)
hKB[["RNA"]] <- as(hKB[["RNA"]], Class = "Assay")

#Convert to homologous genes
#Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/Projects/Mouse_Kidney/Atlas_V1/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% rownames(hKB),]

rownames(hs.ms.table) <- hs.ms.table$Gene.name

counts <- hKB[["RNA"]]$counts
counts <- counts[rownames(counts) %in% rownames(hs.ms.table),]
rownames(counts) <- as.character(hs.ms.table[rownames(counts),]$Mouse.gene.name)

meta <- hKB@meta.data
hKB <- CreateSeuratObject(counts = counts, project = "Human Kidney Atlas V2", min.cells = 3, min.features = 200, 
                          meta.data = meta)
Idents(hKB) <- "v2.subclass.l3"
order <- c("EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","EC-EA","M-EC-PTC","C-EC-PTC","dEC-PTC",
           "angEC-PTC","infEC-PTC","iaEC-PTC",
           "EC-V","EC-PCV","EC-AVR","dEC-AVR","iaEC-AVR",
           "cycEC","EC-LYM"
)
Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

load("AtlasV2_VariableFeatures.rda")
od.genes <- unique(c(pt.od.genes[1:500],dt.od.genes[1:500],imm.od.genes[1:500],
                     str.od.genes[1:500],ec.od.genes[1:500]))
sn.od.genes <- str_to_title(ec.od.genes)[str_to_title(ec.od.genes) %in% VariableFeatures(mKB)]

hKB <- NormalizeData(hKB)
hKB <- ScaleData(hKB, features = sn.od.genes, assay = "RNA")

mKB <- NormalizeData(mKB)
mKB <- ScaleData(mKB, features = sn.od.genes, assay = "RNA")

ave.hKB <- AverageExpression(hKB, features = sn.od.genes, assays = "RNA", slot = "scale.data")
ave.mKB <- AverageExpression(mKB, features = sn.od.genes, assay = "RNA", slot = "scale.data")
to.use <- rownames(ave.hKB$RNA)[rownames(ave.hKB$RNA) %in% rownames(ave.mKB$RNA)]
ave.cor<-cor(as.data.frame(ave.mKB$RNA[to.use,]),as.data.frame(ave.hKB$RNA[to.use,]))
colnames(ave.cor) <- gsub("RNA.", "", colnames(ave.cor))
rownames(ave.cor) <- gsub("RNA.", "", rownames(ave.cor))
write.table(ave.cor, file="mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_Vasculature_corr.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_Vasculature_corr_Plot.pdf',width=14,height=14)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))
dev.off()




###Stroma
#load mouse object
mKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Stroma-Subset_C.Rds")
mKB@meta.data <- mouse.meta[rownames(mKB@meta.data),]
mKB[["RNA"]] <- as(object = mKB[["RNA"]], Class = "Assay")
Idents(mKB) <- "v2.subclass.l3"
order <- c("POD","PEC","PT-S1","PT-S1/2","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT",
           "aPT1-frPT","aPT-S1/S2","aPT-S3","dPT","cycPT","DTL1","DTL2","aDTL",
           "DTL3","ATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B",
           "C-TAL-B","dC-TAL","aTAL1","aTAL2","frTAL","infTAL","cycTAL","DCT1",
           "dDCT1","DCT2","dDCT2","aDCT","CNT","CNT-PC","dCNT-PC","CCD-PC","OMCD-PC",
           "dPC","IMCD","PapE","CCD-IC-A","OMCD-IC-A","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-EA","M-EC-PTC","C-EC-PTC","dEC-PTC","angEC-PTC","infEC-PTC","iaEC-PTC",
           "EC-V","EC-PCV","EC-AVR","P-EC-AVR","dEC-AVR","cycEC","EC-LYM","IM-FIB",
           "OM-FIB","dM-FIB","C/M-FIB","C-FIB","C-FIB-PATH","eaC-FIB","C-FIB-Cxcl14+",
           "C-FIB-Cxcl10+","C-MYOF","cycFIB","pvFIB-Rspo3+","pvFIB-Pi16+","pvFIB",
           "pvMYOF","MC","REN","VSMC","VSMC/P","Ad","B","PL","T","T-REG","Tcm/Naive",
           "CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","resMAC","resMAC-H2hi","resMAC-Cpd+",
           "MAC","moMAC-IA","moMAC-Cxcl10+","moFAM","moMAC","moMAC-Ccr2+","moMAC-Agr1+",
           "MON","cDC2","mDC","cDC1","pDC","cycMAC","cycDC","SC/NEU")
Idents(mKB) <- factor(Idents(mKB), levels = order)

mKB <- subset(mKB, downsample = 500)

#load human object
hKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(hKB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

hKB <- subset(hKB, idents = c("IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
                              "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
                              "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
                              "Ad","SC/NEU"))

hKB <- subset(hKB, downsample = 1000)
hKB[["RNA"]] <- as(hKB[["RNA"]], Class = "Assay")

#Convert to homologous genes
#Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/Projects/Mouse_Kidney/Atlas_V1/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% rownames(hKB),]

rownames(hs.ms.table) <- hs.ms.table$Gene.name

counts <- hKB[["RNA"]]$counts
counts <- counts[rownames(counts) %in% rownames(hs.ms.table),]
rownames(counts) <- as.character(hs.ms.table[rownames(counts),]$Mouse.gene.name)

meta <- hKB@meta.data
hKB <- CreateSeuratObject(counts = counts, project = "Human Kidney Atlas V2", min.cells = 3, min.features = 200, 
                          meta.data = meta)
Idents(hKB) <- "v2.subclass.l3"
Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

load("AtlasV2_VariableFeatures.rda")
od.genes <- unique(c(pt.od.genes[1:500],dt.od.genes[1:500],imm.od.genes[1:500],
                     str.od.genes[1:500],ec.od.genes[1:500]))
sn.od.genes <- str_to_title(str.od.genes)[str_to_title(str.od.genes) %in% VariableFeatures(mKB)]

hKB <- NormalizeData(hKB)
hKB <- ScaleData(hKB, features = sn.od.genes, assay = "RNA")

mKB <- NormalizeData(mKB)
mKB <- ScaleData(mKB, features = sn.od.genes, assay = "RNA")

ave.hKB <- AverageExpression(hKB, features = sn.od.genes, assays = "RNA", slot = "scale.data")
ave.mKB <- AverageExpression(mKB, features = sn.od.genes, assay = "RNA", slot = "scale.data")
to.use <- rownames(ave.hKB$RNA)[rownames(ave.hKB$RNA) %in% rownames(ave.mKB$RNA)]
ave.cor<-cor(as.data.frame(ave.mKB$RNA[to.use,]),as.data.frame(ave.hKB$RNA[to.use,]))
colnames(ave.cor) <- gsub("RNA.", "", colnames(ave.cor))
rownames(ave.cor) <- gsub("RNA.", "", rownames(ave.cor))
write.table(ave.cor, file="mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_Stroma_corr.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_Stroma_corr_Plot.pdf',width=14,height=14)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))
dev.off()



###Immune
#load mouse object
mKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_10-2024_Filtered_Immune-Subset_C.Rds")
mKB@meta.data <- mouse.meta[rownames(mKB@meta.data),]
mKB[["RNA"]] <- as(object = mKB[["RNA"]], Class = "Assay")
Idents(mKB) <- "v2.subclass.l3"
order <- c("POD","PEC","PT-S1","PT-S1/2","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT",
           "aPT1-frPT","aPT-S1/S2","aPT-S3","dPT","cycPT","DTL1","DTL2","aDTL",
           "DTL3","ATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B",
           "C-TAL-B","dC-TAL","aTAL1","aTAL2","frTAL","infTAL","cycTAL","DCT1",
           "dDCT1","DCT2","dDCT2","aDCT","CNT","CNT-PC","dCNT-PC","CCD-PC","OMCD-PC",
           "dPC","IMCD","PapE","CCD-IC-A","OMCD-IC-A","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-EA","M-EC-PTC","C-EC-PTC","dEC-PTC","angEC-PTC","infEC-PTC","iaEC-PTC",
           "EC-V","EC-PCV","EC-AVR","P-EC-AVR","dEC-AVR","cycEC","EC-LYM","IM-FIB",
           "OM-FIB","dM-FIB","C/M-FIB","C-FIB","C-FIB-PATH","eaC-FIB","C-FIB-Cxcl14+",
           "C-FIB-Cxcl10+","C-MYOF","cycFIB","pvFIB-Rspo3+","pvFIB-Pi16+","pvFIB",
           "pvMYOF","MC","REN","VSMC","VSMC/P","Ad","B","PL","T","T-REG","Tcm/Naive",
           "CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","resMAC","resMAC-H2hi","resMAC-Cpd+",
           "MAC","moMAC-IA","moMAC-Cxcl10+","moFAM","moMAC","moMAC-Ccr2+","moMAC-Agr1+",
           "MON","cDC2","mDC","cDC1","pDC","cycMAC","cycDC","SC/NEU")
Idents(mKB) <- factor(Idents(mKB), levels = order)

mKB <- subset(mKB, downsample = 500)

#load human object
hKB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(hKB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = hKB) <- factor(Idents(object = hKB), levels = order)

hKB.IMM <- subset(hKB, idents = c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
                                  "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
                                  "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC"))

hKB.IMM <- subset(hKB.IMM, downsample = 500)
hKB.IMM[["RNA"]] <- as(hKB.IMM[["RNA"]], Class = "Assay")

#Convert to homologous genes
#Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/Projects/Mouse_Kidney/Atlas_V1/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% rownames(hKB.IMM),]

rownames(hs.ms.table) <- hs.ms.table$Gene.name

counts <- hKB.IMM[["RNA"]]$counts
counts <- counts[rownames(counts) %in% rownames(hs.ms.table),]
rownames(counts) <- as.character(hs.ms.table[rownames(counts),]$Mouse.gene.name)

meta <- hKB.IMM@meta.data
hKB.IMM <- CreateSeuratObject(counts = counts, project = "Human Kidney Atlas V2", min.cells = 3, min.features = 200, 
                              meta.data = meta)
Idents(hKB.IMM) <- "v2.subclass.l3"
order <- c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "MON","ncMON","cDC2","mDC","cDC1","pDC","N","cycT","cycMAC")
Idents(object = hKB.IMM) <- factor(Idents(object = hKB.IMM), levels = order)

load("AtlasV2_VariableFeatures.rda")
od.genes <- unique(c(pt.od.genes[1:500],dt.od.genes[1:500],imm.od.genes[1:500],
                     str.od.genes[1:500],ec.od.genes[1:500]))
sn.od.genes <- str_to_title(imm.od.genes)[str_to_title(imm.od.genes) %in% VariableFeatures(mKB)]

hKB.IMM <- NormalizeData(hKB.IMM)
hKB.IMM <- ScaleData(hKB.IMM, features = sn.od.genes, assay = "RNA")

mKB <- NormalizeData(mKB)
mKB <- ScaleData(mKB, features = sn.od.genes, assay = "RNA")

ave.hKB <- AverageExpression(hKB.IMM, features = sn.od.genes, assays = "RNA", slot = "scale.data")
ave.mKB <- AverageExpression(mKB, features = sn.od.genes, assay = "RNA", slot = "scale.data")
to.use <- rownames(ave.hKB$RNA)[rownames(ave.hKB$RNA) %in% rownames(ave.mKB$RNA)]
ave.cor<-cor(as.data.frame(ave.mKB$RNA[to.use,]),as.data.frame(ave.hKB$RNA[to.use,]))
colnames(ave.cor) <- gsub("RNA.", "", colnames(ave.cor))
rownames(ave.cor) <- gsub("RNA.", "", rownames(ave.cor))
write.table(ave.cor, file="mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_Immune_corr.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='mouse_IRI/Mouse-v2-subclasses_vs_Human-v2-subclass_Immune_corr_Plot.pdf',width=14,height=14)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))
dev.off()
