library(Seurat)
library(dplyr)
library(Matrix)
library(pagoda2)
require(parallel)
library(foreach)
library(doParallel)
library(BPCells)
registerDoParallel(cores=4)

options(future.globals.maxSize = 1e9)

setwd("~/data/net/home/Human_Kidney/Atlas_V2")
source("misc/utils.R")

###Append Library IDs to Cell Barcodes
#snRNA-seq Experiments - exp set 1
processed.dir <- "~/data/pod/processed/cellranger/count/sample-"
dir2 <- "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/" 

exp.table <- read.delim("Data_Processing/snRNA_Set1_Exp_Table_03162023.txt", row.names = 1)
exps <- exp.table$Exp
libs <- exp.table$Lib

#Read in counts, append library ID and save
foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exp[[1]]) 1 else 2
  
  lib = exp.table[exp,]$Lib
  dir1 = paste0(processed.dir,exp)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  #counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.mtx"))
}


#snRNA-seq Experiments - exp set 2
processed.dir <- "~/data/pod/processed/cellranger/count/sample-"
dir2 <- "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/" 

exp.table <- read.delim("Data_Processing/snRNA_Set2_Exp_Table_04072023.txt", row.names = 1)
exps <- exp.table$Exp
libs <- exp.table$Lib

#Read in counts, append library ID and save
foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exp[[1]]) 1 else 2
  
  lib = exp.table[exp,]$Lib
  dir1 = paste0(processed.dir,exp)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  #counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.mtx"))
}


#snRNA-seq Experiments - exp set 3 - KB54 which showed >80K cells... needs trimming down
processed.dir <- "~/data/pod/processed/cellranger/count/sample-"
dir2 <- "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/" 

exp <- "KPMP_20211208C"
lib <- "KB54"

#Read in counts, append library ID, trim based on UMI counts and save

dir1 = paste0(processed.dir,exp)
counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
colnames(counts) <- paste(lib, colnames(counts), sep="_")
top.cells <- names(colSums(counts)[colSums(counts) > 400])
counts <- counts[,top.cells]
saveRDS(counts, file = paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.rds"))
writeMM(counts, file= paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.mtx"))




#Get Cellranger QC stats
QCpath = "~/data/pod/processed/cellranger/count/"
exp.list = list.files(path=QCpath, pattern="/*") 
QCpath = paste("~/data/pod/processed/cellranger/count/",exp.list, sep = "")
QC.list = paste(QCpath,"/outs/metrics_summary.csv", sep = "")
QC_files_df <- lapply(QC.list[1:179], function(x) {read.csv(file = x)})
QC.df <- do.call(rbind, QC_files_df)
rownames(QC.df) <- gsub("sample-","",exp.list[1:179])
write.table(QC.df, file="QC_Plots/Cellranger_QC_Table_03-2023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)




#Multiome Experiments - exp set 1
processed.dir <- "~/data/pod/processed/cellranger-arc/"
dir2 <- "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/" 

exp.table <- read.delim("Data_Processing/Multiome_Set1_Exp_Table_0317202.txt", row.names = 1)
exps <- exp.table$Exp
libs <- exp.table$Lib

#Read in counts, append library ID and save
foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exp[[1]]) 1 else 2
  
  lib = exp.table[exp,]$Lib
  dir1 = paste0(processed.dir,exp)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.mtx"))
}


#Multiome Experiments - exp set 2
processed.dir <- "~/data/pod/processed/cellranger-arc/"
dir2 <- "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/" 

exp.table <- read.delim("Data_Processing/Multiome_Set2_Exp_Table_04072023.txt", row.names = 1)
exps <- exp.table$Exp
libs <- exp.table$Lib

#Read in counts, append library ID and save
foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exp[[1]]) 1 else 2
  
  lib = exp.table[exp,]$Lib
  dir1 = paste0(processed.dir,exp)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.mtx"))
}


#Multiome Experiments - exp set 3
processed.dir <- "~/data/pod/processed/cellranger-arc/"
dir2 <- "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/" 

exps <- c("HuBMAP_20230222D")
libs <- c("PA14")

#Read in counts, append library ID and save
foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exp[[1]]) 1 else 2
  
  lib = libs
  dir1 = paste0(processed.dir,exp)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.mtx"))
}


#Multiome Experiments - exp set 4
processed.dir <- "~/data/pod/processed/cellranger-arc/"
dir2 <- "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/" 

exp.table <- read.delim("Data_Processing/Multiome_Exp_Table_04082023.txt", row.names = 1)
exps <- c("KTRC_20220721A","KTRC_20220721B","KTRC_20220721C")
libs <- c("KC131","KM142","KM143")

#Read in counts, append library ID and save
foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  lib = exp.table[exp,]$Lib
  dir1 = paste0(processed.dir,exp)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.mtx"))
}



#Get Cellranger QC stats
QCpath = "~/data/pod/processed/cellranger-arc/"
exp.list = list.files(path=QCpath, pattern="/*") 
QCpath = paste("~/data/pod/processed/cellranger-arc/",exp.list, sep = "")
QC.list = paste(QCpath,"/outs/summary.csv", sep = "")
QC_files_df <- lapply(QC.list, function(x) {read.csv(file = x, row.names = 1)})
QC.df <- do.call(rbind, QC_files_df)
write.table(QC.df, file="QC_Plots/Cellranger-arc_QC_Table_03-2023.txt", sep = "\t", row.names=TRUE, col.names=TRUE)




###DoubletDetection (v3.0) Filter - ran on python  

dir1 = "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/"
dir2 = "~/data/pod/blake_LTS/Atlas_V2/DoubletDetection/"

#Prepare doublet factors - 10X-R Experiments
exp.table <- read.delim("Data_Processing/snRNA_Exp_Table_04082023.txt", row.names = 1)
exps <- exp.table$Exp

lapply(exps, function(x) {
  print(paste("Running for Experiment:", x))
  
  countMatrix <- readRDS(paste0(dir1, x, "_10X-R_UMI_counts_EmptyBC_Filter.rds"))
  
  #Run DoubletDetection on Jupiter Notebook
  doublet.det<-read.delim(paste0(dir2, x, "_10X-R_DoubletDetection_doublets.txt"),sep="\t",header=FALSE)
  doublet.det<-as.character(doublet.det$V1)
  names(doublet.det)<-colnames(countMatrix)
  doublet.det[doublet.det == "1"] <- "doublet"
  doublet.det[doublet.det == "0"] <- "singlet"
  
  saveRDS(doublet.det, file = paste0(dir2, x, "_10X-R_DD_Factor.rds"))
  
})

#Prepare doublet factors - 10X multiome Experiments
exp.table <- read.delim("Data_Processing/Multiome_Exp_Table_04082023.txt", row.names = 1)
exps <- exp.table$Exp

lapply(exps, function(x) {
  print(paste("Running for Experiment:", x))
  
  countMatrix <- readRDS(paste0(dir1, x, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  
  #Run DoubletDetection on Jupiter Notebook
  doublet.det<-read.delim(paste0(dir2, x, "_10X-Dual_DoubletDetection_doublets.txt"),sep="\t",header=FALSE)
  doublet.det<-as.character(doublet.det$V1)
  names(doublet.det)<-colnames(countMatrix)
  doublet.det[doublet.det == "1"] <- "doublet"
  doublet.det[doublet.det == "0"] <- "singlet"
  
  saveRDS(doublet.det, file = paste0(dir2, x, "_10X-Dual_DD_Factor.rds"))
  
})


#Prepare doublet factors - 10X multiome Experiments (2)
exps <- "HuBMAP_20230222D"

lapply(exps, function(x) {
  print(paste("Running for Experiment:", x))
  
  countMatrix <- readRDS(paste0(dir1, x, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  
  #Run DoubletDetection on Jupiter Notebook
  doublet.det<-read.delim(paste0(dir2, x, "_10X-Dual_DoubletDetection_doublets.txt"),sep="\t",header=FALSE)
  doublet.det<-as.character(doublet.det$V1)
  names(doublet.det)<-colnames(countMatrix)
  doublet.det[doublet.det == "1"] <- "doublet"
  doublet.det[doublet.det == "0"] <- "singlet"
  
  saveRDS(doublet.det, file = paste0(dir2, x, "_10X-Dual_DD_Factor.rds"))
  
})


#Prepare doublet factors - 10X multiome Experiments (3)
exps <- c("KTRC_20220721A","KTRC_20220721B","KTRC_20220721C")

lapply(exps, function(x) {
  print(paste("Running for Experiment:", x))
  
  countMatrix <- readRDS(paste0(dir1, x, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  
  #Run DoubletDetection on Jupiter Notebook
  doublet.det<-read.delim(paste0(dir2, x, "_10X-Dual_DoubletDetection_doublets.txt"),sep="\t",header=FALSE)
  doublet.det<-as.character(doublet.det$V1)
  names(doublet.det)<-colnames(countMatrix)
  doublet.det[doublet.det == "1"] <- "doublet"
  doublet.det[doublet.det == "0"] <- "singlet"
  
  saveRDS(doublet.det, file = paste0(dir2, x, "_10X-Dual_DD_Factor.rds"))
  
})





### Combine Samples and Gene/UMI Filter 
##10X-R - set 1 
CNTpath = "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/"
CNT.list = list.files(path=CNTpath, pattern="*10X-R_UMI_counts_EmptyBC_Filter.rds") 
CNT.list = paste(CNTpath,CNT.list, sep = "")
CNT_files_df <- lapply(CNT.list[1:70], function(x) {readRDS(file = x)})
KB1 <- merge.sparse(CNT_files_df)

##10X-R - set 2 
CNT_files_df <- lapply(CNT.list[71:140], function(x) {readRDS(file = x)})
KB2 <- merge.sparse(CNT_files_df)

##10X-Dual - set 3 
CNTpath = "~/data/pod/blake_LTS/Atlas_V2/Raw_Counts/"
CNT.list = list.files(path=CNTpath, pattern="*_10X-Dual_UMI_counts_EmptyBC_Filter.rds") 
CNT.list = paste(CNTpath,CNT.list, sep = "")
CNT_files_df <- lapply(CNT.list[1:55], function(x) {readRDS(file = x)})
KB3 <- merge.sparse(CNT_files_df)

##10X-Dual - set 4 
CNT_files_df <- lapply(CNT.list[56:107], function(x) {readRDS(file = x)})
KB4 <- merge.sparse(CNT_files_df)


##Create Seurat Object and Apply Gene/UMI Filters
##Prepare H5AD files for all genes and merge 
KB1 <- convert_matrix_type(KB1, type = "uint32_t")
KB2 <- convert_matrix_type(KB2, type = "uint32_t")
KB3 <- convert_matrix_type(KB3, type = "uint32_t")
KB4 <- convert_matrix_type(KB4, type = "uint32_t")

# Output BPCells matrices on-disk
write_matrix_dir(
  mat = KB1,
  dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set1_BP",
  overwrite = TRUE
)
write_matrix_dir(
  mat = KB2,
  dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set2_BP",
  overwrite = TRUE
)
write_matrix_dir(
  mat = KB3,
  dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set3_BP",
  overwrite = TRUE
)
write_matrix_dir(
  mat = KB4,
  dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set4_BP",
  overwrite = TRUE
)

# Load in BP matrices
KB1.mat <- open_matrix_dir(dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set1_BP")
KB2.mat <- open_matrix_dir(dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set2_BP")
KB3.mat <- open_matrix_dir(dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set3_BP")
KB4.mat <- open_matrix_dir(dir = "~/data/pod/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2023_Set4_BP")

data.list <- list(KB1.mat, KB2.mat, KB3.mat, KB4.mat)
names(data.list) <- c("snRNA1", "snRNA2", "Multiome1", "Multiome2")
options(Seurat.object.assay.version = "v5")
KB <- CreateSeuratObject(counts = data.list)
KB
#An object of class Seurat 
#36601 features across 1899301 samples within 1 assay 
#Active assay: RNA (36601 features, 0 variable features)
#4 layers present: counts.snRNA1, counts.snRNA2, counts.Multiome1, counts.Multiome2

##Gene Filter >400 and <7500 genes
KB <- subset(x = KB, subset = nFeature_RNA > 400 & nFeature_RNA < 7500)
FeatureScatter(object = KB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
KB
#36601 features across 1681418 samples

KB[["RNA"]]$counts.snRNA1

##Merge layers
Layers(KB[["RNA"]])
KB[["RNA"]] <- JoinLayers(KB[["RNA"]])

##Add in doublet information
DDpath = "~/data/pod/blake_LTS/Atlas_V2/DoubletDetection/"
DD.list = list.files(path=DDpath, pattern="*_DD_Factor.rds") 
DD.list = paste(DDpath,DD.list, sep = "")
DD_files_df <- lapply(DD.list, function(x) {readRDS(file = x)})
DD.df <- unlist(DD_files_df)

KB@meta.data$doubletdetection <- as.character(DD.df[rownames(KB@meta.data)])
table(KB$doubletdetection)
#doublet singlet 
#188216 1492674

KB@meta.data$orig.ident <- sub("\\_.*", "", rownames(KB@meta.data))
table(KB@meta.data$orig.ident)

##Remove Doublets
Idents(KB) <- "doubletdetection"
table(Idents(KB))
KB <- subset(KB, idents = "singlet")
KB
#36080 features across 1492674 samples


#Add in Ribosomal gene proportions
#Ribosomal Protein Gene Database: http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=orglist&org=Homo%20sapiens
ER.GO <- read.delim("misc/human_ribosomal_transcripts.txt",header=FALSE)
ER.GO <- as.character(ER.GO$V1)
ER.GO <- intersect(ER.GO, rownames(KB))
KB[["percent.er"]] <- PercentageFeatureSet(KB, features = ER.GO)

Idents(KB) <- "orig.ident"
VlnPlot(KB, features = c("nFeature_RNA", "nCount_RNA", "percent.er"), 
        ncol = 1, pt.size = -1)

#Add in Mitochondrial gene proportions
KB[["percent.mt"]] <- PercentageFeatureSet(KB, pattern = "^MT-")

VlnPlot(KB, features = c("nFeature_RNA", "percent.er", "percent.mt"), 
        ncol = 1, pt.size = -1)

##Update metadata to include count matrix batches loaded
c("snRNA1", "snRNA2", "Multiome1", "Multiome2")
KB$set <- "snRNA1"
KB@meta.data[rownames(KB@meta.data) %in% colnames(KB2.mat),]$set <- "snRNA2"
KB@meta.data[rownames(KB@meta.data) %in% colnames(KB3.mat),]$set <- "Multiome1"
KB@meta.data[rownames(KB@meta.data) %in% colnames(KB4.mat),]$set <- "Multiome2"
table(KB$set)
#Multiome1 Multiome2    snRNA1    snRNA2 
#338832    327688    403110    423044 

KB[["RNA"]] <- split(KB[["RNA"]], f = KB$set)
Layers(KB[["RNA"]])

saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_Combined_Counts_04-2023_Object.RDS",
  destdir = "~/data/pod/blake_LTS/Atlas_V2/matrices/merged_object"
)

#KB <- readRDS(file = "~/data/pod/blake_LTS/Atlas_V2/matrices/merged_object/Kidney_Atlas_V2_Combined_Counts_04-2023_Object.RDS")


###Remove mitochondrial genes
KB1 <- KB[["RNA"]]$counts.snRNA1
KB2 <- KB[["RNA"]]$counts.snRNA2
KB3 <- KB[["RNA"]]$counts.Multiome1
KB4 <- KB[["RNA"]]$counts.Multiome2

metadata <- KB@meta.data

mt.genes1 <- grep("^MT-", rownames(KB1), value = T)
KB1 <- KB1[!(rownames(KB1) %in% mt.genes1),]
mt.genes2 <- grep("^MT-", rownames(KB2), value = T)
KB2 <- KB2[!(rownames(KB2) %in% mt.genes2),]
mt.genes3 <- grep("^MT-", rownames(KB3), value = T)
KB3 <- KB3[!(rownames(KB3) %in% mt.genes3),]
mt.genes4 <- grep("^MT-", rownames(KB4), value = T)
KB4 <- KB4[!(rownames(KB4) %in% mt.genes4),]

options(Seurat.object.assay.version = "v5")
KB.list <- list(KB1, KB2, KB3, KB4)
names(KB.list) <- c("snRNA1", "snRNA2", "Multiome1", "Multiome2")

options(Seurat.object.assay.version = "v5")
KB <- CreateSeuratObject(counts = KB.list, meta.data = metadata)
KB
#An object of class Seurat 
#36588 features across 1492674 samples within 1 assay 
#Active assay: RNA (36588 features, 0 variable features)
#4 layers present: counts.snRNA1, counts.snRNA2, counts.Multiome1, counts.Multiome2
KB@meta.data

##Remove libraries with less than 100 data sets
to.remove <- names(table(KB$orig.ident)[table(KB$orig.ident) < 100])
Idents(KB) <- "orig.ident"
KB <- subset(KB, idents = to.remove, invert = TRUE)
KB
#36588 features across 1492632 samples

KB[["RNA"]] <- JoinLayers(KB[["RNA"]])
KB[["RNA"]]$counts

saveRDS(
  object = KB,
  file = "Kidney_Atlas_V2_Combined_Counts_04-2023_Object.RDS",
  destdir = "~/data/pod/blake_LTS/Atlas_V2/matrices/merged_object"
)

#KB <- readRDS(file = "~/data/pod/blake_LTS/Atlas_V2/matrices/merged_object/Kidney_Atlas_V2_Combined_Counts_04-2023_Object.RDS")



#Save combined count matrix (filtered and MT transcripts removed)
counts <- KB[["RNA"]]$counts
write_matrix_10x_hdf5(
  counts,
  path = "~/tmp_scratch/Kidney_Atlas_V2_Combined_Counts_04-2023.h5",
  barcodes = colnames(counts),
  feature_ids = rownames(counts),
  feature_names = rownames(counts),
  feature_types = "Gene Expression"
)

saveRDS(KB@meta.data, file = "Kidney_Atlas_V2_Combined_Metadata_04132023.RDS")



###Prepare simplified kidney data objects in scratch directory for fast loading
kidney.data <- open_matrix_10x_hdf5(
  path = "~/tmp_scratch/Kidney_Atlas_V2_Combined_Counts_04-2023.h5"
)

# Write the matrix to a directory
write_matrix_dir(
  mat = kidney.data,
  dir = "~/tmp_scratch/kidney_counts",
  overwrite = TRUE
)

kidney.mat <- open_matrix_dir(dir = "~/tmp_scratch/kidney_counts")
Kid <- CreateSeuratObject(counts = kidney.mat, meta.data = KB@meta.data)
Kid
#An object of class Seurat 
#36588 features across 1492632 samples within 1 assay 
#Active assay: RNA (36588 features, 0 variable features)
#1 layer present: counts

saveRDS(
  object = Kid,
  file = "Kidney_Atlas_V2_04-2023_Object.Rds",
  destdir = "~/tmp_scratch/kidney_object"
)

