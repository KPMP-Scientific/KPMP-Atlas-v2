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

setwd("~/Projects/Human_Kidney/Atlas_V2")
source("misc/utils.R")

###Append Library IDs to Cell Barcodes
#snRNA-seq Experiments
processed.dir <- "~/hsKidAt/processed/cellranger-arc/"
dir2 <- "~/hsKidAt/blake_LTS/Atlas_V2/Raw_Counts/" 

exp.table <- read.delim("Data_Processing/New_exp_data_table_04112024.txt", row.names = 1)
exps <- exp.table[exp.table$assay == "10X snRNA-seq",]$exp
libs <- exp.table[exp.table$assay == "10X snRNA-seq",]$lib

#Read in counts, append library ID and save

foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exp[[1]]) 1 else 2
  
  lib = exp.table[exp,]$lib
  folder = exp.table[exp,]$folder
  dir1 = paste0(processed.dir,folder)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  #counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-R_UMI_counts_EmptyBC_Filter.mtx"))
}


#Get Cellranger QC stats
QCpath = paste("~/hsKidAt/processed/cellranger-arc/",exp.table[exps,]$folder, sep = "")
QC.list = paste(QCpath,"/outs/metrics_summary.csv", sep = "")
QC_files_df <- lapply(QC.list, function(x) {read.csv(file = x)})
QC.df <- do.call(rbind, QC_files_df)
rownames(QC.df) <- exps
write.table(QC.df, file="QC_Plots/Cellranger_QC_Table_04-2024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)




#Multiome Experiments
processed.dir <- "~/hsKidAt/processed/cellranger-arc/"
dir2 <- "~/hsKidAt/blake_LTS/Atlas_V2/Raw_Counts/" 

exp.table <- read.delim("Data_Processing/New_exp_data_table_04112024.txt", row.names = 1)
exps <- exp.table[exp.table$assay == "10X Multiome",]$exp
libs <- exp.table[exp.table$assay == "10X Multiome",]$lib

#Read in counts, append library ID and save
foreach(exp=exps, .final=function(x) setNames(x, exp)) %dopar% {
  i <- if(exp %in% exp[[1]]) 1 else 2
  
  lib = exp.table[exp,]$lib
  folder = exp.table[exp,]$folder
  dir1 = paste0(processed.dir,folder)
  counts <- Read10X_h5(paste0(dir1,"/outs/filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  saveRDS(counts, file = paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.mtx"))
}

#Get Cellranger QC stats
QCpath = paste("~/hsKidAt/processed/cellranger-arc/",exp.table[exps,]$folder, sep = "")
QC.list = paste(QCpath,"/outs/summary.csv", sep = "")
QC_files_df <- lapply(QC.list, function(x) {read.csv(file = x, row.names = 1)})
QC.df <- do.call(rbind, QC_files_df)
write.table(QC.df, file="QC_Plots/Cellranger-arc_QC_Table_04-2024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



###DoubletDetection Filter - ran on python  

dir1 = "~/hsKidAt/blake_LTS/Atlas_V2/Raw_Counts/"
dir2 = "~/hsKidAt/blake_LTS/Atlas_V2/DoubletDetection/"

#Prepare doublet factors - 10X-R Experiments
exp.table <- read.delim("Data_Processing/New_exp_data_table_04112024.txt", row.names = 1)
exps <- exp.table[exp.table$assay == "10X snRNA-seq",]$exp
libs <- exp.table[exp.table$assay == "10X snRNA-seq",]$lib

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
exps <- exp.table[exp.table$assay == "10X Multiome",]$exp
libs <- exp.table[exp.table$assay == "10X Multiome",]$lib

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
CNTpath = "~/hsKidAt/blake_LTS/Atlas_V2/Raw_Counts/"
CNT.list = list.files(path=CNTpath, pattern="*UMI_counts_EmptyBC_Filter.rds") 
CNT.list <- CNT.list[CNT.list %in% c(paste0(exp.table[exp.table$assay == "10X Multiome",]$exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"),
                                     paste0(exp.table[exp.table$assay == "10X snRNA-seq",]$exp, "_10X-R_UMI_counts_EmptyBC_Filter.rds"))]

CNT.list = paste(CNTpath,CNT.list, sep = "")
CNT_files_df <- lapply(CNT.list, function(x) {readRDS(file = x)})
KB <- merge.sparse(CNT_files_df)


##Create Seurat Object and Apply Gene/UMI Filters
##Prepare H5AD files for all genes and merge 
KB <- convert_matrix_type(KB, type = "uint32_t")

# Output BPCells matrices on-disk
write_matrix_dir(
  mat = KB,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2024_BP",
  overwrite = TRUE
)

# Load in BP matrices
KB.mat <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/matrices/Kidney_Atlas_V2_Combined_Counts_04-2024_BP")
options(Seurat.object.assay.version = "v5")
KB <- CreateSeuratObject(counts = KB.mat)
KB
#An object of class Seurat 
#36601 features across 519355 samples within 1 assay

##Gene Filter >400 and <7500 genes
KB <- subset(x = KB, subset = nFeature_RNA > 400 & nFeature_RNA < 7500)
FeatureScatter(object = KB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
KB
#36601 features across 481565

##Add in doublet information
DDpath = "~/hsKidAt/blake_LTS/Atlas_V2/DoubletDetection/"
DD.list = list.files(path=DDpath, pattern="*_DD_Factor.rds") 
DD.list = paste(DDpath,DD.list, sep = "")
DD_files_df <- lapply(DD.list, function(x) {readRDS(file = x)})
DD.df <- unlist(DD_files_df)

KB@meta.data$doubletdetection <- as.character(DD.df[rownames(KB@meta.data)])
table(KB$doubletdetection)
#doublet singlet 
#59395  422050

table(KB@meta.data$orig.ident)

##Remove Doublets
Idents(KB) <- "doubletdetection"
table(Idents(KB))
KB <- subset(KB, idents = "singlet")
KB
#36601 features across 422050 samples


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


###Remove mitochondrial genes
KB1 <- KB[["RNA"]]$counts
metadata <- KB@meta.data

mt.genes1 <- grep("^MT-", rownames(KB1), value = T)
KB1 <- KB1[!(rownames(KB1) %in% mt.genes1),]

options(Seurat.object.assay.version = "v5")
KB <- CreateSeuratObject(counts = KB1, meta.data = metadata)
KB
#An object of class Seurat 
#36588 features across 422050 samples within 1 assay 
KB@meta.data

##Remove libraries with less than 100 data sets
to.remove <- names(table(KB$orig.ident)[table(KB$orig.ident) < 100])
Idents(KB) <- "orig.ident"
KB <- subset(KB, idents = to.remove, invert = TRUE)
KB
#36588 features across 421974 samples

#Save Seurat Object
# Write the matrix to a directory
write_matrix_dir(
  mat = KB[["RNA"]]$counts,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/0424_kidney_count_subset",
  overwrite = TRUE
)

kidney.mat <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/0424_kidney_count_subset")
KB[["RNA"]]$counts <- kidney.mat
KB
#An object of class Seurat 
#36588 features across 421974 samples within 1 assay 
#Active assay: RNA (36588 features, 0 variable features)
#1 layer present: counts

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_Atlas_V2_04-2024-subset_Object.Rds")
