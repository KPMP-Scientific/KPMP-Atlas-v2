library(Seurat)
library(dplyr)
library(Matrix)
library(parallel)
library(pagoda2)
library(foreach)
library(doParallel)
library(BPCells)

registerDoParallel(cores=4)
options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
source("misc/utils.R")

#Experiments
exp.table <- read.delim("mouse_IRI/Mouse_Kidney_Exp_Table_082024.txt", row.names = 1)
exps <- exp.table$exp
libs <- exp.table$lib

exps.use <- exps

#Read in counts, append library ID and save
foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  lib = exp.table[exp,]$lib
  dir1 = exp.table[exp,]$s3_folder
  dir2 <- "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/RNA_RDS/" 
  
  counts <- Read10X_h5(paste0(dir1,"filtered_feature_bc_matrix.h5"), use.names = TRUE, unique.features = TRUE)
  counts <- counts$`Gene Expression`
  colnames(counts) <- paste(lib, colnames(counts), sep="_")
  
  print(dim(counts))
  
  saveRDS(counts, file = paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.rds"))
  writeMM(counts, file= paste0(dir2, exp, "_10X-Dual_UMI_counts_EmptyBC_Filter.mtx"))
  
}


###DoubletDetection Filter
dir1 = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/RNA_RDS/"
dir2 = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/DoubletDetection/"

#Combine Samples
exp.table <- read.delim("mouse_IRI/Mouse_Kidney_Exp_Table_082024.txt", row.names = 1)
exps <- exp.table$exp
exps <- exps

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



###Combining kidney count matrices 
KIDpath = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/RNA_RDS/"
KID_CM.list = list.files(path=KIDpath, pattern="*_10X-Dual_UMI_counts_EmptyBC_Filter.rds") 
KID_CM.list = paste(KIDpath,"/",KID_CM.list, sep = "")

CM_files_df <- lapply(KID_CM.list, function(x) {readRDS(file = x)})
KB <- merge.sparse(CM_files_df)
dim(KB)
#34183 249100

##Create Seurat Object and Apply Gene/UMI Filters
##Prepare H5AD files for all genes and merge 
KB <- convert_matrix_type(KB, type = "uint32_t")

# Output BPCells matrices on-disk
write_matrix_dir(
  mat = KB,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Kidney_Atlas_V1_Combined_Counts_08-2024_BP",
  overwrite = TRUE
)

# Load in BP matrices
KB.mat <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Kidney_Atlas_V1_Combined_Counts_08-2024_BP")
options(Seurat.object.assay.version = "v5")
KB <- CreateSeuratObject(counts = KB.mat)
KB
#34183 features across 249100 samples

#Gene Filter >400 and <5000 genes
KB <- subset(x = KB, subset = nFeature_RNA > 400 & nFeature_RNA < 5000)
FeatureScatter(object = KB, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()
KB
#34183 features across 231944 samples


###Gene Molecule Filter
countMatrix <- KB[["RNA"]]$counts
countMatrix <- as(object = countMatrix, Class = "dgCMatrix")
dim(countMatrix)
KB.gmcf <- gene.vs.molecule.cell.filter(countMatrix,min.cell.size=200)
KB <- subset(KB, cells = colnames(KB.gmcf))
KB
#34183 features across 229115 samples


#Add in Ribosomal gene proportions
#Ribosomal Protein Gene Database: http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=orglist&org=Mus%20musculus
ER.GO <- read.delim("~/Projects/Mouse_Kidney/Atlas_V1/misc/mouse_ribosomal_transcripts.txt",sep="\t",header=FALSE)
ER.GO <- as.character(ER.GO$V1)
ER.GO <- intersect(ER.GO, rownames(KB))
KB[["percent.er"]] <- PercentageFeatureSet(KB, features = ER.GO)
VlnPlot(KB, features = c("nFeature_RNA", "nCount_RNA", "percent.er"), 
        ncol = 1, pt.size = -1)


#Add in Mitochondrial gene proportions
KB[["percent.mt"]] <- PercentageFeatureSet(KB, pattern = "^mt-")

VlnPlot(KB, features = c("nFeature_RNA", "percent.er", "percent.mt"), 
        ncol = 1, pt.size = -1)

save(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Combined_082024_Seurat.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Intermediate_Objects/Mouse_Kidney_Combined_082024_Seurat.rda")


##Add in doublet information
DDpath = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/DoubletDetection/"
DD.list = list.files(path=DDpath, pattern="*_DD_Factor.rds") 
DD.list = paste(DDpath,DD.list, sep = "")
DD_files_df <- lapply(DD.list, function(x) {readRDS(file = x)})
DD.df <- unlist(DD_files_df)

KB@meta.data$doubletdetection <- as.character(DD.df[rownames(KB@meta.data)])
table(KB$doubletdetection)
#doublet singlet 
#19173  209940 

##Remove Doublets
KB <- subset(KB, doubletdetection %in% "singlet")
KB
#34183 features across 209940 samples


###Remove mitochondrial genes
KB.counts <- KB[["RNA"]]$counts
metadata <- KB@meta.data
mt.genes <- grep("^mt-", rownames(KB.counts), value = T)
KB.counts <- KB.counts[!(rownames(KB.counts) %in% mt.genes),]

# Write the matrix to a directory
write_matrix_dir(
  mat = KB.counts,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/kidney_counts_082024",
  overwrite = TRUE
)

kidney.mat <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/kidney_counts_082024")
mmKidAt <- CreateSeuratObject(counts = kidney.mat, meta.data = metadata)
mmKidAt
#An object of class Seurat 
#34170 features across 209940 samples within 1 assay 
#Active assay: RNA (34170 features, 0 variable features)
#1 layer present: counts
saveRDS(mmKidAt, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Atlas_V2_08-2024_Object.Rds")
