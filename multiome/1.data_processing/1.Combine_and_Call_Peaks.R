library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(parallel)
library(foreach)
#library(iterators)
#library(doParallel)

set.seed(1234)
setwd("~/Projects/Human_Kidney/Atlas_V2")

###Set up Parallel processing (only run in terminal)
library(future)
plan("multicore", workers = 10)
plan()
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM



# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

save(annotation, file = "multiome/Gene_Annotations_hg38.rda")
load("multiome/Gene_Annotations_hg38.rda")


###Prepare combined peak set
ATACfolder <- "~/hsKidAt/blake_LTS/Atlas_V2/Raw_ATAC/"
exp.table <- read.delim("Data_Processing/Exp_data_table_04222024_multiome.txt", row.names = 1)
exps <- exp.table[exp.table$assay == "10X Multiome",]$exp
exps.use <- exps

peak.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # read Cellranger peaks
  peaks.exp <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_atac_peaks.bed", sep = ""),
    col.names = c("chr", "start", "end"))
  
  # convert to genomic ranges
  gr.peaks.exp <- makeGRangesFromDataFrame(peaks.exp)
  
  gr.peaks.exp 
}

length(peak.list)
#164

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(peak.list[[1]],peak.list[[2]],peak.list[[3]],peak.list[[4]],
                               peak.list[[5]],peak.list[[6]],peak.list[[7]],peak.list[[8]],peak.list[[9]],peak.list[[10]],
                               peak.list[[11]],peak.list[[12]],peak.list[[13]],peak.list[[14]],peak.list[[15]],peak.list[[16]],
                               peak.list[[17]],peak.list[[18]],peak.list[[19]],peak.list[[20]],peak.list[[21]],peak.list[[22]],
                               peak.list[[23]],peak.list[[24]],peak.list[[25]],peak.list[[26]],peak.list[[27]],peak.list[[28]],
                               peak.list[[29]],peak.list[[30]],peak.list[[31]],peak.list[[32]],peak.list[[33]],peak.list[[34]],
                               peak.list[[35]],peak.list[[36]],peak.list[[37]],peak.list[[38]],peak.list[[39]],peak.list[[40]],
                               peak.list[[41]],peak.list[[42]],peak.list[[43]],peak.list[[44]],peak.list[[45]],peak.list[[46]],
                               peak.list[[47]],peak.list[[48]],peak.list[[49]],peak.list[[50]],peak.list[[51]],peak.list[[52]],
                               peak.list[[53]],peak.list[[54]],peak.list[[55]],peak.list[[56]],peak.list[[57]],peak.list[[58]],
                               peak.list[[59]],peak.list[[60]],peak.list[[61]],peak.list[[62]],peak.list[[63]],peak.list[[64]],
                               peak.list[[65]],peak.list[[66]],peak.list[[67]],peak.list[[68]],peak.list[[69]],peak.list[[70]],
                               peak.list[[71]],peak.list[[72]],peak.list[[73]],peak.list[[74]],peak.list[[75]],peak.list[[76]],
                               peak.list[[77]],peak.list[[78]],peak.list[[79]],peak.list[[80]],peak.list[[81]],peak.list[[82]],
                               peak.list[[83]],peak.list[[84]],peak.list[[85]],peak.list[[86]],peak.list[[87]],peak.list[[88]],
                               peak.list[[89]],peak.list[[90]],peak.list[[91]],peak.list[[92]],peak.list[[93]],peak.list[[94]],
                               peak.list[[95]],peak.list[[96]],peak.list[[97]],peak.list[[98]],peak.list[[99]],peak.list[[100]],
                               peak.list[[101]],peak.list[[102]],peak.list[[103]],peak.list[[104]],peak.list[[105]],
                               peak.list[[106]],peak.list[[107]],peak.list[[108]],peak.list[[109]],peak.list[[110]],
                               peak.list[[111]],peak.list[[112]],peak.list[[113]],peak.list[[114]],peak.list[[115]],
                               peak.list[[116]],peak.list[[117]],peak.list[[118]],peak.list[[119]],peak.list[[120]],
                               peak.list[[121]],peak.list[[122]],peak.list[[123]],peak.list[[124]],peak.list[[125]],
                               peak.list[[126]],peak.list[[127]],peak.list[[128]],peak.list[[129]],peak.list[[130]],
                               peak.list[[131]],peak.list[[132]],peak.list[[133]],peak.list[[134]],peak.list[[135]],
                               peak.list[[136]],peak.list[[137]],peak.list[[138]],peak.list[[139]],peak.list[[140]],
                               peak.list[[141]],peak.list[[142]],peak.list[[143]],peak.list[[144]],peak.list[[145]],
                               peak.list[[146]],peak.list[[147]],peak.list[[148]],peak.list[[149]],peak.list[[150]],
                               peak.list[[151]],peak.list[[152]],peak.list[[153]],peak.list[[154]],peak.list[[155]],
                               peak.list[[156]],peak.list[[157]],peak.list[[158]],peak.list[[159]],peak.list[[160]],
                               peak.list[[161]],peak.list[[162]],peak.list[[163]],peak.list[[164]]))


# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks


###Create Seurat Objects 
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_04192024.rds")

##Object 1
# merge all datasets, adding a cell ID to make sure cell names are unique
exps <- exp.table$exp
exps.use <- exps[1:25]
libs <- exp.table[exps.use,]$lib

sobj.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = paste(ATACfolder, exp.table[exp,]$lib, "_atac_fragments.tsv.gz", sep = ""),
    cells = rownames(md)
  )
  
  ###Quantify Peaks in each dataset
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    process_n = 20000,
    cells = rownames(md)
  )
  
  ###Create objects
  assay <- CreateChromatinAssay(counts, fragments = frags)
  sobj <- CreateSeuratObject(assay, assay = "ATAC")
  
  # add information to identify dataset of origin
  sobj$library <- exp.table[exp,]$lib
  
  sobj 
}

KID.ATAC.1 <- merge(
  x = sobj.list[[1]],
  y = sobj.list[c(2:25)],
  add.cell.ids = libs
)

KID.ATAC.1[["ATAC"]]
#ChromatinAssay data with 506736 features for 255715 cells
#Variable features: 0 
#Genome: 
#  Annotation present: FALSE 
#Motifs present: FALSE 
#Fragment files: 25 

#Trim to RNA QC'ed data sets
to.keep <- rownames(KID.ATAC.1@meta.data)[rownames(KID.ATAC.1@meta.data) %in% rownames(KB@meta.data)]
KID.ATAC.1 <- subset(KID.ATAC.1, cells = to.keep)
#506736 features across 146988 samples

save(KID.ATAC.1, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj1.rda")
rm(KID.ATAC.1)
gc(reset = TRUE)


##Object 2
# merge all datasets, adding a cell ID to make sure cell names are unique
exps <- exp.table$exp
exps.use <- exps[26:52]
libs <- exp.table[exps.use,]$lib

sobj.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = paste(ATACfolder, exp.table[exp,]$lib, "_atac_fragments.tsv.gz", sep = ""),
    cells = rownames(md)
  )
  
  ###Quantify Peaks in each dataset
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    process_n = 20000,
    cells = rownames(md)
  )
  
  ###Create objects
  assay <- CreateChromatinAssay(counts, fragments = frags)
  sobj <- CreateSeuratObject(assay, assay = "ATAC")
  
  # add information to identify dataset of origin
  sobj$library <- exp.table[exp,]$lib
  
  sobj 
}

KID.ATAC.2 <- merge(
  x = sobj.list[[1]],
  y = sobj.list[c(2:length(exps.use))],
  add.cell.ids = libs
)

KID.ATAC.2[["ATAC"]]
#ChromatinAssay data with 506736 features for 223806 cells
#Variable features: 0
#Genome:
#  Annotation present: FALSE
#Motifs present: FALSE
#Fragment files: 27

#Trim to RNA QC'ed data sets
to.keep <- rownames(KID.ATAC.2@meta.data)[rownames(KID.ATAC.2@meta.data) %in% rownames(KB@meta.data)]
KID.ATAC.2 <- subset(KID.ATAC.2, cells = to.keep)
#06736 features across 134460 samples

save(KID.ATAC.2, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj2.rda")
rm(KID.ATAC.2)
gc(reset = TRUE)



##Object 3
# merge all datasets, adding a cell ID to make sure cell names are unique
exps <- exp.table$exp
exps.use <- exps[53:78]
libs <- exp.table[exps.use,]$lib

sobj.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = paste(ATACfolder, exp.table[exp,]$lib, "_atac_fragments.tsv.gz", sep = ""),
    cells = rownames(md)
  )
  
  ###Quantify Peaks in each dataset
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    process_n = 20000,
    cells = rownames(md)
  )
  
  ###Create objects
  assay <- CreateChromatinAssay(counts, fragments = frags)
  sobj <- CreateSeuratObject(assay, assay = "ATAC")
  
  # add information to identify dataset of origin
  sobj$library <- exp.table[exp,]$lib
  
  sobj 
}

KID.ATAC.3 <- merge(
  x = sobj.list[[1]],
  y = sobj.list[c(2:length(exps.use))],
  add.cell.ids = libs
)

KID.ATAC.3[["ATAC"]]
#ChromatinAssay data with 461949 features for 132172 cells
#Variable features: 0 
#Genome: 
#  Annotation present: FALSE 
#Motifs present: FALSE 
#Fragment files: 26

#Trim to RNA QC'ed data sets
to.keep <- rownames(KID.ATAC.3@meta.data)[rownames(KID.ATAC.3@meta.data) %in% rownames(KB@meta.data)]
KID.ATAC.3 <- subset(KID.ATAC.3, cells = to.keep)
#506736 features across 74501 samples

save(KID.ATAC.3, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj3.rda")
rm(KID.ATAC.3, sobj.list)
gc(reset = TRUE)




##Object 4
# merge all datasets, adding a cell ID to make sure cell names are unique
exps <- exp.table$exp
exps.use <- exps[79:105]
libs <- exp.table[exps.use,]$lib

sobj.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = paste(ATACfolder, exp.table[exp,]$lib, "_atac_fragments.tsv.gz", sep = ""),
    cells = rownames(md)
  )
  
  ###Quantify Peaks in each dataset
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    process_n = 20000,
    cells = rownames(md)
  )
  
  ###Create objects
  assay <- CreateChromatinAssay(counts, fragments = frags)
  sobj <- CreateSeuratObject(assay, assay = "ATAC")
  
  # add information to identify dataset of origin
  sobj$library <- exp.table[exp,]$lib
  
  sobj 
}

KID.ATAC.4 <- merge(
  x = sobj.list[[1]],
  y = sobj.list[c(2:length(exps.use))],
  add.cell.ids = libs
)

KID.ATAC.4[["ATAC"]]
#ChromatinAssay data with 461949 features for 229807 cells
#Variable features: 0 
#Genome: 
#  Annotation present: FALSE 
#Motifs present: FALSE 
#Fragment files: 27

#Trim to RNA QC'ed data sets
to.keep <- rownames(KID.ATAC.4@meta.data)[rownames(KID.ATAC.4@meta.data) %in% rownames(KB@meta.data)]
KID.ATAC.4 <- subset(KID.ATAC.4, cells = to.keep)
#506736 features across 146988 samples

save(KID.ATAC.4, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj4.rda")
rm(KID.ATAC.4, sobj.list)
gc(reset = TRUE)


##Object 5
# merge all datasets, adding a cell ID to make sure cell names are unique
exps <- exp.table$exp
exps.use <- exps[106:133]
libs <- exp.table[exps.use,]$lib

sobj.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = paste(ATACfolder, exp.table[exp,]$lib, "_atac_fragments.tsv.gz", sep = ""),
    cells = rownames(md)
  )
  
  ###Quantify Peaks in each dataset
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    process_n = 20000,
    cells = rownames(md)
  )
  
  ###Create objects
  assay <- CreateChromatinAssay(counts, fragments = frags)
  sobj <- CreateSeuratObject(assay, assay = "ATAC")
  
  # add information to identify dataset of origin
  sobj$library <- exp.table[exp,]$lib
  
  sobj 
}

KID.ATAC.5 <- merge(
  x = sobj.list[[1]],
  y = sobj.list[c(2:length(exps.use))],
  add.cell.ids = libs
)

KID.ATAC.5[["ATAC"]]
#ChromatinAssay data with 461949 features for 258507 cells
#Variable features: 0 
#Genome: 
#  Annotation present: FALSE 
#Motifs present: FALSE 
#Fragment files: 28

#Trim to RNA QC'ed data sets
to.keep <- rownames(KID.ATAC.5@meta.data)[rownames(KID.ATAC.5@meta.data) %in% rownames(KB@meta.data)]
KID.ATAC.5 <- subset(KID.ATAC.5, cells = to.keep)
#506736 features across 165770 samples

save(KID.ATAC.5, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj5.rda")
rm(KID.ATAC.5, sobj.list)
gc(reset = TRUE)



##Object 6
# merge all datasets, adding a cell ID to make sure cell names are unique
exps <- exp.table$exp
exps.use <- exps[134:164]
libs <- exp.table[exps.use,]$lib

sobj.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = paste(ATACfolder, exp.table[exp,]$lib, "_atac_fragments.tsv.gz", sep = ""),
    cells = rownames(md)
  )
  
  ###Quantify Peaks in each dataset
  counts <- FeatureMatrix(
    fragments = frags,
    features = combined.peaks,
    process_n = 20000,
    cells = rownames(md)
  )
  
  ###Create objects
  assay <- CreateChromatinAssay(counts, fragments = frags)
  sobj <- CreateSeuratObject(assay, assay = "ATAC")
  
  # add information to identify dataset of origin
  sobj$library <- exp.table[exp,]$lib
  
  sobj 
}

KID.ATAC.6 <- merge(
  x = sobj.list[[1]],
  y = sobj.list[c(2:length(exps.use))],
  add.cell.ids = libs
)

KID.ATAC.6[["ATAC"]]
#ChromatinAssay data with 506736 features for 191025 cells
#Variable features: 0 
#Genome: 
#  Annotation present: FALSE 
#Motifs present: FALSE 
#Fragment files: 31

#Trim to RNA QC'ed data sets
to.keep <- rownames(KID.ATAC.6@meta.data)[rownames(KID.ATAC.6@meta.data) %in% rownames(KB@meta.data)]
KID.ATAC.6 <- subset(KID.ATAC.6, cells = to.keep)
#506736 features across 120987 samples

save(KID.ATAC.6, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj6.rda")
rm(KID.ATAC.6, sobj.list)
gc(reset = TRUE)





###Update metadata
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj1.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj2.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj3.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj4.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj5.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Multiome_ATAC_04232024_Seurat_obj6.rda")

meta <- KB@meta.data
KID.ATAC.1 <- AddMetaData(KID.ATAC.1, meta[rownames(KID.ATAC.1@meta.data),])
KID.ATAC.2 <- AddMetaData(KID.ATAC.2, meta[rownames(KID.ATAC.2@meta.data),])
KID.ATAC.3 <- AddMetaData(KID.ATAC.3, meta[rownames(KID.ATAC.3@meta.data),])
KID.ATAC.4 <- AddMetaData(KID.ATAC.4, meta[rownames(KID.ATAC.4@meta.data),])
KID.ATAC.5 <- AddMetaData(KID.ATAC.5, meta[rownames(KID.ATAC.5@meta.data),])
KID.ATAC.6 <- AddMetaData(KID.ATAC.6, meta[rownames(KID.ATAC.6@meta.data),])



###Create Seurat Objects using combined peak counts

##Prepare object for Proximal-Intermediate Epithelial subclasses
subclasses <- c("POD","PEC","PT","DTL","ATL")
KA1.PT <- subset(KID.ATAC.1, v2.subclass.l1 %in% subclasses)
KA2.PT <- subset(KID.ATAC.2, v2.subclass.l1 %in% subclasses)
KA3.PT <- subset(KID.ATAC.3, v2.subclass.l1 %in% subclasses)
KA4.PT <- subset(KID.ATAC.4, v2.subclass.l1 %in% subclasses)
KA5.PT <- subset(KID.ATAC.5, v2.subclass.l1 %in% subclasses)
KA6.PT <- subset(KID.ATAC.6, v2.subclass.l1 %in% subclasses)

KID.PEpi <- merge(
  x = KA1.PT,
  y = c(KA2.PT, KA3.PT, KA4.PT, KA5.PT, KA6.PT)
)
KID.PEpi
#506736 features across 255833 samples

# add the gene information to the object
save(KID.PEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_04242024.rda")
rm(KID.PEpi)
rm(KA1.PT, KA2.PT, KA3.PT, KA4.PT, KA5.PT, KA6.PT)
gc(reset = TRUE)


##Prepare object for Distal Epithelial subclasses 1
subclasses <- c("TAL","DCT")
KA1.DT <- subset(KID.ATAC.1, v2.subclass.l1 %in% subclasses)
KA2.DT <- subset(KID.ATAC.2, v2.subclass.l1 %in% subclasses)
KA3.DT <- subset(KID.ATAC.3, v2.subclass.l1 %in% subclasses)
KA4.DT <- subset(KID.ATAC.4, v2.subclass.l1 %in% subclasses)
KA5.DT <- subset(KID.ATAC.5, v2.subclass.l1 %in% subclasses)
KA6.DT <- subset(KID.ATAC.6, v2.subclass.l1 %in% subclasses)

KID.D1Epi <- merge(
  x = KA1.DT,
  y = c(KA2.DT, KA3.DT, KA4.DT, KA5.DT, KA6.DT)
)
KID.D1Epi
#506736 features across 246468 samples

# add the gene information to the object
save(KID.D1Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set1_ATAC_Seurat_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set1_ATAC_Seurat_04242024.rda")
rm(KID.D1Epi)
rm(KA1.DT, KA2.DT, KA3.DT, KA4.DT, KA5.DT, KA6.DT)
gc(reset = TRUE)


##Prepare object for Distal Epithelial subclasses 2
subclasses <- c("CNT","PC","IC","PapE")
KA1.DT <- subset(KID.ATAC.1, v2.subclass.l1 %in% subclasses)
KA2.DT <- subset(KID.ATAC.2, v2.subclass.l1 %in% subclasses)
KA3.DT <- subset(KID.ATAC.3, v2.subclass.l1 %in% subclasses)
KA4.DT <- subset(KID.ATAC.4, v2.subclass.l1 %in% subclasses)
KA5.DT <- subset(KID.ATAC.5, v2.subclass.l1 %in% subclasses)
KA6.DT <- subset(KID.ATAC.6, v2.subclass.l1 %in% subclasses)

KID.D2Epi <- merge(
  x = KA1.DT,
  y = c(KA2.DT, KA3.DT, KA4.DT, KA5.DT, KA6.DT)
)
KID.D2Epi
#461949 features across 148308 samples

# add the gene information to the object
save(KID.D2Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set2_ATAC_Seurat_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set2_ATAC_Seurat_04242024.rda")
rm(KID.D2Epi)
rm(KA1.DT, KA2.DT, KA3.DT, KA4.DT, KA5.DT, KA6.DT)
gc(reset = TRUE)


##Prepare object for Non-Epithelial subclasses
subclasses <- c("EC","VSM/P","FIB","Ad","Lymphoid","Myeloid","NEU")

KA1.NE <- subset(KID.ATAC.1, v2.subclass.l1 %in% subclasses)
KA2.NE <- subset(KID.ATAC.2, v2.subclass.l1 %in% subclasses)
KA3.NE <- subset(KID.ATAC.3, v2.subclass.l1 %in% subclasses)
KA4.NE <- subset(KID.ATAC.4, v2.subclass.l1 %in% subclasses)
KA5.NE <- subset(KID.ATAC.5, v2.subclass.l1 %in% subclasses)
KA6.NE <- subset(KID.ATAC.6, v2.subclass.l1 %in% subclasses)

KID.nonEpi <- merge(
  x = KA1.NE,
  y = c(KA2.NE, KA3.NE, KA4.NE, KA5.NE, KA6.NE)
)
KID.nonEpi
#506736 features across 133952 samples

# add the gene information to the object
save(KID.nonEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_04242024.rda")


###Compile QC metrics
ATACfolder <- "~/hsKidAt/blake_LTS/Atlas_V2/Raw_ATAC/"
exp.table <- read.delim("Data_Processing/Exp_data_table_04222024_multiome.txt", row.names = 1)
exps <- exp.table[exp.table$assay == "10X Multiome",]$exp
exps.use <- exps

md.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps.use)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolder, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # add information to identify dataset of origin
  rownames(md) <- paste(exp.table[exp,]$lib, rownames(md), sep = "_")
  
  md$barcode <- rownames(md)
  
  md
  
}

md.list[[1]]

md <- do.call("rbind", lapply(md.list, as.data.frame)) 
rownames(md) <- md$barcode

save(md, file = "multiome/Kidney_QC_Metric_Metadata_04242024.rda")
load("multiome/Kidney_QC_Metric_Metadata_04242024.rda")




###Call Peaks
##PT
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_04242024.rda")
DefaultAssay(KID.PEpi) <- "ATAC"

#cluster peaks
tab <- table(KID.PEpi$v2.clusters)
clusters <- names(tab[tab > 50])
peaks.cl <- CallPeaks(
  object = KID.PEpi,
  group.by = "v2.clusters",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/clusters/",
  cleanup = FALSE
)

save(peaks.cl, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_PT-TL_ATAC_v2_cluster_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_PT-TL_ATAC_v2_cluster_peaks.rda")


#Subclass.l3
tab <- table(KID.PEpi$v2.subclass.l3)
clusters <- names(tab[tab > 50])
peaks.l3 <- CallPeaks(
  object = KID.PEpi,
  group.by = "v2.subclass.l3",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l3/",
  cleanup = FALSE
)

save(peaks.l3, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_PT-TL_ATAC_v2_subclass_l3_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_PT-TL_ATAC_v2_subclass_l3_peaks.rda")

#Subclass.l1
tab <- table(KID.PEpi$v2.subclass.l1)
clusters <- names(tab[tab > 50])
peaks.l1 <- CallPeaks(
  object = KID.PEpi,
  group.by = "v2.subclass.l1",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l1/",
  cleanup = FALSE
)

save(peaks.l1, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_PT-TL_ATAC_v2_subclass_l1_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_PT-TL_ATAC_v2_subclass_l1_peaks.rda")

# Create a unified set of peaks
combined.peaks <- GenomicRanges::reduce(x = c(peaks.cl, peaks.l3, peaks.l1))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#288051 ranges

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
combined.peaks
#287644

save(combined.peaks, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_PT-TL_peaks_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_PT-TL_peaks_04242024.rda")

rm(KID.PEpi, peaks.cl, peaks.l3, peaks.l1)
gc(reset = TRUE)


##DT Set1
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set1_ATAC_Seurat_04242024.rda")
DefaultAssay(KID.D1Epi) <- "ATAC"

#cluster peaks
tab <- table(KID.D1Epi$v2.clusters)
clusters <- names(tab[tab > 50])
peaks.cl <- CallPeaks(
  object = KID.D1Epi,
  group.by = "v2.clusters",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/clusters/",
  cleanup = FALSE
)

save(peaks.cl, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_TAL-DCT_ATAC_v2_cluster_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_TAL-DCT_ATAC_v2_cluster_peaks.rda")


#Subclass.l3
tab <- table(KID.D1Epi$v2.subclass.l3)
clusters <- names(tab[tab > 50])
peaks.l3 <- CallPeaks(
  object = KID.D1Epi,
  group.by = "v2.subclass.l3",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l3/",
  cleanup = FALSE
)

save(peaks.l3, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_TAL-DCT_ATAC_v2_subclass_l3_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_TAL-DCT_ATAC_v2_subclass_l3_peaks.rda")

#Subclass.l1
tab <- table(KID.D1Epi$v2.subclass.l1)
clusters <- names(tab[tab > 50])
peaks.l1 <- CallPeaks(
  object = KID.D1Epi,
  group.by = "v2.subclass.l1",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l1/",
  cleanup = FALSE
)

save(peaks.l1, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_TAL-DCT_ATAC_v2_subclass_l1_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_TAL-DCT_ATAC_v2_subclass_l1_peaks.rda")

# Create a unified set of peaks
combined.peaks <- GenomicRanges::reduce(x = c(peaks.cl, peaks.l3, peaks.l1))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#288051 ranges

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
combined.peaks
#287644

save(combined.peaks, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_TAL-DCT_peaks_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_TAL-DCT_peaks_04242024.rda")

rm(KID.D1Epi, peaks.cl, peaks.l3, peaks.l1)
gc(reset = TRUE)


##DT Set2
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set2_ATAC_Seurat_04242024.rda")
DefaultAssay(KID.D2Epi) <- "ATAC"

#cluster peaks
tab <- table(KID.D2Epi$v2.clusters)
clusters <- names(tab[tab > 50])
peaks.cl <- CallPeaks(
  object = KID.D2Epi,
  group.by = "v2.clusters",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/clusters/",
  cleanup = FALSE
)

save(peaks.cl, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_CD_ATAC_v2_cluster_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_CD_ATAC_v2_cluster_peaks.rda")


#Subclass.l3
tab <- table(KID.D2Epi$v2.subclass.l3)
clusters <- names(tab[tab > 50])
peaks.l3 <- CallPeaks(
  object = KID.D2Epi,
  group.by = "v2.subclass.l3",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l3/",
  cleanup = FALSE
)

save(peaks.l3, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_CD_ATAC_v2_subclass_l3_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_CD_ATAC_v2_subclass_l3_peaks.rda")

#Subclass.l1
tab <- table(KID.D2Epi$v2.subclass.l1)
clusters <- names(tab[tab > 50])
peaks.l1 <- CallPeaks(
  object = KID.D2Epi,
  group.by = "v2.subclass.l1",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l1/",
  cleanup = FALSE
)

save(peaks.l1, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_CD_ATAC_v2_subclass_l1_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_CD_ATAC_v2_subclass_l1_peaks.rda")

# Create a unified set of peaks
combined.peaks <- GenomicRanges::reduce(x = c(peaks.cl, peaks.l3, peaks.l1))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#288051 ranges

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
combined.peaks
#287644

save(combined.peaks, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_CD_peaks_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_CD_peaks_04242024.rda")
rm(KID.D2Epi, peaks.cl, peaks.l3, peaks.l1)
gc(reset = TRUE)


##nonEpi
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_04242024.rda")
DefaultAssay(KID.nonEpi) <- "ATAC"

#cluster peaks
tab <- table(KID.nonEpi$v2.clusters)
clusters <- names(tab[tab > 50])
peaks.cl <- CallPeaks(
  object = KID.nonEpi,
  group.by = "v2.clusters",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/clusters/",
  cleanup = FALSE
)

save(peaks.cl, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_nonEpi_ATAC_v2_cluster_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_nonEpi_ATAC_v2_cluster_peaks.rda")


#Subclass.l3
tab <- table(KID.nonEpi$v2.subclass.l3)
clusters <- names(tab[tab > 50])
peaks.l3 <- CallPeaks(
  object = KID.nonEpi,
  group.by = "v2.subclass.l3",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l3/",
  cleanup = FALSE
)

save(peaks.l3, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_nonEpi_ATAC_v2_subclass_l3_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_nonEpi_ATAC_v2_subclass_l3_peaks.rda")

#Subclass.l1
tab <- table(KID.nonEpi$v2.subclass.l1)
clusters <- names(tab[tab > 50])
peaks.l1 <- CallPeaks(
  object = KID.nonEpi,
  group.by = "v2.subclass.l1",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/bed_files/subclass_l1/",
  cleanup = FALSE
)

save(peaks.l1, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_nonEpi_ATAC_v2_subclass_l1_peaks.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_nonEpi_ATAC_v2_subclass_l1_peaks.rda")

# Create a unified set of peaks
combined.peaks <- GenomicRanges::reduce(x = c(peaks.cl, peaks.l3, peaks.l1))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#288051 ranges

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
combined.peaks
#287644

save(combined.peaks, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_nonEpi_peaks_04242024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_nonEpi_peaks_04242024.rda")

rm(KID.nonEpi, peaks.cl, peaks.l3, peaks.l1)
gc(reset = TRUE)

