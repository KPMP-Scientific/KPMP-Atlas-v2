library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm39)
library(JASPAR2022)
library(TFBSTools)
library(ggplot2)
library(patchwork)
library(parallel)

library(foreach)

set.seed(1234)
setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Set up Parallel processing (only run in terminal, not RStudio)
library(future)
plan("multicore", workers = 30)
plan("multisession")
plan()
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

###Load RNA object
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")

###Get gene annotations for mm39
#S3://altos-lab-reference-data/genomes/Mus_musculus/Ensembl/GRCm39/Index/CellrangerArc/genes/genes.gtf.gz
annotation <- rtracklayer::import('genes.gtf')
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm39"
save(annotation, file = "mouse_IRI/multiome/Gene_Annotations_mm39.rda")
#load("mouse_IRI/multiome/Gene_Annotations_mm39.rda")

###Get mm39 blacklist regions from https://doi.org/10.1093/bioinformatics/btad198
library(excluderanges)
suppressMessages(library(GenomicRanges))
suppressMessages(library(AnnotationHub))
ah <- AnnotationHub()
query_data <- subset(ah, preparerclass == "excluderanges")
mm39.excluderanges <- query_data[["AH107321"]]
mm39.excluderanges <- mm39.excluderanges %>% 
  sort() %>% keepStandardChromosomes(pruning.mode = "tidy")
mm39.excluderanges
seqlevels(mm39.excluderanges) <- sub("chr", "", seqlevels(mm39.excluderanges))
save(mm39.excluderanges, file = "mouse_IRI/multiome/mm39.excluderanges.rda")
#load("mouse_IRI/multiome/mm39.excluderanges.rda")

###Prepare combined peak set
ATACfolders <- "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Raw_ATAC/"
exp.table <- read.delim("mouse_IRI/Mouse_Kidney_Exp_Table_112024.txt", row.names = 1)
exps <- exp.table$exp
libs <- exp.table$lib

exps.use <- exps

peak.list <- foreach(exp=exps.use, .final=function(x) setNames(x, exps)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # read Cellranger peaks
  peaks.exp <- read.table(
    file = paste(ATACfolders, exp.table[exp,]$lib, "_atac_peaks.bed", sep = ""),
    col.names = c("chr", "start", "end"))
  
  # convert to genomic ranges
  gr.peaks.exp <- makeGRangesFromDataFrame(peaks.exp)
  
  gr.peaks.exp 
}

length(peak.list)
#46

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(peak.list[[1]],peak.list[[2]],peak.list[[3]],peak.list[[4]],
                               peak.list[[5]],peak.list[[6]],peak.list[[7]],peak.list[[8]],peak.list[[9]],peak.list[[10]],
                               peak.list[[11]],peak.list[[12]],peak.list[[13]],peak.list[[14]],peak.list[[15]],peak.list[[16]],
                               peak.list[[17]],peak.list[[18]],peak.list[[19]],peak.list[[20]],peak.list[[21]],
                               peak.list[[22]],peak.list[[23]],peak.list[[24]],peak.list[[25]],peak.list[[26]],
                               peak.list[[27]],peak.list[[28]],peak.list[[29]],peak.list[[30]],peak.list[[31]],
                               peak.list[[32]],peak.list[[33]],peak.list[[34]],peak.list[[35]],peak.list[[36]],
                               peak.list[[37]],peak.list[[38]],peak.list[[39]],peak.list[[40]],peak.list[[41]],
                               peak.list[[42]],peak.list[[43]],peak.list[[44]],peak.list[[45]],peak.list[[46]]))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

###Create Seurat Objects (slow, consider parallel computing? https://stuartlab.org/signac/articles/future.html)
sobj.list <- foreach(exp=exps, .final=function(x) setNames(x, exps)) %dopar% {
  i <- if(exp %in% exps[[1]]) 1 else 2
  
  # load metadata
  md <- read.table(
    file = paste(ATACfolders, exp.table[exp,]$lib, "_per_barcode_metrics.csv", sep = ""),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )
  
  # Keep only filtered cells
  md <- md[md$is_cell > 0, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = paste(ATACfolders, exp.table[exp,]$lib, "_atac_fragments.tsv.gz", sep = ""),
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


###Merge Objects
# merge all datasets, adding a cell ID to make sure cell names are unique - set 1
KID.ATAC.1 <- merge(
  x = sobj.list[[1]],
  y = sobj.list[c(2:23)],
  add.cell.ids = libs[1:23]
)

KID.ATAC.1[["ATAC"]]
#ChromatinAssay data with 293659 features for 132416 cells
#Variable features: 0
#Genome:
#  Annotation present: FALSE
#Motifs present: FALSE
#Fragment files: 23

save(KID.ATAC.1, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_Multiome_ATAC_11212024_Seurat_set1.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_Multiome_ATAC_11212024_Seurat_set1.rda")


KID.ATAC.2 <- merge(
  x = sobj.list[[24]],
  y = sobj.list[c(25:46)],
  add.cell.ids = libs[24:46]
)

KID.ATAC.2[["ATAC"]]
#ChromatinAssay data with 293659 features for 116684 cells
#Variable features: 0 
#Genome: 
#  Annotation present: FALSE 
#Motifs present: FALSE 
#Fragment files: 23
save(KID.ATAC.2, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_Multiome_ATAC_11212024_Seurat_set2.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_Multiome_ATAC_11212024_Seurat_set2.rda")

rm(sobj.list)

#Free up unused memory
gc(reset = TRUE)




###Update Metadata and Call peaks
meta <- KB@meta.data

load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_Multiome_ATAC_11212024_Seurat_set1.rda")
KID.ATAC.1 <- subset(KID.ATAC.1, cells = colnames(KID.ATAC.1)[colnames(KID.ATAC.1) %in% rownames(meta)])
KID.ATAC.1 <- AddMetaData(KID.ATAC.1, meta[rownames(KID.ATAC.1@meta.data),])

load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_Multiome_ATAC_11212024_Seurat_set2.rda")
KID.ATAC.2 <- subset(KID.ATAC.2, cells = colnames(KID.ATAC.2)[colnames(KID.ATAC.2) %in% rownames(meta)])
KID.ATAC.2 <- AddMetaData(KID.ATAC.2, meta[rownames(KID.ATAC.2@meta.data),])


###Create Seurat Objects using combined peak counts
##Prepare object for Proximal-Intermediate Epithelial subclasses
subclasses <- c("POD","PEC","PT","DTL","ATL")
KA1.PT <- subset(KID.ATAC.1, v2.subclass.l1 %in% subclasses)
KA2.PT <- subset(KID.ATAC.2, v2.subclass.l1 %in% subclasses)

KID.PEpi <- merge(
  x = KA1.PT,
  y = c(KA2.PT)
)
KID.PEpi
#293659 features across 99814 samples

# add the gene information to the object
save(KID.PEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
rm(KID.PEpi)
rm(KA1.PT, KA2.PT)
gc(reset = TRUE)


##Prepare object for Distal Epithelial subclasses 1
subclasses <- c("TAL","DCT","CNT","PC","IC","PapE")
KA1.DT <- subset(KID.ATAC.1, v2.subclass.l1 %in% subclasses)
KA2.DT <- subset(KID.ATAC.2, v2.subclass.l1 %in% subclasses)

KID.D1Epi <- merge(
  x = KA1.DT,
  y = c(KA2.DT)
)
KID.D1Epi
#293659 features across 65373 samples

# add the gene information to the object
save(KID.D1Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Distal_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Distal_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
rm(KID.D1Epi)
rm(KA1.DT, KA2.DT)
gc(reset = TRUE)


##Prepare object for Non-Epithelial subclasses
subclasses <- c("EC","VSM/P","FIB","Ad","Lymphoid","Myeloid","NEU")

KA1.NE <- subset(KID.ATAC.1, v2.subclass.l1 %in% subclasses)
KA2.NE <- subset(KID.ATAC.2, v2.subclass.l1 %in% subclasses)

KID.nonEpi <- merge(
  x = KA1.NE,
  y = c(KA2.NE)
)
KID.nonEpi
#293659 features across 27429 samples

# add the gene information to the object
save(KID.nonEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_11222024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_11222024.rda")
rm(KID.nonEpi)
rm(KA1.NE, KA2.NE)
gc(reset = TRUE)


###Compile QC metrics
ATACfolder <- "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/Raw_ATAC/"
exp.table <- read.delim("mouse_IRI/Mouse_Kidney_Exp_Table_112024.txt", row.names = 1)
exps <- exp.table$exp
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

save(md, file = "mouse_IRI/multiome/Kidney_QC_Metric_Metadata_11222024.rda")
#load("mouse_IRI/multiome/Kidney_QC_Metric_Metadata_11222024.rda")




###Call Peaks
##PT
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
DefaultAssay(KID.PEpi) <- "ATAC"

#cluster peaks
tab <- table(KID.PEpi$v2.clusters)
clusters <- names(tab[tab > 50])
peaks.cl <- CallPeaks(
  object = KID.PEpi,
  group.by = "v2.clusters",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/clusters/",
  cleanup = FALSE
)

save(peaks.cl, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_PT-TL_ATAC_v2_cluster_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_PT-TL_ATAC_v2_cluster_peaks_11212024.rda")


#Subclass.l3
tab <- table(KID.PEpi$v2.subclass.l3)
clusters <- names(tab[tab > 50])
peaks.l3 <- CallPeaks(
  object = KID.PEpi,
  group.by = "v2.subclass.l3",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/subclass_l3/",
  cleanup = FALSE
)

save(peaks.l3, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_PT-TL_ATAC_v2_subclass_l3_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_PT-TL_ATAC_v2_subclass_l3_peaks_11212024.rda")

#Subclass.l1
tab <- table(KID.PEpi$v2.subclass.l1)
clusters <- names(tab[tab > 50])
peaks.l1 <- CallPeaks(
  object = KID.PEpi,
  group.by = "v2.subclass.l1",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/subclass_l1/",
  cleanup = FALSE
)

save(peaks.l1, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_PT-TL_ATAC_v2_subclass_l1_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_PT-TL_ATAC_v2_subclass_l1_peaks_11212024.rda")

# Create a unified set of peaks
combined.peaks <- GenomicRanges::reduce(x = c(peaks.cl, peaks.l3, peaks.l1))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#351764 ranges

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = mm39.excluderanges, invert = TRUE)
combined.peaks
#343948

save(combined.peaks, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_PT-TL_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_PT-TL_peaks_11212024")

rm(KID.PEpi, peaks.cl, peaks.l3, peaks.l1)
gc(reset = TRUE)


##DT
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Distal_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
DefaultAssay(KID.D1Epi) <- "ATAC"

#cluster peaks
tab <- table(KID.D1Epi$v2.clusters)
clusters <- names(tab[tab > 50])
peaks.cl <- CallPeaks(
  object = KID.D1Epi,
  group.by = "v2.clusters",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/clusters/",
  cleanup = FALSE
)

save(peaks.cl, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_TAL-IC_ATAC_v2_cluster_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_TAL-IC_ATAC_v2_cluster_peaks_11212024.rda")


#Subclass.l3
tab <- table(KID.D1Epi$v2.subclass.l3)
clusters <- names(tab[tab > 50])
peaks.l3 <- CallPeaks(
  object = KID.D1Epi,
  group.by = "v2.subclass.l3",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/subclass_l3/",
  cleanup = FALSE
)

save(peaks.l3, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_TAL-IC_ATAC_v2_subclass_l3_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_TAL-IC_ATAC_v2_subclass_l3_peaks_11212024.rda")

#Subclass.l1
tab <- table(KID.D1Epi$v2.subclass.l1)
clusters <- names(tab[tab > 50])
peaks.l1 <- CallPeaks(
  object = KID.D1Epi,
  group.by = "v2.subclass.l1",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/subclass_l1/",
  cleanup = FALSE
)

save(peaks.l1, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_TAL-IC_ATAC_v2_subclass_l1_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_TAL-IC_ATAC_v2_subclass_l1_peaks_11212024.rda")

# Create a unified set of peaks
combined.peaks <- GenomicRanges::reduce(x = c(peaks.cl, peaks.l3, peaks.l1))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#241862 ranges

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = mm39.excluderanges, invert = TRUE)
combined.peaks
#237498

save(combined.peaks, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_TAL-IC_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_TAL-IC_peaks_11212024.rda")

rm(KID.D1Epi, peaks.cl, peaks.l3, peaks.l1)
gc(reset = TRUE)




##nonEpi
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_11222024.rda")
DefaultAssay(KID.nonEpi) <- "ATAC"

#cluster peaks
tab <- table(KID.nonEpi$v2.clusters)
clusters <- names(tab[tab > 50])
peaks.cl <- CallPeaks(
  object = KID.nonEpi,
  group.by = "v2.clusters",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/clusters/",
  cleanup = FALSE
)

save(peaks.cl, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_nonEpi_ATAC_v2_cluster_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_nonEpi_ATAC_v2_cluster_peaks_11212024.rda")


#Subclass.l3
tab <- table(KID.nonEpi$v2.subclass.l3)
clusters <- names(tab[tab > 50])
peaks.l3 <- CallPeaks(
  object = KID.nonEpi,
  group.by = "v2.subclass.l3",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/subclass_l3/",
  cleanup = FALSE
)

save(peaks.l3, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_nonEpi_ATAC_v2_subclass_l3_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_nonEpi_ATAC_v2_subclass_l3_peaks_11212024.rda")

#Subclass.l1
tab <- table(KID.nonEpi$v2.subclass.l1)
clusters <- names(tab[tab > 50])
peaks.l1 <- CallPeaks(
  object = KID.nonEpi,
  group.by = "v2.subclass.l1",
  idents = clusters,
  macs2.path = "/share/blake/anaconda3/envs/macs/bin/macs3",
  outdir = "/weka/blake/Projects/Human_Kidney/Atlas_V2/mouse_IRI/multiome/peaks/bed_files_112024/subclass_l1/",
  cleanup = FALSE
)

save(peaks.l1, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_nonEpi_ATAC_v2_subclass_l1_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_nonEpi_ATAC_v2_subclass_l1_peaks_11212024.rda")

# Create a unified set of peaks
combined.peaks <- GenomicRanges::reduce(x = c(peaks.cl, peaks.l3, peaks.l1))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
#220774

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = mm39.excluderanges, invert = TRUE)
combined.peaks
#217064

save(combined.peaks, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_nonEpi_peaks_11212024.rda")
#load("mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_nonEpi_peaks_11212024.rda")

rm(KID.nonEpi, peaks.cl, peaks.l3, peaks.l1)
gc(reset = TRUE)

