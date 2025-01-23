library(Signac)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
library(dplyr)
library(JASPAR2022)
library(TFBSTools)

set.seed(1234)
setwd("~/Projects/Human_Kidney/Atlas_V2/")

###Set up Parallel processing (only run in terminal)
library(future)
plan("multicore", workers = 10)
#plan("multisession")
plan()
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM

# get gene annotations for hg38
load("multiome/Gene_Annotations_hg38.rda")

###Load and combine all cell type peaks
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_PT-TL_peaks_04242024.rda")
pt.peaks <- combined.peaks
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_TAL-DCT_peaks_04242024.rda")
dt1.peaks <- combined.peaks
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_CD_peaks_04242024.rda")
dt2.peaks <- combined.peaks
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_nonEpi_peaks_04242024.rda")
nonEpi.peaks <- combined.peaks

combined.peaks <- GenomicRanges::reduce(x = c(pt.peaks, dt1.peaks, dt2.peaks, nonEpi.peaks))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = blacklist_hg38_unified, invert = TRUE)
combined.peaks
#586864

save(combined.peaks, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/peaks/Kidney_Multiome_ATAC_All-Cell-Type_peaks_04272024.rda")

###Update multiome objects with final RNA QC and metadata
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_04192024.rds")
load("multiome/Kidney_QC_Metric_Metadata_04242024.rda")




###New pEpi object
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_04242024.rda")
KID.PEpi
PT.cells <- colnames(KID.PEpi)[colnames(KID.PEpi) %in% colnames(KB)]
KID.PEpi <- subset(KID.PEpi, cells = PT.cells)
KRAC.PEpi <- subset(KB, cells = PT.cells)
pt.counts <- FeatureMatrix(fragments = Fragments(KID.PEpi),
                           cells = PT.cells,
                           features = combined.peaks,
                           process_n = 20000,
                           sep = c('-','-'))
dim(pt.counts)
#586864 255833

KRAC.PEpi[["ATAC"]] <- CreateChromatinAssay(pt.counts, fragments = Fragments(KID.PEpi))
Annotation(KRAC.PEpi[["ATAC"]]) <- annotation

save(KRAC.PEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_Dual_Seurat_all-peaks_04272024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_Dual_Seurat_all-peaks_04272024.rda")
rm(KID.PEpi)
gc(reset = TRUE)

DefaultAssay(KRAC.PEpi) <- "ATAC"
KRAC.PEpi <- AddMetaData(KRAC.PEpi, metadata = md[rownames(KRAC.PEpi@meta.data),])

# compute nucleosome signal score per cell
KRAC.PEpi <- NucleosomeSignal(KRAC.PEpi)
KRAC.PEpi$nucleosome_group <- ifelse(KRAC.PEpi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(KRAC.PEpi$nucleosome_group) #all < 4
pdf(file='multiome/QC_Plots/Kidney_Multiome_Nucleosome-signal_High-low_04272024_PEpi.pdf',width=4,height=4)
FragmentHistogram(object = KRAC.PEpi, group.by = 'nucleosome_group')
dev.off()

#Compute TSS Enrichment score per cell
KRAC.PEpi <- TSSEnrichment(KRAC.PEpi, fast = FALSE)
KRAC.PEpi$high.tss <- ifelse(KRAC.PEpi$TSS.enrichment > 2, 'High', 'Low')

pdf(file='multiome/QC_Plots/Kidney_KPMP_Multiome_TSS_High-low_04272024_PEpi.pdf',width=8,height=4)
TSSPlot(KRAC.PEpi, group.by = 'high.tss') + NoLegend()
dev.off()

#Compute fraction reads in peaks per cell
KRAC.PEpi <- FRiP(KRAC.PEpi, assay = "ATAC", total.fragments = "atac_fragments",
                  col.name = "FRiP", verbose = TRUE)

Idents(KRAC.PEpi) <- "patient"
pdf(file='multiome/QC_Plots/Kidney_Multiome_Violin_stats_pre-ATAC-QC-Filter_04272024_PEpi.pdf',width=12,height=16)
VlnPlot(
  object = KRAC.PEpi,
  features = c("nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "FRiP"),
  pt.size = -1,
  ncol = 1
)
dev.off()

##Filter low quality
KRAC.PEpi <- subset(
  x = KRAC.PEpi,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 70000 &
    FRiP > 0.25 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
KRAC.PEpi
#623452 features across 227547 samples

Idents(KRAC.PEpi) <- "v2.clusters"
Idents(KRAC.PEpi) <- factor(Idents(KRAC.PEpi), levels = paste0("P_",c(1:33)))
table(Idents(KRAC.PEpi))

save(KRAC.PEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_Dual_Seurat_all-peaks_04272024_B.rda")
rm(KRAC.PEpi)
gc(reset = TRUE)





###New dEpi object Set1
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set1_ATAC_Seurat_04242024.rda")
KID.D1Epi
DT1.cells <- colnames(KID.D1Epi)[colnames(KID.D1Epi) %in% colnames(KB)]
KID.D1Epi <- subset(KID.D1Epi, cells = DT1.cells)
KRAC.D1Epi <- subset(KB, cells = DT1.cells)
dt1.counts <- FeatureMatrix(fragments = Fragments(KID.D1Epi),
                           cells = DT1.cells,
                           features = combined.peaks,
                           process_n = 20000,
                           sep = c('-','-'))
dim(dt1.counts)
#517498 153067

KRAC.D1Epi[["ATAC"]] <- CreateChromatinAssay(dt1.counts, fragments = Fragments(KID.D1Epi))
Annotation(KRAC.D1Epi[["ATAC"]]) <- annotation

save(KRAC.D1Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_04272024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_04272024.rda")
rm(KID.D1Epi)
gc(reset = TRUE)

DefaultAssay(KRAC.D1Epi) <- "ATAC"
KRAC.D1Epi <- AddMetaData(KRAC.D1Epi, metadata = md[rownames(KRAC.D1Epi@meta.data),])

# compute nucleosome signal score per cell
KRAC.D1Epi <- NucleosomeSignal(KRAC.D1Epi)
KRAC.D1Epi$nucleosome_group <- ifelse(KRAC.D1Epi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(KRAC.D1Epi$nucleosome_group) #all < 4
pdf(file='multiome/QC_Plots/Kidney_Multiome_Nucleosome-signal_High-low_04272024_dEpi_TAL-DCT.pdf',width=4,height=4)
FragmentHistogram(object = KRAC.D1Epi, group.by = 'nucleosome_group')
dev.off()

#Compute TSS Enrichment score per cell
KRAC.D1Epi <- TSSEnrichment(KRAC.D1Epi, fast = FALSE)
KRAC.D1Epi$high.tss <- ifelse(KRAC.D1Epi$TSS.enrichment > 2, 'High', 'Low')

pdf(file='multiome/QC_Plots/Kidney_KPMP_Multiome_TSS_High-low_04272024_dEpi_TAL-DCT.pdf',width=8,height=4)
TSSPlot(KRAC.D1Epi, group.by = 'high.tss') + NoLegend()
dev.off()

#Compute fraction reads in peaks per cell
KRAC.D1Epi <- FRiP(KRAC.D1Epi, assay = "ATAC", total.fragments = "atac_fragments",
                  col.name = "FRiP", verbose = TRUE)

Idents(KRAC.D1Epi) <- "patient"
pdf(file='multiome/QC_Plots/Kidney_Multiome_Violin_stats_pre-ATAC-QC-Filter_04272024_dEpi_TAL-DCT.pdf',width=12,height=16)
VlnPlot(
  object = KRAC.D1Epi,
  features = c("nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "FRiP"),
  pt.size = -1,
  ncol = 1
)
dev.off()

##Filter low quality
KRAC.D1Epi <- subset(
  x = KRAC.D1Epi,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 70000 &
    FRiP > 0.25 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
KRAC.D1Epi
#499020 features across 171103 samples

Idents(KRAC.D1Epi) <- "v2.clusters"
Idents(KRAC.D1Epi) <- factor(Idents(KRAC.D1Epi), levels = paste0("D_",c(1:23)))
table(Idents(KRAC.D1Epi))

save(KRAC.D1Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_04272024_B.rda")

rm(KRAC.D1Epi)
gc(reset = TRUE)




###New dEpi object Set2
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Distal_Epithelial_Subclasses_Set2_ATAC_Seurat_04242024.rda")
KID.D2Epi
DT2.cells <- colnames(KID.D2Epi)[colnames(KID.D2Epi) %in% colnames(KB)]
KID.D2Epi <- subset(KID.D2Epi, cells = DT2.cells)
KRAC.D2Epi <- subset(KB, cells = DT2.cells)
dt2.counts <- FeatureMatrix(fragments = Fragments(KID.D2Epi),
                            cells = DT2.cells,
                            features = combined.peaks,
                            process_n = 20000,
                            sep = c('-','-'))
dim(dt2.counts)
#586864 148308

KRAC.D2Epi[["ATAC"]] <- CreateChromatinAssay(dt2.counts, fragments = Fragments(KID.D2Epi))
Annotation(KRAC.D2Epi[["ATAC"]]) <- annotation

save(KRAC.D2Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_CD_Dual_Seurat_all-peaks_04272024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_CD_Dual_Seurat_all-peaks_04272024.rda")
rm(KID.D2Epi)
gc(reset = TRUE)

DefaultAssay(KRAC.D2Epi) <- "ATAC"
KRAC.D2Epi <- AddMetaData(KRAC.D2Epi, metadata = md[rownames(KRAC.D2Epi@meta.data),])

# compute nucleosome signal score per cell
KRAC.D2Epi <- NucleosomeSignal(KRAC.D2Epi)
KRAC.D2Epi$nucleosome_group <- ifelse(KRAC.D2Epi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(KRAC.D2Epi$nucleosome_group) #all < 4
pdf(file='multiome/QC_Plots/Kidney_Multiome_Nucleosome-signal_High-low_04272024_dEpi_CD.pdf',width=4,height=4)
FragmentHistogram(object = KRAC.D2Epi, group.by = 'nucleosome_group')
dev.off()

#Compute TSS Enrichment score per cell
KRAC.D2Epi <- TSSEnrichment(KRAC.D2Epi, fast = FALSE)
KRAC.D2Epi$high.tss <- ifelse(KRAC.D2Epi$TSS.enrichment > 2, 'High', 'Low')

pdf(file='multiome/QC_Plots/Kidney_KPMP_Multiome_TSS_High-low_04272024_dEpi_CD.pdf',width=8,height=4)
TSSPlot(KRAC.D2Epi, group.by = 'high.tss') + NoLegend()
dev.off()

#Compute fraction reads in peaks per cell
KRAC.D2Epi <- FRiP(KRAC.D2Epi, assay = "ATAC", total.fragments = "atac_fragments",
                   col.name = "FRiP", verbose = TRUE)

Idents(KRAC.D2Epi) <- "patient"
pdf(file='multiome/QC_Plots/Kidney_Multiome_Violin_stats_pre-ATAC-QC-Filter_04272024_dEpi_CD.pdf',width=12,height=16)
VlnPlot(
  object = KRAC.D2Epi,
  features = c("nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "FRiP"),
  pt.size = -1,
  ncol = 1
)
dev.off()

##Filter low quality
KRAC.D2Epi <- subset(
  x = KRAC.D2Epi,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 70000 &
    FRiP > 0.25 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
KRAC.D2Epi
#623452 features across 131169 samples

Idents(KRAC.D2Epi) <- "v2.clusters"
Idents(KRAC.D2Epi) <- factor(Idents(KRAC.D2Epi), levels = paste0("D_",c(24:48)))
table(Idents(KRAC.D2Epi))

save(KRAC.D2Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_CD_Dual_Seurat_all-peaks_04272024_B.rda")
rm(KRAC.D2Epi)
gc(reset = TRUE)


###New nonEpi object
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_04242024.rda")
KID.nonEpi
nonEpi.cells <- colnames(KID.nonEpi)[colnames(KID.nonEpi) %in% colnames(KB)]
KID.nonEpi <- subset(KID.nonEpi, cells = nonEpi.cells)
KRAC.nonEpi <- subset(KB, cells = nonEpi.cells)
nonEpi.counts <- FeatureMatrix(fragments = Fragments(KID.nonEpi),
                            cells = nonEpi.cells,
                            features = combined.peaks,
                            process_n = 20000,
                            sep = c('-','-'))
dim(nonEpi.counts)
#586864 133952

KRAC.nonEpi[["ATAC"]] <- CreateChromatinAssay(nonEpi.counts, fragments = Fragments(KID.nonEpi))
Annotation(KRAC.nonEpi[["ATAC"]]) <- annotation

save(KRAC.nonEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_04272024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_04272024.rda")
rm(KID.nonEpi)
gc(reset = TRUE)

DefaultAssay(KRAC.nonEpi) <- "ATAC"
KRAC.nonEpi <- AddMetaData(KRAC.nonEpi, metadata = md[rownames(KRAC.nonEpi@meta.data),])

# compute nucleosome signal score per cell
KRAC.nonEpi <- NucleosomeSignal(KRAC.nonEpi)
KRAC.nonEpi$nucleosome_group <- ifelse(KRAC.nonEpi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(KRAC.nonEpi$nucleosome_group) #all < 4
pdf(file='multiome/QC_Plots/Kidney_Multiome_Nucleosome-signal_High-low_04272024_nonEpi.pdf',width=4,height=4)
FragmentHistogram(object = KRAC.nonEpi, group.by = 'nucleosome_group')
dev.off()

#Compute TSS Enrichment score per cell
KRAC.nonEpi <- TSSEnrichment(KRAC.nonEpi, fast = FALSE)
KRAC.nonEpi$high.tss <- ifelse(KRAC.nonEpi$TSS.enrichment > 2, 'High', 'Low')

pdf(file='multiome/QC_Plots/Kidney_KPMP_Multiome_TSS_High-low_04272024_nonEpi.pdf',width=8,height=4)
TSSPlot(KRAC.nonEpi, group.by = 'high.tss') + NoLegend()
dev.off()

#Compute fraction reads in peaks per cell
KRAC.nonEpi <- FRiP(KRAC.nonEpi, assay = "ATAC", total.fragments = "atac_fragments",
                   col.name = "FRiP", verbose = TRUE)

Idents(KRAC.nonEpi) <- "patient"
pdf(file='multiome/QC_Plots/Kidney_Multiome_Violin_stats_pre-ATAC-QC-Filter_04272024_nonEpi.pdf',width=12,height=16)
VlnPlot(
  object = KRAC.nonEpi,
  features = c("nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "FRiP"),
  pt.size = -1,
  ncol = 1
)
dev.off()

##Filter low quality
KRAC.nonEpi <- subset(
  x = KRAC.nonEpi,
  subset = nCount_ATAC > 1000 &
    nCount_ATAC < 70000 &
    FRiP > 0.25 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
KRAC.nonEpi
#623452 features across 114197 samples within 2 assays
#Active assay: ATAC (586864 features, 0 variable features)
#2 layers present: counts, data
#1 other assay present: RNA
#2 dimensional reductions calculated: pca, umap

Idents(KRAC.nonEpi) <- "v2.clusters"
Idents(KRAC.nonEpi) <- factor(Idents(KRAC.nonEpi), levels = c(paste0("E_",c(1:22)),
                                                              paste0("S_",c(1:27)),
                                                              paste0("I_",c(1:28))))
table(Idents(KRAC.nonEpi))

save(KRAC.nonEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_04272024_B.rda")



###Update Metadata
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.PEpi
cols <- colnames(KRAC.PEpi@meta.data)[colnames(KRAC.PEpi@meta.data) %in% colnames(KB@meta.data)]
KRAC.PEpi@meta.data[,cols] <- KB@meta.data[rownames(KRAC.PEpi@meta.data),cols]
head(KRAC.PEpi)
save(KRAC.PEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_PT-TL_Dual_Seurat_all-peaks_05162024.rda")
rm(KRAC.PEpi)

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.D1Epi
cols <- colnames(KRAC.D1Epi@meta.data)[colnames(KRAC.D1Epi@meta.data) %in% colnames(KB@meta.data)]
KRAC.D1Epi@meta.data[,cols] <- KB@meta.data[rownames(KRAC.D1Epi@meta.data),cols]
head(KRAC.D1Epi)
save(KRAC.D1Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_05162024.rda")
rm(KRAC.D1Epi)

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_CD_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.D2Epi
cols <- colnames(KRAC.D2Epi@meta.data)[colnames(KRAC.D2Epi@meta.data) %in% colnames(KB@meta.data)]
KRAC.D2Epi@meta.data[,cols] <- KB@meta.data[rownames(KRAC.D2Epi@meta.data),cols]
head(KRAC.D2Epi)
save(KRAC.D2Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_CNT-IC_Dual_Seurat_all-peaks_05162024.rda")
rm(KRAC.D2Epi)

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.nonEpi
cols <- colnames(KRAC.nonEpi@meta.data)[colnames(KRAC.nonEpi@meta.data) %in% colnames(KB@meta.data)]
KRAC.nonEpi@meta.data[,cols] <- KB@meta.data[rownames(KRAC.nonEpi@meta.data),cols]
head(KRAC.nonEpi)
save(KRAC.nonEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_05162024.rda")
rm(KRAC.nonEpi)

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_v2_Downsampled_RNA-ATAC_Seurat_05012024.rda")
KRAC
cols <- colnames(KRAC@meta.data)[colnames(KRAC@meta.data) %in% colnames(KB@meta.data)]
KRAC@meta.data[,cols] <- KB@meta.data[rownames(KRAC@meta.data),cols]
head(KRAC)
save(KRAC, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_v2_Downsampled_RNA-ATAC_Seurat_05162024.rda")
rm(KRAC)




###Create a down-sampled combined object
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_pEpi_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.PEpi
meta.PEpi <- KRAC.PEpi@meta.data
KRAC.PEpi <- subset(KRAC.PEpi, downsample = 5000) #cap at 5K per v2 cluster

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_TAL-DCT_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.D1Epi
meta.D1Epi <- KRAC.D1Epi@meta.data
KRAC.D1Epi <- subset(KRAC.D1Epi, downsample = 5000) #cap at 5K per v2 cluster

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_dEpi_CD_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.D2Epi
meta.D2Epi <- KRAC.D2Epi@meta.data
KRAC.D2Epi <- subset(KRAC.D2Epi, downsample = 5000) #cap at 5K per v2 cluster

load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_nonEpi_Dual_Seurat_all-peaks_04272024_B.rda")
KRAC.nonEpi
meta.nonEpi <- KRAC.nonEpi@meta.data
KRAC.nonEpi <- subset(KRAC.nonEpi, downsample = 5000) #cap at 5K per v2 cluster

KRAC <- merge(x = KRAC.PEpi,
              y = list(KRAC.D1Epi,KRAC.D2Epi,KRAC.nonEpi))
KRAC
#An object of class Seurat 
#623452 features across 356821 samples within 2 assays 
#Active assay: ATAC (586864 features, 0 variable features)
#2 layers present: counts, data
#1 other assay present: RNA

#Add in embeddings
umap.embeddings <- rbind(Embeddings(object = KRAC.PEpi[["umap"]]),
                         Embeddings(object = KRAC.D1Epi[["umap"]]),
                         Embeddings(object = KRAC.D2Epi[["umap"]]),
                         Embeddings(object = KRAC.nonEpi[["umap"]]))
umap.embeddings <- umap.embeddings[colnames(KRAC),]
KRAC[["umap"]] <- CreateDimReducObject(embeddings = umap.embeddings, key = "umap_", assay = DefaultAssay(KRAC))
pca.embeddings <- rbind(Embeddings(object = KRAC.PEpi[["pca"]]),
                        Embeddings(object = KRAC.D1Epi[["pca"]]),
                        Embeddings(object = KRAC.D2Epi[["pca"]]),
                        Embeddings(object = KRAC.nonEpi[["pca"]]))
pca.embeddings <- pca.embeddings[colnames(KRAC),]
KRAC[["pca"]] <- CreateDimReducObject(embeddings = pca.embeddings, key = "pca_", assay = DefaultAssay(KRAC))
KRAC

DimPlot(KRAC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l1", repel = TRUE) + NoLegend()

save(KRAC, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_v2_Downsampled_RNA-ATAC_Seurat_05012024.rda")



###Create combined metadata
meta <- rbind(meta.PEpi, meta.D1Epi, meta.D2Epi, meta.nonEpi)
dim(meta)
#693211     89

save(meta, file = "multiome/Kidney_ATAC_QC_MetaData_05012024.rda")


###QC Stats
##Final Stats per library
load("multiome/Kidney_ATAC_QC_MetaData_05012024.rda")
stats <- cbind(
  data.frame(meta %>%
               group_by(library) %>%
               summarise_at(vars(nCount_ATAC,TSS.enrichment,FRiP), list(mean))),
  data.frame(meta %>%
               group_by(library) %>%
               tally())
)
rownames(stats) <- stats$library
stats <- stats[,-c(1,5)]
stats  

write.table(stats, file="multiome/QC_Plots/10X_ATAC_post-QC_Stats_05092024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

##Final Stats per cluster
load("multiome/Kidney_ATAC_QC_MetaData_05012024.rda")
stats <- cbind(
  data.frame(meta %>%
               group_by(v2.clusters) %>%
               summarise_at(vars(nCount_ATAC,TSS.enrichment,FRiP), list(mean))),
  data.frame(meta %>%
               group_by(v2.clusters) %>%
               tally())
)
rownames(stats) <- stats$v2.clusters
stats <- stats[,-c(1,5)]
stats  

write.table(stats, file="multiome/QC_Plots/10X_ATAC_post-QC_Stats_by_Clusters_05092024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Joint UMAP Visualization
#Process ATAC
DefaultAssay(KRAC) <- "ATAC"
KRAC <- FindTopFeatures(KRAC, min.cutoff = 5)
KRAC <- RunTFIDF(KRAC)
KRAC <- RunSVD(KRAC)


# build a joint neighbor graph using both assays
KRAC <- FindMultiModalNeighbors(
  object = KRAC,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
KRAC <- RunUMAP(
  object = KRAC,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.key = "atacumap_",
  verbose = TRUE
)

load("color_factors_v2-clusters.robj")
Idents(object = KRAC) <- "v2.subclass.l3"
pdf(file='UMAP_Plots/Multiome_v2_Subclass.l3_umap.pdf',width=10,height=8)
DimPlot(KRAC, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(KRAC))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()
Idents(object = KRAC) <- "v2.subclass.l1"
pdf(file='UMAP_Plots/Multiome_v2_Subclass.l1_umap.pdf',width=10,height=8)
DimPlot(KRAC, reduction = "umap", label = TRUE, raster=FALSE,alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l1"
        ) + NoLegend() + scale_color_manual(values = v2.scl1.cols[levels(Idents(KRAC))], name = "Subclass"
        ) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
dev.off()



save(KRAC, file = "~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_v2_Downsampled_RNA-ATAC_Seurat_05012024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_v2_Downsampled_RNA-ATAC_Seurat_05012024.rda")
