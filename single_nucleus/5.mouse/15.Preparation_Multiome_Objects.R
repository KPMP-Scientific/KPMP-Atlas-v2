library(Signac)
library(Seurat)
library(SeuratDisk)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm39)
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

# get gene annotations for mm39
load("mouse_IRI/multiome/Gene_Annotations_mm39.rda")

#get blacklist regions
load("mouse_IRI/multiome/mm39.excluderanges.rda")

###Load and combine all cell type peaks
load("mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_PT-TL_peaks_11212024.rda")
pt.peaks <- combined.peaks
load("mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_TAL-IC_peaks_11212024.rda")
dt1.peaks <- combined.peaks
load("mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_nonEpi_peaks_11212024.rda")
nonEpi.peaks <- combined.peaks

combined.peaks <- GenomicRanges::reduce(x = c(pt.peaks, dt1.peaks, nonEpi.peaks))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
combined.peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
combined.peaks <- subsetByOverlaps(x = combined.peaks, ranges = mm39.excluderanges, invert = TRUE)
combined.peaks
#504443

save(combined.peaks, file = "mouse_IRI/multiome/peaks/Mouse_Kidney_Multiome_ATAC_All-Cell-Type_peaks_11252024.rda")

###Update multiome objects with final RNA QC and metadata
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
load("mouse_IRI/multiome/Kidney_QC_Metric_Metadata_11222024.rda")




###New pEpi object
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Proximal-Intermediate_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
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
#504443  99814

KRAC.PEpi[["ATAC"]] <- CreateChromatinAssay(pt.counts, fragments = Fragments(KID.PEpi))
Annotation(KRAC.PEpi[["ATAC"]]) <- annotation

save(KRAC.PEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_pEpi_Dual_Seurat_all-peaks_11252024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_pEpi_Dual_Seurat_all-peaks_11252024.rda")
rm(KID.PEpi)
gc(reset = TRUE)

DefaultAssay(KRAC.PEpi) <- "ATAC"
KRAC.PEpi <- AddMetaData(KRAC.PEpi, metadata = md[rownames(KRAC.PEpi@meta.data),])

# compute nucleosome signal score per cell
KRAC.PEpi <- NucleosomeSignal(KRAC.PEpi)
KRAC.PEpi$nucleosome_group <- ifelse(KRAC.PEpi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(KRAC.PEpi$nucleosome_group) 
pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_Multiome_Nucleosome-signal_High-low_11252024_PEpi.pdf',width=4,height=4)
FragmentHistogram(object = KRAC.PEpi, group.by = 'nucleosome_group')
dev.off()

#Compute TSS Enrichment score per cell
KRAC.PEpi <- TSSEnrichment(KRAC.PEpi, fast = FALSE)
KRAC.PEpi$high.tss <- ifelse(KRAC.PEpi$TSS.enrichment > 2, 'High', 'Low')

pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_KPMP_Multiome_TSS_High-low_11252024_PEpi.pdf',width=8,height=4)
TSSPlot(KRAC.PEpi, group.by = 'high.tss') + NoLegend()
dev.off()

#Compute fraction reads in peaks per cell
KRAC.PEpi <- FRiP(KRAC.PEpi, assay = "ATAC", total.fragments = "atac_fragments",
                  col.name = "FRiP", verbose = TRUE)

Idents(KRAC.PEpi) <- "patient"
pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_Multiome_Violin_stats_pre-ATAC-QC-Filter_11252024_PEpi.pdf',width=12,height=16)
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
#535972 features across 76422 samples

Idents(KRAC.PEpi) <- "v2.clusters"
Idents(KRAC.PEpi) <- factor(Idents(KRAC.PEpi), levels = paste0("P_",c(1:31)))
table(Idents(KRAC.PEpi))

save(KRAC.PEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_pEpi_Dual_Seurat_all-peaks_11252024_B.rda")
rm(KRAC.PEpi)
gc(reset = TRUE)





###New dEpi object
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_Distal_Epithelial_Subclasses_ATAC_Seurat_11222024.rda")
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
#504443  65373

KRAC.D1Epi[["ATAC"]] <- CreateChromatinAssay(dt1.counts, fragments = Fragments(KID.D1Epi))
Annotation(KRAC.D1Epi[["ATAC"]]) <- annotation

save(KRAC.D1Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_dEpi_TAL-IC_Dual_Seurat_all-peaks_11252024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_dEpi_TAL-IC_Dual_Seurat_all-peaks_11252024.rda")
rm(KID.D1Epi)
gc(reset = TRUE)

DefaultAssay(KRAC.D1Epi) <- "ATAC"
KRAC.D1Epi <- AddMetaData(KRAC.D1Epi, metadata = md[rownames(KRAC.D1Epi@meta.data),])

# compute nucleosome signal score per cell
KRAC.D1Epi <- NucleosomeSignal(KRAC.D1Epi)
KRAC.D1Epi$nucleosome_group <- ifelse(KRAC.D1Epi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(KRAC.D1Epi$nucleosome_group) #all < 4
pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_Multiome_Nucleosome-signal_High-low_11252024_dEpi_TAL-IC.pdf',width=4,height=4)
FragmentHistogram(object = KRAC.D1Epi, group.by = 'nucleosome_group')
dev.off()

#Compute TSS Enrichment score per cell
KRAC.D1Epi <- TSSEnrichment(KRAC.D1Epi, fast = FALSE)
KRAC.D1Epi$high.tss <- ifelse(KRAC.D1Epi$TSS.enrichment > 2, 'High', 'Low')

pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_KPMP_Multiome_TSS_High-low_11252024_dEpi_TAL-IC.pdf',width=8,height=4)
TSSPlot(KRAC.D1Epi, group.by = 'high.tss') + NoLegend()
dev.off()

#Compute fraction reads in peaks per cell
KRAC.D1Epi <- FRiP(KRAC.D1Epi, assay = "ATAC", total.fragments = "atac_fragments",
                   col.name = "FRiP", verbose = TRUE)

Idents(KRAC.D1Epi) <- "patient"
pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_Multiome_Violin_stats_pre-ATAC-QC-Filter_11252024_dEpi_TAL-IC.pdf',width=12,height=16)
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
#535972 features across 42702 samples

Idents(KRAC.D1Epi) <- "v2.clusters"
Idents(KRAC.D1Epi) <- factor(Idents(KRAC.D1Epi), levels = paste0("D_",c(1:47)))
table(Idents(KRAC.D1Epi))

save(KRAC.D1Epi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_dEpi_TAL-IC_Dual_Seurat_all-peaks_11252024_B.rda")

rm(KRAC.D1Epi)
gc(reset = TRUE)





###New nonEpi object
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Kidney_NonEpithelial_Subclasses_ATAC_Seurat_11222024.rda")
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
#504443  27429

KRAC.nonEpi[["ATAC"]] <- CreateChromatinAssay(nonEpi.counts, fragments = Fragments(KID.nonEpi))
Annotation(KRAC.nonEpi[["ATAC"]]) <- annotation

save(KRAC.nonEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_nonEpi_Dual_Seurat_all-peaks_11252024.rda")
#load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_nonEpi_Dual_Seurat_all-peaks_11252024.rda")
rm(KID.nonEpi)
gc(reset = TRUE)

DefaultAssay(KRAC.nonEpi) <- "ATAC"
KRAC.nonEpi <- AddMetaData(KRAC.nonEpi, metadata = md[rownames(KRAC.nonEpi@meta.data),])

# compute nucleosome signal score per cell
KRAC.nonEpi <- NucleosomeSignal(KRAC.nonEpi)
KRAC.nonEpi$nucleosome_group <- ifelse(KRAC.nonEpi$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
table(KRAC.nonEpi$nucleosome_group) #all < 4
pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_Multiome_Nucleosome-signal_High-low_11252024_nonEpi.pdf',width=4,height=4)
FragmentHistogram(object = KRAC.nonEpi, group.by = 'nucleosome_group')
dev.off()

#Compute TSS Enrichment score per cell
KRAC.nonEpi <- TSSEnrichment(KRAC.nonEpi, fast = FALSE)
KRAC.nonEpi$high.tss <- ifelse(KRAC.nonEpi$TSS.enrichment > 2, 'High', 'Low')

pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_KPMP_Multiome_TSS_High-low_11252024_nonEpi.pdf',width=8,height=4)
TSSPlot(KRAC.nonEpi, group.by = 'high.tss') + NoLegend()
dev.off()

#Compute fraction reads in peaks per cell
KRAC.nonEpi <- FRiP(KRAC.nonEpi, assay = "ATAC", total.fragments = "atac_fragments",
                    col.name = "FRiP", verbose = TRUE)

Idents(KRAC.nonEpi) <- "patient"
pdf(file='mouse_IRI/multiome/QC_Plots/Kidney_Multiome_Violin_stats_pre-ATAC-QC-Filter_11252024_nonEpi.pdf',width=12,height=16)
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
#535972 features across 18908 samples within 2 assays
#Active assay: ATAC (504443 features, 0 variable features)
#2 layers present: counts, data
#1 other assay present: RNA
#2 dimensional reductions calculated: integrated.rpca, umap

Idents(KRAC.nonEpi) <- "v2.clusters"
Idents(KRAC.nonEpi) <- factor(Idents(KRAC.nonEpi), levels = c(paste0("E_",c(1:21)),
                                                              paste0("S_",c(1:25)),
                                                              paste0("I_",c(1:27)),
                                                              paste0("N_",c(1))))
table(Idents(KRAC.nonEpi))

save(KRAC.nonEpi, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_nonEpi_Dual_Seurat_all-peaks_11252024_B.rda")






###Combined object
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_pEpi_Dual_Seurat_all-peaks_11252024_B.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_dEpi_TAL-IC_Dual_Seurat_all-peaks_11252024_B.rda")
load("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_nonEpi_Dual_Seurat_all-peaks_11252024_B.rda")


KRAC <- merge(x = KRAC.PEpi,
              y = list(KRAC.D1Epi,KRAC.nonEpi))
KRAC
#535972 features across 138032 samples within 2 assays 
#Active assay: ATAC (504443 features, 0 variable features)
#2 layers present: counts, data
#1 other assay present: RNA

#Add in embeddings
umap.embeddings <- rbind(Embeddings(object = KRAC.PEpi[["umap"]]),
                         Embeddings(object = KRAC.D1Epi[["umap"]]),
                         Embeddings(object = KRAC.nonEpi[["umap"]]))
umap.embeddings <- umap.embeddings[colnames(KRAC),]
KRAC[["umap"]] <- CreateDimReducObject(embeddings = umap.embeddings, key = "umap_", assay = DefaultAssay(KRAC))
pca.embeddings <- rbind(Embeddings(object = KRAC.PEpi[["integrated.rpca"]]),
                        Embeddings(object = KRAC.D1Epi[["integrated.rpca"]]),
                        Embeddings(object = KRAC.nonEpi[["integrated.rpca"]]))
pca.embeddings <- pca.embeddings[colnames(KRAC),]
KRAC[["integrated.rpca"]] <- CreateDimReducObject(embeddings = pca.embeddings, key = "rpca_", assay = DefaultAssay(KRAC))
KRAC

DimPlot(KRAC, reduction = "umap", pt.size = 0.5, label = TRUE, alpha = 0.1,
        group.by = "v2.subclass.l1", repel = TRUE) + NoLegend()

save(KRAC, file = "~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/multiome/Mouse_Kidney_AtlasV2_Dual_Seurat_11272024.rda")




###QC Stats
##Final Stats per library
meta <-KRAC@meta.data
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

write.table(stats, file="mouse_IRI/multiome/QC_Plots/10X_ATAC_post-QC_Stats_11272024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

##Final Stats per cluster
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

write.table(stats, file="mouse_IRI/multiome/QC_Plots/10X_ATAC_post-QC_Stats_by_Clusters_11272024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
