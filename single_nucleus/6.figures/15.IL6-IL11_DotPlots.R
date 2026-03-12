# IL6 IL11 LIF Expression Plots (TAL / pvFIB) -----------------------------

library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

KB <- readRDS("~/datasets/altos-lab-restricted-project-kpmp/blake_LTS/Atlas_V2/scratch/HKAv2_Seurat_01162026.rds")

KB.TAL <- subset(KB, v2.subclass.l3 %in% c("aTAL1","aTAL2", "C/M-TAL-A","C-TAL-A",
                                           "C/M-TAL-B","C-TAL-B","frTAL"))
Idents(KB.TAL) <- "v2.subclass.l3"
Idents(KB.TAL) <- factor(Idents(KB.TAL), levels = c("aTAL1","aTAL2", "C/M-TAL-A","C-TAL-A",
                                                    "C/M-TAL-B","C-TAL-B","frTAL"))
KB.TAL <- subset(KB.TAL, downsample = 1000)
KB.TAL[["RNA"]] <- as(KB.TAL[["RNA"]], Class = "Assay")
KB.TAL <- NormalizeData(KB.TAL)

KB.FIB <- subset(KB, v2.subclass.l3 %in% c("pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))
Idents(KB.FIB) <- "v2.subclass.l3"
Idents(KB.FIB) <- factor(Idents(KB.FIB), levels = c("pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))
KB.FIB <- subset(KB.FIB, downsample = 1000)
KB.FIB[["RNA"]] <- as(KB.FIB[["RNA"]], Class = "Assay")
KB.FIB <- NormalizeData(KB.FIB)

tal.markers <- c(
  "IL6","IL11","LIF"#,                       #Cytokines
  #"IL6R", "IL6ST", "IL11RA","LIFR"                                        #receptors
)

pdf(file='trajectories/TAL-Lineages_Marker_Dotplot_subclassl3_subset_IL6.pdf',width=5,height=3.2)
DotPlot(KB.TAL, features = tal.markers) + RotatedAxis()
dev.off()


fib.markers <- c(
  #"IL6","IL11","LIF",                       #Cytokines
  "IL6R", #"IL6ST",
  "IL11RA","LIFR"                                        #receptors
)


pdf(file='trajectories/pvFIB-Lineages_Marker_Dotplot_subclassl3_subset_IL6R.pdf',width=5,height=3.2)
DotPlot(KB.FIB, features = fib.markers) + RotatedAxis()
dev.off()


