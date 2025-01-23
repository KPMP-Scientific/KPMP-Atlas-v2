library(Seurat)
library(SeuratData)
library(Matrix)
library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(foreach)
library(BPCells)
library(spacexr)
library("corrplot")
library(tidyr)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")
load("color_factors_v2-clusters.robj")

###Load and combine objects
load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_nonEpi_subset_07212024.Rda")
kss.non.epi <- kss.sub
load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_pEpi_subset_07232024.Rda")
kss.p.epi <- kss.sub
load("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_dEpi_subset_07282024.Rda")
kss.d.epi <- kss.sub
rm(kss.sub)

kss <- merge(x = kss.non.epi, y = c(kss.p.epi,kss.d.epi))

table(kss$maxCellType.l3)

#Update annotations
meta <- kss@meta.data
cl.meta <- read.delim("slide-seq/slide-seq_annotation_table_07292024.txt")
emc <- c("v2.subclass.sp","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$maxCellType.l3, cl.meta$maxCellType.l3)]
}
colnames(meta)
table(meta$v2.subclass.sp)
table(kss@meta.data$maxCellType.l3)
#export l3 prediction assays
l3.predictions <- kss[["l3.predictions"]]$counts

kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_filtered_09-02-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09-2023_Object_filtered.Rds")
DefaultAssay(kss) <- "Spatial"

kss@meta.data <- meta[rownames(kss@meta.data),]

length(colnames(l3.predictions))
kss <- subset(kss, cells = colnames(l3.predictions))
kss[["l3.predictions"]] <- CreateAssayObject(counts = l3.predictions)

order <- c("POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","DTL1","DTL3","ATL",
           "M-TAL","C-TAL","MD","aTAL1","aTAL2","frTAL","DCT","aDCT","CNT","aCNT",
           "PC","IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR",
           "EC-AVR","infEC-AVR","EC-V","EC-PCV","EC-PTC","angEC-PTC","EC-EA",
           "infEC-PTC","EC-LYM","M-FIB","C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-Path",
           "C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+","pvFIB-PI16+",
           "pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
           "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+",
           "resMAC-HLAIIhi","moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON",
           "mDC","cDC1","pDC","N")
kss <- subset(kss, v2.subclass.sp %in% order)
Idents(kss) <- "v2.subclass.sp"
Idents(kss) <- factor(Idents(kss), levels = order)
table(Idents(kss))

##Subset low prediction scores
kss <- subset(kss, maxWeight.l3 > 25)

##Save Filtered subset
saveRDS(kss, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_Filtered_Object_07292024.Rds")
#kss
#An object of class Seurat 
#103153 features across 1623079 samples within 4 assays 
#Active assay: Spatial (51531 features, 0 variable features)
#2 layers present: counts, data
#3 other assays present: l1.predictions, sketch, l3.predictions
##5 dimensional reductions calculated: pca, umap, pca.full, full.umap, full.umap.1
#71 images present



###Subset out spots outside tissue using excluding lines
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_Filtered_Object_07292024.Rds")
source("misc/utils.R")

#structure colors
meta <- kss@meta.data
meta <- meta[!duplicated(meta$v2.subclass.sp),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "v2.structure", sort_label = F, colorset = "varibow")

v2.str.cols <- meta$v2.structure_color; names(v2.str.cols) <- meta$v2.structure_label
Idents(kss) <- "v2.structure"

#Function for selecting spots between defined lines
IsBetweenLines <- function(coord, lineList){
  # given two lines in lineList and giotto object, returns a boolean vector whether
  # the beads fall between the two lines. Line names have to be "topLine" and "botLine".
  # Each line is specified by two points, p1 and p2
  isBet <- rep_len(FALSE, length(rownames(coord)))
  locs.df = coord
  if('topLine' %in% names(lineList)){
    topLine = lineList[['topLine']]
    m = (topLine$p2['y'] - topLine$p1['y']) / (topLine$p2['x'] - topLine$p1['x'])
    # y = m * (x - x1) + y1
    locs.df[['topLine']] <- locs.df$y - (m * (locs.df$x - topLine$p1['x']) + topLine$p1['y'])
    locs.df[['isUnderTopLine']] <- locs.df[['topLine']] <= 0
  } else
    locs.df[['isUnderTopLine']] <- TRUE
  if('botLine' %in% names(lineList)){
    botLine = lineList[['botLine']]
    m = (botLine$p2['y'] - botLine$p1['y']) / (botLine$p2['x'] - botLine$p1['x'])
    # y = m * (x - x1) + y1
    locs.df[['botLine']] <- locs.df$y - (m * (locs.df$x - botLine$p1['x']) + botLine$p1['y'])
    locs.df[['isAboveBotLine']] <- locs.df[['botLine']] >= 0
  }
  else
    locs.df[['isAboveBotLine']] <- TRUE
  isBet <- locs.df[['isAboveBotLine']] & locs.df[['isUnderTopLine']]
  return(isBet)
}

kss$inBoundary <- "false"

#Example for Puck_210113_21
puck <- "Puck_210113_21"
coord <- kss@images[[puck]]@coordinates[,1:2]
p1 <- setNames(c(0,3500),c("x","y"))
p2 <- setNames(c(4200,4000),c("x","y"))
topLine <- list(p1,p2) ; names(topLine) <- c("p1", "p2")
p1 <- setNames(c(0,500),c("x","y"))
p2 <- setNames(c(4400,1400),c("x","y"))
botLine <- list(p1,p2) ; names(botLine) <- c("p1", "p2")
Puck_210113_21 <- list(topLine,botLine) ; names(Puck_210113_21) <- c("topLine","botLine")

cells <- rownames(coord[IsBetweenLines(coord, Puck_210113_21),])
SpatialDimPlot(kss,cells.highlight = cells, images = puck)
SpatialDimPlot(kss, stroke = 0, images = puck, cols = v2.str.cols[levels(Idents(kss))])
kss@meta.data[cells,]$inBoundary <- "true"

#repeated for all pucks, then removed spots that fall outside the boundaries
kss <- subset(kss, inBoundary %in% "true")

##Save Filtered subset
saveRDS(kss, file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_Filtered_Object_07292024_InBoundary.Rds")
kss
#An object of class Seurat 
#103153 features across 1512878 samples within 4 assays 
#Active assay: Spatial (51531 features, 0 variable features)
#2 layers present: counts, data
#3 other assays present: l1.predictions, sketch, l3.predictions
#5 dimensional reductions calculated: pca, umap, pca.full, full.umap, full.umap.1
#71 images present
