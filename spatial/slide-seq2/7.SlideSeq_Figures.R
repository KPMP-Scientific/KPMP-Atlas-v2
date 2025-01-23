library(Seurat)
library(SeuratData)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(BPCells)
library("corrplot")
library(tidyr)
library(Polychrome)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")
load("color_factors_v2-clusters.robj")
source("misc/utils.R")

###Reference Correlation Plots
#Prepare Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(ref) <- "v2.subclass.sp"
Idents(ref) <- factor(Idents(ref), levels = c(
  "POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","aDTL","DTL1","DTL3","ATL",
  "M-TAL","C-TAL","MD","frTAL","aTAL1","aTAL2","DCT","aDCT","CNT","aCNT","PC",
  "IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
  "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM","M-FIB","C/M-FIB",
  "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
  "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
  "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
  "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N","-"))
ref <- subset(ref, idents = c("-"), invert = TRUE)
ref <- subset(ref, downsample = 1000)

load("AtlasV2_VariableFeatures.rda")
od.genes <- unique(c(pt.od.genes[1:500],dt.od.genes[1:500],imm.od.genes[1:500],
                     str.od.genes[1:500],ec.od.genes[1:500]))
sn.od.genes <- od.genes[od.genes %in% rownames(kss)]

ref[["RNA"]] <- as(object = ref[["RNA"]], Class = "Assay")
ref <- NormalizeData(ref)
ref <- ScaleData(ref, features = sn.od.genes, assay = "RNA")

#prepare query
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09042024_Final_Object.Rds")
kss.sub <- subset(kss, downsample = 1000)
table(kss.sub$v2.subclass.sp)
kss.sub[["Spatial"]] <- as(object = kss.sub[["Spatial"]], Class = "Assay")
kss.sub[["Spatial"]] <- NormalizeData(kss.sub[["Spatial"]])
kss.sub[["Spatial"]] <- ScaleData(kss.sub[["Spatial"]], features = sn.od.genes)

ave.ref <- AverageExpression(ref, features = od.genes, assays = "RNA", slot = "scale.data")
ave.kss <- AverageExpression(kss.sub, features = od.genes, assay = "Spatial", slot = "scale.data")

order <- c("POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","DTL1","DTL3","ATL",
           "M-TAL","C-TAL","MD","aTAL1","aTAL2","frTAL","DCT","aDCT","CNT","aCNT",
           "PC","IMCD","IC-A","tPC-IC","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-AVR","infEC-AVR","EC-V","EC-PCV","EC-PTC","angEC-PTC","EC-EA",
           "infEC-PTC","EC-LYM","M-FIB","C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-Path",
           "C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+","pvFIB-PI16+",
           "pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
           "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+",
           "resMAC-HLAIIhi","moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON",
           "mDC","cDC1","pDC","N")

ave.cor<-cor(as.data.frame(ave.kss$Spatial[,order]),as.data.frame(ave.ref$RNA[rownames(ave.kss$Spatial),order]))
colnames(ave.cor) <- gsub("RNA.", "", colnames(ave.cor))
rownames(ave.cor) <- gsub("Spatial.", "", rownames(ave.cor))
write.table(ave.cor, file="slide-seq/QC_Plots/10X_RNA_v2-subclasses_versus_slide-seq_pred-subclasses_corr_all.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

pdf(file='slide-seq/QC_Plots/10X_RNA_v2-subclasses_versus_slide-seq_pred-subclasses_corr_Plot_all.pdf',width=14,height=14)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor, col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))
dev.off()



###Final Post-V2 annotation/alignment Stats
meta <- kss@meta.data
stats <- cbind(
  data.frame(meta %>%
               group_by(library) %>%
               summarise_at(vars(nCount_Spatial,nFeature_Spatial), list(mean))),
  data.frame(meta %>%
               group_by(library) %>%
               tally())
)
rownames(stats) <- stats$library
stats <- stats[,-c(1,4)]
stats  

#Determine number of cell types
df <- t(table(kss$library, kss$v2.subclass.sp))
cols <- vector()
exps <- colnames(df)

for(i in exps){
  cols[i] <- length(df[which(df[,i] > 5),i]) 
}
stats$n_clusters <- cols[rownames(stats)]

write.table(stats, file="slide-seq/QC_Plots/Slide-Seq_post-labeling_Stats_09112024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###UMAP Plot
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09042024_Final_Object.Rds")
Idents(kss) <- "v2.subclass.l1"
DefaultAssay(kss) <- "spatial"
pdf(file='slide-seq/Plots/Slide-seq_v2_Subclass.l1_umap.pdf',width=6,height=6)
DimPlot(kss, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.4,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl1.cols[levels(Idents(kss))], name = "Subclass") 
dev.off()
pdf(file='slide-seq/Plots/Slide-seq_v2_Subclass.l1_umap_full.pdf',width=6,height=6)
DimPlot(kss, reduction = "full.umap", label = TRUE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl1.cols[levels(Idents(kss))], name = "Subclass") 
dev.off()
pdf(file='slide-seq/Plots/Slide-seq_v2_Subclass.l1_umap_full_unlabeled.pdf',width=6,height=6)
DimPlot(kss, reduction = "full.umap", label = FALSE, raster=FALSE, alpha = 0.1,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl1.cols[levels(Idents(kss))], name = "Subclass") 
dev.off()




###Cell type Proportions by region (based on code in Lake et al., 2023)
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09042024_Final_Object.Rds")
maxW.l3.df <- kss@meta.data[,c("library","v2.subclass.sp","region_level2","region_level1")]

ct.puck.count <- maxW.l3.df %>% count(library, v2.subclass.sp) # counting celltypes per puck
ct.puck.count <- ct.puck.count %>% group_by(library) %>% mutate(inPuck.frac = n / sum(n)) # fraction of cell type per puck
ct.puck.count <- ct.puck.count %>% group_by(v2.subclass.sp) %>% mutate(norm.frac = inPuck.frac / max(inPuck.frac)) # normalizing fraction per puck across puck per cell type

exp.tab <- maxW.l3.df %>% select(library,region_level2,region_level1) %>% distinct()
rownames(exp.tab) <- exp.tab$library
ct.puck.count2 <- cbind(ct.puck.count, 
                        exp.tab[ct.puck.count$library,c("region_level2","region_level1")])
region.order <- c("Cortex", "Outer Medulla", "Inner/Outer Medulla", "Inner Medulla")

regMean <- ct.puck.count2 %>% filter(`region_level2` != "Cortex / Medulla") %>% 
  group_by(v2.subclass.sp, `region_level2`) %>%
  summarise(regMeanFrac=mean(inPuck.frac))
regMean <- regMean %>% group_by(v2.subclass.sp) %>%
  mutate(normRegionMean = regMeanFrac / max(regMeanFrac))

regMean.wide <- pivot_wider(regMean[c('region_level2', 'v2.subclass.sp', 'normRegionMean')], 
                            names_from = 'region_level2', values_from = 'normRegionMean')
regMean.wide[is.na(regMean.wide)] <- 0 # filling 0 for unobserved cell types

toPlot.df <- reshape2::melt(regMean.wide, variable.name = 'region', value.name = 'frac.norm')
toPlot.df$region <- factor(toPlot.df$region, levels = region.order)

# Selecting a subset of cell types 
sel.ctypes <- order

toPlot.df <- toPlot.df[toPlot.df$v2.subclass.sp %in% sel.ctypes, ]
toPlot.df$v2.subclass.sp <- factor(toPlot.df$v2.subclass.sp, levels = sel.ctypes)

pdf("slide-seq/QC_Plots/Subclass_level3_Region_CelltypeProp.pdf", width = 4, height = 10)
ggplot(data = toPlot.df, aes(x = region, y = v2.subclass.sp, fill = frac.norm)) + 
  geom_tile() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  scale_fill_continuous(type='viridis')
dev.off()



###Dotplot of marker genes
#predicted labels
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09042024_Final_Object.Rds")
Idents(kss) <- "v2.subclass.sp"
Idents(kss) <- factor(Idents(kss), levels = order)
table(Idents(kss))
kss.sub <- subset(kss, downsample = 1000)
table(kss.sub$v2.subclass.sp)
kss.sub[["Spatial"]] <- as(object = kss.sub[["Spatial"]], Class = "Assay")
kss.sub[["Spatial"]] <- NormalizeData(kss.sub[["Spatial"]])

pt.markers <- c(
  "NPHS1",                         #POD
  "CFH",                            #PEC
  
  "LRP2",                                                #PT
  "PRODH2",                          #S1
  "SLC34A1",                       #S2                                  
  "CDH6",
  "VCAM1",                   #aPT2
  "DLGAP1",                                             #frPTS1/S2
  
  "AQP1",                  #DTL2
  #"SLC14A2",                                  #DTL3
  "AKR1B1",                             #DTL3/ATL
  "PROX1"                         #ATL
)
dt.markers <- c(
  "EGF",                  #TAL
  "CLCNKA",      #M-TAL
  "ENOX1",                        #C-TAL
  "NOS1",                       #MD
  
  "ITGB6",                                 #aTAL
  
  "ADAMTS1",                              #aTAL2
  
  "PROM1",  #frTAL
  
  "SLC12A3",                       #DCT
  "SLC8A1",                                                #DCT2 / CNT
  "GATA3",                                  #PC
  "SLC14A2",                                      #IMCD
  
  "ATP6V0D2",                          #IC
  "SLC26A7",                                      #IC-A                                   
  
  "SLC26A4"               #IC-B
)
ec.markers <- c("EHD3",                             #EC-GC      
                "SULF1",                           #EC-AA
                "PALMD",               #EC-DVR  
                
                "PLVAP",                              #PTC/AVR
                "EDIL3",                               #EC-AVR
                "IFIT1",                             #iaEC-AVR
                
                "VWF",                     #EC-V
                "DOK6",                                                    #EC-PCV
                "NAV3",                          #EC-PTC
                "SLC6A6",                          #EC-EA
                
                'ICAM1',                   #infEC-PTC
                "GBP1",                           #iaEC-PTC
                "MMRN1"                             #EC-LYM
)
str.markers <- c(
  "PDGFRA",                                     #Pan FIB
  "TNC",                         #Pan Medullary FIB
  "ADAMTSL1",                            #C/M-FIB
  "ACTG2",              #IM-pvMYOF
  
  "CCN1",                   #C-FIB (interstitial fib)
  "LUM",                  #C-FIB-OSMRlo
  "OSMR",                     #C-FIB-OSMRhi
  "SULF1","FAP",               #C-MYOF
  "RSPO3",                #pvFIB-RSPO3+
  "C3",                        #pvFIB-PI16+
  "MGAT4C",                                        #pvFIB
  "ITGA8",                                                  #pvMYOF & VSMC
  "MYH11",               #pvMYOF
  "PDGFRB",                                       #Pan VSM markers
  "POSTN",                   #MC
  "REN",                 #REN
  "MCAM",          #VSMC
  "NOTCH3",                                                 #VSMC &VSMC/P
  "RGS5"                             #VSMC/P
  
)
imm.markers <- c(                                            #Broad Immune
                 "BANK1",                     #B Cells
                 "JCHAIN",               #PL Cells
                 "CD96",                    #T
                 "IL7R",                                              #T/MAIT/ILC3
                 "KLRB1",                            #MAIT
                 "TOX2",                         #ILC3
                 "IKZF2",            #T-REG
                 "CD8A",                                              #CD8+
                 "CCL5",                        #CD8+ & NK 
                 "GNLY",               #CD8+ TEM/TEMRA & NK
                 "HBB",                                #Erythrocyte
                 "CPA3",                   #MAST
                 "CD163","CD14",                      #MAC
                 "LYVE1",           #resMAC-LYVE1+
                 "C1QB",                           #resMAC-HLAIIhi
                 "NAMPT",      #moMAC-HBEGF+
                 "CXCL10",    #moMAC-CXCL10+
                 "GPNMB",      #moFAM
                 "C3",             #moMAC-C3+
                 "FCER1A",                          #cDC2
                 "COTL1",                          #ncMON
                 "FCN1",                                             #MON/ncMON
                 "SLCO5A1",           #mDC
                 "WDFY4",                   #cDC1
                 "BCL11A",         #pDC
                 "FCGR3B"                #N
)

pdf("slide-seq/Plots/Slide-seq_Subclass_curated_markers_DotPlot.pdf", width = 20, height = 10)
DotPlot(kss.sub, features = unique(c(pt.markers, dt.markers,ec.markers,str.markers,imm.markers)), dot.scale = 8, scale.max = 50) + RotatedAxis()
dev.off()



###Molecule plots
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09042024_Final_Object.Rds")
Idents(kss) <- "v2.subclass.sp"
Idents(kss) <- factor(Idents(kss), levels = order)

pdf("slide-seq/Image_Plots/Puck_200903_01_feature_NPHS1_Plot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss, images = "Puck_200903_01", 
                   features = c("NPHS1"), ncol = 2, 
                   alpha = c(0.1, 1), stroke = NA, min.cutoff = "q5", max.cutoff = "q95")  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Image_Plots/Puck_200903_01_feature_SLC5A12_Plot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss, images = "Puck_200903_01", 
                   features = c("SLC5A12"), ncol = 2, 
                   alpha = c(0.1, 1), stroke = NA, min.cutoff = "q5", max.cutoff = "q95")  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Image_Plots/Puck_200903_01_feature_SLC12A1_Plot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss, images = "Puck_200903_01", 
                   features = c("SLC12A1"), ncol = 2, 
                   alpha = c(0.1, 1), stroke = NA, min.cutoff = "q5", max.cutoff = "q95")  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Image_Plots/Puck_200903_01_feature_AQP2_Plot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss, images = "Puck_200903_01", 
                   features = c("AQP2","SLC12A1"), ncol = 2, 
                   alpha = c(0.1, 1), stroke = NA, min.cutoff = "q5", max.cutoff = "q95")  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

v2.str.cols <- setNames(c("#A10300","#004CFF","#005300","#720055"),
                        c("POD","PT-S1/2","C-TAL","PC"))


kss.sub <- subset(kss, idents = c("POD","PT-S1/2","C-TAL","PC"))
pdf("slide-seq/Image_Plots/Puck_200903_01_POD-PT-S12-C-TAL-PC_Plot.pdf", width = 6, height = 7)
SpatialDimPlot(kss.sub, images = "Puck_200903_01", cols = v2.str.cols[levels(Idents(kss))],
               alpha = c(1), stroke = NA) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



###Structure Plots
# build-in color palette
Glasbey = glasbey.colors(32)
swatch(Glasbey)
levels(Idents(kss))

v2.str.cols <- setNames(c("#004CFF","#783FC1","#A10300","#00FF00","#DD00FF","#005300","#FFD300",
                          "#720055","#858567","#DC5E93","#00FFBE","#FE8F42"),
                        c("proximal tubules","vessels","renal corpuscle","interstitium - stroma",     
                          "interstitium - immune","Distal tubules" ,"inner intermediate tubules",
                          "outer Collecting tubules","inner Collecting tubules",
                          "interstitium","intermediate tubules","outer intermediate tubules"))


Idents(kss) <- "v2.structure"
table(kss$library)
SpatialDimPlot(kss, images = "Puck_200903_06", cols = v2.str.cols[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) 


fovs <- c("Puck_200903_06","Puck_210113_39","Puck_210113_23",
          "Puck_220122_36","Puck_220122_37","Puck_220122_38","Puck_220122_39",
          "Puck_220122_33","Puck_220122_34","Puck_220122_28")
lapply(fovs, function(x) {
  print(paste("Running for FOV:", x))
  cl.plot <- SpatialDimPlot(kss, images = x, cols = v2.str.cols[levels(Idents(kss))],
                            alpha = c(1), stroke = NA) 
  pdf(file=paste0("slide-seq/Image_Plots/fov_plots/",x,"structure_segmentation_Plot.pdf"),
      width=7,height=6)
  print(cl.plot)
  dev.off()
  
})




###Cell type sets
##Renal corpuscles
Glasbey = glasbey.colors(25)
swatch(Glasbey)
v2.rc.cols <- setNames(c("#0000FF","#FF0000","#00FF00","#DD00FF",
                         "#005300",
                         "#FFD300","#009FFF","#9A4D42","#00FFBE"),
                       c("POD","PEC","MD","EC-GC",
                         "EC-AA",
                         "EC-EA","MC","REN","VSMC"))

Idents(kss) <- "v2.subclass.sp"
na.cols <- setNames(rep("gray10", length(levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.rc.cols)])),
                    levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.rc.cols)])

SpatialDimPlot(kss, images = "Puck_200903_06", cols = c(v2.rc.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()

cl.plot1 <- SpatialDimPlot(kss, images = "Puck_200903_06", cols = c(v2.rc.cols,na.cols)[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + NoLegend() + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
cl.plot2 <- SpatialDimPlot(kss, images = "Puck_200903_06", cols = v2.rc.cols[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pdf(file="slide-seq/Image_Plots/fov_plots/Puck_200903_06_Renal_Corpuscle_segmentation_Plot.pdf",
      width=7,height=6)
cl.plot1
cl.plot2
dev.off()

kss.sub <- subset(kss, idents = c("POD","PEC","MD","EC-GC",
                                  "EC-AA",
                                  "EC-EA","MC","REN","VSMC"))

pdf(file="slide-seq/Image_Plots/fov_plots/Puck_200903_06_Renal_Corpuscle_segmentation_Plot_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_200903_06", cols = v2.rc.cols[levels(Idents(kss))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend() 
dev.off()



##Vasculature
Glasbey = glasbey.colors(35)
swatch(Glasbey)
v2.vsc.cols <- setNames(c("#009FFF","#1F9698","#FE8F42","#FE8F42",
                          "#783FC1","#720055","#004CFF",
                          "#FFFFFF","#FFFFFF","#FFFFFF",
                          "#F1085C","#DD00FF",
                          "#C8FF00","#FFD300","#02AD24","#00FF00"),
                        c("EC-AA","EC-DVR","EC-AVR","infEC-AVR", 
                          "EC-V","EC-LYM","EC-EA",
                          "EC-PTC","angEC-PTC","infEC-PTC",
                          "VSMC","VSMC/P",
                          "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))





na.cols <- setNames(rep("gray10", length(levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.vsc.cols)])),
                    levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.vsc.cols)])
SpatialDimPlot(kss, images = "Puck_200903_06", cols = c(v2.vsc.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()
SpatialDimPlot(kss, images = "Puck_210113_23", cols = c(v2.vsc.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()
SpatialDimPlot(kss, images = "Puck_210113_39", cols = c(v2.vsc.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()


fovs <- c("Puck_200903_06","Puck_210113_23","Puck_210113_39",
          "Puck_220122_36","Puck_220122_37","Puck_220122_38","Puck_220122_39",
          "Puck_220122_33","Puck_220122_34","Puck_220122_28")
lapply(fovs, function(x) {
  print(paste("Running for FOV:", x))
  ln <- length(rownames(kss@meta.data[kss@meta.data$v2.subclass.sp %in% c("EC-AA","EC-DVR","EC-AVR","infEC-AVR", 
                                                                          "EC-V","EC-LYM","EC-EA",
                                                                          "EC-PTC","angEC-PTC","infEC-PTC",
                                                                          "VSMC","VSMC/P",
                                                                          "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"),]))
  if(ln > 10) {
    cl.plot1 <- SpatialDimPlot(kss, images = x, cols = c(v2.vsc.cols,na.cols)[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + NoLegend() + DarkTheme() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    cl.plot2 <- SpatialDimPlot(kss, images = x, cols = v2.vsc.cols[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + DarkTheme() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else {
    print(paste("Not enough cells for FOV:", x))
  }
  pdf(file=paste0("slide-seq/Image_Plots/fov_plots/",x,"Vasculature_segmentation_Plot.pdf"),
      width=7,height=6)
  print(cl.plot1)
  print(cl.plot2)
  dev.off()
})



kss.sub <- subset(kss, idents = c("EC-AA","EC-EA","EC-DVR",
                                  "EC-V",
                                  
                                  "VSMC","VSMC/P",
                                  "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))

pdf(file="slide-seq/Image_Plots/fov_plots/Puck_210113_39_Vasculature_segmentation_Plot_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_210113_39", cols = v2.vsc.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()
pdf(file="slide-seq/Image_Plots/fov_plots/Puck_200903_06_Vasculature_segmentation_Plot_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_200903_06", cols = v2.vsc.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()
pdf(file="slide-seq/Image_Plots/fov_plots/Puck_220122_38_Vasculature_segmentation_Plot_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_220122_38", cols = v2.vsc.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()




##Fibrotic areas
P33 = createPalette(33,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P33)
Glasbey = glasbey.colors(30)
swatch(Glasbey)


kss$v2.plot <- kss$v2.structure
Idents(kss) <- "v2.structure"
kss@meta.data[WhichCells(kss, idents = c("Distal tubules" ,
                         "outer Collecting tubules","outer intermediate tubules")),]$v2.plot <- "other epi"
Idents(kss) <- "v2.subclass.sp"
kss@meta.data[WhichCells(kss, idents = c("PT-S1/2","PT-S3")),]$v2.plot <- "PT"
kss@meta.data[WhichCells(kss, idents = c("frPT","frTAL")),]$v2.plot <- "frEpi"
kss@meta.data[WhichCells(kss, idents = c("M-TAL","C-TAL","MD")),]$v2.plot <- "TAL"
kss@meta.data[WhichCells(kss, idents = c("aTAL1","aTAL2")),]$v2.plot <- "aTAL"
kss@meta.data[WhichCells(kss, idents = c("aPT")),]$v2.plot <- "aPT"
kss@meta.data[WhichCells(kss, idents = c("B","PL","T","MAIT",
                                         "ILC3","T-REG","CD8+ TEM/TRM",
                                         "CD8+ TEM/TEMRA","NK")),]$v2.plot <- "Lymphoid"
kss@meta.data[WhichCells(kss, idents = c("MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
                                         "moMAC-INF","moFAM","moMAC-C3+","cDC2",
                                         "ncMON","MON","mDC","cDC1","pDC","N")),]$v2.plot <- "Myeloid"


Idents(kss) <- "v2.plot"
v2.fn.cols <- setNames(c("#00479E","#1F9698",
                         "#005300","#FE8F42","#DC5E93",
                         "#783FC1","#A10300","#886C00",
                         "#C8FF00","#14F9FF",
                         "gray10"
),
c("PT","aPT",
  "TAL","aTAL","frEpi",
  "vessels","renal corpuscle","interstitium - stroma",
  "Lymphoid","Myeloid",
  "other epi"))


na.cols <- setNames(rep("gray10", length(levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.vsc.cols)])),
                    levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.vsc.cols)])
SpatialDimPlot(kss, images = "Puck_200903_06", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()

SpatialDimPlot(kss, images = "Puck_220122_38", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()
SpatialDimPlot(kss, images = "Puck_220122_39", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()
SpatialDimPlot(kss, images = "Puck_220122_36", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()


fovs <- c("Puck_200903_06","Puck_210113_23","Puck_210113_39",
          "Puck_220122_36","Puck_220122_37","Puck_220122_38","Puck_220122_39",
          "Puck_220122_33","Puck_220122_34","Puck_220122_28")
lapply(fovs, function(x) {
  print(paste("Running for FOV:", x))
  ln <- length(rownames(kss@meta.data[kss@meta.data$v2.plot %in% c("PT","aPT",
                                                                          "TAL","aTAL","frEpi",
                                                                          "vessels","renal corpuscle","interstitium - stroma",
                                                                          "Lymphoid","Myeloid",
                                                                          "other epi"),]))
  if(ln > 10) {
    cl.plot1 <- SpatialDimPlot(kss, images = x, cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + NoLegend() + DarkTheme() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    cl.plot2 <- SpatialDimPlot(kss, images = x, cols = v2.fn.cols[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + DarkTheme() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else {
    print(paste("Not enough cells for FOV:", x))
  }
  pdf(file=paste0("slide-seq/Image_Plots/fov_plots/",x,"Fibrosis_segmentation_Plot.pdf"),
      width=7,height=6)
  print(cl.plot1)
  print(cl.plot2)
  dev.off()
})


kss.sub <- subset(kss, idents = c("aPT",
                                  "aTAL","frEpi",
                                  "renal corpuscle","interstitium - stroma",
                                  "Lymphoid","Myeloid"))

pdf(file="slide-seq/Image_Plots/fov_plots/Puck_220122_37_Fibrosis_segmentation_Plot_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_220122_37", cols = v2.fn.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()




##Fibroblast subtypes

kss$v2.plot <- kss$v2.structure
Idents(kss) <- "v2.structure"
kss@meta.data[WhichCells(kss, idents = c("proximal tubules","Distal tubules" ,
                                         "outer Collecting tubules","outer intermediate tubules")),]$v2.plot <- "epi"
Idents(kss) <- "v2.subclass.sp"
kss@meta.data[WhichCells(kss, idents = c("frPT","frTAL")),]$v2.plot <- "frEpi"
kss@meta.data[WhichCells(kss, idents = c("B","PL","T","MAIT",
                                         "ILC3","T-REG","CD8+ TEM/TRM",
                                         "CD8+ TEM/TEMRA","NK")),]$v2.plot <- "Immune"
kss@meta.data[WhichCells(kss, idents = c("MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
                                         "moMAC-INF","moFAM","moMAC-C3+","cDC2",
                                         "ncMON","MON","mDC","cDC1","pDC","N")),]$v2.plot <- "Immune"
kss@meta.data[WhichCells(kss, idents = c("C-FIB","C-FIB-Path","C-FIB-OSMRlo")),]$v2.plot <- "C-FIB"
kss@meta.data[WhichCells(kss, idents = c("C-FIB-OSMRhi")),]$v2.plot <- "C-FIB-OSMRhi"
kss@meta.data[WhichCells(kss, idents = c("C-MYOF")),]$v2.plot <- "C-MYOF"
kss@meta.data[WhichCells(kss, idents = c("pvFIB-RSPO3+")),]$v2.plot <- "pvFIB-RSPO3+"
kss@meta.data[WhichCells(kss, idents = c("pvFIB-PI16+")),]$v2.plot <- "pvFIB-PI16+"
kss@meta.data[WhichCells(kss, idents = c("pvFIB")),]$v2.plot <- "pvFIB"
kss@meta.data[WhichCells(kss, idents = c("pvMYOF")),]$v2.plot <- "pvMYOF"


Idents(kss) <- "v2.plot"
v2.fn.cols <- setNames(c("gray10","#783FC1","#A10300","white","gray10",
                         "#886C00","#009FFF","#FF0000",
                         "#C8FF00","#FFD300","#02AD24","#00FF00"),
                       c("frEpi","vessels","renal corpuscle","Immune","other epi",
                         "C-FIB","C-FIB-OSMRhi","C-MYOF",
                         "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))


na.cols <- setNames(rep("gray10", length(levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.fn.cols)])),
                    levels(Idents(kss))[!levels(Idents(kss)) %in% names(v2.fn.cols)])
SpatialDimPlot(kss, images = "Puck_200903_06", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()

SpatialDimPlot(kss, images = "Puck_220122_38", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()
SpatialDimPlot(kss, images = "Puck_220122_39", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()
SpatialDimPlot(kss, images = "Puck_220122_36", cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
               alpha = c(0.8), stroke = NA) + DarkTheme()


fovs <- c("Puck_200903_06","Puck_210113_23","Puck_210113_39",
          "Puck_220122_36","Puck_220122_37","Puck_220122_38","Puck_220122_39",
          "Puck_220122_33","Puck_220122_34","Puck_220122_28")
lapply(fovs, function(x) {
  print(paste("Running for FOV:", x))
  ln <- length(rownames(kss@meta.data[kss@meta.data$v2.plot %in% c("frEpi","vessels","renal corpuscle","Immune","other epi",
                                                                   "C-FIB","C-FIB-OSMRhi","C-MYOF",
                                                                   "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"),]))
  if(ln > 10) {
    cl.plot1 <- SpatialDimPlot(kss, images = x, cols = c(v2.fn.cols,na.cols)[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + NoLegend() + DarkTheme() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    cl.plot2 <- SpatialDimPlot(kss, images = x, cols = v2.fn.cols[levels(Idents(kss))],
                               alpha = c(0.8), stroke = NA) + DarkTheme() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else {
    print(paste("Not enough cells for FOV:", x))
  }
  pdf(file=paste0("slide-seq/Image_Plots/fov_plots/",x,"Fibrosis_segmentation_Plot_FIB-subtypes.pdf"),
      width=7,height=6)
  print(cl.plot1)
  print(cl.plot2)
  dev.off()
})


kss.sub <- subset(kss, idents = c("C-FIB-OSMRhi","C-MYOF",
                                  "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))

pdf(file="slide-seq/Image_Plots/fov_plots/Puck_220122_37_Fibrosis_segmentation_Plot_FIB-subtypes_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_220122_37", cols = v2.fn.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()
pdf(file="slide-seq/Image_Plots/fov_plots/Puck_220122_36_Fibrosis_segmentation_Plot_FIB-subtypes_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_220122_36", cols = v2.fn.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()
pdf(file="slide-seq/Image_Plots/fov_plots/Puck_220122_38_Fibrosis_segmentation_Plot_FIB-subtypes_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_220122_38", cols = v2.fn.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()
pdf(file="slide-seq/Image_Plots/fov_plots/Puck_220122_39_Fibrosis_segmentation_Plot_FIB-subtypes_subset.pdf",
    width=7,height=6)
SpatialDimPlot(kss.sub, images = "Puck_220122_39", cols = v2.fn.cols[levels(Idents(kss.sub))]) + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + NoLegend()
dev.off()
