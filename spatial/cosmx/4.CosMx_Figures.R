library(Seurat)
library(SeuratData)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(BPCells)
library(corrplot)
library(tidyr)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")
load("color_factors_v2-clusters.robj")
source("misc/utils.R")
load("AtlasV2_VariableFeatures.rda")

nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated.RDS")
order <- c("POD","PEC","PT-S1/S2","PT-S3","aPT1","aPT2","frPT","DTL/ATL","M-TAL",
           "C-TAL","MD","aTAL","frTAL","dTAL","DCT","CNT","aCNT","C-PC","M-PC",
           "IC-A","IC-B","EC","EC-GC","EC-AEA/DVR","EC-AVR/V","M-EC-PTC",
           "EC-PTC","EC-LYM","FIB","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvMYOF","MC","REN","VSMC","VSMC/P",
           "B","PL","T","CD8+ T / NK","ERY","MAST","resMAC","moMAC-INF","moMAC/DC",
           "ncMON/N","mDC","cyc")
Idents(nano.obj) <- "v2.subclass.l3"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = order)
table(Idents(nano.obj))
write.table(nano.obj@meta.data, file="nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated_Metadata.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Final Post-V2 annotation/alignment Stats
meta <- nano.obj@meta.data
stats <- cbind(
  data.frame(meta %>%
               group_by(fov) %>%
               summarise_at(vars(nCount_Nanostring,nFeature_Nanostring), list(mean))),
  data.frame(meta %>%
               group_by(fov) %>%
               tally())
)
rownames(stats) <- stats$fov
stats <- stats[,-c(1,4)]
stats  

#Determine number of cell types
df <- t(table(nano.obj$fov, nano.obj$v2.subclass.l3))
cols <- vector()
exps <- colnames(df)

for(i in exps){
  cols[i] <- length(df[which(df[,i] > 5),i]) 
}
stats$n_clusters <- cols[rownames(stats)]

write.table(stats, file="nanostring/QC_Plots/CosMx_post-clustering_Stats_09112024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)



###UMAP Plots
Idents(nano.obj) <- "v2.subclass.l3"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = order)

meta <- nano.obj@meta.data
meta <- meta[!duplicated(meta$v2.subclass.l3),]
varibow <- function(n_colors) {
  sats <- rep_len(c(0.25,0.4,0.55,0.7),length.out = n_colors)
  vals <- rep_len(c(0.9,0.7,0.5),length.out = n_colors)
  grDevices::rainbow(n_colors, s = sats, v = vals)
}
meta <- annotate_cat(meta, col = "v2.subclass.l3", sort_label = F, colorset = "varibow")
v2.scl3.cols <- meta$v2.subclass.l3_color; names(v2.scl3.cols) <- meta$v2.subclass.l3_label

pdf(file='nanostring/Plots/Nanostring_v2_Subclass.l3_umap.pdf',width=8,height=7)
DimPlot(nano.obj, reduction = "umap", label = TRUE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(nano.obj))], name = "Subclass") 
dev.off()
pdf(file='nanostring/Plots/Nanostring_v2_Subclass.l3_umap_unlabeled.pdf',width=8,height=7)
DimPlot(nano.obj, reduction = "umap", label = FALSE, raster=FALSE, alpha = 0.3,
        label.size = 4, repel = TRUE) + ggtitle("Subclass.l3"
        ) + NoLegend() + scale_color_manual(values = v2.scl3.cols[levels(Idents(nano.obj))], name = "Subclass") 
dev.off()


pdf("nanostring/Plots/Nanostring_Protein_Stain_FeaturePlots.pdf", width = 9, height = 8)
FeaturePlot(nano.obj, features = c("Mean.CD3","Mean.PanCK","Mean.CD45","Mean.CD298"), alpha = 0.1, min.cutoff = 'q10', max.cutoff = 'q90') 
dev.off()


###Correlation with sn reference atlas
#Prepare Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
var.genes <- VariableFeatures(ref)
Idents(object = ref) <- "v2.subclass.sp"
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

od.genes <- c(pt.od.genes,dt.od.genes,ec.od.genes,str.od.genes,imm.od.genes)
sn.od.genes <- od.genes[od.genes %in% rownames(nano.obj)]
#sn.od.genes <- rownames(nano.obj)[rownames(nano.obj) %in% rownames(ref)]

ref[["RNA"]] <- as(object = ref[["RNA"]], Class = "Assay")
ref <- NormalizeData(ref)
ref <- ScaleData(ref, features = sn.od.genes, assay = "RNA")

#prepare query
nano.obj.sub <- subset(nano.obj, downsample = 5000)
table(nano.obj.sub$v2.subclass.l3)
nano.obj.sub[["Nanostring"]] <- NormalizeData(nano.obj.sub[["Nanostring"]])
nano.obj.sub[["Nanostring"]] <- ScaleData(nano.obj.sub[["Nanostring"]], features = sn.od.genes)

ave.ref <- AverageExpression(ref, features = sn.od.genes, assays = "RNA", slot = "scale.data")
ave.nano.obj <- AverageExpression(nano.obj.sub, features = sn.od.genes, assay = "Nanostring", slot = "scale.data")

ave.cor<-cor(as.data.frame(ave.nano.obj),as.data.frame(ave.ref))
colnames(ave.cor) <- gsub("RNA.", "", colnames(ave.cor))
rownames(ave.cor) <- gsub("Nanostring.", "", rownames(ave.cor))
write.table(ave.cor, file="nanostring/QC_Plots/10X_RNA_v2-subclasses_versus_nanostring_pred-subclasses_corr_all_subclass.l3_08262024.txt", sep = "\t", row.names=TRUE, col.names=TRUE)

order1 <- c("POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","aDTL","DTL1","DTL3","ATL",
            "M-TAL","C-TAL","MD","aTAL1","aTAL2","frTAL","DCT","aDCT","CNT","aCNT","PC",
            "IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR","EC-PTC","EC-AVR",
            "infEC-AVR","EC-V","EC-PCV","angEC-PTC","EC-EA","infEC-PTC","EC-LYM","M-FIB","C/M-FIB",
            "IM-pvMYOF","C-FIB","C-FIB-Path","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+",
            "pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
            "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+","resMAC-HLAIIhi",
            "moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON","mDC","cDC1","pDC","N")

order1 <- gsub("-",".", order1)
order1 <- gsub("/",".", order1)
order <- gsub("-",".", order)
order <- gsub("/",".", order)

col.order <- order1[order1 %in% colnames(ave.cor)]
row.order <- order[order %in% rownames(ave.cor)]

pdf(file='nanostring/QC_Plots/10X_RNA_v2-subclasses_versus_nanostring_pred-subclasses_corr_Plot_all_subclass.l3_08262024.pdf',width=14,height=14)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))
corrplot(ave.cor[row.order,col.order], col = col2(500),  method="color",order = "original",cl.lim=c(-1,1))
dev.off()




###Cell type Proportions by region (based on code in Lake et al., 2023)
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated.RDS")
order <- c("POD","PEC","PT-S1/S2","PT-S3","aPT1","aPT2","frPT","DTL/ATL","M-TAL",
           "C-TAL","MD","aTAL","frTAL","dTAL","DCT","CNT","aCNT","C-PC","M-PC",
           "IC-A","IC-B","EC","EC-GC","EC-AEA/DVR","EC-AVR/V","M-EC-PTC",
           "EC-PTC","EC-LYM","FIB","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvMYOF","MC","REN","VSMC","VSMC/P",
           "B","PL","T","CD8+ T / NK","ERY","MAST","resMAC","moMAC-INF","moMAC/DC",
           "ncMON/N","mDC","cyc")
Idents(nano.obj) <- "v2.subclass.l3"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = order)
table(Idents(nano.obj))
maxW.l3.df <- nano.obj@meta.data[,c("fov","v2.subclass.l3","region_level1")]
maxW.l3.df <- maxW.l3.df[!maxW.l3.df$fov %in% "NA_NA",]
ct.puck.count <- maxW.l3.df %>% count(fov, v2.subclass.l3) # counting celltypes per fov
ct.puck.count <- ct.puck.count %>% group_by(fov) %>% mutate(inPuck.frac = n / sum(n)) # fraction of cell type per fov
ct.puck.count <- ct.puck.count %>% group_by(v2.subclass.l3) %>% mutate(norm.frac = inPuck.frac / max(inPuck.frac)) # normalizing fraction per fov across fov per cell type

exp.tab <- read.delim("nanostring/Kidney_AtlasV2_Nanostring_Experiment_Metadata.txt")
rownames(exp.tab) <- exp.tab$fov
ct.puck.count2 <- cbind(ct.puck.count, 
                        exp.tab[ct.puck.count$fov,c("region_level2","region_level1")])
region.order <- c("Cortex","CMJ","Medulla")

regMean <- ct.puck.count2 %>% 
  group_by(v2.subclass.l3, `region_level1`) %>%
  summarise(regMeanFrac=mean(inPuck.frac))
regMean <- regMean %>% group_by(v2.subclass.l3) %>%
  mutate(normRegionMean = regMeanFrac / max(regMeanFrac))

regMean.wide <- pivot_wider(regMean[c('region_level1', 'v2.subclass.l3', 'normRegionMean')], 
                            names_from = 'region_level1', values_from = 'normRegionMean')
regMean.wide[is.na(regMean.wide)] <- 0 # filling 0 for unobserved cell types

toPlot.df <- reshape2::melt(regMean.wide, variable.name = 'region', value.name = 'frac.norm')
toPlot.df$region <- factor(toPlot.df$region, levels = region.order)

# Selecting a subset of cell types 
sel.ctypes <- c("POD","PEC","PT-S1/S2","PT-S3","aPT1","aPT2","frPT","DTL/ATL","M-TAL",
                "C-TAL","MD","aTAL","frTAL","dTAL","DCT","CNT","aCNT","C-PC","M-PC",
                "IC-A","IC-B","EC","EC-GC","EC-AEA/DVR","EC-AVR/V","M-EC-PTC",
                "EC-PTC","EC-LYM","FIB","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
                "pvFIB-RSPO3+","pvFIB-PI16+","pvMYOF","MC","REN","VSMC","VSMC/P",
                "B","PL","T","CD8+ T / NK","ERY","MAST","resMAC","moMAC-INF","moMAC/DC",
                "ncMON/N","mDC","cyc")

toPlot.df <- toPlot.df[toPlot.df$v2.subclass.l3 %in% sel.ctypes, ]
toPlot.df$v2.subclass.l3 <- factor(toPlot.df$v2.subclass.l3, levels = sel.ctypes)

pdf("nanostring/QC_Plots/Subclass_level3_Region_CelltypeProp_08262024.pdf", width = 4, height = 10)
ggplot(data = toPlot.df, aes(x = region, y = v2.subclass.l3, fill = frac.norm)) + 
  geom_tile() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + 
  scale_fill_continuous(type='viridis')
dev.off()



###Spatial image plots
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated.RDS")

library(Polychrome)

# build-in color palette
Glasbey = glasbey.colors(32)
swatch(Glasbey)
levels(Idents(nano.obj))

v2.str.cols <- setNames(c("#004CFF","#783FC1","#A10300","#00FF00","#DD00FF","#005300",
                          "#720055","#FFD300"),
                        c("proximal tubules","vessels","renal corpuscle","interstitium - stroma",     
                          "interstitium - immune","distal tubules" ,
                          "collecting tubules","intermediate tubules"))


Idents(nano.obj) <- "v2.structure"

pdf("nanostring/Image_Plots/R5446.S2_Structure-level_Plot_08262024.pdf", width = 20, height = 20)
ImageDimPlot(nano.obj, fov = "R5446.S2", axes = FALSE, cols = v2.str.cols[levels(Idents(nano.obj))], border.size = NA, size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
ImageDimPlot(nano.obj, fov = "R5446.S2", cells = WhichCells(nano.obj, 
                                                            idents = c("renal corpuscle", "proximal tubules","interstitium - stroma", "interstitium - immune")),
             cols = v2.str.cols[levels(Idents(nano.obj))], size = 1, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("nanostring/Image_Plots/R5446.S2_Structure-level_Plot_08262024_B.pdf", width = 20, height = 20)
ImageDimPlot(nano.obj, fov = "R5446.S2", cols = v2.str.cols[levels(Idents(nano.obj))], coord.fixed = FALSE,
             size = 1,border.size = 0.1, border.color = "black",dark.background = F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


pdf("nanostring/Image_Plots/R5341.S1_Structure-level_Plot_08262024.pdf", width = 20, height = 20)
ImageDimPlot(nano.obj, fov = "R5341.S1", axes = TRUE, cols = v2.str.cols[levels(Idents(nano.obj))], border.size = NA, size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ImageDimPlot(nano.obj, fov = "R5341.S1", cells = WhichCells(nano.obj, 
                                                            idents = c("renal corpuscle", "proximal tubules","interstitium - stroma", "interstitium - immune")),
             cols = v2.str.cols[levels(Idents(nano.obj))], size = 1, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("nanostring/Image_Plots/R5341.S1_Structure-level_Plot_08262024_B.pdf", width = 20, height = 20)
ImageDimPlot(nano.obj, fov = "R5341.S1", cols = v2.str.cols[levels(Idents(nano.obj))], coord.fixed = FALSE,
             size = 1,border.size = 0.1, border.color = "black",dark.background = F)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


###Dotplot of marker genes
nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated.RDS")
order <- c("POD","PEC","PT-S1/S2","PT-S3","aPT1","aPT2","frPT","DTL/ATL","M-TAL",
           "C-TAL","MD","aTAL","frTAL","dTAL","DCT","CNT","aCNT","C-PC","M-PC",
           "IC-A","IC-B","EC","EC-GC","EC-AEA/DVR","EC-AVR/V","M-EC-PTC",
           "EC-PTC","EC-LYM","FIB","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvMYOF","MC","REN","VSMC","VSMC/P",
           "B","PL","T","CD8+ T / NK","ERY","MAST","resMAC","moMAC-INF","moMAC/DC",
           "ncMON/N","mDC","cyc")
Idents(nano.obj) <- "v2.subclass.l3"
Idents(nano.obj) <- factor(Idents(nano.obj), levels = order)
table(Idents(nano.obj))
nano.obj <- NormalizeData(nano.obj)
nano.obj <- ScaleData(nano.obj, features = rownames(nano.obj))


##Single plot of all markers
markers <- c(
  "NPHS1",                         #POD
  "VCAM1",                            #PEC
  "LRP2",                                                #PT
  "SLC5A12",                                  #S1/S2
  "SLC7A13",                        #PT-S3
  "HAVCR1",
  "IL32","MMP7","CCL2",                   #aPT2
  "PROM1",                                             #frPT
  "TACSTD2","JAG1",                       #TL

  "SLC12A1","EGF",                  #TAL
  "BMP6","NOS1","ITGA3","CD44",                       #aTAL1                                      
  
  "ITGB6",
  "SLC12A3",                       #DCT
  "SLC8A1",                                                #DCT2 / CNT
  "CALB1",                      #CNT
  "CTSD",                                                  #dCNT
  "AQP2",                                 #PC
  "SLC26A7",                                     #IC-A                                   
  "SLC4A9", "SLC26A4",                 #IC-B
  
  "PECAM1","FLT1",                                   #Broad EC
  "EMCN","KDR",                                    #EC-GC
  "AQP1",
  "PALMD",               #EC-DVR  
  "CEACAM1","PLVAP",                             #PTC/AVR
  "FABP5",                         #M-EC-PTC
  "GPM6A","NR2F2",                                           #PTC/AVR
  "VWF","TEK",            #EC-AVR, EC-V
  "CD36","PROX1",                             #EC-LYM
  
  "DCN","PDGFRA",                                     #Pan FIB
  "MEG3",                           #Pan cortical FIB
  "SELENOP","CXCL12",                  #C-FIB-OSMRlo
  "OSMR","IL1R1","CCL2","CCL19",                      #C-FIB-OSMRhi
  "INHBA","COL1A2","COL1A1",               #C-MYOF
  'IGF1',"RSPO3",                #pvFIB-RSPO3+
  "MYH11","ACTA2",                #pvMYOF
  "PDGFRB",                                       #Pan VSM markers
  "GATA3","IL1RL1",                   #MC
  "REN",                 #REN
  "NOTCH3","RGS5",                                 #VSMC/P
  
  "MS4A1",                      #B Cells
  "JCHAIN",               #PL Cells
  "CD3D", "IL7R",                             #T
  "CD8A","CCL5",                       #CD8+ & NK 
  "GNLY",               #CD8+ TEM/TEMRA & NK
  "HBB","HBA1","HBA2",               #ERY
  "CPA3",                   #MAST
  "CD163","CD14",                      #MAC
  "C1QB",                           #resMAC-HLAIIhi
  "CXCL10","CCL2","CCL3","IL1B",    #moMAC-CXCL10+
  "CLEC10A",                         #cDC2
  "COTL1","FCGR3A",                          #ncMON
  "LYZ",                                #MON
  "S100A9","S100A8",                #N
  "LAMP3","CCR7",           #mDC
  "MKI67"                 #cycling
)

pdf("nanostring/Plots/Kidney_v2_Nanostring_Subclassl3_curated_markers_DotPlot.pdf", width = 20, height = 10)
DotPlot(nano.obj, features = unique(markers), dot.scale = 8) + RotatedAxis()
dev.off()



###Molecule plots
pdf("nanostring/Plots/R5341.S1_molecule_NPHS1-SLC5A12-SLC12A1-AQP2_Plot.pdf", width = 20, height = 20)
ImageDimPlot(nano.obj, fov = "R5341.S1", group.by = NA, alpha = 0.4,size = 0.8,
             molecules = c("NPHS1", "SLC5A12","SLC12A1","AQP2"), nmols = 20000, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


new.nano.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/nanostring/Kidney_combined_seurat_obj_08052024_Combined_Annotated_FOVs.RDS")

Idents(new.nano.obj) <- "v2.structure"
DefaultBoundary(new.nano.obj[["R5341.S1.6"]]) <- "segmentation"
DefaultBoundary(new.nano.obj[["R5341.S1.10"]]) <- "segmentation"

pdf("nanostring/Plots/R5341.S1.6_structure_Plot_10112023.pdf", width = 10, height = 8)
ImageDimPlot(new.nano.obj, fov = "R5341.S1.6", cols = v2.str.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("nanostring/Plots/R5341.S1.10_structure_Plot_10112023.pdf", width = 10, height = 8)
ImageDimPlot(new.nano.obj, fov = "R5341.S1.10", cols = v2.str.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.4)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("nanostring/Plots/R5341.S1.6_molecule_NPHS1-SLC5A12-SLC12A1-AQP2_Plot.pdf", width = 10, height = 8)
ImageDimPlot(new.nano.obj, fov = "R5341.S1.6", cols = v2.str.cols[levels(Idents(new.nano.obj))], alpha = 0.3, 
             molecules = c("NPHS1", "SLC5A12","SLC12A1","AQP2"), 
             mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("nanostring/Plots/R5341.S1.10_molecule_NPHS1-SLC5A12-SLC12A1-AQP2_Plot.pdf", width = 10, height = 8)
ImageDimPlot(new.nano.obj, fov = "R5341.S1.10", cols = v2.str.cols[levels(Idents(new.nano.obj))], alpha = 0.3, 
             molecules = c("NPHS1", "SLC5A12","SLC12A1","AQP2"), 
             mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

P71 = createPalette(71,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P71)
names(P71) <- levels(Idents(new.nano.obj))
ImageDimPlot(new.nano.obj, fov = "R5341.S1.6", cols = as.character(P71), alpha = 0.3, 
             molecules = c("NPHS1", "SLC5A12","SLC12A1","AQP2"), 
             mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)



###Structure Plots
Idents(new.nano.obj) <- "v2.structure"
fovs <- unique(new.nano.obj$fov)

lapply(fovs, function(x) {
  print(paste("Running for FOV:", x))
  fov <- gsub("_",".",x)
  DefaultBoundary(new.nano.obj[[fov]]) <- "segmentation"
  cl.plot <- ImageDimPlot(new.nano.obj, fov = fov, cols = v2.str.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
                          size = 0.5,border.size = 0.5, border.color = "black",dark.background = F)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pdf(file=paste0("nanostring/Image_Plots/fov_plots/",x,"structure_segmentation_Plot.pdf"),
      width=10,height=8)
  print(cl.plot)
  dev.off()
  
})


Idents(new.nano.obj) <- "v2.subclass.l3"
order <- c("POD","PEC","PT-S1/S2","PT-S3","aPT1","aPT2","frPT","DTL/ATL","M-TAL",
           "C-TAL","MD","aTAL","frTAL","dTAL","DCT","CNT","aCNT","C-PC","M-PC",
           "IC-A","IC-B","EC","EC-GC","EC-AEA/DVR","EC-AVR/V","M-EC-PTC",
           "EC-PTC","EC-LYM","FIB","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvMYOF","MC","REN","VSMC","VSMC/P",
           "B","PL","T","CD8+ T / NK","ERY","MAST","resMAC","moMAC-INF","moMAC/DC",
           "ncMON/N","mDC","cyc")





###Cell type sets
##Renal corpuscles
Glasbey = glasbey.colors(11)
swatch(Glasbey)
v2.rc.cols <- setNames(c("#0000FF","#FF0000","#00FF00","#FF00B6","#005300",
                         "#009FFF","#9A4D42","#00FFBE"),
                       c("POD","PEC","MD","EC-GC","EC-AEA/DVR",
                         "MC","REN","VSMC"))

Idents(new.nano.obj) <- "v2.subclass.l3"
fovs <- c("R5341.S1.2","R5341.S1.5","R5341.S1.6","R5341.S1.20","R5446.S2.15","R5446.S2.3")
pdf(file="nanostring/Image_Plots/Renal_Corpuscle_segmentation_Plot_R5341.S1.2.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5341.S1.2", cols = v2.rc.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("POD","PEC","MD","EC-GC","EC-AEA/DVR",
                                                         "MC","REN","VSMC")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



##Epithelial Cells
Glasbey = glasbey.colors(25)
swatch(Glasbey)
v2.epi.cols <- setNames(c("#0000FF","#009FFF",
                          "#FF00B6","#B1CC71","#00FF00",
                          "#005300","#FFD300",
                          "#9A4D42","#886C00","#1F9698","#766C95",
                          "#FFFFFF","#FFB79F","#000033",
                          "#FF0000"),
                        c("PT-S1/S2","PT-S3",
                          "DTL/ATL","M-TAL","C-TAL",
                          "DCT","CNT",
                          "C-PC","M-PC","IC-A","IC-B",
                          "EC-AEA/DVR","EC-AVR/V","M-EC-PTC",
                          "VSMC/P"))

fov = c("R5341.S1.6","R5341.S1.2","R5341.S1.16","R5341.S1.10","R5341.S1.11")
pdf(file="nanostring/Image_Plots/Tubules_segmentation_Plot_R5341.S1.11.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5341.S1.11", cols = v2.epi.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("PT-S1/S2","PT-S3",
                                                         "DTL/ATL","M-TAL","C-TAL",
                                                         "DCT","CNT",
                                                         "C-PC","M-PC","IC-A","IC-B",
                                                         "EC-AEA/DVR","EC-AVR/V","M-EC-PTC",
                                                         "VSMC/P")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()



##Fibrotic areas
P33 = createPalette(33,  c("#ff0000", "#00ff00", "#0000ff"))
swatch(P33)
Glasbey = glasbey.colors(30)
swatch(Glasbey)


new.nano.obj$v2.plot <- new.nano.obj$v2.structure
Idents(new.nano.obj) <- "v2.structure"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("distal tubules" ,
                                                           "collecting tubules","intermediate tubules")),]$v2.plot <- "other epi"
Idents(new.nano.obj) <- "v2.subclass.l3"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("PT-S1/S2","PT-S3")),]$v2.plot <- "PT"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("aPT1","aPT2")),]$v2.plot <- "aPT"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("frPT","frTAL")),]$v2.plot <- "frPT/frTAL"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("M-TAL","C-TAL","MD")),]$v2.plot <- "TAL"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("aTAL")),]$v2.plot <- "aTAL"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("B","PL","T","CD8+ T / NK")),]$v2.plot <- "Lymphoid"
new.nano.obj@meta.data[WhichCells(new.nano.obj, idents = c("MAST","resMAC","moMAC-INF",
                                                           "moMAC/DC","ncMON/N","mDC")),]$v2.plot <- "Myeloid"


Idents(new.nano.obj) <- "v2.plot"
v2.fn.cols <- setNames(c("#00479E","#1F9698",
                         "#005300","#FE8F42","#DC5E93",
                         "#783FC1","#A10300","#886C00",
                         "#C8FF00","#14F9FF",
                         "gray10"
                          ),
                        c("PT","aPT",
                          "TAL","aTAL","frPT/frTAL",
                          "vessels","renal corpuscle","interstitium - stroma",
                          "Lymphoid","Myeloid",
                          "other epi"))

fov <- c("R5446.S2.19","R5446.S2.17","R5446.S2.15","R5446.S2.12",
         "R5446.S2.1","R5446.S2.2")


pdf(file="nanostring/Image_Plots/Fibrosis_segmentation_Plot_R5446.S2.2.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5446.S2.2", cols = v2.fn.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("PT","aPT",
                                                         "TAL","aTAL","frPT/frTAL",
                                                         "vessels","renal corpuscle","interstitium - stroma",
                                                         "Lymphoid","Myeloid",
                                                         "other epi")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf(file="nanostring/Image_Plots/Fibrosis_segmentation_Plot_R5446.S2.19.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5446.S2.19", cols = v2.fn.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("PT","aPT",
                                                         "TAL","aTAL","frPT/frTAL",
                                                         "vessels","renal corpuscle","interstitium - stroma",
                                                         "Lymphoid","Myeloid",
                                                         "other epi")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf(file="nanostring/Image_Plots/Fibrosis_segmentation_Plot_R5446.S2.17.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5446.S2.17", cols = v2.fn.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("PT","aPT",
                                                         "TAL","aTAL","frPT/frTAL",
                                                         "vessels","renal corpuscle","interstitium - stroma",
                                                         "Lymphoid","Myeloid",
                                                         "other epi")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf(file="nanostring/Image_Plots/Fibrosis_segmentation_Plot_R5446.S2.15.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5446.S2.15", cols = v2.fn.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("PT","aPT",
                                                         "TAL","aTAL","frPT/frTAL",
                                                         "vessels","renal corpuscle","interstitium - stroma",
                                                         "Lymphoid","Myeloid",
                                                         "other epi")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()




#Plot fibroblast subtypes
Idents(new.nano.obj) <- "v2.subclass.l3"
v2.fn.cols <- setNames(c("#9A4D42","#009FFF","#337ca6","#FF0000",
                         "#FFD300","#FFD300","#00FF00"),
                       c("FIB","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF",
                         "pvFIB-RSPO3+","pvFIB-PI16+","pvMYOF"))
pdf(file="nanostring/Image_Plots/Fibrosis_segmentation_Plot_R5446.S2.17_FIB.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5446.S2.17", cols = v2.fn.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("FIB","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF",
                                                         "pvFIB-RSPO3+","pvFIB-PI16+","pvMYOF")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#Plot immune subtypes
Idents(new.nano.obj) <- "v2.subclass.l3"

v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+",
                         "moFAM","moMAC-C3+","cDC2","ncMON","MON"))

v2.fn.cols <- setNames(c("#00479E","#009FFF",
                         "#FFACFD","#720055",
                         "#886C00","#005300","#00FF00",
                         "#FFD300","#A10300"),
                       c("B","PL",
                         "T","CD8+ T / NK",
                         "MAST","resMAC","moMAC-INF",
                         "moMAC/DC","ncMON/N"))
pdf(file="nanostring/Image_Plots/Fibrosis_segmentation_Plot_R5446.S2.17_IMM.pdf",
    width=8,height=6)
ImageDimPlot(new.nano.obj, fov = "R5446.S2.17", cols = v2.fn.cols[levels(Idents(new.nano.obj))], coord.fixed = FALSE,
             size = 0.5,border.size = 0.1, border.color = "white",dark.background = T,
             cells = WhichCells(new.nano.obj, idents = c("B","PL",
                                                         "T","CD8+ T / NK",
                                                         "MAST","resMAC","moMAC-INF",
                                                         "moMAC/DC","ncMON/N")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
