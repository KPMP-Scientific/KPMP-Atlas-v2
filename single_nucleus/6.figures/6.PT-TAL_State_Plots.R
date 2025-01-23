library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)
library(slingshot)
library(Polychrome)
library(dendextend)
library("RColorBrewer")
library(gplots)
library(tibble)
library(tidyr)
library(pheatmap)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")


###PT
##UMAP Plots - color by subclass.l3
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_slingshot_0424-newData.rda")

v2.fn.cols <- setNames(c("#00479E","#009FFF","#0000FF",
                         "#14F9FF","#00FFBE","#1F9698","#B1CC71","#00FFBE","#DC5E93","#A10300","#1F9698","#1F9698","#1F9698","#A10300"),
                       c("P_4","P_8","P_10",
                         "P_12","P_13","P_14","P_15","P_16","P_17","P_18","P_19","P_20","P_21","P_22"))

pdf(file='trajectories/PT-Lineages_Slingshot_subclassl3_colors.pdf',width=6,height=5)
plot(reducedDims(sce)$UMAP, col = adjustcolor(v2.fn.cols[sce$v2.clusters], alpha.f = 0.5), pch = 19, cex = 0.4)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1,4))
lines(SlingshotDataSet(sce), lwd = 3, col = "darkred", type = 'curves', linInd = c(2,5))
lines(SlingshotDataSet(sce), lwd = 3, col = "gray30", type = 'curves', linInd = c(3))
dev.off()




##Clinical/Pathological Heatmaps (Subclasses)
#Read in sample and pvalue tables
cond.ref.aki.pVal <- readRDS("Condition_Ref-AKI_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ckd.pVal <- readRDS("Condition_Ref-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.aki.ckd.pVal <- readRDS("Condition_AKI-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ati.pVal <- readRDS("Condition_Ref-AKI-ATI_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ati.ain.pVal <- readRDS("Condition_Ref-ATI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ckd.hi.pVal <- readRDS("Condition_Ref-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ckd.lo.ckd.hi.pVal <- readRDS("Condition_CKD-LowRisk-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values.RDS")
adj.cond.HCKD.DKD.pVal <- readRDS("Adj-Condition_H-CKD-DKD_Group_composition_V2_Subclassl3_P-Values.RDS")

eGFR.pVal <- readRDS("eGFR_Group_composition_V2_Subclassl3_P-Values.RDS")
age.pVal <- readRDS("HRT_Age_Group_composition_V2_Subclassl3_P-Values.RDS")
A1C.pVal <- readRDS("A1C_Group_composition_V2_Subclassl3_P-Values.RDS")
Ptn.pVal <- readRDS("Proteinuria_Group_composition_V2_Subclassl3_P-Values.RDS")
Alb.pVal <- readRDS("Albuminuria_Group_composition_V2_Subclassl3_P-Values.RDS")
IF.pVal <- readRDS("Interstitial-Fibrosis_Group_composition_V2_Subclassl3_P-Values.RDS")
TA.pVal <- readRDS("Tubular-Atrophy_Group_composition_V2_Subclassl3_P-Values.RDS")
TI.pVal <- readRDS("Tubular-Injury_Group_composition_V2_Subclassl3_P-Values.RDS")
ALM.pVal <- readRDS("Abnormal_Luminal_Morphology_Group_composition_V2_Subclassl3_P-Values.RDS")
AC.pVal <- readRDS("Acellular_Cast_Group_composition_V2_Subclassl3_P-Values.RDS")
WBC.pVal <- readRDS("Interstitial_WBC_Group_composition_V2_Subclassl3_P-Values.RDS")

#P value plot
row.order <- c("aPT2","aPT1","aPTS1/S2","PTS1","PTS2","frPTS1/S2",
               "aPTS3","PTS3","frPTS3")
p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]

#Clinical Plot
col.order <- c("cond.ref.aki.pVal","cond.ref.ati.pVal","cond.ref.ain.pVal",
               "cond.aki.ckd.pVal",
               "cond.ref.ckd.pVal","cond.ref.ckd.hi.pVal","cond.ckd.lo.ckd.hi.pVal",
               "adj.cond.HCKD.DKD.pVal",
               "eGFR.pVal","age.pVal")

ColSideColors <- c("#B3B336","#B3B336","#B3B336","#8eab9b","#8eab9b","#B3B336","#B3B336","#8eab9b","#B3B336")

pdf(file='trajectories/PT-Lineages_Clinical-Path_Heatmap_subclasses_clin.pdf',width=6,height=5)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(5, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

#Pathology Plot
col.order <- c("IF.pVal","TA.pVal","TI.pVal","WBC.pVal")

pdf(file='trajectories/PT-Lineages_Clinical-Path_Heatmap_subclasses_path.pdf',width=6,height=5)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(20, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

##Marker Gene Plots
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_proximal-tubules_subset_filtered_aPT_0424-newData.rda")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.PT@meta.data$v2.subclass.l3 <- KB@meta.data[rownames(KB.PT@meta.data),]$v2.subclass.l3

Idents(KB.PT) <- "v2.subclass.l3"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = c("aPT2","aPT1","aPT-S1/S2","PT-S1","PT-S2","frPT-S1/S2",
                                                  "aPT-S3","PT-S3","frPT-S3"))

pt.markers <- c(
  "LRP2","CUBN","HNF1B",                                        #PT
  "IGFBP7","SPP1","ITGB8","CDH6","TMEM178B","ALPK2","HAVCR1",
  
  "ITGB3",
  "CST3","CLU","VIM","PIGR",#"APOE",                            #aPT2
  "IL32","SOX4","VCAM1","MMP7","SOX9","CCL2",                   #aPT2
  
  "DCC",                                                        #aPT1
  "GDA","GLIS1",                                                #aPT-S1/S2
  "SLC5A12","SLC22A6","HNF4A",                                  #S1/S2
  "RALYL","PCDH15","PRODH2","SLC22A8",                          #S1
  "SLC34A1","ANKS1B","SLC5A10","SLC5A11",                       #S2                                  
  "DLGAP1","PROM1",                                             #frPTS1/S2
  "APBB1IP","ROBO2","COL23A1","MEG3",                           #frPTS1/S2
  "LSAMP","KCNIP1","NRXN3","WNT2B",                             #frPTS1/S2
  "KCTD16","SPON1",                                             #aPT-S3
  "SLC5A8","GPM6A","SLC22A24","SLC7A13",                        #PT-S3
  "NRG3","FAM189A1","DTNA","KITLG","GRM8"                       #frPT-S3
  
)

pdf(file='trajectories/PT-Lineages_Marker_Dotplot_subclassl3_full.pdf',width=14,height=4)
DotPlot(KB.PT, features = pt.markers) + RotatedAxis()
dev.off()


pt.markers <- c(
  "ITGB8","CDH6","HAVCR1",
  
  "ITGB3",
  "IL32","SOX4","VCAM1",                   #aPT2
  
  "DCC",                                                        #aPT1
  "GDA",                                                #aPT-S1/S2
  "HNF4A",                                  #S1/S2
  "PROM1",                                             #frPTS1/S2
  "ROBO2","MEG3",                           #frPTS1/S2
  "SPON1",                                             #aPT-S3
  "SLC7A13",                       #PT-S3
  "KITLG"                       #frPT-S3
  
)

pdf(file='trajectories/PT-Lineages_Marker_Dotplot_subclassl3_subset.pdf',width=6.5,height=4)
DotPlot(KB.PT, features = pt.markers) + RotatedAxis()
dev.off()


##Plot TF Activities
library(Signac)
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_PT_Dual_Seurat_PT-aPT_Only_05202024.rda")


#TFBS activities - All PT subclass.l3
DefaultAssay(KB.PT) <- "chromvar"
Idents(KB.PT) <- "v2.subclass.l3"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = c("aPT2","aPT1","aPT-S1/S2","PT-S1","PT-S2",
                                                  "frPT-S1/S2",
                                                  "aPT-S3","PT-S3","frPT-S3"))
table(Idents(KB.PT))
rownames(KB.PT@assays$chromvar@data) <- ConvertMotifID(KB.PT,
                                                       id=rownames(KB.PT@assays$chromvar@data),
                                                       assay='ATAC')
rownames(KB.PT@assays$chromvar@data) <- toupper(rownames(KB.PT@assays$chromvar@data))
TFs.plot <- c("REL","CEBPD","STAT1","KLF6","SOX4","SOX9",
              "SP1","HES1","EMX2","ZNF148",
              "MEIS1","NR2F1",
              "HNF4A","PPARD","THRB","RXRA")

pdf(file='trajectories/PT-Lineages_TFBS_Dotplot_subclassl3.pdf',width=6.5,height=4)
DotPlot(KB.PT, features = TFs.plot, dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
dev.off()

TFs.plot <- c("REL","KLF6","SOX9",
              "MEIS1","NR2F1",
              "HNF4A","RXRA")

pdf(file='trajectories/PT-Lineages_TFBS_Dotplot_subclassl3_subset.pdf',width=4.7,height=4)
DotPlot(KB.PT, features = TFs.plot, dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
dev.off()





###TAL
##UMAP Plots - color by subclass.l3
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_distal-tubules_subset_filtered_aTAL_slingshot_0424-newData_TALA.rda")

v2.fn.cols <- setNames(c("#005300","#02AD24","#005300","#02AD24","#009FFF","#14F9FF","#14F9FF",
                         "#DC5E93","#FFD300","#FE8F42"),
                       c("D_5","D_6","D_7","D_8","D_11","D_12","D_13",
                         "D_15","D_16","D_17"))

pdf(file='trajectories/TAL-Lineages_Slingshot_subclassl3_colors.pdf',width=6,height=5)
plot(reducedDims(sce)$UMAP, col = adjustcolor(v2.fn.cols[sce$v2.clusters], alpha.f = 0.5), pch = 19, cex = 0.4)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1,2))
lines(SlingshotDataSet(sce), lwd = 3, col = "darkred", type = 'curves', linInd = c(3))
dev.off()


##Clinical/Pathological Heatmaps (Subclasses)
#Read in sample and pvalue tables
cond.ref.aki.pVal <- readRDS("Condition_Ref-AKI_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ckd.pVal <- readRDS("Condition_Ref-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.aki.ckd.pVal <- readRDS("Condition_AKI-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ati.pVal <- readRDS("Condition_Ref-AKI-ATI_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ati.ain.pVal <- readRDS("Condition_Ref-ATI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ckd.hi.pVal <- readRDS("Condition_Ref-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ckd.lo.ckd.hi.pVal <- readRDS("Condition_CKD-LowRisk-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values.RDS")
adj.cond.HCKD.DKD.pVal <- readRDS("Adj-Condition_H-CKD-DKD_Group_composition_V2_Subclassl3_P-Values.RDS")
eGFR.pVal <- readRDS("eGFR_Group_composition_V2_Subclassl3_P-Values.RDS")
age.pVal <- readRDS("HRT_Age_Group_composition_V2_Subclassl3_P-Values.RDS")
A1C.pVal <- readRDS("A1C_Group_composition_V2_Subclassl3_P-Values.RDS")
Ptn.pVal <- readRDS("Proteinuria_Group_composition_V2_Subclassl3_P-Values.RDS")
Alb.pVal <- readRDS("Albuminuria_Group_composition_V2_Subclassl3_P-Values.RDS")
IF.pVal <- readRDS("Interstitial-Fibrosis_Group_composition_V2_Subclassl3_P-Values.RDS")
TA.pVal <- readRDS("Tubular-Atrophy_Group_composition_V2_Subclassl3_P-Values.RDS")
TI.pVal <- readRDS("Tubular-Injury_Group_composition_V2_Subclassl3_P-Values.RDS")
ALM.pVal <- readRDS("Abnormal_Luminal_Morphology_Group_composition_V2_Subclassl3_P-Values.RDS")
AC.pVal <- readRDS("Acellular_Cast_Group_composition_V2_Subclassl3_P-Values.RDS")
WBC.pVal <- readRDS("Interstitial_WBC_Group_composition_V2_Subclassl3_P-Values.RDS")

#P value plot
row.order <- c("aTAL1","aTAL2", "C/MTALA","CTALA","C/MTALB","CTALB","frTAL")
p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]

#Clinical Plot
col.order <- c("cond.ref.aki.pVal","cond.ref.ati.pVal","cond.ref.ain.pVal",
               "cond.aki.ckd.pVal",
               "cond.ref.ckd.pVal","cond.ref.ckd.hi.pVal","cond.ckd.lo.ckd.hi.pVal",
               "adj.cond.HCKD.DKD.pVal",
               "eGFR.pVal","age.pVal")

ColSideColors <- c("#B3B336","#8eab9b","#8eab9b","#8eab9b","#8eab9b","#B3B336","#B3B336")
p.stats <- p.table[row.order,paste0(col.order,".p")]

pdf(file='trajectories/TAL-Lineages_Clinical-Path_Heatmap_subclasses_clin.pdf',width=5,height=5)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(5, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

#Pathology Plot
col.order <- c("IF.pVal","TA.pVal","TI.pVal","WBC.pVal")

pdf(file='trajectories/TAL-Lineages_Clinical-Path_Heatmap_subclasses_path.pdf',width=5,height=5)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(20, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

##Marker Gene Plots
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.TAL <- subset(KB, v2.subclass.l3 %in% c("aTAL1","aTAL2", "C/M-TAL-A","C-TAL-A",
                                           "C/M-TAL-B","C-TAL-B","frTAL"))
Idents(KB.TAL) <- "v2.subclass.l3"
Idents(KB.TAL) <- factor(Idents(KB.TAL), levels = c("aTAL1","aTAL2", "C/M-TAL-A","C-TAL-A",
                                                  "C/M-TAL-B","C-TAL-B","frTAL"))
KB.TAL <- subset(KB.TAL, downsample = 1000)
KB.TAL[["RNA"]] <- as(KB.TAL[["RNA"]], Class = "Assay")
KB.TAL <- NormalizeData(KB.TAL)

TALA.B.markers <- FindMarkers(KB.TAL, ident.1 = c("C/M-TAL-A","C-TAL-A"), ident.2 = c("C/M-TAL-B","C-TAL-B"), min.pct = 0.25)

tal.markers <- c(
  "CASR","SLC12A1","UMOD","EGF","ESRRB",                       #TAL
  "CREB5","COL4A1","ITGA3","HIVEP3","CYTOR","DIAPH3","SYT16",      #aTAL1
  "SPP1","IGFBP7","ARHGEF38","ITGB6","NRP1","TFPI","ADAMTS1","DNAH5","HAVCR1",   #aTAL2
  
  
  "PHACTR1","SLCO3A1","CXCL12","CNTN1","CABP1","KCNMB2","RGS6",               #C-TAL-A
  "ENOX1","CALCR","RBM20","PDE3A",
  
  "DACH1","LRMDA",              #C-TAL-B
           
  "TENM4","FGF14","PXDNL","GRM8","KCNJ10",                     #C-TAL-B
  "TMEM52B","CLDN16","SCN7A","TMEM207","JAG1",                          
  "COL8A1","LINGO2",
  
  
  "ITGB8","PROM1",                                       
  "ARHGAP26","RNF144B","TMPRSS4","RHEX","CLU"                      #frTAL
)

pdf(file='trajectories/TAL-Lineages_Marker_Dotplot_subclassl3_full.pdf',width=14,height=3.2)
DotPlot(KB.TAL, features = tal.markers) + RotatedAxis()
dev.off()


tal.markers <- c(
  "SLC12A1","EGF",                       #TAL
  "CREB5","ITGA3",      #aTAL1
  "ITGB6","NRP1",   #aTAL2

  "PHACTR1","SLCO3A1",               #C-TAL-A
  "CALCR",
  
  "DACH1","LRMDA",              #C-TAL-B
  
  "SCN7A",                      #C-TAL-B

  "ITGB8","PROM1",                                       
  "TMPRSS4","RHEX"                      #frTAL
)

pdf(file='trajectories/TAL-Lineages_Marker_Dotplot_subclassl3_subset.pdf',width=6.5,height=3.2)
DotPlot(KB.TAL, features = tal.markers) + RotatedAxis()
dev.off()


##Plot TF Activities
library(Signac)
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_PT_Dual_Seurat_TAL-aTAL_Only_05212024.rda")

#TFBS activities - All TAL subclass.l3
DefaultAssay(KB.DT) <- "chromvar"
Idents(KB.DT) <- "v2.subclass.l3"
Idents(KB.DT) <- factor(Idents(KB.DT), levels = c("aTAL1","aTAL2","C/M-TAL-A","C-TAL-A",
                                                  "C/M-TAL-B","C-TAL-B","frTAL"))
table(Idents(KB.DT))
rownames(KB.DT@assays$chromvar@data) <- ConvertMotifID(KB.DT,
                                                       id=rownames(KB.DT@assays$chromvar@data),
                                                       assay='ATAC')
rownames(KB.DT@assays$chromvar@data) <- toupper(rownames(KB.DT@assays$chromvar@data))
TFs.plot <- c("FOS","EGR1","TFAP2B",
              "RUNX1","SOX4","SOX9",
              "MAZ","SP3","SMAD2","EBF1","ZNF148",
              "STAT3","ELF3",
              "GRHL2","NFIA",
              "ESRRA","ESRRB","NR2F2")

pdf(file='trajectories/TAL-Lineages_TFBS_Dotplot_subclassl3.pdf',width=6.5,height=3.2)
DotPlot(KB.DT, features = TFs.plot, dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
dev.off()

TFs.plot <- c("FOS","TFAP2B","SOX9","ELF3",
              
              "NFIA",
              "ESRRB","NR2F2")

pdf(file='trajectories/TAL-Lineages_TFBS_Dotplot_subclassl3_subset.pdf',width=4.7,height=3.2)
DotPlot(KB.DT, features = TFs.plot, dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
dev.off()





###Pathway scores (PT/TAL)
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_09122024.RDS")
KB <- AddMetaData(KB, meta.sub)

# Extract mlm and store it in pathwaysmlm in KB
KB[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(meta.sub[colnames(KB),])))

# Change assay
DefaultAssay(object = KB) <- "pathwaysmlm"

# Scale the data
KB <- ScaleData(KB)
KB@assays$pathwaysmlm@data <- KB@assays$pathwaysmlm@scale.data

KB.PT <- subset(KB, v2.subclass.l3 %in% c("aPT2","aPT1","aPT-S1/S2","PT-S1","PT-S2",
                                          "frPT-S1/S2",
                                          "aPT-S3","PT-S3","frPT-S3",
                                          "aTAL1","aTAL2","C/M-TAL-A","C-TAL-A",
                                          "C/M-TAL-B","C-TAL-B","frTAL"))
Idents(KB.PT) <- "v2.subclass.l3"
Idents(KB.PT) <- factor(Idents(KB.PT), levels = c("aPT2","aPT1","aPT-S1/S2","PT-S1","PT-S2",
                                                  "frPT-S1/S2",
                                                  "aPT-S3","PT-S3","frPT-S3",
                                                  "aTAL1","aTAL2","C/M-TAL-A","C-TAL-A",
                                                  "C/M-TAL-B","C-TAL-B","frTAL"))

# Extract activities from object as a long dataframe
df <- t(as.matrix(KB.PT@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.PT)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

# Transform to wide matrix
top_acts_mat <- df %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot Subset
to.use <- c("TGFb-UCell","SenMayo-UCell","Gly-UCell","MAPK-UCell",
            "NFkB-UCell","Inf-UCell","Aging-UCell","TNFa-UCell","JAK-STAT-UCell",
            "ECM-UCell","EMT-UCell","ProtAge-UCell",
            "FatAcid-UCell","OxPhos-UCell","Plasma-ATI-ARIC-UCell", "Plasma-ATI-BKBC-UCell",
            "CRIC-ARIC-CKD-UCell")
pdf(file='trajectories/PT-TAL-Lineages_Pathway_Gene-set_Scores_subclassl3_subset.pdf',width=7,height=9)
pheatmap(top_acts_mat[,to.use],
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D2", cluster_rows = FALSE) 
dev.off()
