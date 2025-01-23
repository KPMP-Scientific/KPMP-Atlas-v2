library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
library("RColorBrewer")
library(gplots)
library(dendextend)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")
source("misc/utils.R")


load("color_factors.robj")
load("color_factors_v2-clusters.robj")
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB


#load order of clusters based on HC plot
load("Plots/Atlas_V2_Cluster_Dendrogram_05082024_label-order.rda")
order

###Read in pvalue tables
cond.ref.aki.pVal <- readRDS("Condition_Ref-AKI_Group_composition_V2_Clusters_P-Values.RDS")
cond.ref.ckd.pVal <- readRDS("Condition_Ref-CKD_Group_composition_V2_Clusters_P-Values.RDS")
cond.aki.ckd.pVal <- readRDS("Condition_AKI-CKD_Group_composition_V2_Clusters_P-Values.RDS")
cond.ref.ati.pVal <- readRDS("Condition_Ref-AKI-ATI_Group_composition_V2_Clusters_P-Values.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Clusters_P-Values.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Clusters_P-Values.RDS")
cond.ati.ain.pVal <- readRDS("Condition_Ref-ATI-AIN_Group_composition_V2_Clusters_P-Values.RDS")
cond.ref.ckd.hi.pVal <- readRDS("Condition_Ref-CKD-HighRisk_Group_composition_V2_Clusters_P-Values.RDS")
cond.ckd.lo.ckd.hi.pVal <- readRDS("Condition_CKD-LowRisk-CKD-HighRisk_Group_composition_V2_Clusters_P-Values.RDS")
eGFR.pVal <- readRDS("eGFR_Group_composition_V2_Clusters_P-Values.RDS")
age.pVal <- readRDS("HRT_Age_Group_composition_V2_Clusters_P-Values.RDS")
IF.pVal <- readRDS("Interstitial-Fibrosis_Group_composition_V2_Clusters_P-Values.RDS")
TA.pVal <- readRDS("Tubular-Atrophy_Group_composition_V2_Clusters_P-Values.RDS")
TI.pVal <- readRDS("Tubular-Injury_Group_composition_V2_Clusters_P-Values.RDS")
ALM.pVal <- readRDS("Abnormal_Luminal_Morphology_Group_composition_V2_Clusters_P-Values.RDS")
AC.pVal <- readRDS("Acellular_Cast_Group_composition_V2_Clusters_P-Values.RDS")
WBC.pVal <- readRDS("Interstitial_WBC_Group_composition_V2_Clusters_P-Values.RDS")
AS.pVal <- readRDS("Arteriosclerosis_Group_composition_V2_Clusters_P-Values.RDS")
AH.pVal <- readRDS("Arteriolar_Hyalinosis_Group_composition_V2_Clusters_P-Values.RDS")

###P value plot
row.order <- c(paste0("P_",1:33), paste0("D_",1:48),paste0("E_",1:22),paste0("S_",1:27),paste0("I_",1:28),paste0("N_",1))
p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]
write.table(p.table[row.order,], file = "QC_Plots/Composition_Analysis_p-Values.txt", sep = "\t", quote = FALSE)

###Clinical Phenotype Plots
col.order <- c("cond.ref.aki.pVal","cond.ref.ati.pVal","cond.ref.ain.pVal",
               "cond.aki.ckd.pVal",
               "cond.ref.ckd.pVal","cond.ref.ckd.hi.pVal","cond.ckd.lo.ckd.hi.pVal",
               "eGFR.pVal", "age.pVal")
p.stats <- p.table[row.order,paste0(col.order,".p")]
states <- setNames(KB$v2.state.l2,KB$v2.clusters)[order]
ColSideColors <- as.character(state.l2.cols[as.character(states)])
p.stats <- p.table[order,paste0(col.order,".p")]

pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Clinical_ordered.pdf',width=30,height=5.1)
heatmap.2(as.matrix(t(t.stats[order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 


###Pathology Descriptor plots
col.order <- paste0(c("IF","TA","AC","TI","WBC","AS","AH"), ".pVal")
p.stats <- p.table[row.order,paste0(col.order,".p")]
states <- setNames(KB$v2.state.l2,KB$v2.clusters)[order]
ColSideColors <- as.character(state.l2.cols[as.character(states)])
p.stats <- p.table[order,paste0(col.order,".p")]

pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Pathology_ordered.pdf',width=30,height=5)
heatmap.2(as.matrix(t(t.stats[order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 
