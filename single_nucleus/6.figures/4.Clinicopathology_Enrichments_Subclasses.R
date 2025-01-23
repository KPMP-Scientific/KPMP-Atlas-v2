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

#order
order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3",
           "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

###Read in sample and pvalue tables
cond.ref.aki.pVal <- readRDS("Condition_Ref-AKI_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ckd.pVal <- readRDS("Condition_Ref-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.aki.ckd.pVal <- readRDS("Condition_AKI-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ati.pVal <- readRDS("Condition_Ref-AKI-ATI_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ati.ain.pVal <- readRDS("Condition_Ref-ATI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ref.ckd.hi.pVal <- readRDS("Condition_Ref-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ckd.lo.ckd.hi.pVal <- readRDS("Condition_CKD-LowRisk-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values.RDS")

adj.cond.ref.aki.pVal <- readRDS("Adj-Condition_Ref-AKI_Group_composition_V2_Subclassl3_P-Values.RDS")
adj.cond.ref.DKD.pVal <- readRDS("Adj-Condition_Ref-DKD_Group_composition_V2_Subclassl3_P-Values.RDS")
adj.cond.ref.HCKD.pVal <- readRDS("Adj-Condition_Ref-H-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")
adj.cond.HCKD.DKD.pVal <- readRDS("Adj-Condition_H-CKD-DKD_Group_composition_V2_Subclassl3_P-Values.RDS")
adj.cond.AKI.DKD.pVal <- readRDS("Adj-Condition_AKI-DKD_Group_composition_V2_Subclassl3_P-Values.RDS")
adj.cond.AKI.HCKD.pVal <- readRDS("Adj-Condition_AKI-H-CKD_Group_composition_V2_Subclassl3_P-Values.RDS")


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
AS.pVal <- readRDS("Arteriosclerosis_Group_composition_V2_Subclassl3_P-Values.RDS")
AH.pVal <- readRDS("Arteriolar_Hyalinosis_Group_composition_V2_Subclassl3_P-Values.RDS")


###P value plot
row.order <- order
row.order <- gsub("-","",row.order)
row.order <- gsub(" ","",row.order)

p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]
write.table(p.table[row.order,], file = "QC_Plots/Composition_Analysis_p-Values_addit_Adj_Groups.txt", sep = "\t", quote = FALSE)

###Clinical Phenotype Plots
col.order <- c("cond.ref.aki.pVal",
               "cond.aki.ckd.pVal",
               "cond.ref.ckd.pVal",
               "eGFR.pVal", "age.pVal")
p.stats <- p.table[row.order,paste0(col.order,".p")]

states <- setNames(KB$v2.state.l2,KB$v2.subclass.l3)[order]
ColSideColors <- as.character(state.l2.cols[as.character(states)])
p.stats <- p.table[row.order,paste0(col.order,".p")]

pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Clinical_ordered_subclass_level3.pdf',width=30,height=5)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Clinical_ordered_subclass_level3_B.pdf',width=30,height=4)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

#additional adjudication categories
col.order <- c("adj.cond.ref.aki.pVal","cond.ref.ati.pVal","cond.ref.ain.pVal",
               "cond.ref.ckd.hi.pVal","cond.ckd.lo.ckd.hi.pVal",
               "adj.cond.ref.DKD.pVal","adj.cond.AKI.DKD.pVal",
               "adj.cond.ref.HCKD.pVal","adj.cond.AKI.HCKD.pVal","adj.cond.HCKD.DKD.pVal")
p.stats <- p.table[row.order,paste0(col.order,".p")]

states <- setNames(KB$v2.state.l2,KB$v2.subclass.l3)[order]
ColSideColors <- as.character(state.l2.cols[as.character(states)])
p.stats <- p.table[row.order,paste0(col.order,".p")]

pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Clinical_ordered_subclass_level3_Adju_Cat.pdf',width=30,height=6)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 
pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Clinical_ordered_subclass_level3_Adju_Cat_B.pdf',width=30,height=5)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 




###Pathology Descriptor plots
col.order <- paste0(c("IF","TA","AC","TI","WBC","AS","AH"), ".pVal")
p.stats <- p.table[row.order,paste0(col.order,".p")]

states <- setNames(KB$v2.state.l2,KB$v2.subclass.l3)[order]
ColSideColors <- as.character(state.l2.cols[as.character(states)])
p.stats <- p.table[row.order,paste0(col.order,".p")]

pdf(file='QC_Plots/Composition_Analysis_t-statistic_Heatmap_Pathology_ordered_subclass_level3.pdf',width=30,height=5)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 





###Plotting Pathway scores for subclasses
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_09122024.RDS")
KB <- AddMetaData(KB, meta.sub)

# Extract mlm and store it in pathwaysmlm in KB
KB[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(meta.sub[colnames(KB),])))

# Change assay
DefaultAssay(object = KB) <- "pathwaysmlm"

# Scale the data
KB <- ScaleData(KB)
KB@assays$pathwaysmlm@data <- KB@assays$pathwaysmlm@scale.data

Idents(KB) <- "v2.subclass.l3"
Idents(KB) <- factor(Idents(KB), levels = order)

# Extract activities from object as a long dataframe
df <- t(as.matrix(KB@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB)) %>%
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
brewer.pal(11,"RdBu")
my_color = colorRampPalette(c("#0852a6", "white","#c4281f"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

to.use <- c("ProtAge-UCell","ProtAge20-UCell","Plasma-ATI-CHROME-UCell",
            "Plasma-ATI-ARIC-UCell","Plasma-ATI-BKBC-UCell",
            "CKD-Progression-UCell","CRIC-ARIC-CKD-UCell"
)
# Plot All
pdf(file='QC_Plots/Plasma_Gene_Set_Enrichment_Heatmap_ordered_subclass_level3.pdf',width=30,height=5)
pheatmap(t(top_acts_mat[,to.use]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D2", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

