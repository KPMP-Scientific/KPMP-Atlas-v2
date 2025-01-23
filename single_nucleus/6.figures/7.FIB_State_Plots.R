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


##UMAP Plots - color by subclass.l3
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_16-19_slingshot_0424-newData.rda")
v2.fn.cols <- setNames(c("#C8FF00","#FFD300",
                         "#02AD24","#00FF00"),
                       c("S_16","S_17",
                         "S_18","S_19"))

pdf(file='trajectories/pvFIB-Lineages_Slingshot_subclassl3_colors.pdf',width=6,height=5)
plot(reducedDims(sce)$UMAP, col = adjustcolor(v2.fn.cols[sce$v2.clusters], alpha.f = 0.5), pch = 19, cex = 0.4)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()

load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Stroma_subset_filtered_cortex-S_9-14_slingshot_0424-newData.rda")
v2.fn.cols <- setNames(c("#9A4D42","#0000FF","#0000FF","#009FFF","#337ca6","#FF0000"),
                       c("S_9","S_10","S_11","S_12","S_13","S_14"))

pdf(file='trajectories/intFIB-Lineages_Slingshot_subclassl3_colors.pdf',width=6,height=5)
plot(reducedDims(sce)$UMAP, col = adjustcolor(v2.fn.cols[sce$v2.clusters], alpha.f = 0.5), pch = 19, cex = 0.4)
lines(SlingshotDataSet(sce), lwd = 3, col = "gray30", type = 'curves', linInd = c(2,3))
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves', linInd = c(1,4))
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
row.order <- c("CFIB","CFIBPATH","CFIBOSMRlo","CFIBOSMRhi","CMYOF",
               "pvFIBRSPO3+","pvFIBPI16+","pvFIB","pvMYOF")
p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]

#Clinical Plot
col.order <- c("cond.ref.aki.pVal","cond.ref.ati.pVal","cond.ref.ain.pVal",
               "cond.aki.ckd.pVal",
               "cond.ref.ckd.pVal","cond.ref.ckd.hi.pVal","cond.ckd.lo.ckd.hi.pVal",
               "adj.cond.HCKD.DKD.pVal",
               "eGFR.pVal","age.pVal")

ColSideColors <- c("#8eab9b","#B3B336","#B3B336","#B3B336","#B3B336",
                   "#B3B336","#B3B336","#B3B336","#B3B336")
p.stats <- p.table[row.order,paste0(col.order,".p")]

pdf(file='trajectories/FIB-Lineages_Clinical-Path_Heatmap_subclasses_clin.pdf',width=6,height=5)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(5, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

#Pathology Plot
col.order <- c("IF.pVal","TA.pVal","TI.pVal","WBC.pVal")

pdf(file='trajectories/FIB-Lineages_Clinical-Path_Heatmap_subclasses_path.pdf',width=6,height=5)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(18, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 



##Marker Gene Plots
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.FIB <- subset(KB, v2.subclass.l3 %in% c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
                                           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))
Idents(KB.FIB) <- "v2.subclass.l3"
Idents(KB.FIB) <- factor(Idents(KB.FIB), levels = c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
                                                    "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))
KB.FIB <- subset(KB.FIB, downsample = 1000)
KB.FIB[["RNA"]] <- as(KB.FIB[["RNA"]], Class = "Assay")
KB.FIB <- NormalizeData(KB.FIB)

fib.markers <- c(
  "DCN","C7","PDGFRA",                                     #Pan FIB
  "NEGR1","LAMA2","ABCA8","MEG3",                           #Pan cortical FIB
  "CCN1","CCN2","ELL2","SAMHD1","SLC2A3",                   #C-FIB (interstitial fib)
  "GRIN2A","EMID1",                                         #C-FIB-PATH
  "SELENOP","LUM","CXCL12",'GGT5',"ECRG4",                  #C-FIB-OSMRlo
  "OSMR","SOD2","UGCG","IL1R1","CCL2",                      #C-FIB-OSMRhi
  "RELB","CXCL10","CCL19",                                  #C-FIB-OSMRhi
  "SULF1","GLI2","NTM","INHBA","FAP","POSTN",               #C-MYOF
  
  "FLRT2","COL12A1","FGF14",                                #Pan pvFIB
  "PDZRN4",'IGF1','ADAMTS3',"RSPO3","WNT5B",                #pvFIB-RSPO3+
  "C3","EBF2","SFRP2","CD34","PI16",                        #pvFIB-PI16+
  "ITGBL1","PLD5","CNTN5",                                  #pvFIB & pvMYOF
  "MGAT4C","EPHA3",                                         #pvFIB
  "ITGA8",                                                  #pvMYOF & VSMC
  "ADGRL3","MYH11","ACTA2",'KCNMA1',"PCDH7",                #pvMYOF
  "PRUNE2","MYOCD" #"SYNPO2",'MACROD2',                     #pvMYOF
  )

pdf(file='trajectories/FIB-Lineages_Marker_Dotplot_subclassl3_full.pdf',width=14,height=4)
DotPlot(KB.FIB, features = fib.markers) + RotatedAxis()
dev.off()

fib.markers <- c(
  "C7",                                     #Pan FIB
  "MEG3",                           #Pan cortical FIB
  "CCN1",                   #C-FIB (interstitial fib)
  #C-FIB-PATH
  "SELENOP","CXCL12",                  #C-FIB-OSMRlo
  "OSMR","IL1R1","CCL2",                      #C-FIB-OSMRhi
  "SULF1","FAP","POSTN",               #C-MYOF
  "FLRT2",                                #Pan pvFIB
  'IGF1',"RSPO3",                #pvFIB-RSPO3+
  "C3","CD34","PI16",                        #pvFIB-PI16+
  "EPHA3",                                         #pvFIB
  "MYH11","ACTA2"                     #pvMYOF
)

pdf(file='trajectories/FIB-Lineages_Marker_Dotplot_subclassl3_subset.pdf',width=7.5,height=3.5)
DotPlot(KB.FIB, features = fib.markers) + RotatedAxis()
dev.off()



##Plot TF Activities
library(Signac)
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Stroma_Dual_Seurat_Cortex_Only_05172024.rda")

#TFBS activities - subclass.l3
DefaultAssay(KB.STR) <- "chromvar"
Idents(KB.STR) <- "v2.subclass.l3"
Idents(KB.STR) <- factor(Idents(KB.STR), levels = c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
                                                    "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))
table(Idents(KB.STR))
rownames(KB.STR@assays$chromvar@data) <- ConvertMotifID(KB.STR,
                                                       id=rownames(KB.STR@assays$chromvar@data),
                                                       assay='ATAC')
rownames(KB.STR@assays$chromvar@data) <- toupper(rownames(KB.STR@assays$chromvar@data))
TFs.plot <- c("ATF3","RFX2",
              "IRF1","REL","RELB","NFKB1","NFKB2","STAT1","STAT3","SOX4","PLAG1","PLAGL1","EGR3",
              "RUNX1","RUNX2","TCF4","SMAD3",
  "GLI1","GLI2","GLI3","KLF4","KLF5","SOX18","OSR1",
  #"KLF3","KLF8","KLF12",
  "MEF2C")

pdf(file='trajectories/FIB-Lineages_TFBS_Dotplot_subclassl3.pdf',width=8,height=3.5)
DotPlot(KB.STR, features = TFs.plot, dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
dev.off()

TFs.plot <- c("ATF3",
              "IRF1","REL",
              "TCF4","SMAD3","RUNX1","PLAGL1","EGR3",
              "KLF4","GLI1",
              "MEF2C")

pdf(file='trajectories/FIB-Lineages_TFBS_Dotplot_subclassl3_subset.pdf',width=5.7,height=3.5)
DotPlot(KB.STR, features = TFs.plot, dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
dev.off()




##Pathway scores
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

KB.FIB <- subset(KB, v2.subclass.l3 %in% c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
                                           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))
Idents(KB.FIB) <- "v2.subclass.l3"
Idents(KB.FIB) <- factor(Idents(KB.FIB), levels = c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
                                                    "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))


# Extract activities from object as a long dataframe
df <- t(as.matrix(KB.FIB@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.FIB)) %>%
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

# Plot All
pdf(file='trajectories/FIB-Lineages_Pathway_Gene-set_Scores_subclassl3.pdf',width=7,height=5)
pheatmap(top_acts_mat[,!colnames(top_acts_mat) %in% c("Androgen-UCell",
                                                      "Estrogen-UCell","CRIC-ARIC-CKD-UCell","ME48-UCell","Mend-Rand-CKD-UCell","ME15-UCell","X10yr-Risk-CKD-UCell")],
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D2", cluster_rows = FALSE) 
dev.off()

# Plot Subset
to.use <- c("TGFb-UCell","SenMayo-UCell","MAPK-UCell","CKD-Progression-UCell","ProtAge-UCell","ProtAge20-UCell",
            "NFkB-UCell","Inf-UCell","Aging-UCell","TNFa-UCell","JAK-STAT-UCell","EGFR-UCell",
            "ECM-UCell","EMT-UCell","IFNG-UCell","PI3K-UCell","WNT-UCell","Plasma-ATI-ARIC-UCell", "Plasma-ATI-BKBC-UCell",
            "Plasma-ATI-CHROME-UCell")
pdf(file='trajectories/FIB-Lineages_Pathway_Gene-set_Scores_subclassl3_subset.pdf',width=7,height=5)
pheatmap(top_acts_mat[,to.use],
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_rows = FALSE) 
dev.off()
