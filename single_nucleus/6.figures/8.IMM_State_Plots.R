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
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I16-22_slingshot_0424-newData.rda")

v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("I_14","I_15","I_16","I_17",
                         "I_18","I_19","I_20","I_21","I_22"))

pdf(file='trajectories/moMAC-Lineages_Slingshot_subclassl3_colors.pdf',width=6,height=5)
plot(reducedDims(sce)$UMAP, col = adjustcolor(v2.fn.cols[sce$v2.clusters], alpha.f = 0.5), pch = 19, cex = 0.4)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()


###UMAP Plots - color by subclass.l3
load("~/hsKidAt/blake_LTS/Atlas_V2/Intermediate_Objects/Kidney_AtlasV2_Seurat_Immune_subset_filtered_Myeloid_I14-15_slingshot_0424-newData.rda")

pdf(file='trajectories/resMAC-Lineages_Slingshot_subclassl3_colors.pdf',width=6,height=5)
plot(reducedDims(sce)$UMAP, col = adjustcolor(v2.fn.cols[sce$v2.clusters], alpha.f = 0.5), pch = 19, cex = 0.4)
lines(SlingshotDataSet(sce), lwd = 3, col = "black", type = 'curves')
dev.off()



##Clinical/Pathological Heatmaps (Subclasses)
#Read in sample and pvalue tables
cond.ref.aki.pVal <- readRDS("Condition_Ref-AKI_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.ref.ckd.pVal <- readRDS("Condition_Ref-CKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.aki.ckd.pVal <- readRDS("Condition_AKI-CKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.ref.ati.pVal <- readRDS("Condition_Ref-AKI-ATI_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.ref.ain.pVal <- readRDS("Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.ati.ain.pVal <- readRDS("Condition_Ref-ATI-AIN_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.ref.ckd.hi.pVal <- readRDS("Condition_Ref-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
cond.ckd.lo.ckd.hi.pVal <- readRDS("Condition_CKD-LowRisk-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
adj.cond.HCKD.DKD.pVal <- readRDS("Adj-Condition_H-CKD-DKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")

eGFR.pVal <- readRDS("eGFR_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
age.pVal <- readRDS("HRT_Age_Group_composition_V2_Subclassl3_P-Values.RDS")
IF.pVal <- readRDS("Interstitial-Fibrosis_Group_composition_V2_Subclassl3_P-Values.RDS")
TA.pVal <- readRDS("Tubular-Atrophy_Group_composition_V2_Subclassl3_P-Values.RDS")
TI.pVal <- readRDS("Tubular-Injury_Group_composition_V2_Subclassl3_P-Values.RDS")
ALM.pVal <- readRDS("Abnormal_Luminal_Morphology_Group_composition_V2_Subclassl3_P-Values.RDS")
AC.pVal <- readRDS("Acellular_Cast_Group_composition_V2_Subclassl3_P-Values.RDS")
WBC.pVal <- readRDS("Interstitial_WBC_Group_composition_V2_Subclassl3_P-Values.RDS")

#P value plot
row.order <- c("B","PL","NaïveTh","MAIT","ILC3","TREG","CD8+TEM/TRM","CD8+TEM/TEMRA",
               "NK","MAST","resMACLYVE1+","resMACHLAIIhi","MON","moMACHBEGF+",
               "moMACCXCL10+","moFAM","moMACC3+","ncMON","cDC2","cDC1","mDC","pDC","N")
p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]

#Clinical Plot
col.order <- c("cond.ref.aki.pVal","cond.ref.ati.pVal","cond.ref.ain.pVal",
               "cond.aki.ckd.pVal",
               "cond.ref.ckd.pVal","cond.ref.ckd.hi.pVal","cond.ckd.lo.ckd.hi.pVal",
               "adj.cond.HCKD.DKD.pVal",
               "eGFR.pVal","age.pVal")

p.stats <- p.table[row.order,paste0(col.order,".p")]

pdf(file='trajectories/IMM-Lineages_Clinical-Path_Heatmap_subclasses_clin.pdf',width=8,height=5)
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(5, 8), Rowv = NA, Colv = NA,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 


#Pathology Plot
col.order <- c("IF.pVal","TA.pVal","TI.pVal","WBC.pVal")

pdf(file='trajectories/IMM-Lineages_Clinical-Path_Heatmap_subclasses_path.pdf',width=8,height=5)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(14, 8), Rowv = NA, Colv = NA,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 



##Marker Gene Plots
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.IMM <- subset(KB, v2.subclass.l3 %in% c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                           "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                           "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))
Idents(KB.IMM) <- "v2.subclass.l3"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                                    "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                                    "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))
KB.IMM <- subset(KB.IMM, downsample = 1000)
KB.IMM[["RNA"]] <- as(KB.IMM[["RNA"]], Class = "Assay")
KB.IMM <- NormalizeData(KB.IMM)

imm.markers <- c("PTPRC",                                             #Broad Immune
                 "BANK1","MS4A1","CD37","CD79A",                      #B Cells
                 "IGKC","XBP1","MZB1","JCHAIN",#"SDC1",               #PL Cells
                 "CD96","CD247","BCL11B","THEMIS",                    #T
                 "INPP4B","TRAC","CD3D",                              #T
                 "IL7R",                                              #T/MAIT/ILC3
                 "LEF1","CD4","SELL",                                 #Naïve Th
                 "SLC4A10","KLRB1","CCR6",                            #MAIT
                 "PCDH9","TOX2","KIT","RORC",                         #ILC3
                 "IKZF2","RTKN2","IL2RA","CTLA4",#"FOXP3",            #T-REG
                 "CD8A",                                              #CD8+
                 "GZMK",                                              #CD8+ TEM/TRM
                 "CCL5","SAMD3","GZMA","CCL4",                        #CD8+ & NK 
                 "NKG7","KLRD1","GNLY","GZMB","CX3CR1",               #CD8+ TEM/TEMRA & NK
                 "GZMH",                                              #CD8+ TEM/TEMRA
                 "TXK","KLRF1","NCAM1",#"PDGFD",                      #NK
                 
                 
                 "HBB","HBA2","HBA1",                                #Erythrocyte
                 "CPA3","IL18R1","TPSB2","TPSAB1",                   #MAST
                 "CD163","MSR1","CSF1R","CD14",                      #MAC
                 "MRC1","F13A1","STAB1","CD163L1","LYVE1",           #resMAC-LYVE1+
                 "HLA-DPA1","C1QA","C1QB",                           #resMAC-HLAIIhi
                 "HIF1A","NAMPT","PLAUR","ITGAX","HBEGF","OSM",      #moMAC-HBEGF+
                 "PSTPIP2","CXCL10","CXCL9","CCL2","CCL3","IL1B",    #moMAC-CXCL10+
                 "GPNMB","SPP1","APOC1","PLA2G7","CD68","CAPG",      #moFAM
                 "HMOX1","TREM2",                                    #moFAM
                 "C3","KCNQ3","ADGRB3","VASH1","CX3CR1",             #moMAC-C3+
                 "CLEC10A","FCER1A","CD1C",                          #cDC2
                 "TCF7L2","COTL1","FCGR3A",                          #ncMON
                 "FCN1",                                             #MON/ncMON
                 "VCAN","LYZ","CD36",                                #MON
                 "LAMP3","SLCO5A1","CCR7",#"EBI3","CCL19",           #mDC
                 "WDFY4","CADM1","CLEC9A","BATF3",                   #cDC1
                 "BCL11A","CLEC4C","IL3RA","PLD4",#"LILRA4",         #pDC
                 "S100A9","FCGR3B","S100A8","IFITM2",                #N
                 "TOP2A","MKI67",                                    #cycling
                 "NRXN1","GRIK2","CDH19","NCAM2"                     #SC/NEU
                 
)

DotPlot(KB.IMM, features = unique(imm.markers)) + RotatedAxis()


imm.markers <- c("BANK1",                      #B Cells
                 "XBP1",              #PL Cells
                 "CD3D",                              #T
                 "IL7R",                                              #T/MAIT/ILC3
                 "SLC4A10",                            #MAIT
                 "KIT",                         #ILC3
                 "IKZF2",            #T-REG
                 "CD8A",                                              #CD8+
                 "GZMK",                                              #CD8+ TEM/TRM
                 "CCL5",                        #CD8+ & NK 
                 "NKG7",               #CD8+ TEM/TEMRA & NK
                 "KLRF1",                     #NK
                 
                 
                 "CPA3",                   #MAST
                 "MRC1","LYVE1",           #resMAC-LYVE1+
                 "C1QA",                           #resMAC-HLAIIhi
                 "FCN1",                                             #MON/ncMON
                 "HBEGF",      #moMAC-HBEGF+
                 "CXCL10",    #moMAC-CXCL10+
                 "GPNMB",      #moFAM
                 "C3",             #moMAC-C3+
                 "TCF7L2","FCGR3A",                          #ncMON
                 "CLEC10A",                          #cDC2
                 "WDFY4",                   #cDC1
                 "CCR7",           #mDC
                 "IL3RA",         #pDC
                 "FCGR3B"                #N
                 
)

pdf(file='trajectories/IMM-Lineages_Marker_Dotplot_subclassl3_subset.pdf',width=8.5,height=4.5)
DotPlot(KB.IMM, features = unique(imm.markers)) + RotatedAxis()
dev.off()



##UMAP visualization of main immune cell types
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.IMM <- subset(KB, v2.subclass.l3 %in% c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                           "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                           "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))
KB.IMM[["RNA"]] <- as(KB.IMM[["RNA"]], Class = "Assay")

library(pagoda2)
countMatrix <- KB.IMM[["RNA"]]$counts
countMatrix.batch <- as.factor(unlist(lapply(colnames(countMatrix),function(x) unlist(strsplit(x,"_"))[1])))
names(countMatrix.batch) <- colnames(countMatrix)
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10, batch=countMatrix.batch)
p2$adjustVariance(plot = T, gam.k = 10)
#2062 overdispersed genes
p2$calculatePcaReduction(nPcs = 50, n.odgenes = 2000, maxit=1000)

#pagoda2 clusters and PCA values for umap
var.info <- p2$misc$varinfo; var.info <- var.info[order(var.info$lp, decreasing = F),];
sn.od.genes <- rownames(var.info)[1:2000]
VariableFeatures(KB.IMM) <- sn.od.genes

cell.embeddings <- p2$reductions$PCA
KB.IMM[["pca"]] <- CreateDimReducObject(embeddings = cell.embeddings, key = "pca_", assay = DefaultAssay(KB.IMM))
KB.IMM <- RunUMAP(object = KB.IMM, reduction = "pca", dims = 1:30, n.neighbors = 20L, spread = 3,
                  min.dist = 0.01)

load("color_factors_v2-clusters.robj")


v2.fn.cols <- setNames(c("#005300","#1F9698","#00FF00","#C8FF00",
                         "#FFD300","#DD00FF","#B1CC71","#A10300","#FE8F42"),
                       c("resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+",
                         "moFAM","moMAC-C3+","cDC2","ncMON","MON"))
v2.fn.cols <- c(v2.fn.cols, v2.scl3.cols[!names(v2.scl3.cols) %in% names(v2.fn.cols)])


Idents(KB.IMM) <- "v2.subclass.l3"
pdf(file='trajectories/IMM-Lineages_subclassl3_subset_UMAP.pdf',width=7,height=6)
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.6, label = TRUE, alpha = 0.1,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.IMM))],
) + NoLegend()
dev.off()
pdf(file='trajectories/IMM-Lineages_subclassl3_subset_UMAP_Unlabeled.pdf',width=7,height=6)
DimPlot(KB.IMM, reduction = "umap", pt.size = 0.6, label = FALSE, alpha = 0.1,
        repel = TRUE, cols = v2.fn.cols[levels(Idents(KB.IMM))],
) + NoLegend()
dev.off()



##Myeloid Trajectory markers
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB.mye <- subset(KB, v2.subclass.l3 %in% c("resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                           "moMAC-CXCL10+","moFAM","moMAC-C3+"))
Idents(KB.mye) <- "v2.subclass.l3"
Idents(KB.mye) <- factor(Idents(KB.mye), levels = c("resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                                    "moMAC-CXCL10+","moFAM","moMAC-C3+"))
KB.mye <- subset(KB.mye, downsample = 1000)
KB.mye[["RNA"]] <- as(KB.mye[["RNA"]], Class = "Assay")
KB.mye <- NormalizeData(KB.mye)

mye.markers <- c(
  "PDGFC","MERTK","IGF1","PDGFB", "TGFB1", #Reparative
  "HLA-DRA","C1QA","C1QB","CD81","HLA-DQA1","TMEM176B","TMEM176A", #Tissue resident MAC - Zimmerman et al. 2019 https://doi.org/10.1681/ASN.2018090931; Fu et al 2019 https://doi.org/10.1681/ASN.2018090896
  
  "FCN1","CD300E","CFP", #FCN1+ macrophages from Zhang et al., 2021
  "PLAUR",
  "IL1B","CCL3","AREG","CXCL8","OSM", #Inflammatory monocytes from Smillie et al. 2019; #Leukocyte recruitment (CXCL8) and activation (TNF, IL1B)
  #"AREG",#reparative
  "STAT1","GBP1","CXCL10","GBP5","NFKBIA","IRF1","CXCL9","FCGR1A","CCL2","CXCL11",#CXCL10+CCR2+ inflammatory macrophages from Zhang et al., 2021; ##INF-gamma and TNF stimulated genes from Zhang et al., 2021
  "TNF",#Leukocyte recruitment (CXCL8) and activation (TNF, IL1B)
  
  
  "GPNMB","APOE","SPP1","PPARG","LIPA","CD63","PLA2G7","LGALS3","CD9","CAPG","NR1H3","CD68","FABP5","CHIT1","CHI3L1","TREM2", #LAM markers Eraslan et al, SAM markers from Rodríguez-Morales et al., 2023; SPP1+ mac from Ouyang et al., 2023;SPP1hi MERTK+ From Morse et al., 2019 
  
  "MMP2","MMP14","MMP9","TREM2","AXL"
)

pdf(file='trajectories/Myeloid_Trajectories_Marker_Dotplot.pdf',width=12,height=3)
DotPlot(KB.mye, features = unique(mye.markers)) + RotatedAxis()
dev.off()


                                             
##Plot TF Activities
library(Signac)
load("~/hsKidAt/blake_LTS/Atlas_V2/multiome/Kidney_Immune_Dual_Seurat_10012024.rda")


#TFBS activities - All PT subclass.l3
DefaultAssay(KB.IMM) <- "chromvar"
Idents(KB.IMM) <- "v2.subclass.l3"
KB.IMM <- subset(KB.IMM, v2.subclass.l3 %in% c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                               "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                               "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))

Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                                    "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                                    "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))
table(Idents(KB.IMM))
rownames(KB.IMM@assays$chromvar@data) <- ConvertMotifID(KB.IMM,
                                                        id=rownames(KB.IMM@assays$chromvar@data),
                                                        assay='ATAC')
rownames(KB.IMM@assays$chromvar@data) <- toupper(rownames(KB.IMM@assays$chromvar@data))
TFs.plot <- c("POU2F2","TCF4","BCL11B","RUNX1","EOMES",
              "GATA2","MAFB","SPI1","CEBPB")

pdf(file='trajectories/IMM-Lineages_TFBS_Dotplot_subclassl3.pdf',width=6,height=7)
DotPlot(KB.IMM, features = TFs.plot, dot.scale = 6, col.min = 0, cols = c("lightgrey","darkred")) + RotatedAxis()
dev.off()





##Pathway scores
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_09122024.RDS")
KB <- AddMetaData(KB, meta.sub)

# Extract mlm and store it in pathwaysmlm in KB.FIB
KB[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(meta.sub[colnames(KB),])))

# Change assay
DefaultAssay(object = KB) <- "pathwaysmlm"

# Scale the data
KB <- ScaleData(KB)
KB@assays$pathwaysmlm@data <- KB@assays$pathwaysmlm@scale.data

KB.IMM <- subset(KB, v2.subclass.l3 %in% c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                           "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                           "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))

Idents(KB.IMM) <- "v2.subclass.l3"
Idents(KB.IMM) <- factor(Idents(KB.IMM), levels = c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA",
                                                    "NK","MAST","resMAC-LYVE1+","resMAC-HLAIIhi","MON","moMAC-HBEGF+",
                                                    "moMAC-CXCL10+","moFAM","moMAC-C3+","ncMON","cDC2","cDC1","mDC","pDC","N"))


# Extract activities from object as a long dataframe
df <- t(as.matrix(KB.IMM@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.IMM)) %>%
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
pdf(file='trajectories/IMM-Lineages_Pathway_Gene-set_Scores_subclassl3.pdf',width=7,height=5)
pheatmap(top_acts_mat[,!colnames(top_acts_mat) %in% c("Androgen-UCell",
                                                      "Estrogen-UCell","CRIC-ARIC-CKD-UCell","ME48-UCell","Mend-Rand-CKD-UCell","ME15-UCell","X10yr-Risk-CKD-UCell")],
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D2", cluster_rows = FALSE) 
dev.off()

# Plot Subset
to.use <- c("IFNG-UCell",
            "NFkB-UCell",
            "TNFa-UCell",
            "Aging-UCell",
            "JAK-STAT-UCell",
            "Plasma-ATI-CHROME-UCell",
            "SenMayo-UCell",
            "Plasma-ATI-ARIC-UCell",
            "Plasma-ATI-BKBC-UCell",
            "VEGF-UCell",
            "CKD-Progression-UCell",
            "MAPK-UCell",
            "TGFb-UCell",
            "Hypoxia-UCell",
            "Trail-UCell",
            "WNT-UCell",
            "EGFR-UCell",
            "p53-UCell")
pdf(file='trajectories/IMM-Lineages_Pathway_Gene-set_Scores_subclassl3_subset.pdf',width=7,height=7)
pheatmap(top_acts_mat[,to.use],
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_rows = FALSE) 
dev.off()

