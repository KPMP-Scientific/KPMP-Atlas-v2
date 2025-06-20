#####TF-Gene Weight Tables 
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(Matrix)
library(BPCells)
library(UCell)
library(tibble)
library(tidyr)
library(igraph)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

###aPT2-PTS1 trajectory GRN scores
load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-S1-2_l1_Trajectory_GRN_0424-newData.RDA")

netobj <- graph_from_data_frame(df.grn2,directed = TRUE)
V(netobj)$type <- ifelse(V(netobj)$name %in% df.grn2$tf,"TF/Gene","Gene")

# Calculate the outdegree of each node
outdegree <- igraph::degree(netobj, mode = c("out"))
centrality <- igraph::betweenness(netobj, directed = TRUE)
indegree <- igraph::degree(netobj, mode = c("in"))
eigen.centrality <- eigen_centrality(netobj)$vector

df.cor$outdegree <- outdegree[match(df.cor$tfs, names(outdegree))]
df.cor$indegree <- indegree[match(df.cor$tfs, names(indegree))]
df.cor$between.centrality <- centrality[match(df.cor$tfs, names(centrality))]
df.cor$eigen.centrality <- eigen.centrality[match(df.cor$tfs, names(eigen.centrality))]

tfs <- df.cor[df.cor$correlation > 0.6 & 
                  df.cor$between.centrality > 0,]$tfs

df.grn1 <- df.grn %>%
  group_by(tf) %>%
  slice_min(order_by = p_value, n = 500)
df.grn2 <- df.grn1 %>%
  select(c(tf, gene, correlation)) %>%
  rename(weights = correlation)

colnames(df.grn2) <- c("source","target","weight")
df.grn2 <- df.grn2 %>%
  filter(source %in% tfs)
write.table(df.grn2, file = "trajectories/aPT2-PT-S1_TF-Gene_Weight_Table.txt", sep = "\t", quote = FALSE, row.names = FALSE)







#####AKI Progression Secreted DEGs---------------------------------------------------------
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(BPCells)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
KB

###All Subclass.l3 Marker Genes
Idents(object = KB) <- "v2.subclass.l3"
order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)
table(Idents(object = KB))


###Progression DEGs (early, mid, healthy/repaired and failed) and overlap with secreted protein list

#Secreted proteins from human protein atlas
s.prot <- read.delim("sa_location_Secreted.tsv")
table(s.prot$Secretome.location)

#subset to cell types of interest
KB.sub <- subset(KB, idents =  c("PT-S1","PT-S2","PT-S3","aPT2",
                                 "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3",
                                 "C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B",
                                 "aTAL1","aTAL2","frTAL",
                                 "C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF",
                                 "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF",
                                 "resMAC-LYVE1+","resMAC-HLAIIhi",
                                 "moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
                                 "MON"))

KB.sub <- subset(KB.sub, downsample = 2000)
KB.sub[["RNA"]] <- as(KB.sub[["RNA"]], Class = "Assay")
KB.sub <- NormalizeData(KB.sub, assay = "RNA")

KB.sub@meta.data[KB.sub@meta.data$v2.subclass.l3 %in% c("pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"),]$v2.subclass.l1 <- "pvFIB"
KB.sub@meta.data[KB.sub@meta.data$v2.subclass.l3 %in% c("resMAC-LYVE1+","resMAC-HLAIIhi"),]$v2.subclass.l1 <- "resMAC"
KB.sub@meta.data[KB.sub@meta.data$v2.subclass.l3 %in% c("moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
                                                        "MON"),]$v2.subclass.l1 <- "moMAC"

clusters <- unique(KB.sub$v2.subclass.l1)

#all genes
for (i in clusters) {
  KB.sub.ct <- subset(KB.sub, v2.subclass.l1 %in% i)
  ct.markers <- FindAllMarkers(KB.sub.ct, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25)
  write.table(ct.markers, file = paste0("gene-sets/Human_Kidney_AtlasV2_",i,"_l3_trajectory_markers.txt"), sep = "\t", row.names=TRUE, col.names=TRUE)
}

#secreted proteins
for (i in clusters) {
  KB.sub.ct <- subset(KB.sub, v2.subclass.l1 %in% i)
  sec <- s.prot$Gene[s.prot$Gene %in% rownames(KB.sub.ct)]
  ct.markers <- FindAllMarkers(KB.sub.ct, features = sec, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25)
  write.table(ct.markers, file = paste0("gene-sets/Human_Kidney_AtlasV2_",i,"_l3_trajectory_secreted-protein_markers.txt"), sep = "\t", row.names=TRUE, col.names=TRUE)
}






#####AKI Progression Secreted Markers (from TRIBE and BKBC patient cohorts): Expression plots
library(Seurat)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
Idents(object = KB) <- "v2.subclass.l3"
order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT-S3",
           "aPT-S1/S2","aPT1","aPT2","frPT-S1/S2","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
           "frTAL","aTAL1","aTAL2","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)
table(Idents(object = KB))

KB.sub <- subset(KB, downsample = 1000)
KB.sub[["RNA"]] <- as(KB.sub[["RNA"]], Class = "Assay")
KB.sub <- NormalizeData(KB.sub, assay = "RNA")

genes <- c("COL15A1","COL6A3","GDF15","HDGF","TNFRSF1A","COL18A1","NBL1","FSTL3",
           "B4GALT1","TNFRSF1B","PTGDS","CFD","IL1R1","IDS","PTPRS","FBLN1","IGFBP7",
           "B2M","CST3","CD55","GM2A","ST8SIA4","CMPK1","SELENOM","TGFBR3","PXDN")

#Use average scaled expression values
KB <- NormalizeData(KB)
KB <- ScaleData(KB, features = genes)
agg_exp <- AverageExpression(object = KB, features = genes, slot = "scale.data")$RNA
scaled_matrix <- as.matrix(agg_exp)

# Create heatmap using pheatmap
library(pheatmap)
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

pdf(file='gene-sets/AKI_progression_genes_expression_heatmap.pdf',width=16,height=5)
pheatmap(
  scaled_matrix,
  border_color = NA, color=my_color, breaks = my_breaks,
  clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = TRUE
)
dev.off()

##Cortical and cortico-medullary junction subtypes
CM.sub <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT-S3",
            "aPT-S1/S2","aPT1","aPT2","frPT-S1/S2","frPT-S3","cycPT","dPT",
            "C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
            "frTAL","aTAL1","aTAL2","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
            "CNT-PC","dCNT-PC","CCD-PC","CCD-IC-A","dCCD-IC-A",
            "IC-B","C/M-FIB","C-FIB","C-FIB-PATH","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","dFIB",
            "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF",
            "B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK",
            "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
            "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC")


# Plot
pdf(file='gene-sets/AKI_progression_genes_expression_heatmap_CM.pdf',width=13,height=4)
pheatmap(scaled_matrix[,CM.sub],
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "complete", cluster_cols = FALSE, cluster_rows = TRUE) 
dev.off()


#Plot PT subtypes
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
Idents(object = KB) <- "v2.subclass.l3"

KB.PT <- subset(KB, idents = c("aPT2","aPT1","aPT-S1/S2","PT-S1","PT-S2","frPT-S1/S2",
                               "aPT-S3","PT-S3","frPT-S3"))
KB.PT[["RNA"]] <- as(KB.PT[["RNA"]], Class = "Assay")
KB.PT <- NormalizeData(KB.PT)
Idents(KB.PT) <- factor(Idents(KB.PT), levels = c("aPT2","aPT1","aPT-S1/S2","PT-S1","PT-S2","frPT-S1/S2",
                                                  "aPT-S3","PT-S3","frPT-S3"))


pdf(file='gene-sets/AKI_progression_genes_Violin_CM_sub_B.pdf',width=6,height=6)
VlnPlot(KB.PT, features = c("SOX4", "PXDN"), ncol = 1, 
        alpha = 0.01)
dev.off()





###Compare PXDN/SOX4 between AKI Groups
load("color_factors_v2-clusters.robj")
condl1.cols
KB.PT <- subset(KB.PT, condition_level1 %in% c("AKI"))

recovered = c('30-10044','30-10125','30-10868','30-10929','30-11084','32-10003','32-10074','32-10333','33-10331','33-10376','34-10050','34-10187','34-10240','34-10393')
progressed = c('30-10018','30-10034','30-10631','30-11051','30-11080','30-11081','32-10205','32-2','34-10184','34-10209','34-10331','34-10579')

KB.PT <- subset(KB.PT, patient %in% c(recovered,progressed))
KB.PT$AKI_Group <- "Recovered"
KB.PT$AKI_Group[KB.PT$patient %in% progressed] <- "Progressed"
KB.PT$AKI_Group <- factor(KB.PT$AKI_Group, levels = c("Recovered","Progressed"))

pdf(file='gene-sets/AKI_progression_genes_Violin_CM_sub_B_AKI-group_Split.pdf',width=6,height=6)
P1 <- VlnPlot(KB.PT, features = c("SOX4"), ncol = 1, split.by = "AKI_Group",
              split.plot = TRUE ,alpha = 0.1, cols = c("#3b4cc0","#b40426"))
P2 <- VlnPlot(KB.PT, features = c("PXDN"), ncol = 1, split.by = "AKI_Group",
              split.plot = TRUE ,alpha = 0.1, cols = c("#3b4cc0","#b40426"))
P1 / P2
dev.off()





#####Composition Analyses - AKI Recovery and Progression Sub-Groups by subclass.l3
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
Idents(object = KB) <- "v2.subclass.l3"
order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT-S3",
           "aPT-S1/S2","aPT1","aPT2","frPT-S1/S2","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
           "frTAL","aTAL1","aTAL2","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)
table(Idents(object = KB))

recovered = c('30-10044','30-10125','30-10868','30-10929','30-11084','32-10003','32-10074','32-10333','33-10331','33-10376','34-10050','34-10187','34-10240','34-10393')
progressed = c('30-10018','30-10034','30-10631','30-11051','30-11080','30-11081','32-10205','32-2','34-10184','34-10209','34-10331','34-10579')

KB <- subset(KB, patient %in% c(recovered, progressed))
KB$aki_group <- "recovered"
KB$aki_group[KB$patient %in% progressed] <- "progressed"


df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$aki_group))[select.donors, ]

df_comp_relative$cond <- "recovered"
df_comp_relative$cond[df_cond$progressed != 0] <- "progressed"

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("recovered", "progressed"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub("/","",colnames(df_comp_relative))


###Statistical significance of enrichment
#T-Tests
x <- "recovered"
y <- "progressed"
clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_AKI-Recovery-Progression_Group_composition_V2_Subclassl3_P-Values.RDS")


###P value plots
library(gplots)
library(dendextend)
library("RColorBrewer")
row.order <- gsub("-","",order)
row.order <- gsub(" ","",row.order)
row.order <- gsub("/","",row.order)

cond.ref.ati.pVal <- readRDS("Condition_Ref-AKI-ATI_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.ati.ain.pVal <- readRDS("Condition_Ref-ATI-AIN_Group_composition_V2_Subclassl3_P-Values.RDS")
cond.aki.recovery.pVal <- readRDS("Condition_AKI-Recovery-Progression_Group_composition_V2_Subclassl3_P-Values.RDS")

p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]

col.order <- c("cond.ref.ati.pVal","cond.ati.ain.pVal","cond.aki.recovery.pVal")


pdf(file='trajectories/AKI-Progression_Group_subclass_Enrichments_Heatmap.pdf',width=24,height=5)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(5, 8), Rowv = NA, Colv = NA,
          #ColSideColors = ColSideColors,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT-S3",
            "aPT-S1/S2","aPT1","aPT2","frPT-S1/S2","frPT-S3","cycPT","dPT",
            "C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
            "frTAL","aTAL1","aTAL2","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
            "CNT-PC","dCNT-PC","CCD-PC","CCD-IC-A","dCCD-IC-A",
            "IC-B","C/M-FIB","C-FIB","C-FIB-PATH","C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","dFIB",
            "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF",
            "B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK",
            "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
            "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC")
row.order <- gsub("-","",order)
row.order <- gsub(" ","",row.order)
row.order <- gsub("/","",row.order)


pdf(file='trajectories/AKI-Progression_Group_subclass_Enrichments_Heatmap_sub.pdf',width=24,height=5)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(5, 8), Rowv = NA, Colv = NA,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

pdf(file='trajectories/AKI-Progression_Group_subclass_Enrichments_Heatmap_sub_B.pdf',width=24,height=3)
p.stats <- p.table[row.order,paste0(col.order,".p")]
heatmap.2(as.matrix(t(t.stats[row.order,paste0(col.order,".t")])),col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(5, 8), Rowv = NA, Colv = NA,
          cellnote = ifelse(as.matrix(t(p.stats)) < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

###
