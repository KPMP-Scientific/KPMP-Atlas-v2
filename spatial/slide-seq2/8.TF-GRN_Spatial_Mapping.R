library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(igraph)
library(BPCells)
library(UCell)

library(tibble)
library(tidyr)
library(pheatmap)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")
source("misc/utils.R")

kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09042024_Final_Object.Rds")

files <- c("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S16-19_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S9-S12_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S9-S14_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-S1-2_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-S1-2_l1_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aTAL-frTAL_Trajectory_GRN_0424-newData_TALA.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aTAL-C-TAL-A_Trajectory_GRN_0424-newData_TALA.RDA")
names <- setNames(c("pvMYOF","CFIBOSMRhi","CMYOF","frPTS1S2","PTS1S2","frTAL","CMTALA"),files)

marker.list <- lapply(files, function(x) {
  print(paste("Running for file:", x))
  
  load(x)
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
  
  genes <- df.cor[df.cor$correlation > 0.6 & 
                    df.cor$between.centrality > 0,]$tfs
  
  df.grn1 <- df.grn %>%
    group_by(tf) %>%
    slice_min(order_by = p_value, n = 200)
  df.grn2 <- df.grn1 %>%
    subset(correlation > 0.6) %>%
    select(c(tf, gene, correlation)) %>%
    rename(weights = correlation)
  
  markers <- lapply(genes, function(g) {
    gene_list = c(df.grn2[df.grn2$tf == g & df.grn2$weights > 0,]$gene, paste0(df.grn2[df.grn2$tf == g & df.grn2$weights < 0,]$gene,"-"))
    return(gene_list)
  })
  names(markers) <- paste0(as.character(names[x]),".",genes)
  return(markers)
  
})


markers <- do.call(c,marker.list)

###Add scores to main objects
table(kss$region_level1)
kss <- subset(kss, region_level1 %in% c("Cortex"))
Idents(kss) <- "library"
celltype <- levels(Idents(kss))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(kss, library %in% ct)
  s.obj[["Spatial"]] <- as(s.obj[["Spatial"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub <- meta[,grep("_UCell", colnames(meta))]

saveRDS(meta.sub, file = "gene-sets/UCell_Trajectory_GRN_Scores_Slide-Seq_09132024.RDS")
meta.sub <- readRDS("gene-sets/UCell_Trajectory_GRN_Scores_Slide-Seq_09132024.RDS")




#Add myeloid GRNs
files <- c("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_MON-moFAM_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_MON-CXCL10_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_resMAC_Trajectory_GRN_0424-newData.RDA")
names <- setNames(c("moFAM","MONCXCL10","resMAC"),files)

marker.list <- lapply(files, function(x) {
  print(paste("Running for file:", x))
  
  load(x)
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
  
  genes <- df.cor[df.cor$correlation > 0.55 & 
                    df.cor$between.centrality > 0,]$tfs
  
  df.grn1 <- df.grn %>%
    group_by(tf) %>%
    slice_min(order_by = p_value, n = 200)
  df.grn2 <- df.grn1 %>%
    subset(correlation > 0.6) %>%
    select(c(tf, gene, correlation)) %>%
    rename(weights = correlation)
  
  markers <- lapply(genes, function(g) {
    gene_list = c(df.grn2[df.grn2$tf == g & df.grn2$weights > 0,]$gene, paste0(df.grn2[df.grn2$tf == g & df.grn2$weights < 0,]$gene,"-"))
    return(gene_list)
  })
  names(markers) <- paste0(as.character(names[x]),".",genes)
  return(markers)
  
})


markers <- do.call(c,marker.list)

###Add scores to main objects
table(kss$region_level1)
kss <- subset(kss, region_level1 %in% c("Cortex"))
Idents(kss) <- "library"
celltype <- levels(Idents(kss))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(kss, library %in% ct)
  s.obj[["Spatial"]] <- as(s.obj[["Spatial"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub2 <- meta[,grep("_UCell", colnames(meta))]
meta.sub <- readRDS("gene-sets/UCell_Trajectory_GRN_Scores_Slide-Seq_09132024.RDS")
meta.sub <- cbind(meta.sub, meta.sub2)

saveRDS(meta.sub, file = "gene-sets/UCell_Trajectory_GRN_Scores_Slide-Seq_10022024.RDS")
meta.sub <- readRDS("gene-sets/UCell_Trajectory_GRN_Scores_Slide-Seq_10022024.RDS")







###Plot scores by subclass
kss <- AddMetaData(kss, metadata = meta.sub)

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

Idents(kss) <- "v2.subclass.sp"
Idents(kss) <- factor(Idents(kss), levels = order)


# Extract mlm and store it in pathwaysmlm in KB.PT
kss[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(kss@meta.data[,grep("_UCell", colnames(kss@meta.data))])))

# Change assay
DefaultAssay(object = kss) <- "pathwaysmlm"

# Scale the data
kss <- ScaleData(kss)
kss@assays$pathwaysmlm@data <- kss@assays$pathwaysmlm@scale.data





##Proximal Tubules
kss.pt <- subset(kss, v2.subclass.sp %in% c("PT-S1/2","PT-S3","aPT","frPT"))

# Extract activities from object as a long dataframe
Idents(kss.pt) <- factor(Idents(kss.pt), levels = c("PT-S1/2","PT-S3","aPT","frPT"))
df <- t(as.matrix(kss.pt@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(kss.pt)) %>%
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

# Plot
subset <- c("PTS1S2.NR2F1-UCell","PTS1S2.HNF1A-UCell",
            "PTS1S2.THRB-UCell","PTS1S2.HNF4A-UCell",
            "frPTS1S2.ATF3-UCell","frPTS1S2.JUN-UCell",
            "frPTS1S2.NFKB1-UCell","frPTS1S2.NFKB2-UCell","frPTS1S2.ARNT2-UCell",
            "frPTS1S2.RELB-UCell","PTS1S2.SOX9-UCell","PTS1S2.SOX4-UCell")

pdf(file='slide-seq/Plots/PTS1S2_GRN_Gene-set_Scores_Heatmap_subset.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[c("PT-S1/2","aPT","frPT"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

pdf("slide-seq/Plots/Puck_220122_38_PTS1S2_NFKB2-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.pt, images = "Puck_220122_38", 
                   features = c("frPTS1S2.NFKB2-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_PTS1S2_SOX9-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.pt, images = "Puck_220122_38", 
                   features = c("PTS1S2.SOX9-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("slide-seq/Plots/Puck_220122_36_PTS1S2_NFKB2-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.pt, images = "Puck_220122_36", 
                   features = c("frPTS1S2.NFKB2-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_PTS1S2_SOX9-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.pt, images = "Puck_220122_36", 
                   features = c("PTS1S2.SOX9-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()






##TAL Tubules
kss.tal <- subset(kss, v2.subclass.sp %in% c("C-TAL","aTAL1","aTAL2","frTAL"))
kss.tal$v2.subclass.sp <- as.character(kss.tal$v2.subclass.sp)
kss.tal$v2.subclass.sp[kss.tal$v2.subclass.sp == "aTAL1"] <- "aTAL"
kss.tal$v2.subclass.sp[kss.tal$v2.subclass.sp == "aTAL2"] <- "aTAL"

# Extract activities from object as a long dataframe
Idents(kss.tal) <- "v2.subclass.sp"
Idents(kss.tal) <- factor(Idents(kss.tal), levels = c("C-TAL","aTAL","frTAL"))
df <- t(as.matrix(kss.tal@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(kss.tal)) %>%
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

# Plot
subset <- c("CMTALA.MITF-UCell", "CMTALA.NR3C1-UCell", "CMTALA.TEF-UCell", 
            "CMTALA.ESRRB-UCell", "CMTALA.RORA-UCell",  
            "frTAL.GLIS3-UCell","frTAL.SOX4-UCell",
            "frTAL.EHF-UCell", "frTAL.ELF3-UCell",
            "frTAL.KLF13-UCell", "frTAL.EGR1-UCell", "frTAL.KLF5-UCell", 
            "frTAL.KLF2-UCell")

pdf(file='slide-seq/Plots/TAL_GRN_Gene-set_Scores_Heatmap_subset.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[c("C-TAL","aTAL","frTAL"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

pdf("slide-seq/Plots/Puck_220122_38_TAL_ELF3-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.tal, images = "Puck_220122_38", 
                   features = c("frTAL.ELF3-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_TAL_SOX4-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.tal, images = "Puck_220122_38", 
                   features = c("frTAL.SOX4-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("slide-seq/Plots/Puck_220122_36_TAL_ELF3-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.tal, images = "Puck_220122_36", 
                   features = c("frTAL.ELF3-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_TAL_SOX4-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.tal, images = "Puck_220122_36", 
                   features = c("frTAL.SOX4-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()





##Interstitial Fibroblasts 
kss.fib <- subset(kss, v2.subclass.sp %in% c("C-FIB","C-FIB-Path","C-FIB-OSMRlo",
                                             "C-FIB-OSMRhi","C-MYOF"))
kss.fib$v2.subclass.sp <- as.character(kss.fib$v2.subclass.sp)
kss.fib$v2.subclass.sp[kss.fib$v2.subclass.sp == "C-FIB-Path"] <- "C-FIB"
kss.fib$v2.subclass.sp[kss.fib$v2.subclass.sp == "C-FIB-OSMRlo"] <- "C-FIB"

# Extract activities from object as a long dataframe
Idents(kss.fib) <- "v2.subclass.sp"
Idents(kss.fib) <- factor(Idents(kss.fib), levels = c("C-FIB","C-FIB-OSMRhi","C-MYOF"))
df <- t(as.matrix(kss.fib@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(kss.fib)) %>%
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

# Plot
subset <- c("CMYOF.ATF3-UCell",
            "CMYOF.CREM-UCell","CMYOF.FOXK1-UCell","CFIBOSMRhi.NR2F2-UCell",
            "CFIBOSMRhi.BACH1-UCell","CFIBOSMRhi.IRF2-UCell","CFIBOSMRhi.NFKB1-UCell",
            "CFIBOSMRhi.STAT1-UCell","CFIBOSMRhi.RELB-UCell","CFIBOSMRhi.IRF1-UCell",
            "CMYOF.LEF1-UCell","CMYOF.SOX4-UCell","CMYOF.STAT1-UCell","CMYOF.GLI2-UCell",
            "CMYOF.RUNX1-UCell","CMYOF.RUNX2-UCell","CMYOF.KLF3-UCell","CMYOF.SMAD3-UCell"
            )

pdf(file='slide-seq/Plots/Int_FIB_OSMRhi_CMYOF_GRN_Gene-set_Scores_Heatmap_subset.pdf',width=4,height=6)
pheatmap(t(top_acts_mat[c("C-FIB","C-FIB-OSMRhi","C-MYOF"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

pdf("slide-seq/Plots/Puck_220122_38_Int_FIB_ATF3-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_38", 
                   features = c("CMYOF.ATF3-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_Int_FIB_IRF2-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_38", 
                   features = c("CFIBOSMRhi.IRF2-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_Int_FIB_RUNX1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_38", 
                   features = c("CMYOF.RUNX1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("slide-seq/Plots/Puck_220122_36_Int_FIB_ATF3-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_36", 
                   features = c("CMYOF.ATF3-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_Int_FIB_IRF2-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_36", 
                   features = c("CFIBOSMRhi.IRF2-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_Int_FIB_RUNX1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_36", 
                   features = c("CMYOF.RUNX1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()





##Perivascular Fibroblasts 
kss.fib <- subset(kss, v2.subclass.sp %in% c("pvFIB-RSPO3+","pvFIB-PI16+",
                                             "pvFIB","pvMYOF"))

# Extract activities from object as a long dataframe
Idents(kss.fib) <- factor(Idents(kss.fib), levels = c("pvFIB-RSPO3+","pvFIB-PI16+",
                                                      "pvFIB","pvMYOF"))
df <- t(as.matrix(kss.fib@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(kss.fib)) %>%
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

# Plot
subset <- c("pvMYOF.KLF9-UCell","pvMYOF.KLF5-UCell","pvMYOF.STAT3-UCell","pvMYOF.GLI2-UCell",
            "pvMYOF.NFE2L1-UCell","pvMYOF.CEBPB-UCell","pvMYOF.CEBPD-UCell",
            "pvMYOF.JUND-UCell","pvMYOF.FOS-UCell","pvMYOF.JDP2-UCell",
            "pvMYOF.RUNX1-UCell","pvMYOF.ELF1-UCell","pvMYOF.ERG-UCell","pvMYOF.ETS1-UCell",
            "pvMYOF.MEF2A-UCell","pvMYOF.MEF2C-UCell")

pdf(file='slide-seq/Plots/pvFIB_GRN_Gene-set_Scores_Heatmap_subset.pdf',width=4,height=6)
pheatmap(t(top_acts_mat[,subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

pdf("slide-seq/Plots/Puck_220122_38_pvFIB_KLF5-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_38", 
                   features = c("pvMYOF.KLF5-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_pvFIB_ELF1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_38", 
                   features = c("pvMYOF.ELF1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_pvFIB_MEF2C-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_38", 
                   features = c("pvMYOF.MEF2C-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("slide-seq/Plots/Puck_220122_36_pvFIB_KLF5-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_36", 
                   features = c("pvMYOF.KLF5-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_pvFIB_ELF1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_36", 
                   features = c("pvMYOF.ELF1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_pvFIB_MEF2C-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.fib, images = "Puck_220122_36", 
                   features = c("pvMYOF.MEF2C-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()





##Myeloid Cells
kss.mye <- subset(kss, v2.subclass.sp %in% c("resMAC-LYVE1+","resMAC-HLAIIhi",
                                             "MON","moMAC-INF","moFAM","moMAC-C3+"))

# Extract activities from object as a long dataframe
Idents(kss.mye) <- factor(Idents(kss.mye), levels = c("resMAC-LYVE1+","resMAC-HLAIIhi",
                                                      "MON","moMAC-INF","moFAM","moMAC-C3+"))
df <- t(as.matrix(kss.mye@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(kss.mye)) %>%
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

# Plot
subset <- c("resMAC.CEBPB-UCell",
            "resMAC.XBP1-UCell",
            "resMAC.HES1-UCell",
            "resMAC.IRF1-UCell",
            "resMAC.KLF13-UCell",
            "resMAC.PBX1-UCell",
            "resMAC.NR2F2-UCell",
            "resMAC.POU3F3-UCell",
            "resMAC.SOX4-UCell",
            "resMAC.TEAD1-UCell",
            
            "MONCXCL10.FOSL2-UCell",
            "MONCXCL10.JUND-UCell",
            "MONCXCL10.CEBPB-UCell",
            "moFAM.ATF3-UCell",
            "moFAM.CREM-UCell",
            "moFAM.HIF1A-UCell",
            "moFAM.JUND-UCell",
            "MONCXCL10.IRF1-UCell",
            "MONCXCL10.NFATC2-UCell",
            "MONCXCL10.IRF8-UCell",
            "MONCXCL10.STAT1-UCell",
            "MONCXCL10.NFKB1-UCell",
            "MONCXCL10.RELB-UCell",
            "moFAM.EGR2-UCell",
            "moFAM.TFEB-UCell",
            "moFAM.MITF-UCell",
            "moFAM.USF2-UCell",
            
            "moFAM.MAF-UCell",
            "moFAM.NFE2L1-UCell")

pdf(file='slide-seq/Plots/Myeloid_GRN_Gene-set_Scores_Heatmap_subset.pdf',width=5,height=6)
pheatmap(t(top_acts_mat[,subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

pdf("slide-seq/Plots/Puck_220122_38_resMAC_TEAD1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.mye, images = "Puck_220122_38", 
                   features = c("resMAC.TEAD1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_MONCXCL10_IRF1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.mye, images = "Puck_220122_38", 
                   features = c("MONCXCL10.IRF1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_38_moFAM_HIF1A-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.mye, images = "Puck_220122_38", 
                   features = c("moFAM.HIF1A-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

pdf("slide-seq/Plots/Puck_220122_36_resMAC_TEAD1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.mye, images = "Puck_220122_36", 
                   features = c("resMAC.TEAD1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_MONCXCL10_IRF1-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.mye, images = "Puck_220122_36", 
                   features = c("MONCXCL10.IRF1-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf("slide-seq/Plots/Puck_220122_36_moFAM_HIF1A-UCell_GRN_SpatialPlot.pdf", width = 6, height = 7)
SpatialFeaturePlot(kss.mye, images = "Puck_220122_36", 
                   features = c("moFAM.HIF1A-UCell"), ncol = 1, 
                   alpha = c(0.4, 1), stroke = NA)  + DarkTheme() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
