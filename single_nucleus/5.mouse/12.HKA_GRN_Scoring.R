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

###Load Mouse Object
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
table(KB$assay)
###Load GRNs
files <- c("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S16-19_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S9-S12_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_S9-S14_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-S1-2_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-S1-2_l1_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aTAL-frTAL_Trajectory_GRN_0424-newData_TALA.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aTAL-C-TAL-A_Trajectory_GRN_0424-newData_TALA.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-frPTS3_l5_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-PTS3_l4_Trajectory_GRN_0424-newData.RDA")
names <- setNames(c("pvMYOF","CFIBOSMRhi","CMYOF","frPTS1S2","PTS1S2","frTAL","CMTALA","frPTS3","PTS3"),files)

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


#Convert to homologous genes
#Download orthologues based on the blog https://www.ensembl.info/2009/01/21/how-to-get-all-the-orthologous-genes-between-two-species/
hs.ms.table <- read.table("~/Projects/Mouse_Kidney/Atlas_V1/Intermediate_Objects/Human_GRCh38.p14_Mouse_Orthologues_mart_export.txt",sep="\t",header=TRUE)
hs.ms.table <- hs.ms.table[hs.ms.table$Mouse.homology.type %in% "ortholog_one2one",]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Mouse.gene.name),]
hs.ms.table <- hs.ms.table[!duplicated(hs.ms.table$Gene.name),]
hs.ms.table <- hs.ms.table[hs.ms.table$Gene.name %in% unique(unlist(markers)),]

rownames(hs.ms.table) <- hs.ms.table$Gene.name

#convert markers from human to mouse
mouse.markers <- lapply(markers, function(x) {
  print(paste("Running for GRN:", x))
  mouse_genes <- hs.ms.table$Mouse.gene.name[match(x, hs.ms.table$Gene.name)]
  return(mouse_genes[!is.na(mouse_genes)])
  
})



###Add scores to main objects
Idents(KB) <- "library"
celltype <- levels(Idents(KB))
KB[["RNA"]]$counts <- as(KB[["RNA"]]$counts, Class = "dgCMatrix")
KB <- AddModuleScore_UCell(KB, features = mouse.markers)
meta <- KB@meta.data
meta.sub <- meta[,grep("_UCell", colnames(meta))]

saveRDS(meta.sub, file = "mouse_IRI/gene_sets/UCell_HKA_GRN_Scores_11172024.RDS")
meta.sub <- readRDS("mouse_IRI/gene_sets/UCell_HKA_GRN_Scores_11172024.RDS")


###GRNs for MAC trajectories (loosened TF selection criteria)
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

#convert markers from human to mouse
mouse.markers <- lapply(markers, function(x) {
  print(paste("Running for GRN:", x))
  mouse_genes <- hs.ms.table$Mouse.gene.name[match(x, hs.ms.table$Gene.name)]
  return(mouse_genes[!is.na(mouse_genes)])
  
})

###Add scores to main objects
Idents(KB) <- "library"
celltype <- levels(Idents(KB))
KB[["RNA"]]$counts <- as(KB[["RNA"]]$counts, Class = "dgCMatrix")
KB <- AddModuleScore_UCell(KB, features = mouse.markers)
meta <- KB@meta.data
meta.sub <- meta[,grep("_UCell", colnames(meta))]

saveRDS(meta.sub, file = "mouse_IRI/gene_sets/UCell_HKA_GRN_Scores_11172024.RDS")
#meta.sub <- readRDS("mouse_IRI/gene_sets/UCell_HKA_GRN_Scores_11172024.RDS")



###Visualizations
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


###Proximal Tubules
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")

meta.sub <- readRDS("mouse_IRI/gene_sets/UCell_HKA_GRN_Scores_11172024.RDS")
KB <- AddMetaData(object = KB, metadata = meta.sub)

# Extract mlm and store it in pathwaysmlm 
KB[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(KB@meta.data[,grep("PT", colnames(KB@meta.data))])))

# Change assay
DefaultAssay(object = KB) <- "pathwaysmlm"

# Scale the data
KB <- ScaleData(KB)
KB@assays$pathwaysmlm@data <- KB@assays$pathwaysmlm@scale.data

#subset to PT
table(KB$v2.subclass.l2)

#prepare comparison groups
control <- rownames(KB@meta.data[KB@meta.data$v2.subclass.l2 %in% c("PT-S1","PT-S1/2","PT-S2","PT-S3") & KB@meta.data$condition_level3 %in% c("Sham"),])
injured <- rownames(KB@meta.data[KB@meta.data$v2.subclass.l2 %in% c("aPT") & !KB@meta.data$condition_level3 %in% c("Ref","Sham"),])
KB <- subset(KB, cells = c(control, injured))

# Extract activities from object as a long dataframe
Idents(KB) <- "condition_level3"
table(KB$condition_level3)
Idents(KB) <- factor(Idents(KB), levels = c("Sham","4hrs","12hrs",
                                                    "2d","2wk","4-6wk","20-24wk"))
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
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# Plot
subset2 <- c("PTS1S2.HNF4A-UCell","PTS1S2.THRB-UCell","PTS1S2.HNF1A-UCell","PTS1S2.NR2F1-UCell",
             "PTS1S2.MAF-UCell",
             
             "PTS1S2.ATF4-UCell",
             "frPTS1S2.KLF6-UCell",
             
             "PTS1S2.STAT1-UCell",
             "PTS1S2.IRF1-UCell",
             "PTS1S2.SOX9-UCell",
             "PTS1S2.HIF1A-UCell","PTS1S2.HES1-UCell",
             "frPTS1S2.JUN-UCell",
             "frPTS1S2.RELB-UCell","frPTS1S2.NFKB1-UCell",
             "frPTS1S2.ARNT2-UCell",
             "PTS1S2.BHLHE40-UCell",
             "PTS1S2.MITF-UCell","PTS1S2.SOX4-UCell","PTS1S2.MYC-UCell"
)

pdf(file='mouse_IRI/gene_sets/PTS1S2_GRN_Gene-set_Scores_Timecourse_Heatmap_subset2.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[,subset2]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()








###TAL Tubules
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")

meta.sub <- readRDS("mouse_IRI/gene_sets/UCell_HKA_GRN_Scores_11172024.RDS")
KB <- AddMetaData(object = KB, metadata = meta.sub)

# Extract mlm and store it in pathwaysmlm 
KB[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(KB@meta.data[,grep("TAL", colnames(KB@meta.data))])))

# Change assay
DefaultAssay(object = KB) <- "pathwaysmlm"

# Scale the data
KB <- ScaleData(KB)
KB@assays$pathwaysmlm@data <- KB@assays$pathwaysmlm@scale.data

#subset to TAL
table(KB$v2.subclass.l2)

#prepare comparison groups
control <- rownames(KB@meta.data[KB@meta.data$v2.subclass.l3 %in% c("C/M-TAL-A","C-TAL-A") & KB@meta.data$condition_level3 %in% c("Sham"),])
injured <- rownames(KB@meta.data[KB@meta.data$v2.subclass.l2 %in% c("aTAL","frTAL") & !KB@meta.data$condition_level3 %in% c("Ref","Sham"),])
KB <- subset(KB, cells = c(control, injured))

# Extract activities from object as a long dataframe
Idents(KB) <- "condition_level3"
table(KB$condition_level3)
Idents(KB) <- factor(Idents(KB), levels = c("Sham","4hrs","12hrs",
                                            "2d","2wk","4-6wk","20-24wk"))
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
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))


# Plot
subset2 <- c("CMTALA.MITF-UCell", "CMTALA.ESRRB-UCell",
             "CMTALA.RORA-UCell", "CMTALA.NR3C1-UCell", "CMTALA.TEF-UCell",
             
             "CMTALA.BACH1-UCell", "CMTALA.KLF6-UCell", 
             "CMTALA.SOX4-UCell", "CMTALA.SMAD3-UCell",
             "CMTALA.PPARG-UCell",
             "CMTALA.CREM-UCell", "CMTALA.TCF7L1-UCell", 
             "CMTALA.SP1-UCell", 
             "CMTALA.NFKB1-UCell", "CMTALA.SOX9-UCell", 
             
             "CMTALA.ELK4-UCell", "CMTALA.ELF1-UCell", 
             "CMTALA.IRF1-UCell",
             
             "CMTALA.STAT1-UCell", 
             "CMTALA.ELF2-UCell", "CMTALA.CEBPD-UCell", 
             
             "frTAL.EGR1-UCell","frTAL.ELF3-UCell", "frTAL.KLF13-UCell")

pdf(file='mouse_IRI/gene_sets/TALA_GRN_Gene-set_Scores_Timecourse_Heatmap_subset2.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[,subset2]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()







##Interstitial Fibroblasts 
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")

meta.sub <- readRDS("mouse_IRI/gene_sets/UCell_HKA_GRN_Scores_11172024.RDS")
KB <- AddMetaData(object = KB, metadata = meta.sub)

# Extract mlm and store it in pathwaysmlm 
KB[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(KB@meta.data)))

# Change assay
DefaultAssay(object = KB) <- "pathwaysmlm"

# Scale the data
KB <- ScaleData(KB)
KB@assays$pathwaysmlm@data <- KB@assays$pathwaysmlm@scale.data

#Use only Sham and IRI groups
KB <- subset(KB, condition_level3 %in% c("Ref"), invert = TRUE)

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
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

subset <- c("CFIBOSMRhi.ATF3-UCell",
            "CFIBOSMRhi.CREM-UCell",
            "CFIBOSMRhi.RFX2-UCell",
            
            "CFIBOSMRhi.NR2F1-UCell",
            "CMYOF.BHLHE40-UCell",
            "CFIBOSMRhi.NR2F2-UCell",
            "CFIBOSMRhi.HIF1A-UCell",
            
            "CFIBOSMRhi.NFKB1-UCell",
            "CFIBOSMRhi.NFKB2-UCell",
            "CFIBOSMRhi.IRF2-UCell",
            "CFIBOSMRhi.BACH1-UCell",
            "CFIBOSMRhi.IRF1-UCell",
            
            "CFIBOSMRhi.KLF10-UCell",
            "CFIBOSMRhi.SOX4-UCell",
            "CFIBOSMRhi.STAT1-UCell",
            "CFIBOSMRhi.FOSL2-UCell",
            "CFIBOSMRhi.RELB-UCell",
            
            "CMYOF.RUNX2-UCell",
            "CMYOF.REST-UCell",
            "CMYOF.KLF3-UCell",
            "CMYOF.RUNX1-UCell",
            "CMYOF.HOXA10-UCell",
            "CMYOF.PRRX1-UCell",
            "CMYOF.BACH1-UCell",
            "CMYOF.BACH2-UCell",
            "CMYOF.SMAD3-UCell",
            "CMYOF.LEF1-UCell",
            "CMYOF.GLI2-UCell",
            "CMYOF.SOX4-UCell",
            "CMYOF.STAT1-UCell",
            "CMYOF.TBX2-UCell",
            "CMYOF.FOXP1-UCell",
            "CMYOF.FOXK1-UCell",
            "CMYOF.ZEB1-UCell"
)

pdf(file='mouse_IRI/gene_sets/intFIB_GRN_Gene-set_Scores_Timecourse_Heatmap_subset.pdf',width=4,height=6)
pheatmap(t(top_acts_mat[c("C-FIB","C-FIB-PATH","eaC-FIB","C-FIB-Cxcl14+","C-FIB-Cxcl10+","C-MYOF"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()





##Perivascular Fibroblasts 

subset <- c("pvMYOF.GLI2-UCell",
            "pvMYOF.SOX4-UCell",
            "pvMYOF.NFKB1-UCell",
            "pvMYOF.ZBTB7C-UCell",
            "pvMYOF.NR2F1-UCell",
            "pvMYOF.NR4A2-UCell",
            "pvMYOF.ESR1-UCell",
            
            
            "pvMYOF.KLF2-UCell",
            "pvMYOF.KLF3-UCell",
            "pvMYOF.KLF5-UCell",
            "pvMYOF.CREB3L1-UCell",
            "pvMYOF.KLF13-UCell",
            "pvMYOF.KLF10-UCell",
            "pvMYOF.KLF9-UCell",
            
            "pvMYOF.STAT3-UCell",
            "pvMYOF.RUNX1-UCell",
            "pvMYOF.RUNX2-UCell",
            "pvMYOF.NFIL3-UCell",
            "pvMYOF.CEBPB-UCell",
            "pvMYOF.CEBPD-UCell",
            "pvMYOF.NFE2L1-UCell",
            "pvMYOF.FOS-UCell",
            "pvMYOF.JUND-UCell",
            "pvMYOF.JDP2-UCell",
            "pvMYOF.JUN-UCell",
            
            
            
            "pvMYOF.ELF1-UCell",
            "pvMYOF.ERG-UCell",
            "pvMYOF.ETS1-UCell",
            "pvMYOF.MEIS1-UCell",
            "pvMYOF.MEIS2-UCell",
            
            
            "pvMYOF.MEF2A-UCell",
            "pvMYOF.MEF2C-UCell"
)

pdf(file='mouse_IRI/gene_sets/pvFIB_GRN_Gene-set_Scores_Timecourse_Heatmap_subset.pdf',width=4,height=6)
pheatmap(t(top_acts_mat[c("pvFIB-Rspo3+","pvFIB-Pi16+",
                          "pvFIB","pvMYOF"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()


