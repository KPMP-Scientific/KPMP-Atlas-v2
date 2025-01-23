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

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")

####Add gene set scores to full snRNA object
###Cytokine pathways from Jiang et al., BioRxIV
markers <- read.delim("~/hsKidAt/blake_LTS/Atlas_V2/public/JiangBioRxIV/JiangBioRxIV_ST3_B.txt")

###Gene sets for signaling pathways from PROGENy
library(progeny)

model <- progeny::model_human_full
model <- model[model$gene %in% rownames(KB),]
model <- model %>%
  group_by(pathway) %>%
  slice_min(order_by = p.value, n = 100)

table(model$pathway)

###Msigdbr
library(msigdbr)
msigdbr_genesets <- msigdbr(species = "human")
msigdbr_genesets$gs_id

#HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
EMT_geneset <- msigdbr_genesets[msigdbr_genesets$gs_id == "M5930",]$gene_symbol

#GOBP_EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION & REACTOME_COLLAGEN_FORMATION
ECM_geneset <- unique(msigdbr_genesets[msigdbr_genesets$gs_id %in% c("M40442","M631"),]$gene_symbol)

#HALLMARK_INFLAMMATORY_RESPONSE; KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION
Inf_geneset <- unique(msigdbr_genesets[msigdbr_genesets$gs_id %in% c("M5932","M9809"),]$gene_symbol)

#DEMAGALHAES_AGING_UP
age_geneset <- msigdbr_genesets[msigdbr_genesets$gs_id %in% "M2144",]$gene_symbol

SenMayo <- read.delim("~/hsKidAt/blake_LTS/Atlas_V2/public/Saul2022/Saul_2022_SenMayo_human.txt")


###Combined gene sets
markers <- list(
  "Androgen" = c(as.character(model[model$pathway == "Androgen" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "Androgen" & model$weight < 0, ]$gene), "-")),
  "EGFR" = c(as.character(model[model$pathway == "EGFR" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "EGFR" & model$weight < 0, ]$gene), "-")),
  "Estrogen" = c(as.character(model[model$pathway == "Estrogen" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "Estrogen" & model$weight < 0, ]$gene), "-")),
  "Hypoxia" = c(as.character(model[model$pathway == "Hypoxia" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "Hypoxia" & model$weight < 0, ]$gene), "-")),
  "JAK-STAT" = c(as.character(model[model$pathway == "JAK-STAT" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "JAK-STAT" & model$weight < 0, ]$gene), "-")), 
  "MAPK" = c(as.character(model[model$pathway == "MAPK" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "MAPK" & model$weight < 0, ]$gene), "-")),
  "NFkB" = c(as.character(model[model$pathway == "NFkB" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "NFkB" & model$weight < 0, ]$gene), "-")),
  "p53" = c(as.character(model[model$pathway == "p53" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "p53" & model$weight < 0, ]$gene), "-")),
  "PI3K" = c(as.character(model[model$pathway == "PI3K" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "PI3K" & model$weight < 0, ]$gene), "-")),
  "TGFb" = c(as.character(model[model$pathway == "TGFb" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "TGFb" & model$weight < 0, ]$gene), "-")),
  "TNFa" = c(as.character(model[model$pathway == "TNFa" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "TNFa" & model$weight < 0, ]$gene), "-")),
  "WNT" = c(as.character(model[model$pathway == "WNT" & model$weight > 0,  ]$gene), paste0(as.character(model[model$pathway == "WNT" & model$weight < 0,  ]$gene), "-")),
  "Trail" = c(as.character(model[model$pathway == "Trail" & model$weight > 0, ]$gene), paste0( as.character(model[model$pathway == "Trail" & model$weight < 0, ]$gene), "-")), 
  "VEGF" = c(as.character(model[model$pathway == "VEGF" & model$weight > 0, ]$gene), paste0(as.character(model[model$pathway == "VEGF" & model$weight < 0, ]$gene), "-")),
  "IFNB" = markers$IFNB_program1,
  "IFNG" = markers$IFNG_program1,
  "TNFA" = markers$TNFA_program1,
  "TGFB1" = markers$TGFB1_program1,
  "EMT" = EMT_geneset,
  "ECM" = ECM_geneset,
  "Inf" = Inf_geneset,
  "SenMayo" = SenMayo$Gene,
  "Aging" = age_geneset
  )

saveRDS(markers, file = "gene-sets/Pathway_Gene_Sets_GRP1.RDS")
markers <- readRDS("gene-sets/Pathway_Gene_Sets_GRP1.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/Pathway_Gene_Sets_GRP1.csv', row.names = FALSE, na = '')

Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub <- meta[,grep("_UCell", colnames(meta))]
saveRDS(meta.sub, file = "gene-sets/UCell_Gene_Set_Scores_snRNA_05212024.RDS")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_05212024.RDS")



###Protein Biomarker gene sets
bio.m <- read.delim("KPMP_ProteinBiomarker_GeneSets_05222024.txt")

markers <- list(
  "CRIC_ARIC_CKD" = bio.m$CRIC_ARIC_CKD,
  "Mend_Rand_CKD" = bio.m$Mend_Rand_CKD,
  "X10yr_Risk_CKD" = bio.m$X10yr_Risk_CKD,
  "CKD_Progression" = c(bio.m$CRIC_ARIC_CKD, bio.m$Mend_Rand_CKD, bio.m$X10yr_Risk_CKD),
  "ME15" = bio.m$ME15,
  "ME48" = bio.m$ME48
)

saveRDS(markers, file = "gene-sets/Pathway_Gene_Sets_GRP2.RDS")
markers <- readRDS("gene-sets/Pathway_Gene_Sets_GRP2.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/Pathway_Gene_Sets_GRP2.csv', row.names = FALSE, na = '')


Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub.2 <- meta[,grep("_UCell", colnames(meta))]

meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_05212024.RDS")
meta.sub <- cbind(meta.sub, meta.sub.2)

saveRDS(meta.sub, file = "gene-sets/UCell_Gene_Set_Scores_snRNA_05222024.RDS")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_05222024.RDS")



###Metabolic pathways

#HALLMARK_OXIDATIVE_PHOSPHORYLATION
OxPhos_geneset <- msigdbr_genesets[msigdbr_genesets$gs_id == "M5936",]$gene_symbol

#HALLMARK_FATTY_ACID_METABOLISM
FatAcid_geneset <- msigdbr_genesets[msigdbr_genesets$gs_id == "M5935",]$gene_symbol

#HALLMARK_GLYCOLYSIS
Gly_geneset <- msigdbr_genesets[msigdbr_genesets$gs_id == "M5937",]$gene_symbol

###Combined gene sets
markers <- list(
  "OxPhos" = OxPhos_geneset,
  "FatAcid" = FatAcid_geneset,
  "Gly" = Gly_geneset
  )
saveRDS(markers, file = "gene-sets/Pathway_Gene_Sets_GRP3.RDS")
markers <- readRDS("gene-sets/Pathway_Gene_Sets_GRP3.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/Pathway_Gene_Sets_GRP3.csv', row.names = FALSE, na = '')

Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta

meta.sub.2 <- meta[,grep("_UCell", colnames(meta))]

meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_05222024.RDS")
meta.sub <- cbind(meta.sub, meta.sub.2)

saveRDS(meta.sub, file = "gene-sets/UCell_Gene_Set_Scores_snRNA_05232024.RDS")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_05232024.RDS")


###Plasma ATI gene sets
ATI_geneset <- read.delim("gene-sets/Schmidt_2024_ATI_Plasma_Markers.txt")

markers <- list(
  "Plasma_ATI_BKBC" = ATI_geneset$Plasma_ATI_BKBC,
  "Plasma_ATI_ARIC" = ATI_geneset$Plasma_ATI_ARIC,
  "Plasma_ATI_CHROME" = ATI_geneset$Plasma_ATI_CHROME
  )

saveRDS(markers, file = "gene-sets/Pathway_Gene_Sets_GRP4.RDS")
markers <- readRDS("gene-sets/Pathway_Gene_Sets_GRP4.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/Pathway_Gene_Sets_GRP4.csv', row.names = FALSE, na = '')


Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub.2 <- meta[,grep("_UCell", colnames(meta))]

meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_05232024.RDS")
meta.sub <- cbind(meta.sub, meta.sub.2)

saveRDS(meta.sub, file = "gene-sets/UCell_Gene_Set_Scores_snRNA_09102024.RDS")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_09102024.RDS")



###ProtAge gene sets
ProtAge_geneset <- read.delim("gene-sets/Argentieri2024_ProtAge_Gene_Set.txt")

markers <- list(
  "ProtAge20" = ProtAge_geneset$ProtAge20,
  "ProtAge" = ProtAge_geneset$ProtAge
)

saveRDS(markers, file = "gene-sets/Pathway_Gene_Sets_GRP5.RDS")
markers <- readRDS("gene-sets/Pathway_Gene_Sets_GRP5.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/Pathway_Gene_Sets_GRP5.csv', row.names = FALSE, na = '')


Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub.2 <- meta[,grep("_UCell", colnames(meta))]

meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_09102024.RDS")
meta.sub <- cbind(meta.sub, meta.sub.2)

saveRDS(meta.sub, file = "gene-sets/UCell_Gene_Set_Scores_snRNA_09122024.RDS")
meta.sub <- readRDS("gene-sets/UCell_Gene_Set_Scores_snRNA_09122024.RDS")




###GRNs
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

saveRDS(markers, file = "gene-sets/GRN_Gene_Sets_GRP1.RDS")
markers <- readRDS("gene-sets/GRN_Gene_Sets_GRP1.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/GRN_Gene_Sets_GRP1.csv', row.names = FALSE, na = '')


###Add scores to main objects
Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub <- meta[,grep("_UCell", colnames(meta))]

saveRDS(meta.sub, file = "gene-sets/UCell_GRN_Scores_snRNA_09162024.RDS")
meta.sub <- readRDS("gene-sets/UCell_GRN_Scores_snRNA_09162024.RDS")





###GRNs for PT-S3 and MAC trajectories
files <- c("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-frPTS3_l5_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-PTS3_l4_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_MON-moFAM_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_MON-CXCL10_Trajectory_GRN_0424-newData.RDA",
           "~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_resMAC_Trajectory_GRN_0424-newData.RDA")
names <- setNames(c("frPTS3","PTS3","moFAM","MONCXCL10","resMAC"),files)

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

saveRDS(markers, file = "gene-sets/GRN_Gene_Sets_GRP2.RDS")
markers <- readRDS("gene-sets/GRN_Gene_Sets_GRP2.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/GRN_Gene_Sets_GRP2.csv', row.names = FALSE, na = '')

###Add scores to main objects
Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub2 <- meta[,grep("_UCell", colnames(meta))]
meta.sub <- readRDS("gene-sets/UCell_GRN_Scores_snRNA_09162024.RDS")
meta.sub <- cbind(meta.sub, meta.sub2)

saveRDS(meta.sub, file = "gene-sets/UCell_GRN_Scores_snRNA_09182024.RDS")
meta.sub <- readRDS("gene-sets/UCell_GRN_Scores_snRNA_09182024.RDS")







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

saveRDS(markers, file = "gene-sets/GRN_Gene_Sets_GRP3.RDS")
markers <- readRDS("gene-sets/GRN_Gene_Sets_GRP3.RDS")
max_length <- max(sapply(markers, length))
padded_markers <- lapply(markers, function(v) {
  length(v) <- max_length  # Set the length of each vector to the maximum length
  v  # Return the modified vector
})
markers_df <- do.call(cbind, padded_markers)
write.csv(markers_df, 'gene-sets/GRN_Gene_Sets_GRP3.csv', row.names = FALSE, na = '')

###Add scores to main objects
Idents(KB) <- "library"
celltype <- levels(Idents(KB))

meta <- do.call(rbind, lapply(celltype, function(ct) {
  s.obj <- subset(KB, library %in% ct)
  s.obj[["RNA"]] <- as(s.obj[["RNA"]], Class = "Assay")
  s.obj <- AddModuleScore_UCell(s.obj, features = markers)
  meta <- s.obj@meta.data
  return(meta)
}))
meta
meta.sub2 <- meta[,grep("_UCell", colnames(meta))]
meta.sub <- readRDS("gene-sets/UCell_GRN_Scores_snRNA_09182024.RDS")
meta.sub <- meta.sub[,!colnames(meta.sub) %in% colnames(meta.sub)[c(grep("resMAC", colnames(meta.sub)),grep("MONCXCL10", colnames(meta.sub)),grep("moFAM", colnames(meta.sub)))]]
meta.sub <- cbind(meta.sub, meta.sub2)

saveRDS(meta.sub, file = "gene-sets/UCell_GRN_Scores_snRNA_10012024.RDS")
meta.sub <- readRDS("gene-sets/UCell_GRN_Scores_snRNA_10012024.RDS")
