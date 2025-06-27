library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
library(anndata)
options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")
source("misc/utils.R")


load("color_factors_v2-clusters.robj")
load("color_factors.robj")



###By Clusters
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB <- subset(KB, tissue_type %in% "Biopsy")
KB
table(KB$region_level2)
table(KB$condition_level1)

#sample info
Ref <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "RT-UCS",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS)
sample.tab
#Ref      41
#AKI      38
#CKD      73
#DMR      11
#NHT       2
#RT.UCS    0


##Differential composition analysis
KB <- subset(KB, condition_level1 %in% c("HRT","AKI","CKD"))
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$condition_level1))[select.donors, ]

df_comp_relative$cond <- "HRT"
df_comp_relative$cond[df_cond$AKI != 0] <- "AKI"
df_comp_relative$cond[df_cond$CKD != 0] <- "CKD"

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("HRT","AKI","CKD"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Condition_Group_composition_V2_Clusters_biopsy.RDS")
df_comp_relative <- readRDS("Condition_Group_composition_V2_Clusters_biopsy.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "HRT"
y <- "AKI"

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
saveRDS(p.table, file = "Condition_Ref-AKI_Group_composition_V2_Clusters_P-Values_biopsy.RDS")

x <- "HRT"
y <- "CKD"

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
saveRDS(p.table, file = "Condition_Ref-CKD_Group_composition_V2_Clusters_P-Values_biopsy.RDS")

x <- "AKI"
y <- "CKD"

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
saveRDS(p.table, file = "Condition_AKI-CKD_Group_composition_V2_Clusters_P-Values_biopsy.RDS")


##Comp of AIN/ATN
table(KB$adj_cm_ati.ain)
AIN <- unique(KB@meta.data[KB@meta.data$adj_cm_ati.ain %in% c("ADJ_CM_AIN","ADJ_CM_ATI_POSSIBLY,ADJ_CM_AIN","ADJ_CM_ATI,ADJ_CM_AIN"),]$patient)
#n=11
ATI <- unique(KB@meta.data[KB@meta.data$adj_cm_ati.ain %in% c("ADJ_CM_ATI","ADJ_CM_ATI,ADJ_CM_AIN_POSSIBLY"),]$patient)
#n=40

#T-Tests
x <- "HRT"
y <- "AKI"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% ATI & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-AKI-ATI_Group_composition_V2_Clusters_P-Values_biopsy.RDS")

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% AIN & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-AKI-AIN_Group_composition_V2_Clusters_P-Values_biopsy.RDS")

x <- "AKI"
y <- "AKI"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[rownames(df_comp_relative) %in% ATI & df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% AIN & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-ATI-AIN_Group_composition_V2_Clusters_P-Values_biopsy.RDS")


##Comp of CKD Risk
table(KB$ckd_stageC)
CKD.hi <- unique(KB@meta.data[KB@meta.data$ckd_stageC %in% c("high risk","very high risk"),]$patient)
#n=26
CKD.lo <- unique(KB@meta.data[KB@meta.data$ckd_stageC %in% c("low risk","moderately increased risk"),]$patient)
#n=9

#T-Tests
x <- "HRT"
y <- "CKD"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% CKD.hi & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-CKD-HighRisk_Group_composition_V2_Clusters_P-Values_biopsy.RDS")


x <- "CKD"
y <- "CKD"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[rownames(df_comp_relative) %in% CKD.lo & df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% CKD.hi & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_CKD-LowRisk-CKD-HighRisk_Group_composition_V2_Clusters_P-Values_biopsy.RDS")




###Cell type proportions for eGFR groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB <- subset(KB, tissue_type %in% "Biopsy")

#Stage 1 with normal or high GFR (GFR > 90 mL/min)
#Stage 2 Mild CKD (GFR = 60-89 mL/min)
#Stage 3A Moderate CKD (GFR = 45-59 mL/min)
#Stage 3B Moderate CKD (GFR = 30-44 mL/min)
#Stage 4 Severe CKD (GFR = 15-29 mL/min)
#Stage 5 End Stage CKD (GFR <15 mL/min
table(KB$baseline_eGFR_binned)

KB$eGFR_group <- "NA"
KB@meta.data[KB@meta.data$baseline_eGFR_binned %in% c("0-9 ml/min/1.73m2",
                                                      "20-29 ml/min/1.73m2"),]$eGFR_group <- "<30"
KB@meta.data[KB@meta.data$baseline_eGFR_binned %in% c("30-39 ml/min/1.73m2",
                                                      "40-49 ml/min/1.73m2","50-59 ml/min/1.73m2"),]$eGFR_group <- "30-60"
KB@meta.data[KB@meta.data$baseline_eGFR_binned %in% c(">60 ml/min/1.73m2","60-69 ml/min/1.73m2",
                                                      "70-79 ml/min/1.73m2","80-89 ml/min/1.73m2",
                                                      "90-99 ml/min/1.73m2",
                                                      "100-109 ml/min/1.73m2","110-119 ml/min/1.73m2",
                                                      "120-129 ml/min/1.73m2","130-139 ml/min/1.73m2",
                                                      "150-159 ml/min/1.73m2","170-179 ml/min/1.73m2"
),]$eGFR_group <- ">60"

KB <- subset(KB, eGFR_group %in% c("<30","30-60",">60"))
table(KB$region_level2)

#sample info
Ref <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB@meta.data[KB@meta.data$eGFR_group == ">60",]$patient))
group2 <- length(unique(KB@meta.data[KB@meta.data$eGFR_group %in% c("<30","30-60"),]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref      16
#AKI      34
#CKD      64
#DMR      11
#NHT       0
#RT.UCS    0
#group1   69
#group2   56

##Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$eGFR_group))[select.donors, ]

df_comp_relative$cond <- "eGFR >60"
df_comp_relative$cond[df_cond$'<30' != 0] <- "eGFR <60"
df_comp_relative$cond[df_cond$'30-60' != 0] <- "eGFR <60"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("eGFR >60", "eGFR <60"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "eGFR_Group_composition_V2_Clusters_biopsy.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "eGFR >60"
y <- "eGFR <60"

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
saveRDS(p.table, file = "eGFR_Group_composition_V2_Clusters_P-Values_biopsy.RDS")




###Cell type proportions for age groups (HRT only)
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB <- subset(KB, condition_level1 == "HRT")
table(KB$age_binned)

KB$age_group <- "NA"
KB@meta.data[KB@meta.data$age_binned %in% c("10-19", "20-29", "30-39", "40-49"),]$age_group <- "<50"
KB@meta.data[KB@meta.data$age_binned %in% c("50-59", "60-69","70-79","90-99"),]$age_group <- ">=50"

KB <- subset(KB, age_group %in% c("<50",">=50"))
table(KB$region_level2)

#sample info
Ref <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB@meta.data[KB@meta.data$age_group == "<50",]$patient))
group2 <- length(unique(KB@meta.data[KB@meta.data$age_group %in% c(">=50"),]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref      72
#group1   34
#group2   38


##Differential composition analysis - Subclass level
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$age_group))[select.donors, ]

df_comp_relative$cond <- "<50"
df_comp_relative$cond[df_cond$'>=50' != 0] <- ">=50"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<50", ">=50"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "HRT_Age_Group_composition_V2_clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<50"
y <- ">=50"

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
saveRDS(p.table, file = "HRT_Age_Group_composition_V2_Clusters_P-Values.RDS")





### By subclass.l3 
###By Condition type
#Load data and remove ref samples with high altered states
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB
table(KB$region_level2)
table(KB$condition_level1)
KB <- subset(KB, tissue_type %in% "Biopsy")

#sample info
Ref <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "RT-UCS",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS)
#saveRDS(sample.tab, file = "Condition_Group_composition_sample_table.RDS")
sample.tab
#Ref      41
#AKI      38
#CKD      73
#DMR      11
#NHT       2
#RT.UCS    0


##Differential composition analysis
KB <- subset(KB, condition_level1 %in% c("HRT","AKI","CKD"))
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$condition_level1))[select.donors, ]

df_comp_relative$cond <- "HRT"
df_comp_relative$cond[df_cond$AKI != 0] <- "AKI"
df_comp_relative$cond[df_cond$CKD != 0] <- "CKD"

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("HRT","AKI","CKD"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Condition_Group_composition_V2_Subclass.l3_biopsy.RDS")
df_comp_relative <- readRDS("Condition_Group_composition_V2_Subclass.l3_biopsy.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "HRT"
y <- "AKI"

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
saveRDS(p.table, file = "Condition_Ref-AKI_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")

x <- "HRT"
y <- "CKD"

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
saveRDS(p.table, file = "Condition_Ref-CKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")

x <- "AKI"
y <- "CKD"

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
saveRDS(p.table, file = "Condition_AKI-CKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")


##Comp of AIN/ATN
table(KB$adj_cm_ati.ain)
AIN <- unique(KB@meta.data[KB@meta.data$adj_cm_ati.ain %in% c("ADJ_CM_AIN","ADJ_CM_ATI_POSSIBLY,ADJ_CM_AIN","ADJ_CM_ATI,ADJ_CM_AIN"),]$patient)
#n=11
ATI <- unique(KB@meta.data[KB@meta.data$adj_cm_ati.ain %in% c("ADJ_CM_ATI","ADJ_CM_ATI,ADJ_CM_AIN_POSSIBLY"),]$patient)
#n=40

#T-Tests
x <- "HRT"
y <- "AKI"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% ATI & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-AKI-ATI_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% AIN & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-AKI-AIN_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")

x <- "AKI"
y <- "AKI"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[rownames(df_comp_relative) %in% ATI & df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% AIN & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-ATI-AIN_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")


##Comp of CKD Risk
table(KB$ckd_stageC)
CKD.hi <- unique(KB@meta.data[KB@meta.data$ckd_stageC %in% c("high risk","very high risk"),]$patient)
#n=26
CKD.lo <- unique(KB@meta.data[KB@meta.data$ckd_stageC %in% c("low risk","moderately increased risk"),]$patient)
#n=9

#T-Tests
x <- "HRT"
y <- "CKD"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% CKD.hi & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_Ref-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")


x <- "CKD"
y <- "CKD"

p.table <- do.call(rbind, lapply(clusters, function(ct) {
  print(paste("Running for Cluster:", ct))
  data.x <- df_comp_relative[rownames(df_comp_relative) %in% CKD.lo & df_comp_relative$cond %in% x,ct]
  data.y <- df_comp_relative[rownames(df_comp_relative) %in% CKD.hi & df_comp_relative$cond %in% y,ct]
  p <- t.test(x = data.x, y = data.y, alternative = c("two.sided", "less", "greater"), mu = 0, 
              paired = FALSE, var.equal = FALSE, conf.level = 0.95)
  p.df <- data.frame(p = p$p.value, t = p$statistic)
  p.df
}))

rownames(p.table) <- clusters
saveRDS(p.table, file = "Condition_CKD-LowRisk-CKD-HighRisk_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")




###Cell type proportions for eGFR groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB <- subset(KB, tissue_type %in% "Biopsy")

#Stage 1 with normal or high GFR (GFR > 90 mL/min)
#Stage 2 Mild CKD (GFR = 60-89 mL/min)
#Stage 3A Moderate CKD (GFR = 45-59 mL/min)
#Stage 3B Moderate CKD (GFR = 30-44 mL/min)
#Stage 4 Severe CKD (GFR = 15-29 mL/min)
#Stage 5 End Stage CKD (GFR <15 mL/min
table(KB$baseline_eGFR_binned)

KB$eGFR_group <- "NA"
KB@meta.data[KB@meta.data$baseline_eGFR_binned %in% c("0-9 ml/min/1.73m2",
                                                      "20-29 ml/min/1.73m2"),]$eGFR_group <- "<30"
KB@meta.data[KB@meta.data$baseline_eGFR_binned %in% c("30-39 ml/min/1.73m2",
                                                      "40-49 ml/min/1.73m2","50-59 ml/min/1.73m2"),]$eGFR_group <- "30-60"
KB@meta.data[KB@meta.data$baseline_eGFR_binned %in% c(">60 ml/min/1.73m2","60-69 ml/min/1.73m2",
                                                      "70-79 ml/min/1.73m2","80-89 ml/min/1.73m2",
                                                      "90-99 ml/min/1.73m2",
                                                      "100-109 ml/min/1.73m2","110-119 ml/min/1.73m2",
                                                      "120-129 ml/min/1.73m2","130-139 ml/min/1.73m2",
                                                      "150-159 ml/min/1.73m2","170-179 ml/min/1.73m2"
),]$eGFR_group <- ">60"

KB <- subset(KB, eGFR_group %in% c("<30","30-60",">60"))
table(KB$region_level2)

#sample info
Ref <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB@meta.data[KB@meta.data$eGFR_group == ">60",]$patient))
group2 <- length(unique(KB@meta.data[KB@meta.data$eGFR_group %in% c("<30","30-60"),]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
#saveRDS(sample.tab, file = "eGFR_Group_composition_sample_table.RDS")
sample.tab
#Ref      16
#AKI      34
#CKD      64
#DMR      11
#NHT       0
#RT.UCS    0
#group1   69
#group2   56

##Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$eGFR_group))[select.donors, ]

df_comp_relative$cond <- "eGFR >60"
df_comp_relative$cond[df_cond$'<30' != 0] <- "eGFR <60"
df_comp_relative$cond[df_cond$'30-60' != 0] <- "eGFR <60"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("eGFR >60", "eGFR <60"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "eGFR_Group_composition_V2_Subclassl3_biopsy.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "eGFR >60"
y <- "eGFR <60"

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
saveRDS(p.table, file = "eGFR_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")





###Cell type proportions for age groups (HRT only)
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
#KB <- subset(KB, tissue_type %in% "Biopsy")
KB <- subset(KB, condition_level1 == "HRT")
table(KB$age_binned)

KB$age_group <- "NA"
KB@meta.data[KB@meta.data$age_binned %in% c("10-19", "20-29", "30-39", "40-49"),]$age_group <- "<50"
KB@meta.data[KB@meta.data$age_binned %in% c("50-59", "60-69","70-79","90-99"),]$age_group <- ">=50"

KB <- subset(KB, age_group %in% c("<50",">=50"))
table(KB$region_level2)

#sample info
Ref <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB@meta.data[KB@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB@meta.data[KB@meta.data$age_group == "<50",]$patient))
group2 <- length(unique(KB@meta.data[KB@meta.data$age_group %in% c(">=50"),]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref      72
#group1   34
#group2   38


##Differential composition analysis - Subclass level
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$age_group))[select.donors, ]

df_comp_relative$cond <- "<50"
df_comp_relative$cond[df_cond$'>=50' != 0] <- ">=50"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<50", ">=50"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "HRT_Age_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<50"
y <- ">=50"

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
saveRDS(p.table, file = "HRT_Age_Group_composition_V2_Subclassl3_P-Values.RDS")











###Additional Clinical Phenotypes by subclass.l3 ----------------------------------------

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB <- subset(KB, tissue_type %in% "Biopsy")
KB
table(KB$region_level2)
table(KB$condition_level1)


##Differential composition analysis
KB <- subset(KB, condition_level1 %in% c("HRT","AKI","CKD"))
clin.data <- read.csv("kpmp-core-clinical-data-2024-08-05.csv")
rownames(clin.data) <- clin.data$study_id
clin.data <- clin.data[rownames(clin.data) %in% KB@meta.data$patient,]

meta <- KB@meta.data
emc <- c("adj_primary_catC","mh_diabetes_type")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- clin.data[,emc[i]][match(meta$patient, clin.data$study_id)]
}
colnames(meta)
KB@meta.data <- meta

KB@meta.data$condition <- KB@meta.data$adj_primary_catC
KB@meta.data[KB@meta.data$condition_level1 %in% "HRT",]$condition <- "HRT"
KB@meta.data[KB@meta.data$condition %in% "Diabetic Kidney Disease",]$condition <- "DKD"
KB@meta.data[KB@meta.data$condition %in% "Hypertensive Kidney Disease",]$condition <- "H-CKD"
KB@meta.data[KB@meta.data$condition %in% "Acute Interstitial Nephritis",]$condition <- "AKI"
KB@meta.data[KB@meta.data$condition %in% "Acute Tubular Injury",]$condition <- "AKI"

unique(KB$condition)
table(KB$condition)

KB <- subset(KB, condition %in% c("HRT","AKI","DKD","H-CKD"))


df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$condition))[select.donors, ]

df_comp_relative$cond <- "HRT"
df_comp_relative$cond[df_cond$AKI != 0] <- "AKI"
df_comp_relative$cond[df_cond$DKD != 0] <- "DKD"
df_comp_relative$cond[df_cond$'H-CKD' != 0] <- "H-CKD"
table(df_comp_relative$cond)
#AKI   DKD H-CKD   HRT 
#28    31    24    41 

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("HRT","AKI","DKD","H-CKD"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Adj-Condition_Group_composition_V2_Subclass.l3_biopsy.RDS")
df_comp_relative <- readRDS("Adj-Condition_Group_composition_V2_Subclass.l3_biopsy.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "HRT"
y <- "AKI"

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
saveRDS(p.table, file = "Adj-Condition_Ref-AKI_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")

x <- "HRT"
y <- "DKD"

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
saveRDS(p.table, file = "Adj-Condition_Ref-DKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")


x <- "HRT"
y <- "H-CKD"

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
saveRDS(p.table, file = "Adj-Condition_Ref-H-CKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")


x <- "H-CKD"
y <- "DKD"

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
saveRDS(p.table, file = "Adj-Condition_H-CKD-DKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")



x <- "AKI"
y <- "DKD"

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
saveRDS(p.table, file = "Adj-Condition_AKI-DKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")


x <- "AKI"
y <- "H-CKD"

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
saveRDS(p.table, file = "Adj-Condition_AKI-H-CKD_Group_composition_V2_Subclassl3_P-Values_biopsy.RDS")
