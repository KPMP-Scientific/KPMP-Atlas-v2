library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(pagoda2)
library(igraph)
library(BPCells)
options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2")
options(Seurat.object.assay.version = "v5")



###By Cluster
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

#Update metadata
path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta


##Interstitial fibrosis groups
table(KB.sub$percent_int_fibrosis)
unique(KB.sub$percent_int_fibrosis)

KB.sub$IF_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_int_fibrosis %in% c("0","1","1 to 4","2","3","5","10","15","20"),]$IF_group <- "<=20"
KB.sub@meta.data[KB.sub@meta.data$percent_int_fibrosis %in% c("25","30","35","40","45","50","60","70","80"),]$IF_group <- ">20"

KB.sub <- subset(KB.sub, IF_group %in% c("<=20",">20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$IF_group == "<=20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$IF_group == ">20",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Interstitial_Fibrosis_Group_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$IF_group))[select.donors, ]

df_comp_relative$cond <- "<=20"
df_comp_relative$cond[df_cond$'>20' != 0] <- ">20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<=20",">20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Interstitial-Fibrosis_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<=20"
y <- ">20"

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
saveRDS(p.table, file = "Interstitial-Fibrosis_Group_composition_V2_Clusters_P-Values.RDS")


##Tubular Atrophy Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Tubular Atrophy groups
table(KB.sub$percent_tub_atrophy)
unique(KB.sub$percent_tub_atrophy)
table(KB.sub$percent_tub_atrophy, KB.sub$condition_level1)

KB.sub$TA_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_atrophy %in% c("0","1","1 to 4","2","3","5","10","15"),]$TA_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_atrophy %in% c("20","25","30","35","40","45","50","60","70","80"),]$TA_group <- ">=20"

KB.sub <- subset(KB.sub, TA_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TA_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TA_group == ">=20",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Tubular-Atrophy_Group_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$TA_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Tubular-Atrophy_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Tubular-Atrophy_Group_composition_V2_Clusters_P-Values.RDS")




##Tubular Injury Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Tubular Injury groups
table(KB.sub$percent_tub_injury)
unique(KB.sub$percent_tub_injury)
table(KB.sub$percent_tub_injury, KB.sub$condition_level1)

KB.sub$TI_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_injury %in% c("0","1","1 to 4","1 TO 4","2","5","10","15"),]$TI_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_injury %in% c("20","25","30","35","40","50","60","70","75","80","90"),]$TI_group <- ">=20"

KB.sub <- subset(KB.sub, TI_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TI_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TI_group == ">=20",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Tubular-Injury_Group_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$TI_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Tubular-Injury_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Tubular-Injury_Group_composition_V2_Clusters_P-Values.RDS")



##Abnormal Luminal Morphology Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Abnormal Luminal Morphology groups
table(KB.sub$percent_ab_lum_morph)
unique(KB.sub$percent_ab_lum_morph)
table(KB.sub$percent_ab_lum_morph, KB.sub$condition_level1)

KB.sub$ALM_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_ab_lum_morph %in% c("0","1","1 to 4","2","5","10","15"),]$ALM_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_ab_lum_morph %in% c("20","25","30","40","50","60","70","80"),]$ALM_group <- ">=20"

KB.sub <- subset(KB.sub, ALM_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$ALM_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$ALM_group == ">=20",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Abnormal_Luminal_Morphology_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$ALM_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Abnormal_Luminal_Morphology_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Abnormal_Luminal_Morphology_Group_composition_V2_Clusters_P-Values.RDS")


##Acellular Cast Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Acellular Cast groups
table(KB.sub$percent_acellular_casts)
unique(KB.sub$percent_acellular_casts)
table(KB.sub$percent_acellular_casts, KB.sub$condition_level1)

KB.sub$AC_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_acellular_casts %in% c("0","0-5","1","1-4%","1 to 4","1 TO 4","1 to 4 ","2","4","5","10","10, ","15"),]$AC_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_acellular_casts %in% c("20","25","30"),]$AC_group <- ">=20"

KB.sub <- subset(KB.sub, AC_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AC_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AC_group == ">=20",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Acellular_Cast_Morphology_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$AC_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Acellular_Cast_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Acellular_Cast_Group_composition_V2_Clusters_P-Values.RDS")


##Interstitial Monuclear WBC Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Interstitial Monuclear WBC groups
table(KB.sub$percent_int_WBCs)
unique(KB.sub$percent_int_WBCs)
table(KB.sub$percent_int_WBCs, KB.sub$condition_level1)

KB.sub$WBC_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_int_WBCs %in% c("0","2","1-4%","5","10"),]$WBC_group <- "<=10"
KB.sub@meta.data[KB.sub@meta.data$percent_int_WBCs %in% c("15","20","25","30","40","45","50","55","70","80","90"),]$WBC_group <- ">10"

KB.sub <- subset(KB.sub, WBC_group %in% c("<=10",">10"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$WBC_group == "<=10",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$WBC_group == ">10",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Interstitial_WBC_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$WBC_group))[select.donors, ]

df_comp_relative$cond <- "<=10"
df_comp_relative$cond[df_cond$'>10' != 0] <- ">10"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<=10",">10"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Interstitial_WBC_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<=10"
y <- ">10"

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
saveRDS(p.table, file = "Interstitial_WBC_Group_composition_V2_Clusters_P-Values.RDS")



##Arteriosclerosis Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Arteriosclerosis groups
table(KB.sub$arteriosclerosis)
unique(KB.sub$arteriosclerosis)
table(KB.sub$arteriosclerosis, KB.sub$condition_level1)

KB.sub$AS_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$arteriosclerosis %in% c("0","1"),]$AS_group <- "0-1"
KB.sub@meta.data[KB.sub@meta.data$arteriosclerosis %in% c("2","3"),]$AS_group <- "2-3"

KB.sub <- subset(KB.sub, AS_group %in% c("0-1","2-3"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AS_group == "0-1",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AS_group == "2-3",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Arteriosclerosis_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$AS_group))[select.donors, ]

df_comp_relative$cond <- "0-1"
df_comp_relative$cond[df_cond$'2-3' != 0] <- "2-3"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("0-1","2-3"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Arteriosclerosis_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "0-1"
y <- "2-3"

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
saveRDS(p.table, file = "Arteriosclerosis_Group_composition_V2_Clusters_P-Values.RDS")


##Arteriolar Hyalinosis Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Arteriolar Hyalinosis groups
table(KB.sub$arteriolar_hyalinosis)
unique(KB.sub$arteriolar_hyalinosis)
table(KB.sub$arteriolar_hyalinosis, KB.sub$condition_level1)

KB.sub$AH_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$arteriolar_hyalinosis %in% c("0","1"),]$AH_group <- "0-1"
KB.sub@meta.data[KB.sub@meta.data$arteriolar_hyalinosis %in% c("2","3"),]$AH_group <- "2-3"

KB.sub <- subset(KB.sub, AH_group %in% c("0-1","2-3"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
HRT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AH_group == "0-1",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AH_group == "2-3",]$patient))
sample.tab <- rbind(HRT,AKI,CKD,DMR,group1,group2)
saveRDS(sample.tab, file = "Arteriolar_Hyalinosis_composition_sample_table.RDS")

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$AH_group))[select.donors, ]

df_comp_relative$cond <- "0-1"
df_comp_relative$cond[df_cond$'2-3' != 0] <- "2-3"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("0-1","2-3"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Arteriolar_Hyalinosis_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "0-1"
y <- "2-3"

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
saveRDS(p.table, file = "Arteriolar_Hyalinosis_Group_composition_V2_Clusters_P-Values.RDS")





###By Subclasses
##Interstitial fibrosis groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Interstitial fibrosis groups
table(KB.sub$percent_int_fibrosis)
unique(KB.sub$percent_int_fibrosis)

KB.sub$IF_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_int_fibrosis %in% c("0","1","1 to 4","2","3","5","10","15","20"),]$IF_group <- "<=20"
KB.sub@meta.data[KB.sub@meta.data$percent_int_fibrosis %in% c("25","30","35","40","45","50","60","70","80"),]$IF_group <- ">20"

KB.sub <- subset(KB.sub, IF_group %in% c("<=20",">20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$IF_group == "<=20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$IF_group == ">20",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      25
#CKD      51
#DMR       5
#NHT       0
#RT.UCS    0
#group1   48
#group2   33


##Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$IF_group))[select.donors, ]

df_comp_relative$cond <- "<=20"
df_comp_relative$cond[df_cond$'>20' != 0] <- ">20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<=20",">20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Interstitial-Fibrosis_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<=20"
y <- ">20"

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
saveRDS(p.table, file = "Interstitial-Fibrosis_Group_composition_V2_Subclassl3_P-Values.RDS")


##Tubular Atrophy Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Tubular Atrophy groups
table(KB.sub$percent_tub_atrophy)
unique(KB.sub$percent_tub_atrophy)
table(KB.sub$percent_tub_atrophy, KB.sub$condition_level1)

KB.sub$TA_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_atrophy %in% c("0","1","1 to 4","2","3","5","10","15"),]$TA_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_atrophy %in% c("20","25","30","35","40","45","50","60","70","80"),]$TA_group <- ">=20"

KB.sub <- subset(KB.sub, TA_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TA_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TA_group == ">=20",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      25
#CKD      51
#DMR       5
#NHT       0
#RT.UCS    0
#group1   45
#group2   36

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$TA_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Tubular-Atrophy_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Tubular-Atrophy_Group_composition_V2_Subclassl3_P-Values.RDS")



##Tubular Injury Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Tubular Injury groups
table(KB.sub$percent_tub_injury)
unique(KB.sub$percent_tub_injury)
table(KB.sub$percent_tub_injury, KB.sub$condition_level1)

KB.sub$TI_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_injury %in% c("0","1","1 to 4","1 TO 4","2","5","10","15"),]$TI_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_tub_injury %in% c("20","25","30","35","40","50","60","70","75","80","90"),]$TI_group <- ">=20"

KB.sub <- subset(KB.sub, TI_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TI_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$TI_group == ">=20",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      25
#CKD      51
#DMR       5
#NHT       0
#RT.UCS    0
#group1   38
#group2   43

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$TI_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Tubular-Injury_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Tubular-Injury_Group_composition_V2_Subclassl3_P-Values.RDS")


##Abnormal Luminal Morphology Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Abnormal Luminal Morphology groups
table(KB.sub$percent_ab_lum_morph)
unique(KB.sub$percent_ab_lum_morph)
table(KB.sub$percent_ab_lum_morph, KB.sub$condition_level1)

KB.sub$ALM_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_ab_lum_morph %in% c("0","1","1 to 4","2","5","10","15"),]$ALM_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_ab_lum_morph %in% c("20","25","30","40","50","60","70","80"),]$ALM_group <- ">=20"

KB.sub <- subset(KB.sub, ALM_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$ALM_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$ALM_group == ">=20",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      25
#CKD      51
#DMR       5
#NHT       0
#RT.UCS    0
#group1   52
#group2   29

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$ALM_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Abnormal_Luminal_Morphology_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Abnormal_Luminal_Morphology_Group_composition_V2_Subclassl3_P-Values.RDS")



##Pathology score groupings - Acellular Cast Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
KB

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Acellular Cast groups
table(KB.sub$percent_acellular_casts)
unique(KB.sub$percent_acellular_casts)
table(KB.sub$percent_acellular_casts, KB.sub$condition_level1)

KB.sub$AC_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_acellular_casts %in% c("0","0-5","1","1-4%","1 to 4","1 TO 4","1 to 4 ","2","4","5","10","10, ","15"),]$AC_group <- "<20"
KB.sub@meta.data[KB.sub@meta.data$percent_acellular_casts %in% c("20","25","30"),]$AC_group <- ">=20"

KB.sub <- subset(KB.sub, AC_group %in% c("<20",">=20"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AC_group == "<20",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AC_group == ">=20",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      25
#CKD      50
#DMR       5
#NHT       0
#RT.UCS    0
#group1   72
#group2    8

##Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$AC_group))[select.donors, ]

df_comp_relative$cond <- "<20"
df_comp_relative$cond[df_cond$'>=20' != 0] <- ">=20"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<20",">=20"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Acellular_Cast_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<20"
y <- ">=20"

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
saveRDS(p.table, file = "Acellular_Cast_Group_composition_V2_Subclassl3_P-Values.RDS")


##Interstitial Monuclear WBC Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Interstitial Monuclear WBC groups
table(KB.sub$percent_int_WBCs)
unique(KB.sub$percent_int_WBCs)
table(KB.sub$percent_int_WBCs, KB.sub$condition_level1)

KB.sub$WBC_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$percent_int_WBCs %in% c("0","2","1-4%","5","10"),]$WBC_group <- "<=10"
KB.sub@meta.data[KB.sub@meta.data$percent_int_WBCs %in% c("15","20","25","30","40","45","50","55","70","80","90"),]$WBC_group <- ">10"

KB.sub <- subset(KB.sub, WBC_group %in% c("<=10",">10"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$WBC_group == "<=10",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$WBC_group == ">10",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      25
#CKD      51
#DMR       5
#NHT       0
#RT.UCS    0
#group1   45
#group2   36

#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$WBC_group))[select.donors, ]

df_comp_relative$cond <- "<=10"
df_comp_relative$cond[df_cond$'>10' != 0] <- ">10"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("<=10",">10"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Interstitial_WBC_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "<=10"
y <- ">10"

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
saveRDS(p.table, file = "Interstitial_WBC_Group_composition_V2_Subclassl3_P-Values.RDS")



##Arteriosclerosis Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Arteriosclerosis groups
table(KB.sub$arteriosclerosis)
unique(KB.sub$arteriosclerosis)
table(KB.sub$arteriosclerosis, KB.sub$condition_level1)

KB.sub$AS_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$arteriosclerosis %in% c("0","1"),]$AS_group <- "0-1"
KB.sub@meta.data[KB.sub@meta.data$arteriosclerosis %in% c("2","3"),]$AS_group <- "2-3"

KB.sub <- subset(KB.sub, AS_group %in% c("0-1","2-3"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AS_group == "0-1",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AS_group == "2-3",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      24
#CKD      50
#DMR       5
#NHT       0
#RT.UCS    0
#group1   45
#group2   34


#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$AS_group))[select.donors, ]

df_comp_relative$cond <- "0-1"
df_comp_relative$cond[df_cond$'2-3' != 0] <- "2-3"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("0-1","2-3"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Arteriosclerosis_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]
#T-Tests
x <- "0-1"
y <- "2-3"

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
saveRDS(p.table, file = "Arteriosclerosis_Group_composition_V2_Subclassl3_P-Values.RDS")


##Arteriolar Hyalinosis Groups
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")

path.meta <- read.delim("KPMP_TIV_Descriptor_Scores_8-29-2024_subset.txt")
path.meta <- path.meta[path.meta$percent_cortex >= 25, ]
colnames(path.meta)
KB.sub <- subset(KB, patient %in% path.meta$patient)

meta <- KB.sub@meta.data
emc <- c("percent_int_fibrosis","percent_tub_atrophy","percent_tub_injury","percent_ab_lum_morph",
         "percent_acellular_casts","percent_int_WBCs","arteriosclerosis","arteriolar_hyalinosis")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- path.meta[,emc[i]][match(meta$patient, path.meta$patient)]
}
colnames(meta)
KB.sub@meta.data <- meta

#Cell type proportions for Arteriolar Hyalinosis groups
table(KB.sub$arteriolar_hyalinosis)
unique(KB.sub$arteriolar_hyalinosis)
table(KB.sub$arteriolar_hyalinosis, KB.sub$condition_level1)

KB.sub$AH_group <- "NA"
KB.sub@meta.data[KB.sub@meta.data$arteriolar_hyalinosis %in% c("0","1"),]$AH_group <- "0-1"
KB.sub@meta.data[KB.sub@meta.data$arteriolar_hyalinosis %in% c("2","3"),]$AH_group <- "2-3"

KB.sub <- subset(KB.sub, AH_group %in% c("0-1","2-3"))
table(KB.sub$region_level2)
#Low coverage of Medulla, remove 
KB.sub <- subset(KB.sub, region_level2 %in% c("M"), invert = TRUE)

#sample info
Ref <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "HRT",]$patient)) 
AKI <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "AKI",]$patient))
CKD <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "CKD",]$patient))
DMR <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "DM-R",]$patient))
NHT <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "NHT",]$patient))
RT.UCS <- length(unique(KB.sub@meta.data[KB.sub@meta.data$condition_level1 == "RT-UCS",]$patient))
group1 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AH_group == "0-1",]$patient))
group2 <- length(unique(KB.sub@meta.data[KB.sub@meta.data$AH_group == "2-3",]$patient))
sample.tab <- rbind(Ref,AKI,CKD,DMR,NHT,RT.UCS,group1,group2)
sample.tab
#Ref       0
#AKI      25
#CKD      51
#DMR       5
#NHT       0
#RT.UCS    0
#group1   49
#group2   32


#Differential composition analysis - Cluster level
df_comp <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB.sub$patient, KB.sub$AH_group))[select.donors, ]

df_comp_relative$cond <- "0-1"
df_comp_relative$cond[df_cond$'2-3' != 0] <- "2-3"
df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("0-1","2-3"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "Arteriolar_Hyalinosis_Group_composition_V2_Subclassl3.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]
#T-Tests
x <- "0-1"
y <- "2-3"

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
saveRDS(p.table, file = "Arteriolar_Hyalinosis_Group_composition_V2_Subclassl3_P-Values.RDS")
