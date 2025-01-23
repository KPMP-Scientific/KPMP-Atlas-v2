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




###By Clusters
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = c(paste0("P_",1:31),
                                            paste0("D_",1:47),
                                            paste0("E_",1:21),
                                            paste0("S_",1:25),
                                            paste0("I_",1:27),
                                            paste0("N_",1)))
table(Idents(KB))
KB

#Differential composition analysis
table(KB$condition_level3)
KB <- subset(KB, condition_level3 %in% c("Sham","4hrs","12hrs","2d","2wk","4-6wk","20-24wk"))
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$condition_level3))[select.donors, ]

df_comp_relative$cond <- "Sham"
df_comp_relative$cond[df_cond$`4hrs` != 0] <- "4hrs"
df_comp_relative$cond[df_cond$`12hrs` != 0] <- "12hrs"
df_comp_relative$cond[df_cond$`2d` != 0] <- "2d"
df_comp_relative$cond[df_cond$`2wk` != 0] <- "2wk"
df_comp_relative$cond[df_cond$`4-6wk` != 0] <- "4-6wk"
df_comp_relative$cond[df_cond$`20-24wk` != 0] <- "20-24wk"

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("Sham","4hrs","12hrs","2d","2wk","4-6wk","20-24wk"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "mouse_IRI/cond_plots/Condition_L3_Group_composition_V2_Clusters.RDS")
df_comp_relative <- readRDS("mouse_IRI/cond_plots/Condition_L3_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "Sham"
y <- "4hrs"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-4hrs_composition_V2_Clusters_P-Values.RDS")

x <- "Sham"
y <- "12hrs"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-12hrs_composition_V2_Clusters_P-Values.RDS")


x <- "Sham"
y <- "2d"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-2d_composition_V2_Clusters_P-Values.RDS")


x <- "Sham"
y <- "2wk"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-2wk_composition_V2_Clusters_P-Values.RDS")


x <- "Sham"
y <- "4-6wk"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-4-6wk_composition_V2_Clusters_P-Values.RDS")


x <- "Sham"
y <- "20-24wk"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-20-24wk_composition_V2_Clusters_P-Values.RDS")


##By Age Group
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
Idents(KB) <- "v2.clusters"
Idents(KB) <- factor(Idents(KB), levels = c(paste0("P_",1:31),
                                            paste0("D_",1:47),
                                            paste0("E_",1:21),
                                            paste0("S_",1:25),
                                            paste0("I_",1:27),
                                            paste0("N_",1)))
table(Idents(KB))
KB

#Differential composition analysis
table(KB$condition_level3)
KB <- subset(KB, condition_level3 %in% c("Ref"))
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.clusters))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$age_group))[select.donors, ]

df_comp_relative$cond <- "mature"
df_comp_relative$cond[df_cond$`middle-aged` != 0] <- "middle-aged"
df_comp_relative$cond[df_cond$`old` != 0] <- "old"

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("mature","middle-aged","old"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "mouse_IRI/cond_plots/Age_Group_composition_V2_Clusters.RDS")
df_comp_relative <- readRDS("mouse_IRI/cond_plots/Age_Group_composition_V2_Clusters.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "mature"
y <- "old"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Age_Group_composition_V2_Clusters_P-Values.RDS")


##Visualizations
IRI.4hrs.pVal <- readRDS("mouse_IRI/cond_plots/Sham-4hrs_composition_V2_Clusters_P-Values.RDS")
IRI.12hrs.pVal <- readRDS("mouse_IRI/cond_plots/Sham-12hrs_composition_V2_Clusters_P-Values.RDS")
IRI.2d.pVal <- readRDS("mouse_IRI/cond_plots/Sham-2d_composition_V2_Clusters_P-Values.RDS")
IRI.2wk.pVal <- readRDS("mouse_IRI/cond_plots/Sham-2wk_composition_V2_Clusters_P-Values.RDS")
IRI.4.6wk.pVal <- readRDS("mouse_IRI/cond_plots/Sham-4-6wk_composition_V2_Clusters_P-Values.RDS")
IRI.20.24wk.pVal <- readRDS("mouse_IRI/cond_plots/Sham-20-24wk_composition_V2_Clusters_P-Values.RDS")
age.pVal <- readRDS("mouse_IRI/cond_plots/Age_Group_composition_V2_Clusters_P-Values.RDS")

#P value plot
row.order <- c(paste0("P_",1:31),
               paste0("D_",1:47),
               paste0("E_",1:21),
               paste0("S_",1:25),
               paste0("I_",1:27),
               paste0("N_",1))
col.order <- c("IRI.4hrs.pVal","IRI.12hrs.pVal","IRI.2d.pVal",
               "IRI.2wk.pVal",
               "IRI.4.6wk.pVal","IRI.20.24wk.pVal","age.pVal")
p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]
t.stats <- as.matrix(t(t.stats[row.order,paste0(col.order,".t")]))

# Cap the values between -10 and 10
t.stats <- pmax(pmin(t.stats, 10), -10)
write.table(p.table[row.order,], file = "mouse_IRI/cond_plots/V2_Clusters_Composition_Analysis_p-Values.txt", sep = "\t", quote = FALSE)

p.stats <- as.matrix(t(p.table[row.order,paste0(col.order,".p")]))
pdf(file='mouse_IRI/cond_plots/V2_Cluster_Composition_Analysis_t-statistic_Heatmap.pdf',width=30,height=5.1)
heatmap.2(t.stats,col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          cellnote = ifelse(p.stats < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 






###By Subclasses
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
Idents(KB) <- "v2.subclass.l3"
order <- c("POD","PEC","PT-S1","PT-S1/2","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT",
           "aPT1-frPT","aPT-S1/S2","aPT-S3","dPT","cycPT","DTL1","DTL2","aDTL",
           "DTL3","ATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B",
           "C-TAL-B","dC-TAL","aTAL1","aTAL2","frTAL","infTAL","cycTAL","DCT1",
           "dDCT1","DCT2","dDCT2","aDCT","CNT","CNT-PC","dCNT-PC","CCD-PC","OMCD-PC",
           "dPC","IMCD","PapE","CCD-IC-A","OMCD-IC-A","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-EA","M-EC-PTC","C-EC-PTC","dEC-PTC","angEC-PTC","infEC-PTC","iaEC-PTC",
           "EC-V","EC-PCV","EC-AVR","P-EC-AVR","dEC-AVR","cycEC","EC-LYM","IM-FIB",
           "OM-FIB","dM-FIB","C/M-FIB","C-FIB","C-FIB-PATH","eaC-FIB","C-FIB-Cxcl14+",
           "C-FIB-Cxcl10+","C-MYOF","cycFIB","pvFIB-Rspo3+","pvFIB-Pi16+","pvFIB",
           "pvMYOF","MC","REN","VSMC","VSMC/P","Ad","B","PL","T","T-REG","Tcm/Naive",
           "CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","resMAC","resMAC-H2hi","resMAC-Cpd+",
           "MAC","moMAC-IA","moMAC-Cxcl10+","moFAM","moMAC","moMAC-Ccr2+","moMAC-Agr1+",
           "MON","cDC2","mDC","cDC1","pDC","cycMAC","cycDC","SC/NEU")
Idents(KB) <- factor(Idents(KB), levels = order)
table(Idents(KB))
KB

#Differential composition analysis
table(KB$condition_level3)
KB <- subset(KB, condition_level3 %in% c("Sham","4hrs","12hrs","2d","2wk","4-6wk","20-24wk"))
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$condition_level3))[select.donors, ]

df_comp_relative$cond <- "Sham"
df_comp_relative$cond[df_cond$`4hrs` != 0] <- "4hrs"
df_comp_relative$cond[df_cond$`12hrs` != 0] <- "12hrs"
df_comp_relative$cond[df_cond$`2d` != 0] <- "2d"
df_comp_relative$cond[df_cond$`2wk` != 0] <- "2wk"
df_comp_relative$cond[df_cond$`4-6wk` != 0] <- "4-6wk"
df_comp_relative$cond[df_cond$`20-24wk` != 0] <- "20-24wk"

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("Sham","4hrs","12hrs","2d","2wk","4-6wk","20-24wk"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "mouse_IRI/cond_plots/Condition_L3_Group_composition_V2_Subclasses.RDS")
df_comp_relative <- readRDS("mouse_IRI/cond_plots/Condition_L3_Group_composition_V2_Subclasses.RDS")

clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "Sham"
y <- "4hrs"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-4hrs_composition_V2_Subclasses_P-Values.RDS")

x <- "Sham"
y <- "12hrs"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-12hrs_composition_V2_Subclasses_P-Values.RDS")


x <- "Sham"
y <- "2d"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-2d_composition_V2_Subclasses_P-Values.RDS")


x <- "Sham"
y <- "2wk"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-2wk_composition_V2_Subclasses_P-Values.RDS")


x <- "Sham"
y <- "4-6wk"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-4-6wk_composition_V2_Subclasses_P-Values.RDS")


x <- "Sham"
y <- "20-24wk"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Sham-20-24wk_composition_V2_Subclasses_P-Values.RDS")


##By Age Group
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/mouse_IRI/seuratV5/Mouse_Kidney_Integrated_Atlas_V2_Seurat_11152024.Rds")
Idents(KB) <- "v2.subclass.l3"
Idents(KB) <- factor(Idents(KB), levels = order)
table(Idents(KB))
KB

#Differential composition analysis
table(KB$condition_level3)
KB <- subset(KB, condition_level3 %in% c("Ref"))
df_comp <- as.data.frame.matrix(table(KB$patient, KB$v2.subclass.l3))
select.donors <- rownames(df_comp)[rowSums(df_comp) > 50]
df_comp <- df_comp[select.donors, ]
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")

df_cond <- as.data.frame.matrix(table(KB$patient, KB$age_group))[select.donors, ]

df_comp_relative$cond <- "mature"
df_comp_relative$cond[df_cond$`middle-aged` != 0] <- "middle-aged"
df_comp_relative$cond[df_cond$`old` != 0] <- "old"

df_comp_relative$cond <- factor(df_comp_relative$cond, levels = c("mature","middle-aged","old"))
colnames(df_comp_relative) <- gsub("-","",colnames(df_comp_relative))
colnames(df_comp_relative) <- gsub(" ","",colnames(df_comp_relative))
saveRDS(df_comp_relative, file = "mouse_IRI/cond_plots/Age_Group_composition_V2_Subclasses.RDS")
df_comp_relative <- readRDS("mouse_IRI/cond_plots/Age_Group_composition_V2_Subclasses.RDS")

#Plot each cluster
clusters <- colnames(df_comp_relative)[!colnames(df_comp_relative) %in% "cond"]

#T-Tests
x <- "mature"
y <- "old"

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
saveRDS(p.table, file = "mouse_IRI/cond_plots/Age_Group_composition_V2_Subclasses_P-Values.RDS")

##Visualizations
IRI.4hrs.pVal <- readRDS("mouse_IRI/cond_plots/Sham-4hrs_composition_V2_Subclasses_P-Values.RDS")
IRI.12hrs.pVal <- readRDS("mouse_IRI/cond_plots/Sham-12hrs_composition_V2_Subclasses_P-Values.RDS")
IRI.2d.pVal <- readRDS("mouse_IRI/cond_plots/Sham-2d_composition_V2_Subclasses_P-Values.RDS")
IRI.2wk.pVal <- readRDS("mouse_IRI/cond_plots/Sham-2wk_composition_V2_Subclasses_P-Values.RDS")
IRI.4.6wk.pVal <- readRDS("mouse_IRI/cond_plots/Sham-4-6wk_composition_V2_Subclasses_P-Values.RDS")
IRI.20.24wk.pVal <- readRDS("mouse_IRI/cond_plots/Sham-20-24wk_composition_V2_Subclasses_P-Values.RDS")
age.pVal <- readRDS("mouse_IRI/cond_plots/Age_Group_composition_V2_Subclasses_P-Values.RDS")

#P value plot
row.order <- c(gsub(" ","",order))
row.order <- c(gsub("-","",row.order))

col.order <- c("IRI.4hrs.pVal","IRI.12hrs.pVal","IRI.2d.pVal",
               "IRI.2wk.pVal",
               "IRI.4.6wk.pVal","IRI.20.24wk.pVal","age.pVal")
p.table <- do.call(cbind, mget(ls(pattern = "\\.pVal")))
t.stats <- p.table[row.order,grep(".t", colnames(p.table))]
t.stats <- as.matrix(t(t.stats[row.order,paste0(col.order,".t")]))

# Cap the values between -10 and 10
t.stats <- pmax(pmin(t.stats, 10), -10)
write.table(p.table[row.order,], file = "mouse_IRI/cond_plots/V2_Subclasses_Composition_Analysis_p-Values.txt", sep = "\t", quote = FALSE)

p.stats <- as.matrix(t(p.table[row.order,paste0(col.order,".p")]))

load("mouse_color_factors_v2-clusters.robj")
states <- setNames(KB$v2.state.l2,KB$v2.subclass.l3)[order]
ColSideColors <- as.character(state.l2.cols[as.character(states)])

pdf(file='mouse_IRI/cond_plots/V2_Subclass_Composition_Analysis_t-statistic_Heatmap.pdf',width=30,height=5.1)
heatmap.2(t.stats,col=brewer.pal(11,"RdBu"),scale="none", trace="none", 
          density.info= "none", cexRow = 1, margins = c(8, 8), Rowv = NA, Colv = NA,
          ColSideColors = ColSideColors,
          cellnote = ifelse(p.stats < 0.05, "*", ""),
          notecol = 'black', notecex = 2)
dev.off() 

