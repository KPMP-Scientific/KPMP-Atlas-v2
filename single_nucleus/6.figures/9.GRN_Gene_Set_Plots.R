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


###snRNA Object
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
meta.sub <- readRDS("gene-sets/UCell_GRN_Scores_snRNA_10012024.RDS")
KB <- AddMetaData(KB, metadata = meta.sub)

# Extract mlm and store it in pathwaysmlm in KB.PT
KB[['pathwaysmlm']] <- Seurat::CreateAssayObject(as.matrix(t(KB@meta.data[,grep("_UCell", colnames(KB@meta.data))])))

# Change assay
DefaultAssay(object = KB) <- "pathwaysmlm"

# Scale the data
KB <- ScaleData(KB)
KB@assays$pathwaysmlm@data <- KB@assays$pathwaysmlm@scale.data



##Proximal Tubules
#snRNA
KB.pt <- subset(KB, v2.subclass.l3 %in% c("PT-S1","PT-S2","PT-S3",
                                           "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3"))

# Extract activities from object as a long dataframe
Idents(KB.pt) <- "v2.subclass.l3"
Idents(KB.pt) <- factor(Idents(KB.pt), levels = c("PT-S1","PT-S2","PT-S3",
                                                    "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3"))
df <- t(as.matrix(KB.pt@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.pt)) %>%
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
pdf(file='Plots/GRN/PTS1S2_GRN_Gene-set_Scores_snRNA_Heatmap.pdf',width=4,height=8)
pheatmap(t(top_acts_mat[c("PT-S1","PT-S2",
                          "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2"),grep("PTS1", colnames(top_acts_mat))]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE) 
dev.off()


subset <- c("PTS1S2.HNF4A-UCell","PTS1S2.THRB-UCell","PTS1S2.HNF1A-UCell","PTS1S2.NR2F1-UCell",
            "PTS1S2.MAF-UCell",
            
            "PTS1S2.ATF4-UCell","frPTS1S2.EGR1-UCell","frPTS1S2.ETS1-UCell",
            "frPTS1S2.GRHL2-UCell","frPTS1S2.RUNX2-UCell",
            "frPTS1S2.CEBPD-UCell","PTS1S2.HOXB7-UCell",
            "frPTS1S2.KLF6-UCell",
            
            "PTS1S2.STAT1-UCell",
            "PTS1S2.ELF3-UCell","PTS1S2.IRF1-UCell",
            "PTS1S2.SOX9-UCell","frPTS1S2.VEZF1-UCell","PTS1S2.POU3F3-UCell",
            "PTS1S2.FOXC1-UCell","PTS1S2.HIF1A-UCell","PTS1S2.HES1-UCell",
            "frPTS1S2.JUN-UCell",
            
            
            
            "frPTS1S2.RELB-UCell","frPTS1S2.NFKB1-UCell","frPTS1S2.NFKB2-UCell",
            
            "PTS1S2.BACH2-UCell",
            "PTS1S2.JUND-UCell","frPTS1S2.FOSL2-UCell","frPTS1S2.JUNB-UCell",
            "frPTS1S2.ARNT2-UCell","frPTS1S2.ATF3-UCell","PTS1S2.BHLHE40-UCell",
            "PTS1S2.MITF-UCell","PTS1S2.SOX4-UCell","PTS1S2.MYC-UCell"
            )

pdf(file='Plots/GRN/PTS1S2_GRN_Gene-set_Scores_snRNA_Heatmap_subset.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[c("PT-S1","PT-S2",
                          "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

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

pdf(file='Plots/GRN/PTS1S2_GRN_Gene-set_Scores_snRNA_Heatmap_subset2.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[c("PT-S1","PT-S2",
                          "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2"),subset2]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()






##Proximal Tubules (S3)
#snRNA
KB.pt <- subset(KB, v2.subclass.l3 %in% c("PT-S1","PT-S2","PT-S3",
                                          "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3"))

# Extract activities from object as a long dataframe
Idents(KB.pt) <- "v2.subclass.l3"
Idents(KB.pt) <- factor(Idents(KB.pt), levels = c("PT-S1","PT-S2","PT-S3",
                                                  "aPT2","aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3"))
df <- t(as.matrix(KB.pt@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.pt)) %>%
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
pdf(file='Plots/GRN/PTS3_GRN_Gene-set_Scores_snRNA_Heatmap.pdf',width=4,height=8)
pheatmap(t(top_acts_mat[c("PT-S3",
                          "aPT2","aPT1","aPT-S3","frPT-S3"),grep("PTS3", colnames(top_acts_mat))]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE) 
dev.off()


subset <- c("PTS3.NR1I3-UCell",
            "PTS3.TEF-UCell",
            "PTS3.THRB-UCell",
            "PTS3.HNF4A-UCell",
            "PTS3.HNF4G-UCell",
            
            "frPTS3.NR3C2-UCell",
            "PTS3.NFIX-UCell",
            "frPTS3.EMX2-UCell",
            "frPTS3.CEBPD-UCell",
            "frPTS3.HOXB7-UCell",
            "PTS3.MYC-UCell",
            "PTS3.MEIS1-UCell",
            "PTS3.RUNX2-UCell",
            "PTS3.KLF6-UCell",
            "PTS3.RUNX1-UCell",
            
            "PTS3.ATF4-UCell",
            "PTS3.ETS1-UCell",
            "PTS3.RELB-UCell",
            "PTS3.TBX2-UCell",
            "PTS3.SREBF2-UCell",
            "PTS3.HES1-UCell",
            "PTS3.HIF1A-UCell",
            "PTS3.ARNT2-UCell",
            "PTS3.BHLHE40-UCell",
            "PTS3.BHLHE41-UCell",
            "PTS3.EGR1-UCell",
            "PTS3.IRF1-UCell",
            "PTS3.MSANTD3-UCell",
            "PTS3.ELF3-UCell",
            "PTS3.STAT1-UCell",
            "PTS3.VEZF1-UCell",
            "PTS3.ATF3-UCell",
            "PTS3.SOX4-UCell",
            "PTS3.JUN-UCell",
            "PTS3.SOX9-UCell",
            "PTS3.MITF-UCell",
            "PTS3.NFKB1-UCell",
            
            "PTS3.FOSL2-UCell",
            "PTS3.JUNB-UCell",
            "PTS3.JUND-UCell",
            "PTS3.BACH2-UCell",
            "PTS3.NFKB2-UCell",
            "PTS3.PBX3-UCell",
            "PTS3.ZEB1-UCell",
            "PTS3.TGIF1-UCell"
            )

pdf(file='Plots/GRN/PTS3_GRN_Gene-set_Scores_snRNA_Heatmap_subset.pdf',width=4,height=7)
pheatmap(t(top_acts_mat[c("PT-S3",
                          "aPT2","aPT1","aPT-S3","frPT-S3"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()


subset2 <- c("PTS3.NR1I3-UCell",
            "PTS3.TEF-UCell",
            "PTS3.THRB-UCell",
            "PTS3.HNF4A-UCell",
            "PTS3.HNF4G-UCell",
            
            "frPTS3.NR3C2-UCell",
            "PTS3.NFIX-UCell",
            "frPTS3.CEBPD-UCell",
            "PTS3.RUNX2-UCell",
            "PTS3.KLF6-UCell",
            
            "PTS3.ATF4-UCell",
            "PTS3.HES1-UCell",
            "PTS3.ARNT2-UCell",
            "PTS3.BHLHE40-UCell",
            "PTS3.IRF1-UCell",
            "PTS3.STAT1-UCell",
            "PTS3.SOX4-UCell",
            "PTS3.SOX9-UCell",
            "PTS3.MITF-UCell",
            "PTS3.NFKB1-UCell",
            
            "PTS3.JUND-UCell",
            "PTS3.ZEB1-UCell",
            "PTS3.TGIF1-UCell"
)

pdf(file='Plots/GRN/PTS3_GRN_Gene-set_Scores_snRNA_Heatmap_subset2.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[c("PT-S3",
                          "aPT2","aPT1","aPT-S3","frPT-S3"),subset2]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()






##TAL Tubules
KB.tal <- subset(KB, v2.subclass.l3 %in% c("M-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B",
                                           "C-TAL-B","aTAL1","aTAL2","frTAL"))

# Extract activities from object as a long dataframe
Idents(KB.tal) <- "v2.subclass.l3"
Idents(KB.tal) <- factor(Idents(KB.tal), levels = c("M-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B",
                                                    "C-TAL-B","aTAL1","aTAL2","frTAL"))
df <- t(as.matrix(KB.tal@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.tal)) %>%
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
pdf(file='Plots/GRN/CM-TAL_GRN_Gene-set_Scores_snRNA_Heatmap.pdf',width=4,height=9)
pheatmap(t(top_acts_mat[c("C/M-TAL-A","C-TAL-A","C/M-TAL-B",
                          "C-TAL-B","aTAL1","aTAL2","frTAL"),grep("TAL", colnames(top_acts_mat))]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE) 
dev.off()


subset <- c("CMTALA.MITF-UCell", "CMTALA.ESRRB-UCell",
            "CMTALA.RORA-UCell", "CMTALA.NR3C1-UCell", "CMTALA.TEF-UCell",
            
            "CMTALA.FOSL2-UCell", "CMTALA.BACH2-UCell", "CMTALA.JUNB-UCell", "CMTALA.CUX1-UCell",
            "CMTALA.HOXB4-UCell", "CMTALA.BACH1-UCell", "CMTALA.KLF6-UCell", "CMTALA.GLIS3-UCell",
            "CMTALA.MAFF-UCell", "CMTALA.RELB-UCell", "CMTALA.PKNOX1-UCell",
            "CMTALA.SOX4-UCell", "CMTALA.MSANTD3-UCell", "CMTALA.ZEB1-UCell", "CMTALA.SMAD3-UCell",
            "frTAL.JUN-UCell", "CMTALA.NFAT5-UCell", "CMTALA.TEAD1-UCell", "CMTALA.TEAD4-UCell",
            "CMTALA.E2F3-UCell", "CMTALA.HMBOX1-UCell", "CMTALA.KLF3-UCell", "CMTALA.PPARG-UCell",
            "CMTALA.ATF7-UCell", "CMTALA.CREM-UCell", "CMTALA.TCF7L1-UCell", "CMTALA.TGIF1-UCell",
            "CMTALA.HNF1B-UCell", "CMTALA.HOXB3-UCell",
            "CMTALA.RFX7-UCell", "CMTALA.ETS1-UCell", "CMTALA.SP1-UCell", "CMTALA.HIF1A-UCell",
            "CMTALA.NFKB1-UCell", "CMTALA.SOX9-UCell", "CMTALA.BCL6-UCell", "CMTALA.ELK3-UCell",
            "CMTALA.ETS2-UCell", 
            
            "CMTALA.ELK4-UCell", "CMTALA.ETV5-UCell",
            "CMTALA.ETV1-UCell", "CMTALA.ETV6-UCell", "CMTALA.ELF1-UCell", 
            "CMTALA.IRF1-UCell","CMTALA.IRF3-UCell",
            
            "CMTALA.STAT1-UCell", "CMTALA.STAT3-UCell",
            "CMTALA.ELF2-UCell", "CMTALA.RFX2-UCell","CMTALA.CEBPB-UCell", "CMTALA.CEBPD-UCell", 
            
            "frTAL.EGR1-UCell", "frTAL.HOXB9-UCell",
            "frTAL.HOXC9-UCell", "frTAL.EHF-UCell", "frTAL.ELF3-UCell", "frTAL.KLF13-UCell",
            "frTAL.KLF2-UCell", "frTAL.KLF5-UCell")

pdf(file='Plots/GRN/CM-TAL_GRN_Gene-set_Scores_snRNA_Heatmap_subset.pdf',width=4,height=8)
pheatmap(t(top_acts_mat[c("C/M-TAL-A","C-TAL-A","C/M-TAL-B",
                          "C-TAL-B","aTAL1","aTAL2","frTAL"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()

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

pdf(file='Plots/GRN/CM-TAL_GRN_Gene-set_Scores_snRNA_Heatmap_subset2.pdf',width=4,height=5)
pheatmap(t(top_acts_mat[c("C/M-TAL-A","C-TAL-A","C/M-TAL-B",
                          "C-TAL-B","aTAL1","aTAL2","frTAL"),subset2]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()







##Interstitial Fibroblasts 
KB.int <- subset(KB, v2.subclass.l3 %in% c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo",
                                          "C-FIB-OSMRhi","C-MYOF"))

# Extract activities from object as a long dataframe
Idents(KB.int) <- "v2.subclass.l3"
Idents(KB.int) <- factor(Idents(KB.int), levels = c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo",
                                                    "C-FIB-OSMRhi","C-MYOF"))
df <- t(as.matrix(KB.int@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.int)) %>%
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
pdf(file='Plots/GRN/Int-FIB_GRN_Gene-set_Scores_snRNA_Heatmap.pdf',width=4,height=8)
pheatmap(t(top_acts_mat[c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo",
                          "C-FIB-OSMRhi","C-MYOF"),c(grep("CFIBOSMRhi", colnames(top_acts_mat)),grep("CMYOF", colnames(top_acts_mat)))]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE) 
dev.off()


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

pdf(file='Plots/GRN/Int-FIB_GRN_Gene-set_Scores_snRNA_Heatmap_subset.pdf',width=4,height=6)
pheatmap(t(top_acts_mat[c("C-FIB","C-FIB-PATH","C-FIB-OSMRlo",
                          "C-FIB-OSMRhi","C-MYOF"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()







##Perivascular Fibroblasts 
KB.pv <- subset(KB, v2.subclass.l3 %in% c("pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))

# Extract activities from object as a long dataframe
Idents(KB.pv) <- "v2.subclass.l3"
Idents(KB.pv) <- factor(Idents(KB.pv), levels = c("pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF"))
df <- t(as.matrix(KB.pv@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.pv)) %>%
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
pdf(file='Plots/GRN/pvFIB_GRN_Gene-set_Scores_snRNA_Heatmap.pdf',width=4,height=8)
pheatmap(t(top_acts_mat[c("pvFIB-RSPO3+","pvFIB-PI16+",
                          "pvFIB","pvMYOF"),grep("pv", colnames(top_acts_mat))]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE) 
dev.off()


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

pdf(file='Plots/GRN/pvFIB_GRN_Gene-set_Scores_snRNA_Heatmap_subset.pdf',width=4,height=6)
pheatmap(t(top_acts_mat[c("pvFIB-RSPO3+","pvFIB-PI16+",
                          "pvFIB","pvMYOF"),subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()






###Myeloid Cells
KB.mye <- subset(KB, v2.subclass.l3 %in% c("resMAC-LYVE1+","resMAC-HLAIIhi",
                                           "MON","moMAC-HBEGF+","moMAC-CXCL10+",
                                             "moFAM","moMAC-C3+"))

# Extract activities from object as a long dataframe
Idents(KB.mye) <- "v2.subclass.l3"
Idents(KB.mye) <- factor(Idents(KB.mye), levels = c("resMAC-LYVE1+","resMAC-HLAIIhi",
                                                  "MON","moMAC-HBEGF+","moMAC-CXCL10+",
                                                  "moFAM","moMAC-C3+"))
df <- t(as.matrix(KB.mye@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(KB.mye)) %>%
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
pdf(file='Plots/GRN/Myeloid_GRN_Gene-set_Scores_snRNA_Heatmap.pdf',width=4,height=8)
pheatmap(t(top_acts_mat[,c(grep("resMAC", colnames(top_acts_mat)),grep("MONCXCL10", colnames(top_acts_mat)),grep("moFAM", colnames(top_acts_mat)))]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE) 
dev.off()

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

pdf(file='Plots/GRN/Myeloid_GRN_Gene-set_Scores_snRNA_Heatmap_subset.pdf',width=5,height=6)
pheatmap(t(top_acts_mat[,subset]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = FALSE, cluster_rows = FALSE) 
dev.off()



