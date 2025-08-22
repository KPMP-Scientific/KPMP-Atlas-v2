library(Seurat)
library(SeuratDisk)
library(patchwork)
library(Matrix)
library(ggplot2)
library(igraph)
library(BPCells)
library(decoupleR)
library(tibble)
library(tidyr)
library(pheatmap)
library(dplyr)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")
source("misc/utils.R")

# Read in the raw counts data
raw_counts <- read.table('in_vitro/ILMN_2401_Eadon_totalRNAseq15_Jan2025_STAR_featureCounts_raw-counts.txt', header=TRUE, row.names=1)

library('biomaRt')
library(DESeq2)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genelist <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  values = rownames(raw_counts), mart = mart)

raw_counts$stableID <-  gsub("\\..*", "", rownames(raw_counts))
raw_counts$annot.symbol <- genelist[match(raw_counts$stableID,genelist$ensembl_gene_id),2]
raw_counts$annot.symbol[is.na(raw_counts$annot.symbol)] <- raw_counts$stableID[is.na(raw_counts$annot.symbol)]

#keep only counts associated with gene symbols
genes <- unique(genelist$hgnc_symbol[!is.na(genelist$hgnc_symbol)])
genes <- genes[!genes %in% ""]
raw_counts <- raw_counts[raw_counts$annot.symbol %in% genes,]
raw_counts <- raw_counts[!duplicated(raw_counts$annot.symbol),]
rownames(raw_counts) <- raw_counts$annot.symbol
colnames <- c("Neg.siRNA_S23", "Neg.siRNA_S24", "Neg.siRNA_S25", "No.siRNA_S36", "No.siRNA_S37", "No.siRNA_S38", "siR.SOX.4plusSOX.9_S32",
              "siR.SOX.4plusSOX.9_S33", "siR.SOX.4plusSOX.9_S34", "siR.SOX.4_S26", "siR.SOX.4_S27", "siR.SOX.4_S28", "siR.SOX.9_S29", "siR.SOX.9_S30",
              "siR.SOX.9_S31")
raw_counts <- raw_counts[,colnames]

# Create a sample information dataframe (modify according to your experimental design)
sample_info <- data.frame(
  row.names = colnames(raw_counts),
  condition = c("Neg.siRNA", "Neg.siRNA", "Neg.siRNA", "No.siRNA", "No.siRNA", "No.siRNA",
                "siR.SOX.4plusSOX.9", "siR.SOX.4plusSOX.9", "siR.SOX.4plusSOX.9",
                "siR.SOX.4", "siR.SOX.4", "siR.SOX.4", "siR.SOX.9", "siR.SOX.9",
                "siR.SOX.9")
)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = sample_info, design = ~ condition)

# Normalize the data
dds <- DESeq(dds)
counts <- counts(dds, normalized=TRUE)


###GRN DEGs
###Full frPT trajectory GRN scores
load("~/hsKidAt/blake_LTS/Atlas_V2/scratch/regulatory/scMEGA_aPT-S1-2_Trajectory_GRN_0424-newData.RDA")

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
  subset(correlation > 0.6 | correlation < -0.6) %>%  
  dplyr::select(c(tf, gene, correlation, p_value)) %>%
  dplyr::rename(source = tf, target = gene, weight = correlation)

df.grn2 <- df.grn2[df.grn2$source %in% genes,]
sample_acts <- run_mlm(mat=counts, net=df.grn2, .source='source', .target='target',
                       .mor='weight', minsize = 5)

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

# Plot
pdf(file='in_vitro/frPT-TF-Gene-set_Scores_bulk-siRNA_04-2025_Heatmap_C.pdf',width=4,height=5)
pheatmap(t(sample_acts_mat[c("Neg.siRNA_S23", "Neg.siRNA_S24", "Neg.siRNA_S25", 
                             "siR.SOX.4_S26", "siR.SOX.4_S27", "siR.SOX.4_S28"),
                           ]),
         border_color = NA, color=my_color, breaks = my_breaks,
         clustering_method = "ward.D", cluster_cols = TRUE, cluster_rows = TRUE) 
dev.off()

pdf(file='in_vitro/frPT-TF-Gene-set_Scores_bulk-siRNA_04-2025_Heatmap_D.pdf',width=4,height=5)
pheatmap(t(sample_acts_mat[c("Neg.siRNA_S23", "Neg.siRNA_S24", "Neg.siRNA_S25", 
                             "siR.SOX.4_S26", "siR.SOX.4_S27", "siR.SOX.4_S28"),
                           c("SOX4","SOX9","ETS1","GRHL2","FOSL2","HNF4G","THRB","KLF6","NFKB1",
                             "ATF3","POU3F3","RELB","HNF1A","HNF4A",
                             "EGR1","ELF3","MAF","JUNB","FOXC1")]),
border_color = NA, color=my_color, breaks = my_breaks,
clustering_method = "ward.D", cluster_cols = TRUE, cluster_rows = TRUE) 
dev.off()

write.table(t(sample_acts_mat[c("Neg.siRNA_S23", "Neg.siRNA_S24", "Neg.siRNA_S25", 
                                "siR.SOX.4_S26", "siR.SOX.4_S27", "siR.SOX.4_S28"),
]), file = "in_vitro/frPT-TF-Gene-set_Scores_bulk-siRNA_04-2025.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


