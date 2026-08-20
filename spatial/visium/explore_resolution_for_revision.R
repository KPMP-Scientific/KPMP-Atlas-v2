library(Seurat)
library(ggplot2)

reduction_niches <- readRDS('E:/atlasV2_data/new_samples/reduction_niches_new_group.RDS')
spatial <- readRDS('E:/atlasV2_data/new_samples/spatial_niches_new_group_FTU.RDS')

annot <- rio::import('E:/atlasV2_data/Supplementary Tables Version 1 05-2024.xlsx',
                     sheet = "Table S12 - sn markers", range = 'B5:N133')
rownames(annot) <- annot$v2.subclass.l3

#spatial <- UpdateSeuratObject(spatial)
unique(spatial$new_group)

assay <- 'predv2subclassl3'
DefaultAssay(spatial) <- assay

for (res in c(0.05, 0.1, 0.2, 0.4, 0.5, 0.8, 1, 1.5, 2)){ 
pdf(paste0('../figures/new_samples/niches_by_new_group_FTU_res_',res,'.pdf'),width = 20,height = 15)
for (group in unique(spatial$new_group)) {
  # group  <-  'PT'
  print(c(group,res))
  cells <- rownames(spatial@meta.data)[spatial$new_group == group]
  #sptgrp <- subset(spatial, idents=group)
  sptgrp <- subset(spatial, cells = cells)
  group <- stringr::str_replace_all(group,' ','.')
  group <- stringr::str_replace_all(group,'/','.')
  VariableFeatures(sptgrp) <- rownames(sptgrp@assays[[assay]])
  sptgrp <- ScaleData(sptgrp)
  sptgrp <- RunPCA(sptgrp, verbose = FALSE,reduction.name = paste0('pca_',group,'_',res)) #should I be setting approx = FALSE to run a full PCA instead
  sptgrp <- FindNeighbors(sptgrp, dims = 1:50,reduction = paste0('pca_',group,'_',res))
  sptgrp <- FindClusters(sptgrp, verbose = FALSE,resolution = res,cluster.name = paste0(group,'_',assay,'_res',res))
  sptgrp <- RunUMAP(sptgrp, dims = 1:50, reduction = paste0('pca_',group,'_',res), verbose = FALSE,umap.name = paste0('umap_',group,'_',res))
  reduction_niches[[paste0(group,'_',res)]] <- sptgrp@reductions[[paste0('umap_',group,'_',res)]]
  plot(DotPlot(sptgrp,features = unique(annot$v2.subclass.l3),group.by = paste0(group,'_',assay,'_res',res))+
         ggtitle(paste0(group,'_res',res))+
         theme(legend.position = 'left',
               axis.text.x = element_text(angle=90,vjust = .5)))
  spatial@meta.data[rownames(sptgrp@meta.data),paste0('niches_new_group_res_',res)] <- as.numeric(as.character(sptgrp$seurat_clusters))
}
dev.off()
}
saveRDS(reduction_niches, 'E:/atlasV2_data/new_samples/reduction_niches_new_group_many_res.RDS')
saveRDS(spatial, 'E:/atlasV2_data/new_samples/spatial_niches_new_group_FTU_many_res.RDS')


cells <- rownames(spatial@meta.data)[spatial$new_group == 'PT']
pts <- subset(spatial, cells = cells)
colnames(pts@meta.data)

pdf('../figures/new_samples/dotplot_niches_by_new_group_FTU_manyres.pdf',width = 20,height = 10)
for (res in c(0.05, 0.1, 0.2, 0.4, 0.5, 0.8, 1, 1.5, 2)){ 
  plot(DotPlot(pts,features=annot$v2.subclass.l3,group.by=paste0('niches_new_group_res_',res))+
         theme(axis.text.x = element_text(angle=90,vjust=0.5)))
}
dev.off()

library(clustree)
unique(annot$v2.subclass.l1)

dotplot_average <- function(x) {
  log1p(mean(expm1(x), na.rm = TRUE))
}

pdf('../figures/new_samples/clustree_lymphoid_niches_by_new_group_FTU_manyres.pdf',width = 15,height = 8)
for (ct in annot[annot$v2.subclass.l1 %in% c("Lymphoid"),'v2.subclass.l3']){
p <- clustree(
  pts,
  prefix = "niches_new_group_res_",
  node_colour = ct,
  node_colour_aggr = "dotplot_average"
)
plot(p)
}
dev.off()

pdf('../figures/new_samples/clustree_myeloid_niches_by_new_group_FTU_manyres.pdf',width = 15,height = 8)
for (ct in annot[annot$v2.subclass.l1 %in% c("Myeloid"),'v2.subclass.l3']){
  p <- clustree(
    pts,
    prefix = "niches_new_group_res_",
    node_colour = ct,
    node_colour_aggr = "dotplot_average"
  )
  plot(p)
}
dev.off()

pdf('../figures/new_samples/clustree_fib_niches_by_new_group_FTU_manyres.pdf',width = 15,height = 8)
for (ct in annot[annot$v2.subclass.l1 %in% c("FIB"),'v2.subclass.l3']){
  p <- clustree(
    pts,
    prefix = "niches_new_group_res_",
    node_colour = ct,
    node_colour_aggr = "dotplot_average"
  )
  plot(p)
}
dev.off()

pdf('../figures/new_samples/clustree_vsm-p_niches_by_new_group_FTU_manyres.pdf',width = 15,height = 8)
for (ct in annot[annot$v2.subclass.l1 %in% c("VSM/P"),'v2.subclass.l3']){
  p <- clustree(
    pts,
    prefix = "niches_new_group_res_",
    node_colour = ct,
    node_colour_aggr = "dotplot_average"
  )
  plot(p)
}
dev.off()

p <- clustree(
  pts,
  prefix = "niches_new_group_res_",
  node_colour = ct,
  node_colour_aggr = "mean"
)
plot(p)