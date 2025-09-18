library(Seurat)
library(BPCells)


args = commandArgs(trailingOnly=TRUE)
sample <- args[1]

print(paste0('../../spatial/',sample,'/outs/spatial'))
list.files(paste0('../../spatial/',sample,'/outs/spatial'))
spatial <- Load10X_Spatial(paste0('../../spatial/',sample,'/outs/'))

kernel <- c()
edge <- c()
coords <- read.csv(paste0('../positions/',sample,'_positions.csv'),)
rownames(coords) <- coords$barcode
coords <- coords[colnames(spatial),]

for (sp in 1:nrow(coords)){
  spot <- rownames(coords)[sp]
  i <- coords[sp,"array_row"]
  j <- coords[sp,"array_col"]
  neigh <- rownames(coords[coords$array_row==i-1 & coords$array_col == j-1,])
  neigh <- c(neigh,rownames(coords[coords$array_row==i-1 & coords$array_col == j+1,]))
  neigh <- c(neigh,rownames(coords[coords$array_row==i & coords$array_col == j-2,]))
  neigh <- c(neigh,rownames(coords[coords$array_row==i & coords$array_col == j+2,]))
  neigh <- c(neigh,rownames(coords[coords$array_row==i+1 & coords$array_col == j-1,]))
  neigh <- c(neigh,rownames(coords[coords$array_row==i+1 & coords$array_col == j+1,]))
  
  if (length(neigh) == 6){
    kernel <- c(kernel,spot)
  } else {
    edge <- c(edge,spot)
  }
}

spatial@meta.data[kernel,'positions'] <- 'Kernel'
spatial@meta.data[edge,'positions'] <- 'Edge'
print(SpatialDimPlot(spatial,group.by = 'positions',crop = F))
Idents(spatial) <- spatial$positions
spatial <- subset(spatial,idents = 'Kernel')
cells <- rownames(spatial@meta.data[spatial@meta.data$nCount_Spatial > 200,])
spatial <- subset(spatial,cells=cells)


spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, dims = 1:30)


atlas2 <- readRDS('../atlasv2_subset_l3_071224.RDS')

Idents(atlas2) <- atlas2@meta.data[["v2.subclass.l3"]]
anchors <- FindTransferAnchors(reference = atlas2, query = spatial, normalization.method = "SCT",query.assay='SCT',recompute.residuals = FALSE)
predictions.assay <- TransferData(anchorset = anchors, refdata = atlas2@meta.data[["v2.subclass.l3"]], prediction.assay = TRUE, 
                                  weight.reduction = spatial[["pca"]], dims = 1:30)
print(sample)
print(head(predictions.assay[,1:5]))
spatial[["predv2subclassl3"]] <- predictions.assay
print(names(spatial@assays))

atlas2 <- readRDS('../atlasv2_subset_sp_071224.RDS')

Idents(atlas2) <- atlas2@meta.data[["v2.subclass.sp"]]
anchors <- FindTransferAnchors(reference = atlas2, query = spatial, normalization.method = "SCT",query.assay='SCT',recompute.residuals = FALSE)
predictions.assay <- TransferData(anchorset = anchors, refdata = atlas2@meta.data[["v2.subclass.sp"]], prediction.assay = TRUE, 
                                  weight.reduction = spatial[["pca"]], dims = 1:30)
print(sample)
print(head(predictions.assay[,1:5]))
spatial[["predv2subclass.sp"]] <- predictions.assay
print(names(spatial@assays))

saveRDS(spatial,paste0(sample,'_atlas2_seurat5.RDS'))

DefaultAssay(spatial) <- 'predv2subclassl3'
pdf(paste0('feature_subclass_',sample,'_atlas2_seurat5.pdf'))
for (ct in unique(atlas2$v2.subclass.l3)){
  plot(SpatialFeaturePlot(spatial,ct,crop = F,pt.size.factor = 1))
}
dev.off()
DefaultAssay(spatial) <- 'predv2subclass.sp'
pdf(paste0('feature_clusters_',sample,'_atlas2_seurat5.pdf'))
for (ct in stringr::str_replace(sort(unique(atlas2$v2.subclass.sp)),'_','-')){
  plot(SpatialFeaturePlot(spatial,ct,crop = F,pt.size.factor = 1))
}
dev.off()
