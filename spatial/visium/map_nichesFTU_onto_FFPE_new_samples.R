library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 5000 * 1024^2)

spatial <- readRDS('../data/new_samples/spatial_niches_FTU.RDS')
spatial$fullniche <- paste0(spatial$group, '_', spatial$niches_group)
DefaultAssay(spatial) <- 'Spatial'
JoinLayers(spatial)
#spatial[['SpatialV3']] <- as(object = spatial[["Spatial"]], Class = "Assay")
#spatial <- SCTransform(spatial,assay = 'SpatialV3',new.assay.name='newSCT')
spatial <- SCTransform(spatial,assay = 'Spatial')

smp <- 'V42D20-364_XY04_2240610'
ffpe <- Load10X_Spatial(paste0()'../../spatial_samples/',smp,'/outs'))
ffpe@meta.data$orig.ident <- smp
ffpe <- SCTransform(ffpe,assay = 'Spatial')
ffpe <- RunPCA(ffpe)


Idents(spatial) <- spatial$fullniche
DefaultAssay(spatial) <- 'newSCT'


anchors <- FindTransferAnchors(reference = spatial, query = ffpe, normalization.method = "SCT",query.assay='SCT',recompute.residuals = FALSE)
predictions.assay <- TransferData(anchorset = anchors, refdata = spatial$fullniche, prediction.assay = TRUE, 
                                  weight.reduction = ffpe[["pca"]], dims = 1:30)
ffpe[["predniches"]] <- predictions.assay
DefaultAssay(ffpe) <- 'predniches'

saveRDS(ffpe, paste0('../data/ffpe_niches/',smp,'.RDS'))