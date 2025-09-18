library(Seurat)
library(SeuratWrappers)
library(batchelor)

setwd('R')

list_sample <- list.files('../map_v2_l3_sp_seurat5/','*RDS')
list_sample <- stringr::str_remove(list_sample,'_atlas2_seurat5.RDS')

samples <- data.frame(list(samples=list_sample,
                           iu_id=stringr::str_split(list_sample,'_',simplify = T)[,3]))

spatial_list <- list()
pdf('remove_spots.pdf')
for (r in 1:nrow(samples)){
  smp <- samples[r,"samples"]
  spatial <- readRDS(paste0('../map_v2_l3_sp_seurat5/',smp,'_atlas2_seurat5.RDS'))
  spatial@meta.data$orig.ident <- smp
  spatial@meta.data$subject <- samples[r,"iu_id"]
  spatial@meta.data$condition <- samples[r,"condition"]
  
  
  print(SpatialDimPlot(spatial,group.by = 'positions',crop = F))
  spatial_list[[smp]] <- spatial
}
dev.off()

length(spatial_list)

stmerged <- merge(spatial_list[[1]],spatial_list[2:153])
names(stmerged@images) <- stringr::str_replace_all(samples$sample,'-','.')

DefaultAssay(stmerged) <- 'Spatial'
saveRDS(stmerged,'../data/new_samples/stmerged_081224.RDS')
