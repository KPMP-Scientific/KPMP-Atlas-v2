library(Seurat)
library(ggplot2)
library(patchwork)

spatial <- readRDS('../data/new_samples/spatial_niches_new_group_FTU.RDS')
reduction_niches <- readRDS('../data/new_samples/reduction_niches_new_group.RDS')
names(reduction_niches)

annot <- rio::import_list('../data/Supplementary Tables Version 1 05-2024.xlsx')
annot <- annot[["Table S - sn markers"]]
rownames(annot) <- annot$v2.subclass.l3
head(annot)
unique(annot$v2.class)

data <- reshape2::melt(as.matrix(spatial@assays$predv2subclassl3@data))
data$altered <- ifelse(annot[match(data$Var1,rownames(annot)),'v2.state.l1'] == 'altered',TRUE,FALSE)
data$stromal <- ifelse(annot[match(data$Var1,rownames(annot)),'v2.class'] == 'stroma cells',TRUE,FALSE)
data$immune <- ifelse(annot[match(data$Var1,rownames(annot)),'v2.class'] == 'immune cells',TRUE,FALSE)
head(data)
tempdf <- aggregate(value ~ Var2 + altered, data = data, FUN = sum)
spatial$altered_prop <- tempdf[tempdf$altered == TRUE, 'value'][match(rownames(spatial@meta.data), tempdf$Var2)]
tempdf <- aggregate(value ~ Var2 + immune, data = data, FUN = sum)
spatial$immune_prop <- tempdf[tempdf$immune == TRUE, 'value'][match(rownames(spatial@meta.data), tempdf$Var2)]
tempdf <- aggregate(value ~ Var2 + stromal, data = data, FUN = sum)
spatial$stromal_prop <- tempdf[tempdf$stromal == TRUE, 'value'][match(rownames(spatial@meta.data), tempdf$Var2)]

condition_meta <- read.csv('../data/new_samples/condition_meta.csv')
condition_meta[condition_meta$Patient == 'IU-REF',c("condition_level3","condition_level2","condition_level1")] <- 'RT-UCS'
condition_meta <- unique(condition_meta[,c("Patient","subject","spot","condition_level3","condition_level2","condition_level1")])
rownames(condition_meta) <- condition_meta$spot

spatial@meta.data[,c("condition_level3","condition_level2","condition_level1")] <- condition_meta[rownames(spatial@meta.data),c("condition_level3","condition_level2","condition_level1")]

groups <- c(
  "Glomerular", "PT", "DTL", "TAL", "ATL", "DCT_CNT",
  "CD", "endothelial cells", "stroma cells", "immune cells"
)

grpspatial <- list()
for (group in groups) {
  #group <- groups[1]
  print(group)
  sptgrp <- subset(spatial, idents=group)
  sptgrp@reductions[['umap']] <- reduction_niches[[stringr::str_replace(group, " ", ".")]]
  grpspatial[[group]] <- sptgrp
}
saveRDS(grpspatial,'../data/new_samples/spatial_group_objects_niches_new_group_FTU.RDS')
grpspatial <- readRDS('../data/new_samples/spatial_group_objects_niches_new_group_FTU.RDS')

colnames(grpspatial[[1]]@meta.data)
lapply(grpspatial,function(x) {length(unique(x@meta.data$niches_new_group))})
ncol

umaplist <- list()
features <- c('str_to_altered','imm_to_altered')
for (ft in features){
    for (group in groups[c(1:4,6:10)]) {
        #group <- groups[1]
        print(group)
        sptgrp <- grpspatial[[group]]
        sptgrp$imm_to_altered <- sptgrp$immune_prop / (sptgrp$altered_prop + 1e-6)
        sptgrp$str_to_altered <- sptgrp$stromal_prop / (sptgrp$altered_prop + 1e-6)

        p <- FeaturePlot(sptgrp, reduction = 'umap', features = ft,order=T) +
                scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1), oob = scales::squish) +
                ggtitle(group) +
                theme_void() +
                theme(
                    plot.title.position = "plot",
                    plot.title = element_text(hjust = 0.5, vjust = 1, size = 14)
                )
        umaplist[[paste0(group,'_',ft)]] <- p
    }

}
wrap_plots(umaplist, ncol = 9, guides = "collect") & theme(legend.position = "right")
ggsave(paste0('../figures/niches_umap_all_row_ratio_capped.pdf'), width = 36.2, height = 8)
