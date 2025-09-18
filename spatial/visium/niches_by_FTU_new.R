library(Seurat)
library(ggplot2)
library(randomcoloR)

spatial <- readRDS('../data/new_samples/stmerged_081224.RDS')

annot <- rio::import_list('../data/Supplementary Tables Version 1 05-2024.xlsx')
annot <- annot[["Table S - sn markers"]]
rownames(annot) <- annot$v2.subclass.l3

scores <- as.data.frame(reshape2::melt(spatial@assays$predv2subclassl3@data))
scores <- scores[scores$Var1 != 'max',]
colnames(scores) <- c('subclass.l3','spot','proportion')
scores$subclass.l3 <- as.character(scores$subclass.l3)
scores$spot <- as.character(scores$spot)
scores$subclass.l1 <- annot[scores$subclass.l3,"v2.subclass.l1"]
scores$class <- annot[scores$subclass.l3,"v2.class"]

classdf <- aggregate(proportion ~ class + spot, data=scores, FUN=sum)
classdf <- reshape2::dcast(classdf, spot ~ class, value.var='proportion')
rownames(classdf) <- classdf$spot
classdf <- classdf[,-1]

class <- apply(classdf, 1, function(x) names(x)[which.max(x)])
classdf$max <- apply(classdf, 1, max)
classdf$class <- as.character(class)

subclassdf <- aggregate(proportion ~ subclass.l1 + spot, data=scores, FUN=sum)
subclassdf <- reshape2::dcast(subclassdf, spot ~ subclass.l1, value.var='proportion')
rownames(subclassdf) <- subclassdf$spot
subclassdf <- subclassdf[,-1]

class <- apply(subclassdf, 1, function(x) names(x)[which.max(x)])
subclassdf$max <- apply(subclassdf, 1, max)
subclassdf$subclass.l1 <- as.character(class)

classdf$group <- ifelse(classdf$class != 'epithelial cells' & classdf$max > 0.2, classdf$class, 'others')
classdf$group <- ifelse(classdf$class == 'epithelial cells' & subclassdf$max > 0.2, subclassdf$subclass.l1, classdf$group)

smallgrp <- names(table(classdf$group)[table(classdf$group) < 100])
classdf$group <- ifelse(classdf$group %in% smallgrp, 'others', classdf$group)

spatial$group <- classdf[rownames(spatial@meta.data),'group']
spatial$group <- ifelse(spatial$group == 'EC','endothelial cells',spatial$group)
spatial$group <- ifelse(spatial$group == 'FIB','stroma cells',spatial$group)
spatial$group <- ifelse(spatial$group == 'Lymphoid','immune cells',spatial$group)
spatial$group <- ifelse(spatial$group == 'Myeloid','immune cells',spatial$group)


spatial$new_group <- ifelse(spatial$group %in% c('POD','PEC'),'Glomerular',
                        ifelse(spatial$group %in% c('DCT','CNT'),'DCT_CNT',
                        ifelse(spatial$group %in% c('PC','IC'),'CD',spatial$group)))
spatial$new_group <- ifelse(spatial$new_group %in% 
                        c('Glomerular',"PT", "DTL", "TAL","DCT_CNT","CD","endothelial cells","stroma cells","immune cells"),
                        spatial$new_group,'other')

unique(spatial$new_group)
spatial$niches_new_group <- spatial$niches_group
Idents(spatial) <- 'new_group'
assay <- 'predv2subclassl3'
DefaultAssay(spatial) <- assay
res <- 2
pdf('../figures/new_samples/niches_by_new_group_FTU.pdf',width = 20,height = 15)
for (group in c('Glomerular','DCT_CNT','CD')) {
    # group  <-  'PT'
    print(group)
    sptgrp <- subset(spatial, idents=group)
    group <- stringr::str_replace_all(group,' ','.')
    group <- stringr::str_replace_all(group,'/','.')
    VariableFeatures(sptgrp) <- rownames(sptgrp@assays[[assay]])
    sptgrp <- ScaleData(sptgrp)
    sptgrp <- RunPCA(sptgrp, verbose = FALSE,reduction.name = paste0('pca_',group)) #should I be setting approx = FALSE to run a full PCA instead
    sptgrp <- FindNeighbors(sptgrp, dims = 1:50,reduction = paste0('pca_',group))
    sptgrp <- FindClusters(sptgrp, verbose = FALSE,resolution = res,cluster.name = paste0(group,'_',assay,'_res',res))
    sptgrp <- RunUMAP(sptgrp, dims = 1:50, reduction = paste0('pca_',group), verbose = FALSE,umap.name = paste0('umap_',group))
    reduction_niches[[group]] <- sptgrp@reductions[['umap']]
    plot(DotPlot(sptgrp,features = unique(annot$v2.subclass.l3),group.by = paste0(group,'_',assay,'_res',res))+
         ggtitle(paste0(group,'_res',res))+
         theme(legend.position = 'left',
         axis.text.x = element_text(angle=90,vjust = .5)))
    spatial@meta.data[rownames(sptgrp@meta.data),'niches_new_group'] <- as.numeric(as.character(sptgrp$seurat_clusters))
}
dev.off()

saveRDS(reduction_niches, '../data/new_samples/reduction_niches_new_group.RDS')

pdf('../figures/new_samples/niches_by_FTU_order.pdf',width = 20,height = 15)
for (group in names(sort(table(spatial$group),decreasing = T))) {
    sptgrp <- subset(spatial, idents=group)
    d <- DotPlot(sptgrp,features = unique(annot$v2.subclass.l3),group.by = 'niches_group')+
         ggtitle(paste0(group,'_res',res))+
         theme(legend.position = 'left',
         axis.text.x = element_text(angle=90,vjust = .5))

    mat <- reshape2::dcast(d$data,features.plot~id,value.var = 'avg.exp.scaled')
    rownames(mat) <- mat$features.plot
    mat <- mat[,-1]

    ser <- seriation::seriate(mat,margin=2L,method='PCA_angle')
    mat2 <- seriation::permute(mat,ser)
    
    plot(DotPlot(sptgrp,features = unique(annot$v2.subclass.l3),group.by = 'niches_group')+
            ggtitle(paste0(group,'_res',res))+
            theme(legend.position = 'left',
            axis.text.x = element_text(angle=90,vjust = .5))+
            scale_y_discrete(limits = colnames(mat2)))

}
dev.off()

annot$new_group <- ifelse(annot$v2.subclass.l1 %in% c('POD','PEC'),'Glomerular',
                        ifelse(annot$v2.subclass.l1 %in% c('DCT','CNT'),'DCT_CNT',
                        ifelse(annot$v2.subclass.l1 %in% c('PC','IC'),'CD',annot$v2.subclass.l1)))
annot$new_group <- ifelse(annot$v2.class != c('epithelial cells'),annot$v2.class,annot$new_group)
annot$new_group <- ifelse(annot$new_group %in% 
                        c('Glomerular',"PT", "DTL", "TAL","DCT_CNT","CD","endothelial cells","stroma cells","immune cells"),
                        annot$new_group,'other')
table(annot[,c('new_group','v2.state.l1')])
Idents(spatial) <- 'new_group'
spatial$state <- 'none'
for (grp in unique(spatial$new_group)[c(1:3,5:10,13)]){ #exclude others
    #grp <- 'PT'
    subannot <- annot[annot$new_group == grp,]
    subspatial <- subset(spatial, idents=grp)
    data <- subspatial@assays$predv2subclassl3@data[rownames(subannot),]
    data <- as.data.frame(reshape2::melt(as.matrix(data)))
    colnames(data) <- c('subclass.l3','spot','proportion')
    data$v2.state.l1 <- subannot[match(data$subclass.l3,subannot$v2.subclass.l3),'v2.state.l1']
    data <- aggregate(proportion ~ v2.state.l1 + spot, data=data, FUN=sum)
    data <- reshape2::dcast(data, spot ~ v2.state.l1, value.var='proportion')
    rownames(data) <- data$spot
    data$state <- apply(data[,-1], 1, function(x) names(x)[which.max(x)])
    spatial@meta.data[rownames(data),'state'] <- data$state
}
table(spatial@meta.data[,c('new_group','state')])
saveRDS(spatial, '../data/new_samples/spatial_niches_new_group_FTU.RDS')


