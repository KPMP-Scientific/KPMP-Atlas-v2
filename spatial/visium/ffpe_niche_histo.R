library(Seurat)
library(ggplot2)
library(patchwork)

list_smp <- list.files('../map_ffpe_niche', pattern = '.RDS')
list_smp <- stringr::str_remove(list_smp, '_ffpe_niches.RDS')

list_spatial <- list()
for (smp in list_smp){
    #smp <- list_smp[1]
    spatial <- readRDS(paste0('../map_ffpe_niche/',smp,'_ffpe_niches.RDS'))
    DefaultAssay(spatial) <- 'predniches'
    histo <- read.csv(paste0('../ffpe_transfer/loupe_annotation/new/',smp,'_ME_Histo.csv'))
    rownames(histo) <- histo[[1]]
    spatial$Histo_annotation <- histo[rownames(spatial@meta.data),2]

    kernel <- c()
    edge <- c()
    coords <- read.csv(paste0('../positions/',smp,'_positions.csv'),)
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


    list_spatial[[smp]] <- spatial
}

stmerged <- merge(list_spatial[[1]],list_spatial[2:length(list_spatial)])
names(stmerged@images) <- stringr::str_replace_all(list_smp,'-','.')

stmerged$Histo_annotation <- ifelse(stmerged$Histo_annotation == 'MED_TI', 'Med_TI',stmerged$Histo_annotation)
stmerged$Histo_annotation <- ifelse(stmerged$Histo_annotation == 'SCL_GLOM', 'Scl_glom',stmerged$Histo_annotation)
stmerged$Histo_annotation <- ifelse(stmerged$Histo_annotation == 'GLOM', 'Glom',stmerged$Histo_annotation)
stmerged$Histo_annotation <- ifelse(stmerged$Histo_annotation == 'FIB1', 'FIB',stmerged$Histo_annotation)
stmerged$Histo_annotation <- ifelse(stmerged$Histo_annotation == 'IMM1', 'IMM',stmerged$Histo_annotation)
stmerged$Histo_annotation <- ifelse(stmerged$Histo_annotation == 'IGNORE', 'Ignore',stmerged$Histo_annotation)
table(stmerged$Histo_annotation)
stmerged$Histo_annotation <- ifelse(stmerged$positions == 'Edge', 'Ignore',stmerged$Histo_annotation)
table(stmerged$Histo_annotation)

# saveRDS(stmerged,'../data/ffpe_niche_histo.RDS')
# stmerged <- readRDS('../data/ffpe_niche_histo.RDS')  

head(stmerged@meta.data)
for (smp in unique(stmerged$orig.ident)){
    #smp <- unique(stmerged$orig.ident)[1]
    annot <- stmerged@meta.data[stmerged$orig.ident == smp,'Histo_annotation',drop = FALSE]
    annot$Barcode <- stringr::str_split(rownames(annot),'_',simplify = TRUE)[,1]
    write.csv(annot,paste0('../ffpe_transfer/loupe_annotation/new/',smp,'_final_Histo.csv'),
              row.names = FALSE,quote = FALSE)
}    

stmerged <- subset(stmerged,subset = Histo_annotation != 'Ignore')


data <- as.data.frame(stmerged@assays$predniches@data)
data <- data[1:(nrow(data)-1),]
data <- reshape2::melt(as.matrix(data),id.vars = 'l1')
data$l1 <- stringr::str_split(data$Var1,'-',simplify = TRUE)[,1]
head(data)
data <- aggregate(data$value,by = list(data$l1,data$Var2),FUN = sum)
data <- reshape2::dcast(data,Group.2~Group.1)
data[1:5,1:5]
rownames(data) <- data[,1]
data <- data[,-1]
main <- colnames(data)[max.col(data)]
main <- as.data.frame(cbind(rownames(data),main))
main$Histo_annotation <- stmerged$Histo_annotation[main[,1]]
head(main)
ptspots <- main[main$main == 'PT',1]
maintab <- as.data.frame(table(main[,2:3]))
head(maintab)
maintab <- maintab[maintab$Histo_annotation != 'Ignore',]


labels <- c('Vessel','Medullary TI','Cortical TI - IMM','Cortical TI - FIB',
            'Cortical TI - Healthy','Sclerotic glom', 'Healthy glom')
nhisto <- length(unique(maintab$Histo_annotation))
forest_objs <- list()
fplot_list <- list()
pdf('../figures/forest_histo.pdf',width = 3.5,height = 4)
for (group in unique(maintab$main)) {
  #group <- 'PT'
  print(group)
  forest <- as.data.frame(matrix(0,ncol=3,nrow=nhisto,dimnames = list(unique(maintab$Histo_annotation),c('low','or','high'))))
  for (n in 1:nhisto){
    print(n)
    h <- unique(maintab$Histo_annotation)[n]

    test <- fisher.test(matrix(c(maintab[maintab$main == group & maintab$Histo_annotation == h,'Freq'],
                                sum(maintab[maintab$main != group & maintab$Histo_annotation == h,'Freq']),
                                sum(maintab[maintab$main == group & maintab$Histo_annotation != h,'Freq']),
                                sum(maintab[maintab$main != group & maintab$Histo_annotation != h,'Freq'])),ncol = 2))
    forest[n,] <- c(log(test$conf.int[1]),log(test$estimate),log(test$conf.int[2]))
    }

  forest$Histo_annotation <- rownames(forest)
  forest$color <- ifelse((forest$low >0 | forest$high < 0) &
                                  is.finite(forest$low >0),'blue','#BBBBBB')
  forest_objs[[group]] <- forest

  fplot <- ggplot(data=forest, aes(x=Histo_annotation, y=or, ymin=low, ymax=high)) + #,color=niches
                geom_pointrange(size = .7,linewidth = .7,shape = 2) + 
                geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
                #  coord_flip() +  # flip coordinates (puts labels on y axis)
                xlab('') + ylab("OR (95% CI)")+
                # scale_color_manual(values = forest[,'color'],
                #                 name = element_blank())+#,
                #                 #guide = guide_legend(reverse = TRUE),
                #                 #breaks = as.character(rev(ptniches)))+
                scale_x_discrete(limits=unique(maintab$Histo_annotation)[c(7,5,4,2,1,6,3)],
                                    labels=labels)+
                coord_flip() +
                theme_classic()+
                ggtitle(group)+
                theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
                    axis.text.x = element_text(size = 10),
                    axis.text.y = element_text(size = 10),panel.grid.major.y = element_line(linetype = 2),
                    legend.position = 'none')#,
                    #plot.margin = unit(c(5.5,5.5,110,5.5),'pt'))
                    
    plot(fplot)
    fplot_list[[group]] <- fplot$data
}
dev.off()
rio::export(fplot_list,'R_final/data_source/SD10.xlsx')


grns <- rio::import('../data/Supplementary Tables Version 1 05-2024.xlsx',
                    sheet= 'Table S21 - Traj GRNs',range='B5:LA191')


tfscoremeta <- list()
for (trajtf in colnames(grns)){
  geneset <- grns[trajtf]
  geneset <- geneset[!is.na(geneset)]
  geneset <- geneset[geneset %in% rownames(stmerged@assays$SCT)]

  p <- AddModuleScore(stmerged,features = list(geneset=geneset),assay='SCT')
  tfscoremeta[[trajtf]] <- p@meta.data
}


ptniches <- c(22,10,11,4,30,5,38,28,48,16,14,41)

tfs <- c('NR2F1','HNF4A','HNF1A','SOX4','SOX9','ELF3','KLF6','HIF1A','BACH2',
         'ATF3','NFKB1','RUNX2','RELB','JUNB','FOSL2')

heatdf <- c()
for (tf in tfs){
    heatdf  <- cbind(heatdf,tfscoremeta[[paste0('PTS1S2.',tf)]]$Cluster1)
}
head(heatdf)
colnames(heatdf) <- tfs
rownames(heatdf) <- rownames(tfscoremeta[[1]])
heatdf <- reshape2::melt(heatdf)
head(heatdf)
heatdf$Histo_annotation <- stmerged$Histo_annotation[heatdf$Var1]
heatdf <- aggregate(heatdf$value,by = list(heatdf$Histo_annotation,heatdf$Var2),FUN = mean)

min(heatdf$x)
max(heatdf$x)
pdf(paste0('../figures/heatmap_histo_selected_HRactivity.pdf'),width=5,height = 1.5)
hm <- ggplot(heatdf,aes(y=Group.1,x=Group.2,fill=x))+
    geom_tile(color='white',size=.7)+
    scale_fill_gradient2(low='#67001f',mid='white',high='#053061',
                         midpoint=0,limits=c(-0.2,0.3),oob=scales::squish)+
    theme_linedraw()+
    RotatedAxis()+
    labs(y='Histology',x='Transcription factor',fill='HR Activity')+
    scale_y_discrete(limits = c('IMM','FIB','Cort_TI'),
                    labels=c('Cortical TI - IMM','Cortical TI - FIB','Cortical TI - Healthy'))+
    scale_x_discrete(limits=tfs)
plot(hm)
dev.off()
write.csv(hm$data,'R_final/data_source/Fig_5b_TF_heatmap.csv',row.names = F)


origmetadata <- stmerged@meta.data
stmerged@meta.data <- origmetadata
Idents(stmerged) <- 'Histo_annotation'
CTIstmerged <- subset(stmerged,idents = c('Cort_TI','IMM','FIB'))

CTIstmerged@meta.data$SOX4 <- tfscoremeta[['PTS1S2.SOX4']][rownames(CTIstmerged@meta.data),'Cluster1']
CTIstmerged@meta.data$HNF4A <- tfscoremeta[['PTS1S2.HNF4A']][rownames(CTIstmerged@meta.data),'Cluster1']


plot_score <- function(spatial,img, img_grob,cell_type,sz,scl){
  if (ncol(spatial) > 1){
    spatial_coord <- GetTissueCoordinates(object = spatial@images[[1]],scale=NULL,cols = c("imagerow", "imagecol"))
    spatial_coord$x <- spatial_coord$x * scl
    spatial_coord$y <- spatial_coord$y * scl
    spatial_coord$data <- spatial@meta.data[,cell_type]

    feature_plt <- suppressMessages(
      ggplot2::ggplot() +
        ggplot2::annotation_custom(
          grob = img_grob,
          xmin = 0,
          xmax = ncol(img),
          ymin = 0,
          ymax = -nrow(img)) +
        ggplot2::geom_point(
          data = spatial_coord,
          ggplot2::aes(x = y,
                       y = x,
                       color= data),
                       size = sz)+
        ggplot2::scale_color_gradientn(colours = viridis::viridis(256,option='H'))+
        ggplot2::scale_y_reverse() +
        ggplot2::ylim(nrow(img), 0) +
        ggplot2::xlim(0, ncol(img)) +
        cowplot::theme_half_open(11, rel_small = 1) +
        ggplot2::theme_void()+
        ggplot2::coord_fixed(ratio = 1,
                              xlim = NULL,
                              ylim = NULL,
                              expand = TRUE,
                              clip = "on"))
    plot(feature_plt)
  }
}

plot_category <- function(spatial,img, img_grob,cell_type,sz,scl){  
  if (ncol(spatial) > 1){
    spatial_coord <- GetTissueCoordinates(object = spatial@images[[1]],scale=NULL,cols = c("imagerow", "imagecol"))
    spatial_coord$x <- spatial_coord$x * scl
    spatial_coord$y <- spatial_coord$y * scl
    spatial_coord$data <- spatial@meta.data[,cell_type]
    
    feature_plt <- suppressMessages(
      ggplot2::ggplot() +
        ggplot2::annotation_custom(
          grob = img_grob,
          xmin = 0,
          xmax = ncol(img),
          ymin = 0,
          ymax = -nrow(img)) +
        ggplot2::geom_point(
          data = spatial_coord,
          ggplot2::aes(x = y,
                       y = x,
                       color= data),
                       size = sz)+
        ggplot2::scale_y_reverse() +
        ggplot2::ylim(nrow(img), 0) +
        ggplot2::xlim(0, ncol(img)) +
        cowplot::theme_half_open(11, rel_small = 1) +
        ggplot2::theme_void()+
        ggplot2::coord_fixed(ratio = 1,
                              xlim = NULL,
                              ylim = NULL,
                              expand = TRUE,
                              clip = "on"))#+
    plot(feature_plt)
  }
}


pdf('../figures/FRactivity_histo.pdf',width = 18,height = 6)
for (img in names(CTIstmerged@images)){
    a <- SpatialFeaturePlot(CTIstmerged,'SOX4',crop=T,images=img)+
        ggtitle(paste0(tf,' - ',img))# +
        #theme(legend.position = 'none')
    b <- SpatialFeaturePlot(CTIstmerged,'HNF4A',crop=T,images=img)+
        ggtitle(paste0(tf,' - ',img))# +
        #theme(legend.position = 'none')
    c <- SpatialDimPlot(CTIstmerged,crop=T,images=img)+
        ggtitle(paste0('Histology - ',img))
    plot(a + b + c)
}
dev.off()

smp <- 'V42N07-395_XY01_235142'
Idents(CTIstmerged) <- 'orig.ident'
CTIsample <- subset(CTIstmerged,idents = smp)
img <- tiff::readTIFF(paste0('../ffpe_transfer/',smp,'.tif'))
img_grob <- grid::rasterGrob(img,
                            interpolate = FALSE,
                            width = grid::unit(1, "npc"),
                            height = grid::unit(1, "npc"))
for (tf in c('SOX4','HNF4A')){
    pdf(paste0('../figures/FRactivity_histo_',smp,'_',tf,'.pdf'),width = 6,height = 6)
    plot(plot_score(CTIsample,img,img_grob,tf,sz=1,scl=1))
    dev.off()
}
tf <- 'Histo_annotation'
pdf(paste0('../figures/FRactivity_histo_',smp,'_',tf,'.pdf'),width = 6,height = 6)
plot(plot_category(CTIsample,img,img_grob,tf,sz=1,scl=1))
dev.off()

biomarkers <- rev(c(
  'TGFBR3','PTGDS','PXDN','IL1R1','COL18A1','PTPRS','TNFRSF1A','B4GALT1',
  'HDGF','CMPK1','IDS','NBL1','FSTL3','SELENOM','GDF15','IGFBP7','FBLN1',
  'CFD','COL6A3','COL15A1','TNFRSF1B','ST8SIA4','GM2A','CD55',
  "CST3", "B2M"
))

unique(stmerged$Histo_annotation)
stmerged$Histo_annotation <- factor(stmerged$Histo_annotation,
                                    levels = c('VESSEL','Med_TI','IMM','FIB','Cort_TI',
                                               'Scl_glom','Glom'))
labels <- c('Vessel','Medullary TI','Cortical TI - IMM','Cortical TI - FIB',
            'Cortical TI - Healthy','Sclerotic glom', 'Healthy glom')

library(dplyr)

biodata <- reshape2::melt(as.matrix(stmerged@assays$SCT@data[biomarkers,]))
head(biodata)
biodata$Histo_annotation <- stmerged$Histo_annotation[biodata$Var2]
biodata <- biodata[biodata$Histo_annotation != 'Ignore',]
biodata <- aggregate(biodata$value,by = list(biodata$Histo_annotation,biodata$Var1),FUN = mean)
head(biodata)
biodata <- biodata %>%
  group_by(Group.2) %>%
  mutate(zscore = scale(x)) %>%
  ungroup()
biodata$Group.1 <- factor(biodata$Group.1,
                           levels = rev(c('VESSEL','Med_TI','IMM','FIB','Cort_TI',
                                      'Scl_glom','Glom')))

hm <- ggplot(biodata,aes(y=Group.2,x=Group.1,fill=zscore))+
    geom_tile(color='white',size=.7)+
    scale_fill_gradient2(low='#053061',mid='white',high='#67001f',
                                        midpoint=0,oob=scales::squish)+#limits=c(fc_min,fc_max),
    theme_linedraw()+
    theme(axis.text.x = element_text(angle=45,vjust = .5,hjust=1),
          axis.text.y = element_text(angle = 0, hjust = 1))+
    labs(y='',x='',fill='Expression')+
    scale_y_discrete(limits = rev(biomarkers))+
    theme(legend.position = 'bottom')
hm
ggsave('../figures/heatmap_biomarkers_histo_blue2red.pdf',width = 2.3,height = 4.5)
write.csv(hm$data,'R_final/data_source/Fig_6e_Visium_heatmap.csv',row.names = F)


