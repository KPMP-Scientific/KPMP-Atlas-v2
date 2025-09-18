library(Seurat)
library(ggplot2)

plot_forest <- function(comparison,niches){
    forest_plot <- ggplot(data=forest[[comparison]][rev(niches)+1,], aes(x=niches, y=or, ymin=low, ymax=high,color=niches)) + #,color=niches
                        geom_pointrange(size = .7,linewidth = 1,shape = 2) + 
                        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
                        #  coord_flip() +  # flip coordinates (puts labels on y axis)
                        xlab('') + ylab("OR (95% CI)")+
                        scale_color_manual(values = forest[[comparison]][rev(niches)+1,'color'],
                                        name = element_blank(),
                                        #guide = guide_legend(reverse = TRUE),
                                        breaks = as.character(rev(niches)))+
                        scale_x_discrete(limits=as.character(rev(niches)))+
                        coord_flip() +
                        theme_classic()+
                        theme(panel.border = element_rect(color = "black", fill=NA, linewidth = 1.5),
                            axis.text.x = element_text(size = 10),
                            axis.text.y = element_text(size = 10),panel.grid.major.y = element_line(linetype = 2),
                            legend.position = 'none')#,
                            #plot.margin = unit(c(5.5,5.5,110,5.5),'pt'))
    return(forest_plot)
}

plot_heat <- function(heatdf,activity,niches,tflist,min,max){
    heatplot <- ggplot(heatdf,aes_string(y='niches_group',x='TF',fill=activity))+
                    geom_tile(color='white',size=.7)+
                    scale_fill_gradient2(low='#67001f',mid='white',high='#053061',
                                        midpoint=0,limits=c(min,max),oob=scales::squish)+
                    theme_linedraw()+
                    theme(axis.text.x = element_text(angle=45,vjust = .5,hjust=1),
                        axis.text.y = element_text(angle = 0, hjust = 1))+
                    labs(x='Niche',y='Transcription factor',fill='Activity')+
                    scale_y_discrete(limits = as.character(rev(niches)))+
                    scale_x_discrete(limits=tflist)
    return(heatplot)
}


spatial <- readRDS('../data/new_samples/spatial_niches_new_group_FTU.RDS')



annot <- rio::import('../data/Supplementary Tables Version 1 05-2024.xlsx',
                    sheet= 'Table S12 - sn markers',range='B5:N133')
rownames(annot) <- annot$v2.subclass.l3
head(annot)


for(grp in unique(spatial$new_group)){
    cells <- rownames(spatial@meta.data[spatial$new_group == grp,])
    sptgrp <- subset(spatial,cells = cells)

    plt_height <- ceiling(2*max(as.numeric(sptgrp$niches_new_group))/10) + 1
    pdf(paste0('../figures/new_samples_abridged/dotplot_new_group_',grp,'_niches.pdf'),width=25,height = plt_height)
        p <- DotPlot(sptgrp,features = unique(annot$v2.subclass.l3),group.by = 'niches_new_group',assay = 'predv2subclassl3')+
                RotatedAxis() +
                guides(
                    size = guide_legend(title = "% Spots"),
                    color = guide_colorbar(title = "Avg. Prop.")
                ) +
                # theme(legend.position = 'left',
                #     axis.text.x = element_text(angle=45,vjust = .5,hjust=1))+
                scale_y_discrete(limits = as.character(rev(0:max(as.numeric(sptgrp$niches_new_group)))))
        plot(p)
    dev.off()
}

forest <- readRDS('../data/new_samples/forest_niches_objs_FTU_updated.rds')

tfscorelist <- readRDS('../data/new_samples/transcription_factor_activity_FTU_niches.rds')


cells <- rownames(spatial@meta.data[spatial$new_group == 'PT',])
sptgrp <- subset(spatial,cells = cells)

ptniches <- c(22,10,11,4,30,5,38,28,48,3,34,14,41) #check 22, 31, 48 mapping
ptcts <- c('aPT2','aPT1','aPT-S1/S2','PT-S1','PT-S2','frPT-S1/S2','PT-S3','frPT-S3','cycPT',
           "C-FIB", "C/M-FIB", "C-FIB-PATH", "C-FIB-OSMRhi", "C-FIB-OSMRlo",
           "C-MYOF",  "pvFIB", "pvMYOF",
           "B", "PL", "Naïve Th", "MAIT",
           "CD8+ TEM/TRM",
           "resMAC-LYVE1+", "resMAC-HLAIIhi",
           "ncMON", "MON", "pDC")
pdf(paste0('../figures/new_samples_abridged/dotplot_PT_niches.pdf'),width=7.5,height = 3.7)
DotPlot(sptgrp,features = ptcts,group.by = 'niches_group',assay = 'predv2subclassl3')+ #unique(annot$v2.subclass.l3)
    theme(legend.position = 'left',
            axis.text.x = element_text(angle=45,vjust = .5,hjust=1))+
    scale_y_discrete(limits = as.character(rev(ptniches)))
dev.off()

pdf(paste0('../figures/new_samples_abridged/dotplot_PT_allniches.pdf'),width=7.5,height = 12)
DotPlot(sptgrp,features = ptcts,group.by = 'niches_group',assay = 'predv2subclassl3')+ 
    theme(legend.position = 'left',
            axis.text.x = element_text(angle=45,vjust = .5,hjust=1))+
    scale_y_discrete(limits = as.character(rev(0:max(as.numeric(sptgrp$niches_new_group)))))
dev.off()


DefaultAssay(sptgrp) <- 'predv2subclassl3'
pdf(paste0('../figures/new_samples_abridged/dotplot_PT_niches_allCT.pdf'),width=25,height = 5)
DotPlot(sptgrp,features = unique(annot$v2.subclass.l3),group.by = 'niches_group',assay = 'predv2subclassl3')+ #unique(annot$v2.subclass.l3)
    theme(legend.position = 'left',
            axis.text.x = element_text(angle=45,vjust = .5,hjust=1))+
    scale_y_discrete(limits = as.character(rev(0:max(as.numeric(sptgrp$niches_new_group)))))
dev.off()


tfs <- c('SOX4','SOX9','HNF4A','RUNX2','ATF3','POU3F3','HES1','ETS1','EGR1',
         'ELF3','NFKB1','KLF6','JUNB','FOSL2')
heatdf <- tfscorelist[['frPTS1S2']]
heatdf <- heatdf[heatdf$group == 'PT',c('niches_group',paste0('frPTS1S2.',tfs))]
heatdf <- reshape2::melt(heatdf)
colnames(heatdf) <- c('niches_group','TF','FR_activity')
heatdf$TF <- stringr::str_replace(heatdf$TF,'frPTS1S2.','')

heatdf2 <- tfscorelist[['PTS1S2']]
heatdf2 <- heatdf2[heatdf2$group == 'PT',c('niches_group',paste0('PTS1S2.',tfs))]
heatdf2 <- reshape2::melt(heatdf2)
colnames(heatdf2) <- c('niches_group','TF','HR_activity')
heatdf2$TF <- stringr::str_replace(heatdf2$TF,'PTS1S2.','')

heatdf <- merge(heatdf,heatdf2,by = c('niches_group','TF'))

heatdf$delta <- heatdf$FR_activity - heatdf$HR_activity
heatdf$FC <- heatdf$FR_activity/heatdf$HR_activity
heatdf$logFC <- log2(heatdf$FC)
min(heatdf$FR_activity)
max(heatdf$FR_activity)
pdf(paste0('../figures/new_samples_abridged/heatmap_PT_niches_FRactivity.pdf'),width=4.5,height = 3.7)
hm <- ggplot(heatdf,aes(y=niches_group,x=TF,fill=FR_activity))+
    geom_tile(color='white',size=.7)+
    scale_fill_gradient2(low='#67001f',mid='white',high='#053061',
                         midpoint=0,limits=c(-0.1,0.3),oob=scales::squish)+
    theme_linedraw()+
    theme(axis.text.x = element_text(angle=45,vjust = .5,hjust=1),
          axis.text.y = element_text(angle = 0, hjust = 1))+
    labs(y='Niche',x='Transcription factor',fill='Activity')+
    scale_y_discrete(limits = as.character(rev(ptniches)))+
    scale_x_discrete(limits=tfs)
plot(hm)
dev.off()
write.csv(hm$data,'R_final/data_source/Fig_4b_PTheatmap.csv',row.names = F,quote = F)


akip <- plot_forest('aki_ckdhigh_PT',ptniches)
pdf(paste0('../figures/new_samples_abridged/forest_PT_niches_AKI_updated.pdf'),width=1.5,height = 3.7)
akip
dev.off()
write.csv(akip$data,'R_final/data_source/Fig_4b_AKIforest.csv',row.names = F,quote = F)

ckdp <- plot_forest('ckdlow_ckdhigh_PT',ptniches)
pdf(paste0('../figures/new_samples_abridged/forest_PT_allniches_CKD.pdf'),width=1.5,height = 12)
ckdp
dev.off()
write.csv(ckdp$data,'R_final/data_source/Fig_4b_CKDforest.csv',row.names = F,quote = F)

 

#Niches for TAL
cells <- rownames(spatial@meta.data[spatial$group == 'TAL',])
sptgrp <- subset(spatial,cells = cells)

talniches <- c(21,20,3,18,41,4,30,32,17,14,22,16,26)
talcts <- c('aTAL1','aTAL2','C/M-TAL-A','C-TAL-A','C/M-TAL-B','C-TAL-B','frTAL','cycTAL',
            "IM-FIB", "OM-FIB", "dFIB", "dOM-FIB",
           "C/M-FIB", "C-MYOF", 
           "B", "PL", "Naïve Th", "MAIT","ncMON", "mDC","cDC1")
pdf(paste0('../figures/new_samples_abridged/dotplot_TAL_niches.pdf'),width=6.7,height = 3.6)
DotPlot(sptgrp,features = talcts,group.by = 'niches_group',assay = 'predv2subclassl3')+ #unique(annot$v2.subclass.l3)
    theme(legend.position = 'left',
            axis.text.x = element_text(angle=45,vjust = .5,hjust=1))+
    scale_y_discrete(limits = as.character(rev(talniches)))
dev.off()

tfs <- c('SOX4','ESRRB','JUN','BACH1','ELF3','KLF5','KLF13','FOXP2')

heatdf <- tfscorelist[['frTAL']]
heatdf <- heatdf[heatdf$group == 'TAL',c('niches_group',paste0('frTAL.',tfs))]
heatdf <- heatdf[heatdf$niches_group %in% talniches,]
heatdf <- reshape2::melt(heatdf)
colnames(heatdf) <- c('niches_group','TF','FR_activity')
heatdf$TF <- stringr::str_replace(heatdf$TF,'frTAL.','')

min(heatdf$FR_activity)
max(heatdf$FR_activity)
pdf(paste0('../figures/new_samples_abridged/heatmap_TAL_niches_FRactivity.pdf'),width=3.3,height = 3.6)
hm <- ggplot(heatdf,aes(y=niches_group,x=TF,fill=FR_activity))+
    geom_tile(color='white',size=.7)+
    scale_fill_gradient2(low='#67001f',mid='white',high='#053061',
                         midpoint=0,limits=c(-0.2,0.8),oob=scales::squish)+
    theme_linedraw()+
    theme(axis.text.x = element_text(angle=45,vjust = .5,hjust=1),
          axis.text.y = element_text(angle = 0, hjust = 1))+
    labs(x='Niche',y='Transcription factor',fill='Activity')+
    scale_y_discrete(limits = as.character(rev(talniches)))+
    scale_x_discrete(limits=tfs)
plot(hm)
dev.off()
write.csv(hm$data,'R_final/data_source/Fig_ED8a_TALheatmap.csv',row.names = F,quote = F)


heatdf <- heatdf[heatdf$TF != 'SOX4',]
min(heatdf$FR_activity)
max(heatdf$FR_activity)
pdf(paste0('../figures/new_samples_abridged/heatmap_TAL_niches_FRactivity_noSOX4.pdf'),width=3.1,height = 3.6)
ggplot(heatdf,aes(y=niches_group,x=TF,fill=FR_activity))+
    geom_tile(color='white',size=.7)+
    scale_fill_gradient2(low='#67001f',mid='white',high='#053061',
                         midpoint=0,limits=c(-0.1,0.3),oob=scales::squish)+
    theme_linedraw()+
    theme(axis.text.x = element_text(angle=45,vjust = .5,hjust=1),
          axis.text.y = element_text(angle = 0, hjust = 1))+
    labs(x='Niche',y='Transcription factor',fill='Activity')+
    scale_y_discrete(limits = as.character(rev(talniches)))+
    scale_x_discrete(limits=tfs[tfs!= 'SOX4'])
dev.off()


akip <- plot_forest('aki_ckdhigh_TAL',talniches)
pdf(paste0('../figures/new_samples_abridged/forest_TAL_niches_AKI_updated.pdf'),width=1.5,height = 3.6)
akip
dev.off()
write.csv(akip$data,'R_final/data_source/Fig_ED8a_TAL_AKIforest.csv',row.names = F,quote = F)

ckdp <- plot_forest('ckdlow_ckdhigh_TAL',talniches)
pdf(paste0('../figures/new_samples_abridged/forest_TAL_niches_CKD_updated.pdf'),width=1.5,height = 3.6)
ckdp
dev.off()
write.csv(ckdp$data,'R_final/data_source/Fig_ED8a_TAL_CKDforest.csv',row.names = F,quote = F)

#Niches for FIB
cells <- rownames(spatial@meta.data[spatial$group == 'stroma cells',])
sptgrp <- subset(spatial,cells = cells)

fibniches1 <- c(18,41,17,0,49,13,46,6,7,44,43,30,1,19)
fibcts1 <- c("C-FIB", "C-FIB-PATH", "C-FIB-OSMRhi", "C-FIB-OSMRlo",
           "C-MYOF", "pvFIB-RSPO3+", "pvFIB-PI16+", "pvFIB", "pvMYOF",
           "B", "PL","NK", "ERY", "MAST",
           "resMAC-LYVE1+", "resMAC-HLAIIhi", "moMAC-HBEGF+", "moMAC-CXCL10+",
           "moMAC-C3+", "cDC2", "mDC",
           "cDC1", "pDC", "cycT", 
            'aPT2','aPT1','aPT-S1/S2','PT-S1','PT-S2','frPT-S1/S2','PT-S3','frPT-S3','cycPT')

pdf(paste0('../figures/new_samples_abridged/dotplot_FIB_final_niches.pdf'),width=9,height = 3.7)
DotPlot(sptgrp,features = fibcts1,group.by = 'niches_group',assay = 'predv2subclassl3')+ #unique(annot$v2.subclass.l3)
    theme(legend.position = 'left',
            axis.text.x = element_text(angle=45,vjust = .5,hjust=1))+
    scale_y_discrete(limits = as.character(rev(fibniches1)))
dev.off()


akip <- plot_forest('aki_ckdhigh_stroma cells',fibniches1)
pdf(paste0('../figures/new_samples_abridged/forest_FIB_final_niches_AKI_updated.pdf'),width=1.5,height = 3.7)
akip
dev.off()
write.csv(akip$data,'R_final/data_source/Fig_ED8b_FIB_AKIforest.csv',row.names = F,quote = F)

ckdp <- plot_forest('ckdlow_ckdhigh_stroma cells',fibniches1)
pdf(paste0('../figures/new_samples_abridged/forest_FIB_final_niches_CKD_updated.pdf'),width=1.5,height = 3.7)
ckdp
dev.off()
write.csv(ckdp$data,'R_final/data_source/Fig_ED8b_FIB_CKDforest.csv',row.names = F,quote = F)


tflist1 <- colnames(tfscorelist[['CFIBOSMRhi']])
tflist1 <- tflist1[grep('CFIBOSMRhi',tflist1)]
tflist1 <- stringr::str_split(tflist1,'\\.',simplify = T)[,2]
tflist1 <- tflist1[c(4:6,1:3,7:11)]

heatdf <- tfscorelist[['CFIBOSMRhi']]
heatdf <- heatdf[heatdf$group == 'stroma cells',c('niches_group',paste0('CFIBOSMRhi.',tflist1))]
heatdf <- heatdf[heatdf$niches_group %in% rev(fibniches1),]
heatdf <- reshape2::melt(heatdf)
colnames(heatdf) <- c('niches_group','TF','OSMR_activity')
heatdf$TF <- stringr::str_replace(heatdf$TF,'CFIBOSMRhi.','')
min(heatdf$OSMR_activity)
max(heatdf$OSMR_activity)
heatplot <- plot_heat(heatdf,'OSMR_activity',fibniches1,tflist1,-0.1,0.9)
pdf(paste0('../figures/new_samples_abridged/heatmap_FIB_final_niches_OSMRactivity.pdf'),width=3.7,height = 3.7)
heatplot
dev.off()
write.csv(heatplot$data,'R_final/data_source/Fig_ED8b_FIB_heatmap.csv',row.names = F,quote = F)

heatdf <- heatdf[!heatdf$TF %in% c("NR2F1","NR2F2"),]
min(heatdf$OSMR_activity)
max(heatdf$OSMR_activity)
heatplot <- plot_heat(heatdf,'OSMR_activity',fibniches1,tflist1[!tflist1 %in% c("NR2F1","NR2F2")],-0.1,0.2)
pdf(paste0('../figures/new_samples_abridged/heatmap_FIB_final_niches_OSMRactivity.pdf'),width=4,height = 3.2)
heatplot
dev.off()



fibniches3 <- c(11,40,50,39,47,45,32,2,34)
fibcts3 <- c("C-FIB", "C-FIB-PATH", "C/M-FIB", "C-FIB-OSMRhi", "C-FIB-OSMRlo",
           "C-MYOF", "pvFIB-RSPO3+", "pvFIB-PI16+", "pvFIB", "pvMYOF",
           'VSMC', 'M-VSMC/P', 'VSMC/P', 'dVSMC')

pdf(paste0('../figures/new_samples_abridged/dotplot_FIB_tr3_niches.pdf'),width=6,height = 3.3)
DotPlot(sptgrp,features = fibcts3,group.by = 'niches_group',assay = 'predv2subclassl3')+ #unique(annot$v2.subclass.l3)
    theme(legend.position = 'left',
            axis.text.x = element_text(angle=45,vjust = .5,hjust=1))+
    scale_y_discrete(limits = as.character(rev(fibniches3)))
dev.off()

akip <- plot_forest('aki_ckdhigh_stroma cells',fibniches3)
pdf(paste0('../figures/new_samples_abridged/forest_FIB_tr3_niches_AKI_updated.pdf'),width=1.9,height = 3.3)
akip
dev.off()
write.csv(akip$data,'R_final/data_source/Fig_ED9i_FIB_AKIforest.csv',row.names = F,quote = F)

ckdp <- plot_forest('ckdlow_ckdhigh_stroma cells',fibniches3)
pdf(paste0('../figures/new_samples_abridged/forest_FIB_tr3_niches_CKD_updated.pdf'),width=1.9,height = 3.3)
ckdp
dev.off()
write.csv(ckdp$data,'R_final/data_source/Fig_ED8a_TAL_CKDforest.csv',row.names = F,quote = F)

tflist3 <- colnames(tfscorelist[['pvMYOF']])
tflist3 <- tflist3[grep('pvMYOF',tflist3)]
tflist3 <- stringr::str_split(tflist3,'\\.',simplify = T)[,2]
heatdf3 <- tfscorelist[['pvMYOF']]
heatdf3 <- heatdf3[heatdf3$group == 'stroma cells',c('niches_group',paste0('pvMYOF.',tflist3))]
heatdf3 <- heatdf3[heatdf3$niches_group %in% rev(fibniches3),]
heatdf3 <- reshape2::melt(heatdf3)
colnames(heatdf3) <- c('niches_group','TF','pvMYOF_activity')
heatdf3$TF <- stringr::str_replace(heatdf3$TF,'pvMYOF.','')

min(heatdf3$pvMYOF_activity)
max(heatdf3$pvMYOF_activity)
heatplot <- plot_heat(heatdf3,'pvMYOF_activity',fibniches3,tflist3,-0.1,0.8)
pdf(paste0('../figures/new_samples_abridged/heatmap_FIB_tr3_niches_pvMYOFactivity.pdf'),width=9,height = 3.3)
heatplot
dev.off()




