library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)


spatial <- readRDS('../data/new_samples/spatial_niches_FTU.RDS')

supptab <- rio::import_list('../data/new_samples/Supplementary Tables Version new-samples.xlsx')
annot <- supptab[["Table S - sn markers"]]
rownames(annot) <- annot$v2.subclass.l3

grns <- supptab[["Table S - Traj GRNs"]]
colnames(grns) <- grns[1,]
grns <- grns[2:nrow(grns),]



tfscoremeta <- list()
for (trajtf in colnames(grns)){
  geneset <- grns[trajtf]
  geneset <- geneset[!is.na(geneset)]
  geneset <- geneset[geneset %in% rownames(spatial@assays$SCT)]

  p <- AddModuleScore(spatial,features = list(geneset=geneset),assay='SCT')
  tfscoremeta[[trajtf]] <- p@meta.data
}
sox4act <- tfscoremeta[['frPTS1S2.SOX4']][,c('orig.ident','Cluster1')]
colnames(sox4act) <- c('Sample','frPTS1S2.SOX4')
sox4act$PTS1S2.SOX4 <- tfscoremeta[['PTS1S2.SOX4']][rownames(sox4act),'Cluster1']

sox4_means <- sox4act %>%
  group_by(Sample) %>%
  summarise(
    mean_frPTS1S2_SOX4 = mean(frPTS1S2.SOX4, na.rm = TRUE),
    mean_PTS1S2_SOX4 = mean(PTS1S2.SOX4, na.rm = TRUE)
  )

trajgroup <- data.frame('trajectory'=unique(stringr::str_split(colnames(grns),'\\.',simplify=T)[,1]),
                        'group'=c('stroma cells','stroma cells','stroma cells','PT','PT','TAL','TAL','PT','PT','immune cells','immune cells','immune cells'))
rownames(trajgroup) <- trajgroup$trajectory

tfscorelist <- list()
for (traj in trajgroup$trajectory){
  tfs4traj <- colnames(grns)[grep(paste0('^',traj),colnames(grns))]

  combined_scores <- data.frame()

  for (trajtf in tfs4traj) {
    current_scores <- aggregate(tfscoremeta[[trajtf]][,'Cluster1'],
                              by=list('group'=tfscoremeta[[trajtf]]$group,'niches_group'=tfscoremeta[[trajtf]]$niches_group),FUN=mean)
    
    if (nrow(combined_scores) == 0) {
      combined_scores <- current_scores
    } else {
      combined_scores <- full_join(combined_scores, current_scores, by = c("group", "niches_group"))
    }
    colnames(combined_scores)[ncol(combined_scores)] <- trajtf
  }
  tfscorelist[[traj]] <- combined_scores
}
saveRDS(tfscorelist,'../data/new_samples/transcription_factor_activity_FTU_niches.rds')
