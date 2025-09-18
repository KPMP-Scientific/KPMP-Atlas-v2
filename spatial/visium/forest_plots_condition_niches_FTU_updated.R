library(Seurat)
library(ggplot2)

spatial <- readRDS('../data/new_samples/spatial_niches_FTU.RDS')
spatial$barcode <- row.names(spatial@meta.data)
meta <- rio::import('../data/Supplementary Tables Version 1 05-2024.xlsx',sheet=3,range='C6:K171')

meta <- merge(spatial@meta.data,meta,by.x='orig.ident',by.y='Specimen',all.x=TRUE)
head(meta)
row.names(meta) <- meta$barcode
meta$barcode <- NULL
spatial@meta.data <- meta

clinical <- rio::import('../data/Atlas v2 Participant Summary.xlsx')
unique(clinical$`Baseline eGFR (ml/min/1.73m2) (Binned)`)
clinical <- clinical[clinical$`Baseline eGFR (ml/min/1.73m2) (Binned)` != 'Unknown',]
clinical[clinical$`Baseline eGFR (ml/min/1.73m2) (Binned)` == '>60 ml/min/1.73m2','Baseline eGFR (ml/min/1.73m2) (Binned)'] <- '60-xx ml/min/1.73m2'
clinical$`simple egfr` <- as.numeric(stringr::str_split(clinical$`Baseline eGFR (ml/min/1.73m2) (Binned)`, "-",simplify = TRUE)[,1])
clinical$`eGFR group` <- cut(clinical$`simple egfr`, breaks = c(-Inf, 59, Inf), labels = c("Low", "High"))
head(clinical)
unique(clinical$patient)
table(clinical$`eGFR group`,clinical$`simple egfr`)

spatial$barcode <- row.names(spatial@meta.data)
meta <- merge(spatial@meta.data,clinical[,c('patient','Baseline eGFR (ml/min/1.73m2) (Binned)','eGFR group')],by.x='Patient',by.y='patient',all.x=TRUE)
head(meta)
row.names(meta) <- meta$barcode
meta$barcode <- NULL
table(unique(meta[,c('Patient','eGFR group','condition_level1')])[,c('eGFR group','condition_level1')])
table(meta[meta$condition_level1 == 'HRT',c('Patient','eGFR group')])
table(meta[,c('eGFR group','condition_level1')])
spatial@meta.data <- meta


write.csv(meta,'../data/new_samples/condition_meta_updated.csv',quote=F,row.names=F)

meta$forest_group <- ifelse(meta$condition_level1 == 'AKI', 'AKI',
                            ifelse(meta$condition_level1 == 'CKD' & meta$`eGFR group` == 'Low', 'CKD_Low',
                                   ifelse(meta$condition_level1 %in% c('CKD', 'DM-R') & meta$`eGFR group` == 'High', 'CKD_High', NA)))

forest_objs <- list()
for (group in names(sort(table(meta$group),decreasing = T))) {
  #group <- 'PT'
  nichetab <- table(meta[meta$group == group,c("forest_group","niches_group")])
  numniches <- max(as.numeric(colnames(nichetab)))
 
  forest_aki_ckdhigh <- as.data.frame(matrix(0,ncol=3,nrow=numniches+1,dimnames = list(as.character(0:numniches),c('low','or','high'))))
  forest_ckdlow_ckdhigh <- as.data.frame(matrix(0,ncol=3,nrow=numniches+1,dimnames = list(as.character(0:numniches),c('low','or','high'))))
  forest_aki_ckdlow <- as.data.frame(matrix(0,ncol=3,nrow=numniches+1,dimnames = list(as.character(0:numniches),c('low','or','high'))))
  for (n in as.character(0:numniches)){
    test <- fisher.test(matrix(c(nichetab['AKI',n],
                                nichetab['CKD_High',n],
                                sum(nichetab['AKI',])-nichetab['AKI',n],
                                sum(nichetab['CKD_High',])-nichetab['CKD_High',n]),ncol = 2))
    forest_aki_ckdhigh[n,] <- c(log(test$conf.int[1]),log(test$estimate),log(test$conf.int[2]))
    test <- fisher.test(matrix(c(nichetab['CKD_Low',n],
                                nichetab['CKD_High',n],
                                sum(nichetab['CKD_Low',])-nichetab['CKD_Low',n],
                                sum(nichetab['CKD_High',])-nichetab['CKD_High',n]),ncol = 2))
    forest_ckdlow_ckdhigh[n,] <- c(log(test$conf.int[1]),log(test$estimate),log(test$conf.int[2]))
    test <- fisher.test(matrix(c(nichetab['AKI',n],
                                nichetab['CKD_Low',n],
                                sum(nichetab['AKI',])-nichetab['AKI',n],
                                sum(nichetab['CKD_Low',])-nichetab['CKD_Low',n]),ncol = 2))
    forest_aki_ckdlow[n,] <- c(log(test$conf.int[1]),log(test$estimate),log(test$conf.int[2]))
    
  }

  forest_aki_ckdhigh$niches <- rownames(forest_aki_ckdhigh)
  forest_ckdlow_ckdhigh$niches <- rownames(forest_ckdlow_ckdhigh)
  forest_aki_ckdlow$niches <- rownames(forest_aki_ckdlow)

  forest_aki_ckdhigh$color <- ifelse((forest_aki_ckdhigh$low >0 | forest_aki_ckdhigh$high < 0) &
                                  is.finite(forest_aki_ckdhigh$low >0),'blue','#BBBBBB')
  forest_ckdlow_ckdhigh$color <- ifelse((forest_ckdlow_ckdhigh$low >0 | forest_ckdlow_ckdhigh$high < 0) &
                                  is.finite(forest_ckdlow_ckdhigh$low >0),'green','#BBBBBB')
  forest_aki_ckdlow$color <- ifelse((forest_aki_ckdlow$low >0 | forest_aki_ckdlow$high < 0) &
                                  is.finite(forest_aki_ckdlow$low >0),'purple','#BBBBBB')
  forest_objs[[paste0('aki_ckdhigh_',group)]] <- forest_aki_ckdhigh
  forest_objs[[paste0('ckdlow_ckdhigh_',group)]] <- forest_ckdlow_ckdhigh
  forest_objs[[paste0('aki_ckdlow_',group)]] <- forest_aki_ckdlow
}

saveRDS(forest_objs,'../data/new_samples/forest_niches_objs_FTU_updated.rds')
