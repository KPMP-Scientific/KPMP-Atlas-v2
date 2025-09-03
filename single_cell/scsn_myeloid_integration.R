library(Seurat)
library(leiden)
library(data.table)
## single cell data . The single cell object can be downloaded from https://cellxgene.cziscience.com/collections
kpmp = readRDS("./KPMP_PREMIERE_SC.RDS")
myl = subset(kpmp, cells = row.names(kpmp@meta.data)[kpmp$subclass.level.l1 %in% 'Myeloid'])
mylmeta = myl@meta.data
scmyldata = as(myl[["RNA"]]$counts, "CsparseMatrix")
scmyl = CreateSeuratObject(scmyldata, min.cells = 3)
scmyl = AddMetaData(scmyl, metadata = mylmeta)
scmyl$tech = 'single cell'
dim(scmyl@meta.data)
head(scmyl@meta.data)
## single nucleus data. Here the counts and metadata were obtained from local directories. The rawcounts ands meta data can be obtained from https://cellxgene.cziscience.com/collections/
sn.mat <- open_matrix_dir(dir = "./sn/v2_object_filtered-selected/v2_counts_filtered/")
sn = readRDS("./v2_object_filtered-selected/Kidney_Atlas_V2_07-2023_Object_Filtered.Rds")
myl = subset(sn, cells = row.names(sn@meta.data)[sn$v2.subclass.l1 %in% 'Myeloid'])
mylmeta = myl@meta.data
sn <- CreateSeuratObject(counts = sn.mat)
snmyl = subset(sn, cells = row.names(myl@meta.data))
snmyl = AddMetaData(snmyl, metadata = mylmeta)
snmyl$tech = 'single nucleus'
sndim(sctb@meta.data)

## integration & annotation . Manual checking on marker list was performed by external researchers 
scsn = merge(snmyl, scmyl)
scsn[['RNA']]$counts.1 <- convert_matrix_type(scsn[['RNA']]$counts.1, "double")
scsn[['RNA']]$counts.2 <- convert_matrix_type(scsn[['RNA']]$counts.2, "double")
scsn[['RNA']] = JoinLayers(scsn[['RNA']])
scsn <- NormalizeData(scsn)
scsn <- FindVariableFeatures(scsn)
scsn <- ScaleData(scsn, vars.to.regress = c('tech'))
scsn <- RunPCA(scsn)
scsn[["RNA"]] <- split(scsn[["RNA"]], f = scsn$tech)
scsn <- IntegrateLayers(object = scsn, method = CCAIntegration,orig.reduction = "pca", new.reduction = "cca",verbose = TRUE)
scsn <- FindNeighbors(scsn, reduction = "cca", dims = 1:30)
scsn <- FindClusters(scsn, algorithm = 4,method = 'igraph' , resolution =0.75, cluster.name = "cca.clusters")
scsn <- RunUMAP(scsn, reduction = "cca", dims = 1:30, reduction.name = "umap.cca")
scsn[['RNA']] = JoinLayers(scsn[['RNA']])

#### Annotating clusters
DotPlot(scsn, features = c('CXCL10'))
scsn = RenameIdents(scsn,'17' = 'moMAC-CXCL10+')
DotPlot(scsn, features = 'GYPC')
DotPlot(scsn, features = 'FCGR3A')
scsn = RenameIdents(scsn,'3' = 'ncMON')
scsn = RenameIdents(scsn,'6' = 'ncMON'))
DotPlot(scsn, features = 'FCN1')
scsn = RenameIdents(scsn,'7' = 'cMON')
scsn = RenameIdents(scsn,'11' = 'cMON')
DotPlot(scsn, features = 'C3')
scsn = RenameIdents(scsn,'8' = 'moMAC-C3+')
DotPlot(scsn, features = 'LYVE1')
scsn = RenameIdents(scsn,'4' = 'resMAC-LYVE1+')
scsn = RenameIdents(scsn,'5' = 'resMAC-HLAIIhi')
DotPlot(scsn, features = 'CPA3')
scsn = RenameIdents(scsn,'26' = 'MAST')
DotPlot(scsn, features = 'HBEGF')
DotPlot(scsn, features = 'LUCAT1'))
scsn = RenameIdents(scsn,'23' = 'moMAC-HBEGF+')
DotPlot(scsn, features = 'LUCAT1')
scsn = RenameIdents(scsn,'19' = 'moMAC-LUCAT1+'))
DotPlot(scsn, features = 'CENPF')
scsn = RenameIdents(scsn,'21' = 'cycMAC')
DotPlot(scsn, features = c('LAMP3'))
scsn = RenameIdents(scsn,'24' = 'mDC')
DotPlot(scsn, features = c('APOC1'))
scsn = RenameIdents(scsn,'10' = 'moMAC-FAM')
DotPlot(scsn, features = c('MZB1'))
scsn = RenameIdents(scsn,'25' = 'pDC')
DotPlot(scsn, features = c('CLEC10A'))
scsn = RenameIdents(scsn,'1' = 'cDC2'))
DotPlot(scsn, features = c('CLEC9A'))
scsn = RenameIdents(scsn,'18' = 'cDC1')
DotPlot(scsn, features = c('LYVE1')))
scsn = RenameIdents(scsn,'13' = 'resMAC-LYVE1+')
DimPlot(scsn, label = T, label.size = 6, repel = T)
DotPlot(scsn, features = c('LYVE1'))
DimPlot(scsn, label = T, label.size = 6, repel = T)
table(scsn$tech, Idents(scsn))
DotPlot(scsn, features = c('CLEC4A'))
FindMarkers(scsn, ident.1 = '22', only.pos = T)
scsn = RenameIdents(scsn,'22' = 'MYL_NonSpecific')
FindMarkers(scsn, ident.1 = '20', only.pos = T)
scsn = RenameIdents(scsn,'20' = 'MYL-EC_doub')
FindMarkers(scsn, ident.1 = '16', only.pos = T)
scsn = RenameIdents(scsn,'16' = 'MYL-PT_doub')
FindMarkers(scsn, ident.1 = '9', only.pos = T)
DotPlot(scsn, features = c('FCN1','AOAH','CSTA','CD36','UTRN','GPCA'))
DotPlot(scsn, features = c('FCN1','AOAH','CSTA','CD36','UTRN','GYPC'))
DotPlot(scsn, features = c('FCN1','AOAH','CSTA','CD36','UTRN','FCGR3A'))
scsn = RenameIdents(scsn,'14' = 'cMON')
scsn = RenameIdents(scsn,'9' = 'cMON')
FindMarkers(scsn, ident.1 = '2', only.pos = T)
scsn = RenameIdents(scsn,'2' = 'MYL-TAL_doub')
FindMarkers(scsn, ident.1 = '12', only.pos = T)
scsn = RenameIdents(scsn,'12' = 'cDC2')
FindMarkers(scsn, ident.1 = '15', only.pos = T)
scsn = RenameIdents(scsn,'15' = 'MYL-T_doub')
scsn$v2.subclass.level2 = Idents(scsn)
scsn$v2.cluster = ''
sort(table(scsn$subclass.level2))
scsn$v2.cluster[scsn$subclass.level2 %in% "cMON"] = 'MYL1'
scsn$v2.cluster[scsn$subclass.level2 %in% "ncMON"] = 'MYL2'
scsn$v2.cluster[scsn$subclass.level2 %in% "cDC2"] = 'MYL3'
scsn$v2.cluster[scsn$subclass.level2 %in% "resMAC-LYVE1+"] = 'MYL4'
scsn$v2.cluster[scsn$subclass.level2 %in% "MYL-TAL_doub"] = 'MYL5'
scsn$v2.cluster[scsn$subclass.level2 %in% "resMAC-HLAIIhi"] = 'MYL6'
scsn$v2.cluster[scsn$subclass.level2 %in% "moMAC-C3+"] = 'MYL7'
scsn$v2.cluster[scsn$subclass.level2 %in% "moMAC-FAM"] = 'MYL8'
scsn$v2.cluster[scsn$subclass.level2 %in% "MYL-T_doub"] = 'MYL9'
scsn$v2.cluster[scsn$subclass.level2 %in% "MYL-PT_doub"] = 'MYL10'
scsn$v2.cluster[scsn$subclass.level2 %in% "moMAC-CXCL10+"] = 'MYL11'
scsn$v2.cluster[scsn$subclass.level2 %in% "cDC1"] = 'MYL12'
scsn$v2.cluster[scsn$subclass.level2 %in% "moMAC-LUCAT1+"] = 'MYL13'
scsn$v2.cluster[scsn$subclass.level2 %in% "MYL-EC_doub"] = 'MYL14'
scsn$v2.cluster[scsn$subclass.level2 %in% "cycMAC"] = 'MYL15'
scsn$v2.cluster[scsn$subclass.level2 %in% "MYL_NonSpecific"] = 'MYL16'
scsn$v2.cluster[scsn$subclass.level2 %in% "moMAC-HBEGF+"] = 'MYL17'
scsn$v2.cluster[scsn$subclass.level2 %in% "mDC"] = 'MYL18'
scsn$v2.cluster[scsn$subclass.level2 %in% "pDC"] = 'MYL19'
scsn$v2.cluster[scsn$subclass.level2 %in% "MAST"] = 'MYL20'
scsn1 = subset(scsn, cells = row.names(scsn@meta.data)[!scsn$subclass.level2 %like% 'doub'])
save(scsn1, file = "scsn_Myl.integrated_052025.Robj")

## figure codes

pastel_palette <- c("#F6A9A9", "#FFCCB3", "#FFF0B3", "#CDEAC0", "#B3E0E6","#C2B3E6", "#E6B3E0", "#E6B3CC", "#FFB3E0", "#FFB3CC","#E6B3FF", "#CCB3FF", "#B3CCFF", "#B3E0FF", "#B3FFF0","#B3FFCC", "#CCFFB3", "#E0FFB3", "#F0FFB3", "#FFCDE6")
DimPlot(scsn1, label = F, repel = T, alpha = 0.75, split.by = 'tech',cols = pastel_palette[c(9,4)], pt.size = 0.6)
DimPlot(scsn1, label = T, repel = T, alpha = 0.75, split.by = 'tech',cols = pastel_palette, pt.size = 0.6)+theme(legend.position = 'none')
DotPlot(scsn1, features = c('CLEC9A','CLEC10A','CPA3','LILRA4','MKI67','CXCL10','C3','APOC1','FCN1','CCR7','GYPC'),   dot.scale = 12, cols = c('lightblue','darkred'))+theme(text = element_text(size=17))+ theme(axis.text.x = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"))+ theme(axis.title.y  = element_blank(),axis.text.y = element_text(color = "grey20", size = 17,  hjust = .5, vjust = .5, face = "plain"))


