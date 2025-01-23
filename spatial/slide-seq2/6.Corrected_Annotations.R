library(Seurat)
library(SeuratData)
library(Matrix)
library(ggplot2)
library(dplyr)
library(BPCells)

options(future.globals.maxSize = 1e9)

setwd("~/Projects/Human_Kidney/Atlas_V2/")
options(Seurat.object.assay.version = "v5")

###Add in experiment metadata
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_Filtered_Object_07292024_InBoundary.Rds")
meta <- kss@meta.data
exp.meta <- read.delim("slide-seq/slide-seq_experiment_table.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level2","condition_level1","condition",
         "region_level2","region_level1",
         "age_binned","sex","race","tissue_type","atlas_version")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
kss@meta.data <- meta
colnames(kss@meta.data)

order <- c("library","cell_ID","nCount_Spatial","nFeature_Spatial","maxWeight.l1.ct","leverage.score","pagoda_k100_infomap","pagoda_k100_infomap_full",
           "pagoda_k100_infomap_full.score","source","assay","experiment","patient","specimen","condition_level2","condition_level1","condition","region_level1",
           "region_level2","age_binned","sex","race","tissue_type","atlas_version","maxWeight.l3","maxCellType.l3","v2.subclass.sp",
           "v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
kss@meta.data <- kss@meta.data[,order]

#Fix count location
counts <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_counts_09-02-2023")
kss[["Spatial"]]$counts <- counts[,colnames(kss)]
kss <- NormalizeData(kss)


###Correct annotations
order <- c("POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","DTL1","DTL3","ATL",
           "M-TAL","C-TAL","MD","aTAL1","aTAL2","frTAL","DCT","aDCT","CNT","aCNT",
           "PC","IMCD","IC-A","tPC-IC","IC-B","EC-GC","aEC-GC","EC-AA","EC-DVR",
           "EC-AVR","infEC-AVR","EC-V","EC-PCV","EC-PTC","angEC-PTC","EC-EA",
           "infEC-PTC","EC-LYM","M-FIB","C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-Path",
           "C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+","pvFIB-PI16+",
           "pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
           "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+",
           "resMAC-HLAIIhi","moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON",
           "mDC","cDC1","pDC","N")
Idents(kss) <- "v2.subclass.sp"
Idents(kss) <- factor(Idents(kss), levels = order)
table(Idents(kss))

meta <- kss@meta.data


##POD-TL predictions and Marker Genes
kss.sub <- subset(kss, idents = c("POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","DTL1","DTL3","ATL"))
kss.sub <- subset(kss.sub, downsample = 1000)
kss.sub <- NormalizeData(kss.sub)

pt.markers <- c(
  "NPHS1","NPHS2","PODXL",                         #POD
  "ALDH1A2","CFH",                            #PEC
  
  "LRP2","CUBN",                                                #PT
  "SLC5A12","SLC22A6",                                  #S1/S2
  "PRODH2",                          #S1
  "SLC34A1",                       #S2                                  
  "SLC5A8",                         #PT-S3
  
  "ITGB8","CDH6",
  "IL32","SOX4","VCAM1",                   #aPT2
  "DCC",                                                        #aPT1
  "GDA",                                                #aPT-S1/S2
  "DLGAP1","PROM1",                                             #frPTS1/S2
  "APBB1IP","ROBO2","MEG3",                           #frPTS1/S2
  "LSAMP",                             #frPTS1/S2
  
  "CRYAB","TACSTD2",                       #TL
  "AQP1", "UNC5D","LRRC4C",                 #DTL2
  "SIM2",                                                       #DTL1/3/ATL
  "JAG1","SMAD9","ID1",                                #DTL1
  "SLC14A2","SMOC2",                                  #DTL3
  "ABCA4","BCL6","AKR1B1","SH3GL3",                             #DTL3/ATL
  "CLCNKA","PROX1"                         #ATL
)

DotPlot(kss.sub, features = unique(pt.markers), dot.scale = 8, scale.max = 50) + RotatedAxis()

#Clean up PECs
PEC <- rownames(meta[meta$v2.subclass.sp == "PEC",])
kss <- SetIdent(kss, cells = PEC, value = "NA")
PEC.keep <- intersect(PEC, WhichCells(object = kss, expression = CFH > 0))
kss <- SetIdent(kss, cells = PEC.keep, value = "PEC")



##TAL-IC predictions and Marker Genes
kss.sub <- subset(kss, idents = c("M-TAL","C-TAL","MD","aTAL1","aTAL2","frTAL","DCT","aDCT","CNT","aCNT",
                                  "PC","IMCD","IC-A","tPC-IC","IC-B"))
kss.sub <- subset(kss.sub, downsample = 1000)
kss.sub <- NormalizeData(kss.sub)

dt.markers <- c(
  "SLC12A1","UMOD","EGF",                  #TAL
  "ANK2","CLCNKA",      #M-TAL
  "ENOX1","TMEM52B","CLDN16",                        #C-TAL
  "BBOX1","NOS1","ROBO2","TMPRSS4",                       #MD
  
  "CREB5","COL4A1","ITGA3","HIVEP3","CYTOR",              #aTAL1
  "CD44",                       #aTAL1                                      
  
  "ITGB6","NRP1","TFPI",                                  #aTAL
  
  "HIF1A","ADAMTS1",                              #aTAL2
  
  "ITGB8","PROM1","ARHGAP26","RNF144B","TMPRSS4","RHEX",  #frTAL
  
  "SLC12A3","CNNM2","KLHL3","TRPM6",                       #DCT
  "ACSL4",                                      #frDCT
  "FGF13","IGF2BP2","FAM155A","NRG1",                      #aDCT 
  
  "SLC8A1",                                                #DCT2 / CNT
  "HSD11B2","CALB1",                      #CNT
  "RAPGEF5","DLGAP1","BIRC3",                              #aCNT
  
  "GATA3","AQP2","PAPPA",                                  #PC
  "SCNN1G","SCNN1B",
  "GREB1L",                                                #OMCD-PC/IMCD
  "SLC14A2","HS3ST5",                                      #IMCD
  
  "ATP6V0D2", "ATP6V1C2", "CLNK",                          #IC
  "SLC26A7", "SLC4A1",                                     #IC-A                                   
  
  "SLC4A9", "SLC26A4", "INSRR", "TLDC2"                #IC-B
)

DotPlot(kss.sub, features = unique(dt.markers), dot.scale = 8, scale.max = 50) + RotatedAxis()


#Clean up MD
TAL <- rownames(meta[meta$v2.subclass.sp == "MD",])
MD <- intersect(TAL, WhichCells(object = kss, expression = BBOX1 > 0))
kss <- SetIdent(kss, cells = TAL, value = "NA")
kss <- SetIdent(kss, cells = MD, value = "MD")



##EC predictions and Marker Genes
kss.sub <- subset(kss, idents = c("EC-GC","aEC-GC","EC-AA","EC-DVR",
                                  "EC-AVR","infEC-AVR","EC-V","EC-PCV","EC-PTC","angEC-PTC","EC-EA",
                                  "infEC-PTC","EC-LYM"))
kss.sub <- subset(kss.sub, downsample = 1000)
kss.sub <- NormalizeData(kss.sub)

ec.markers <- c("PECAM1","PTPRB","FLT1",                                   #Broad EC
                "EMCN","HECW2","ITGA8",                                    #EC-GC
                "EHD3","SOST",                              #EC-GC      
                'PTCHD4',"ZMAT3",                                          #aEC-GC
                
                "PDE3A","SULF1","NKAIN2","NOS1",                           #EC-AA
                "AQP1","ADAMTS6","MCTP1","PALMD","SLC14A1","ITIH5",               #EC-DVR  
                
                "CEACAM1","PLVAP","DNASE1L3",                              #PTC/AVR
                "GPM6A","NR2F2",                                           #PTC/AVR
                "ZNF385D","RANBP3L","EDIL3",                               #EC-AVR
                "MX2","RSAD2","ISG15","IFIT1",                             #iaEC-AVR
                
                "SLCO2A1",                                                
                "VWF","RYR3","ADGRG6",                     #EC-V
                "ADAMTSL1",                                       #EC-V/EC-PCV
                "DOK6",                                                    #EC-PCV
                "NAV3","OSMR",                                    #C-EC-PTC
                "AFAP1L1","MYO1B","LAMA4","NETO2",                 #angEC-PTC
                "SLC6A6","FUT8","ATP13A3","AFF3",                          #EC-EA

                'ICAM1',"TNFAIP3",'CCL2','SELE',                   #infEC-PTC
                "VCAM1",                                                   #inf/iaEC-PTC
                "CXCL10","GBP1",                           #iaEC-PTC
                "MMRN1","CD36","TBX1","PROX1"                             #EC-LYM
)

DotPlot(kss.sub, features = unique(ec.markers), dot.scale = 8, scale.max = 50) + RotatedAxis()

table(Idents(kss))
table(kss$v2.subclass.sp,kss$condition_level1)

#re-label aEC-GC
aECGC <- rownames(meta[meta$v2.subclass.sp %in% c("aEC-GC"),])
kss <- SetIdent(kss, cells = aECGC, value = "EC-GC")

#infEC-AVR (too highly represented)
#angEC-PTC (too highly represented)
#infEC-PTC (too highly represented)


##Fibroblast predictions and Marker Genes
kss.sub <- subset(kss, idents = c("M-FIB","C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-Path",
                                  "C-FIB-OSMRhi","C-FIB-OSMRlo","C-MYOF","pvFIB-RSPO3+","pvFIB-PI16+",
                                  "pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P"))
kss.sub <- subset(kss.sub, downsample = 1000)
kss.sub <- NormalizeData(kss.sub)

str.markers <- c(
  "DCN","C7","PDGFRA",                                     #Pan FIB
  "TNC",                         #Pan Medullary FIB
  "CA8",                                  #IM/OM-FIB
  "FREM1","KCNQ3",                      #OM & C/M-FIB
  "ADAMTSL1",                            #C/M-FIB
  "ACTG2",              #IM-pvMYOF
  
  "NEGR1","LAMA2","MEG3",                           #Pan cortical FIB
  "CCN1","CCN2","ELL2",                   #C-FIB (interstitial fib)
  "GRIN2A","EMID1",                                         #C-FIB-PATH
  "SELENOP","LUM","CXCL12",'GGT5',"ECRG4",                  #C-FIB-OSMRlo
  "OSMR","SOD2","UGCG","IL1R1","CCL2",                      #C-FIB-OSMRhi
  "COL1A2","COL3A1","COL1A1",                              #dFIB & MYOF
  "SULF1","GLI2","NTM","INHBA","FAP",               #C-MYOF
  "SPARC","BGN",                                      #dFIB
  
  "FLRT2","COL12A1","FGF14",                                #Pan pvFIB
  "PDZRN4",'IGF1',"RSPO3",                #pvFIB-RSPO3+
  "C3","SFRP2","PI16",                        #pvFIB-PI16+
  "ITGBL1","PLD5","CNTN5",                                  #pvFIB & pvMYOF
  "MGAT4C","EPHA3",                                         #pvFIB
  "ITGA8",                                                  #pvMYOF & VSMC
  "MYH11","ACTA2",                #pvMYOF
  "PRUNE2","MYOCD",#"SYNPO2",'MACROD2',                     #pvMYOF
  
  "PDGFRB","SLCO3A1",                                       #Pan VSM markers
  "ROBO1","PIEZO2",                                #MC & REN
  "GATA3","POSTN","IL1RL1",                   #MC
  "ROBO2","REN","SLCO2A1",                 #REN
  "MYH11","NTRK3",'RGS6',"MCAM",          #VSMC
  "NOTCH3",                                                 #VSMC &VSMC/P
  "RGS5", "ADGRB3","SLC38A11","C2CD2"                             #VSMC/P
  
  )

DotPlot(kss.sub, features = unique(str.markers), dot.scale = 8, scale.max = 50) + RotatedAxis()


#C-FIB-OSMRlo looks over-represented, likely mixture with C-FIB-PATH


##Immune Cell predictions and Marker Genes
kss.sub <- subset(kss, idents = c("B","PL","T","MAIT","ILC3",
                                  "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+",
                                  "resMAC-HLAIIhi","moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON",
                                  "mDC","cDC1","pDC","N"))
kss.sub <- subset(kss.sub, downsample = 1000)
kss.sub <- NormalizeData(kss.sub)

imm.markers <- c("PTPRC",                                             #Broad Immune
                 "BANK1","MS4A1","CD37",                     #B Cells
                 "XBP1","MZB1","JCHAIN",               #PL Cells
                 "CD96","THEMIS",                    #T
                 "INPP4B","CD3D",                              #T
                 "IL7R",                                              #T/MAIT/ILC3
                 "LEF1","CD4",                                 #Naïve Th
                 "SLC4A10","KLRB1",                            #MAIT
                 "PCDH9","TOX2","KIT",                         #ILC3
                 "IKZF2","RTKN2","IL2RA",            #T-REG
                 "CD8A",                                              #CD8+
                 "GZMK",                                              #CD8+ TEM/TRM
                 "CCL5","SAMD3","GZMA","CCL4",                        #CD8+ & NK 
                 "KLRD1","GNLY","GZMB",               #CD8+ TEM/TEMRA & NK
                 "TXK","KLRF1","NCAM1",                      #NK
                 
                 
                 "HBB","HBA2","HBA1",                                #Erythrocyte
                 "CPA3","IL18R1","TPSB2","TPSAB1",                   #MAST
                 "CD163","MSR1","CSF1R","CD14",                      #MAC
                 "MRC1","F13A1","STAB1","LYVE1",           #resMAC-LYVE1+
                 "HLA-DPA1","C1QA","C1QB",                           #resMAC-HLAIIhi
                 "HIF1A","NAMPT","PLAUR",      #moMAC-HBEGF+
                 "CXCL10","CXCL9","CCL2","IL1B",    #moMAC-CXCL10+
                 "GPNMB","SPP1","CD68","CAPG",      #moFAM
                 "TREM2",                                    #moFAM
                 "C3","KCNQ3","ADGRB3","VASH1",             #moMAC-C3+
                 "CLEC10A","FCER1A","CD1C",                          #cDC2
                 "TCF7L2","COTL1","FCGR3A",                          #ncMON
                 "FCN1",                                             #MON/ncMON
                 "VCAN","LYZ","CD36",                                #MON
                 "SLCO5A1","CCR7",           #mDC
                 "WDFY4","CADM1","CLEC9A",                   #cDC1
                 "BCL11A",         #pDC
                 "S100A9","FCGR3B","S100A8","IFITM2"                #N
                 )

DotPlot(kss.sub, features = unique(imm.markers), dot.scale = 8, scale.max = 50) + RotatedAxis()


#Res Mac
lyve1 <- rownames(meta[meta$v2.subclass.sp %in% c("resMAC-LYVE1+"),])
MHC <- intersect(lyve1, WhichCells(object = kss, expression = C1QA > 0))
mhc <- rownames(meta[meta$v2.subclass.sp %in% c("resMAC-HLAIIhi"),])
kss <- SetIdent(kss, cells = mhc, value = "NA")
mhc.real <- intersect(mhc, WhichCells(object = kss, expression = CD163 > 0 | STAB1 > 0))
mhc.real <- intersect(mhc.real, WhichCells(object = kss, expression = C1QA > 0))
kss <- SetIdent(kss, cells = c(MHC,mhc.real), value = "resMAC-HLAIIhi")




###Finalize annotations
kss$v2.subclass.sp <- Idents(kss)

order <- c("POD","PEC","PT-S1/2","PT-S3","aPT","frPT","DTL2","DTL1","DTL3","ATL",
           "M-TAL","C-TAL","MD","aTAL1","aTAL2","frTAL","DCT","aDCT","CNT","aCNT",
           "PC","IMCD","IC-A","tPC-IC","IC-B","EC-GC","EC-AA","EC-DVR",
           "EC-AVR","infEC-AVR","EC-V","EC-PCV","EC-PTC","angEC-PTC","EC-EA",
           "infEC-PTC","EC-LYM","M-FIB","C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-Path",
           "C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","pvFIB-RSPO3+","pvFIB-PI16+",
           "pvFIB","pvMYOF","MC","REN","VSMC","VSMC/P","B","PL","T","MAIT","ILC3",
           "T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY","MAST","resMAC-LYVE1+",
           "resMAC-HLAIIhi","moMAC-INF","moFAM","moMAC-C3+","cDC2","ncMON","MON",
           "mDC","cDC1","pDC","N")
Idents(kss) <- "v2.subclass.sp"
Idents(kss) <- factor(Idents(kss), levels = order)
table(Idents(kss))

kss <- subset(kss, idents = order)



###Update object and metadata
#Update annotations
meta <- kss@meta.data
cl.meta <- read.delim("slide-seq/Slide-Seq_Annotation_Table.txt")
emc <- c("v2.subclass.full","v2.subclass.l1","v2.state.l2","v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$v2.subclass.sp, cl.meta$v2.subclass.sp)]
}
colnames(meta)

kss@meta.data <- meta

DimPlot(kss, label = T, label.size = 3, reduction = "full.umap.1", raster=FALSE,
        group.by = "v2.subclass.sp", alpha = 0.4) + NoLegend()

kss[["full.umap.1"]] <- NULL
saveRDS(kss,file = "~/hsKidAt/blake_LTS/Atlas_V2/spatial/slide-seq/Spatial_object_09-22-2023/Kidney_Slide-seq_Spatial_Atlas_V2_09042024_Final_Object.Rds")
kss
#An object of class Seurat 
#103153 features across 1431823 samples within 4 assays 
#Active assay: Spatial (51531 features, 0 variable features)
#2 layers present: counts, data
#3 other assays present: l1.predictions, sketch, l3.predictions
#4 dimensional reductions calculated: pca, umap, pca.full, full.umap
#71 images present

