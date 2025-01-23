library(Seurat)
library(SeuratDisk)
library(patchwork)
library(dplyr)
library(Matrix)
library(ggplot2)
library(BPCells)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

setwd("~/Projects/Human_Kidney/Atlas_V2/")


#### PT-TL Marker Plot
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(KB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)

KB.PT <- subset(KB, idents = c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
                               "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
                               "dDTL","ATL","dATL"))

KB.PT <- subset(KB.PT, downsample = 1000)
KB.PT[["RNA"]] <- as(KB.PT[["RNA"]], Class = "Assay")
KB.PT <- NormalizeData(KB.PT)

pt.markers <- c(
  "NPHS1","NPHS2","ST6GALNAC3","PODXL",                         #POD
  "PTGDS","BST2","TPPP3","CHI3L1",                              #dPOD
  "ALDH1A2","CFH","FAM155A","CLDN1",                            #PEC
  
  "LRP2","CUBN",                                                #PT
  "SLC5A12","SLC22A6","HNF4A",                                  #S1/S2
  "RALYL","PCDH15","PRODH2","SLC22A8",                          #S1
  "SLC34A1","ANKS1B","SLC5A10","SLC5A11",                       #S2                                  
  "SLC5A8","GPM6A","SLC22A24","SLC7A13",                        #PT-S3
  
  "IGFBP7","SPP1","ITGB8","CDH6","TMEM178B","ALPK2","HAVCR1",
  "ITGB3",
  "CST3","CLU","VIM","PIGR",#"APOE",                            #aPT2
  "IL32","SOX4","VCAM1","MMP7","SOX9","CCL2",                   #aPT2
  "DCC",                                                        #aPT1
  "GDA","GLIS1",                                                #aPT-S1/S2
  "DLGAP1","PROM1",                                             #frPTS1/S2
  "APBB1IP","ROBO2","COL23A1","MEG3",                           #frPTS1/S2
  "LSAMP","KCNIP1","NRXN3","WNT2B",                             #frPTS1/S2
  "KCTD16","SPON1",                                             #aPT-S3
  "NRG3","FAM189A1","DTNA","KITLG","GRM8",                      #frPT-S3
  "TOP2A","MKI67",                                              #cycling
  "EGR1","FOS",                                                 #dPT
  
  "PAX8-AS1","SLC44A5","CRYAB","TACSTD2",                       #TL
  "ABCA13",                                                     #DTL
  "AQP1", "UNC5D","LRRC4C","DENND2A","SLCO1A2",                 #DTL2
  "IRX3", "SERPINA1",                                           #aDTL2
  "SIM2",                                                       #DTL1/3/ATL
  "ADGRL3","JAG1","SMAD9","ID1",                                #DTL1
  "SLC14A2","FAM169A","SMOC2",                                  #DTL3
  "ABCA4","BCL6","AKR1B1","SH3GL3",                             #DTL3/ATL
  "CLCNKA","PROX1","CACNA1C","COLEC10",                         #ATL
  "SOD3","HSPA2","CSRP2",                                       #dATL
  "CST3","APOE","GATM","ALDOB","CLU","DEFB1",                   #Injury
  "SPP1","IGFBP7","CLU"
)

pdf(file='Plots/Proximal-TL_Subclassl3_Marker_Dotplot.pdf',width=24,height=6)
DotPlot(KB.PT, features = unique(pt.markers)) + RotatedAxis()
dev.off()







### TAL-IC Marker Plot
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(KB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)

KB.DT <- subset(KB, idents = c("M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","MD","C/M-TAL-B","C-TAL-B","dC-TAL",
                               "frTAL","aTAL1","aTAL2","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
                               "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
                               "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B"))

KB.DT <- subset(KB.DT, downsample = 1000)
KB.DT[["RNA"]] <- as(KB.DT[["RNA"]], Class = "Assay")
KB.DT <- NormalizeData(KB.DT)

dt.markers <- c(
  "CASR","SLC12A1","UMOD","EGF","ESRRB",                  #TAL
  "ANK2","CLCNKA","KCTD16","BMP6","RIPOR3","CLDN14",      #M-TAL
  "PHACTR1","SLCO3A1","CXCL12","CNTN1","CABP1",           #TAL-A
  "KCNMB2","RGS6",                                        #C/M-TAL-A
  "ENOX1","CALCR","RBM20","PDE3A",                        #C-TAL-A
  "DACH1","LRMDA",                                        #TAL-B
  
  "TENM4","FGF14","PXDNL","GRM8",                         #C/M-TAL-B
  "KCNJ10","TMEM52B","CLDN16","TMEM207","JAG1",           #TAL-B
  "SCN7A","COL8A1","LINGO2",                              #C-TAL-B                          
  "BBOX1","NOS1","ROBO2","TMPRSS4",                       #MD
  
  "CREB5","COL4A1","ITGA3","HIVEP3","CYTOR",              #aTAL1
  "CD44","DIAPH3","SYT16","HAVCR1",                       #aTAL1                                      

  "ITGB6","NRP1","TFPI",                                  #aTAL
  
  "HIF1A","ADAMTS1","DNAH5",                              #aTAL2
  
  "ITGB8","PROM1","ARHGAP26","RNF144B","TMPRSS4","RHEX",  #frTAL
  "EZH2","CENPP","MKI67",#"TOP2A",                        #cycling

  "SLC12A3","CNNM2","KLHL3","TRPM6",                       #DCT
  "GRAMD1B","ADAMTS17","ZNF385D",                          #DCT1
  
  "UNC5C","HMCN1","FSTL4","TRPV5",                         #DCT2
  "CACNA1D", "ACSL4",                                      #frDCT
  "FGF13","IGF2BP2","FAM155A","NRG1",                      #aDCT 

  "SLC8A1",                                                #DCT2 / CNT
  "HSD11B2","CALB1","ANGPT1","SNTG1",                      #CNT
  "CTSD",                                                  #dCNT
  "RAPGEF5","DLGAP1","BIRC3",                              #aCNT
  
  "GATA3","AQP2","PAPPA",                                  #PC
  "SCNN1G","SCNN1B",
  "SGCD","STC1",                                           #CNT-PC
  "FGF12","PTCSC3","CNTNAP3B",                             #CCD/OMCD PC
  "MCTP1","CPAMD8","GABRA2",                               #OMCD-PC
  "GREB1L",                                                #OMCD-PC/IMCD
  "SLC14A2","HS3ST5","TENM3",#"TGFBR3","AKR1C1","FAT3",    #IMCD
  "AOC1","TFPI2","LCN2",                                   #dIMCD
  "ADIRF","SNCG","DHRS2","GPX2","TP63",#"FXYD3",           #PapE                              #PapE
  
  "ATP6V0D2", "ATP6V1C2", "CLNK",                          #IC
  "SLC26A7", "SLC4A1",                                     #IC-A                                   
  "HS6ST3","NRXN3", "NXPH2", "LEF1",                       #CCD-IC-A
   
  "FAM184B","ADTRP","AQP6","STAP1",                        #OMCD-IC-A
  
  "SLC4A9", "SLC26A4", "INSRR", "TLDC2",                #IC-B

  "CKB","COX8A",#"PEBP1","UQCRB",                         #Injury
  "CST3","DEFB1","SPP1","IGFBP7","CLU"
)

pdf(file='Plots/Distal_TAL-IC_Subclassl3_Marker_Dotplot.pdf',width=26,height=7)
DotPlot(KB.DT, features = unique(dt.markers)) + RotatedAxis()
dev.off()






#### EC Marker Plot
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(KB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)

KB.EC <- subset(KB, idents = c("EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
                               "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
                               "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC"))

KB.EC <- subset(KB.EC, downsample = 1000)
KB.EC[["RNA"]] <- as(KB.EC[["RNA"]], Class = "Assay")
KB.EC <- NormalizeData(KB.EC)

ec.markers <- c("PECAM1","PTPRB","FLT1",                                   #Broad EC
                "EMCN","HECW2","ITGA8",                                    #EC-GC
                "EHD3","RXFP1", "MGP","SOST",                              #EC-GC      
                'PTCHD4',"ZMAT3",                                          #aEC-GC
                "AQP1","FILIP1","H19","ESM1","SLC45A4",                    #EC-GC-FILIP1+
                
                "PDE3A","SULF1","NKAIN2","NOS1",                           #EC-AA
                "ADAMTS6","MCTP1","PALMD","SLC14A1","ITIH5",               #EC-DVR  
                #"LYPD6B","EPHA3","STXBP6","CP",
                
                "CEACAM1","PLVAP","DNASE1L3",                              #PTC/AVR
                "COL15A1","FABP5","ALPL", #"CD36",                         #M-EC-PTC
                "GPM6A","NR2F2",                                           #PTC/AVR
                "ZNF385D","RANBP3L","EDIL3",                               #EC-AVR
                "MX2","RSAD2","ISG15","IFIT1",                             #iaEC-AVR
                
                "SLCO2A1",                                                
                "VWF","RYR3","ADGRG6","CPE","TRABD2B",                     #EC-V
                "ADAMTSL1","CMKLR1",                                       #EC-V/EC-PCV
                "DOK6",                                                    #EC-PCV
                "NAV3","OSMR","SYCP2L",                                    #C-EC-PTC
                "AFAP1L1","USP31","MYO1B","LAMA4","NETO2",                 #angEC-PTC
                "SLC6A6","FUT8","ATP13A3","AFF3",                          #EC-EA
                "IFITM3","HLA-DRA","CAVIN2","CCL14","CA4",                 #dEC-PTC
                
                'ICAM1',"TNFAIP3",'CCL2','SELE',"CXCL2",                   #infEC-PTC
                "VCAM1",                                                   #inf/iaEC-PTC
                "CXCL10","GBP1","CXCL11","CTSS",                           #iaEC-PTC
                "MMRN1","CD36","TBX1","PROX1",                             #EC-LYM
                "TOP2A","MKI67","CENPF"                                    #cycling
)

pdf(file='Plots/Vascular_EC_Subclassl3_Marker_Dotplot.pdf',width=20,height=5.2)
DotPlot(KB.EC, features = unique(ec.markers)) + RotatedAxis()
dev.off()







### STR Marker Plot 

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(KB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)

KB.STR <- subset(KB, idents = c("IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
                               "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
                               "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
                               "Ad"))

KB.STR <- subset(KB.STR, downsample = 1000)
KB.STR[["RNA"]] <- as(KB.STR[["RNA"]], Class = "Assay")
KB.STR <- NormalizeData(KB.STR)

str.markers <- c(
  "DCN","C7","PDGFRA",                                     #Pan FIB
  "SYT1","TNC","IGFBP2","RARRES2",                         #Pan Medullary FIB
  "CA8","HPSE2","GABRG3",                                  #IM/OM-FIB
  "MCTP2","SNAP25","BDKRB2",                               #IM-FIB
  "KIF26B","FREM1","KCNQ3","LOXHD1",                       #OM & C/M-FIB
  "SPP1","TIMP1",                                          #dM-FIB
  "COL1A2","COL3A1","COL1A1",                              #dFIB & MYOF
  "ADAMTSL1","KCNIP1","ADARB2",                            #C/M-FIB
  "RYR2","ZNF536","SEMA3D","ACTG2",# "PAMR1",              #IM-pvMYOF
  
  "NEGR1","LAMA2","ABCA8","MEG3",                           #Pan cortical FIB
  "CCN1","CCN2","ELL2","SAMHD1","SLC2A3",                   #C-FIB (interstitial fib)
  "GRIN2A","EMID1",                                         #C-FIB-PATH
  "SELENOP","LUM","CXCL12",'GGT5',"ECRG4",                  #C-FIB-OSMRlo
  "OSMR","SOD2","UGCG","IL1R1","CCL2",                      #C-FIB-OSMRhi
  "RELB","CXCL10","CCL19",                                  #C-FIB-OSMRhi
  "SULF1","GLI2","NTM","INHBA","FAP","POSTN",               #C-MYOF
  "SPARC","BGN","VIM",                                      #dFIB

  "FLRT2","COL12A1","FGF14",                                #Pan pvFIB
  "PDZRN4",'IGF1','ADAMTS3',"RSPO3","WNT5B",                #pvFIB-RSPO3+
  "C3","EBF2","SFRP2","CD34","PI16",                        #pvFIB-PI16+
  "ITGBL1","PLD5","CNTN5",                                  #pvFIB & pvMYOF
  "MGAT4C","EPHA3",                                         #pvFIB
  "ITGA8",                                                  #pvMYOF & VSMC
  "ADGRL3","MYH11","ACTA2",'KCNMA1',"PCDH7",                #pvMYOF
  "PRUNE2","MYOCD",#"SYNPO2",'MACROD2',                     #pvMYOF

  
  "PDGFRB","SLCO3A1",                                       #Pan VSM markers
  "SAMD12","ROBO1","PIEZO2",                                #MC & REN
  "DAAM2","GATA3","DCC","POSTN","IL1RL1",                   #MC
  "ROBO2","REN","KCNQ5","SLCO2A1","SPOCK3",                 #REN
  "MYH11","NTRK3",'RGS6','KCNMA1',"ADRA1A","MCAM",          #VSMC
  "NOTCH3",                                                 #VSMC &VSMC/P
  "RGS5",'PLCB4',"FRMD3",                                   #VSMC/P
  "ADGRB3","SLC38A11","C2CD2","DNAH11","SLC6A1",            #M-VSMC/P
  'PDE1C',"STEAP4","RXFP1",                                 #VSMC/P (cortical)
  'FLNA', 'TAGLN',"ACTA2","VIM","ACTB",                     #dVSMC
  "CD36","PLIN1","LPL","ADIPOQ","FABP4"                     #Ad
)


pdf(file='Plots/Stromal_Subclassl3_Marker_Dotplot.pdf',width=26,height=5.6)
DotPlot(KB.STR, features = unique(str.markers)) + RotatedAxis()
dev.off()







### IMM Marker Plot 

KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_05162024.rds")
Idents(KB) <- "v2.subclass.l3"

order <- c("POD","dPOD","PEC","PT-S1","dPT-S1","PT-S2","dPT-S2","PT-S3","dPT-S3","aPT2",
           "aPT1","aPT-S1/S2","frPT-S1/S2","aPT-S3","frPT-S3","cycPT","dPT","DTL2","aDTL2","DTL1","DTL3",
           "dDTL","ATL","dATL","M-TAL","dM-TAL","C/M-TAL-A","C-TAL-A","C/M-TAL-B","C-TAL-B","MD","dC-TAL",
           "aTAL1","aTAL2","frTAL","cycTAL","DCT1","DCT2","dDCT","frDCT","aDCT","CNT","dCNT","aCNT",
           "CNT-PC","dCNT-PC","CCD-PC","OMCD-PC","dOMCD-PC","IMCD","dIMCD","PapE","CCD-IC-A","dCCD-IC-A",
           "OMCD-IC-A","dOMCD-IC-A","tOMCD-PC-IC","IC-B","EC-GC","dEC-GC","aEC-GC","EC-GC-FILIP1+","EC-AA",
           "EC-DVR","M-EC-PTC","EC-AVR","dEC-AVR","iaEC-AVR","EC-V","EC-PCV","C-EC-PTC","angEC-PTC",
           "EC-EA","dEC-PTC","infEC-PTC","iaEC-PTC","EC-LYM","cycEC","IM-FIB","OM-FIB","dIM-FIB","dOM-FIB",
           "C/M-FIB","IM-pvMYOF","C-FIB","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C-MYOF","dFIB",
           "pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","MC","REN","VSMC","M-VSMC/P","VSMC/P","dVSMC",
           "Ad","B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
           "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
           "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU")

Idents(object = KB) <- factor(Idents(object = KB), levels = order)

KB.IMM <- subset(KB, idents = c("B","PL","Naïve Th","MAIT","ILC3","T-REG","CD8+ TEM/TRM","CD8+ TEM/TEMRA","NK","ERY",
                                "MAST","resMAC-LYVE1+","resMAC-HLAIIhi","moMAC-HBEGF+","moMAC-CXCL10+","moFAM","moMAC-C3+",
                                "cDC2","ncMON","MON","mDC","cDC1","pDC","N","cycT","cycMAC","SC/NEU"))

KB.IMM <- subset(KB.IMM, downsample = 1000)
KB.IMM[["RNA"]] <- as(KB.IMM[["RNA"]], Class = "Assay")
KB.IMM <- NormalizeData(KB.IMM)

imm.markers <- c("PTPRC",                                             #Broad Immune
                 "BANK1","MS4A1","CD37","CD79A",                      #B Cells
                 "IGKC","XBP1","MZB1","JCHAIN",#"SDC1",               #PL Cells
                 "CD96","CD247","BCL11B","THEMIS",                    #T
                 "INPP4B","TRAC","CD3D",                              #T
                 "IL7R",                                              #T/MAIT/ILC3
                 "LEF1","CD4","SELL",                                 #Naïve Th
                 "SLC4A10","KLRB1","CCR6",                            #MAIT
                 "PCDH9","TOX2","KIT","RORC",                         #ILC3
                 "IKZF2","RTKN2","IL2RA","CTLA4",#"FOXP3",            #T-REG
                 "CD8A",                                              #CD8+
                 "GZMK",                                              #CD8+ TEM/TRM
                 "CCL5","SAMD3","GZMA","CCL4",                        #CD8+ & NK 
                 "NKG7","KLRD1","GNLY","GZMB","CX3CR1",               #CD8+ TEM/TEMRA & NK
                 "GZMH",                                              #CD8+ TEM/TEMRA
                 "TXK","KLRF1","NCAM1",#"PDGFD",                      #NK
                   
                 
                 "HBB","HBA2","HBA1",                                #Erythrocyte
                 "CPA3","IL18R1","TPSB2","TPSAB1",                   #MAST
                 "CD163","MSR1","CSF1R","CD14",                      #MAC
                 "MRC1","F13A1","STAB1","CD163L1","LYVE1",           #resMAC-LYVE1+
                 "HLA-DPA1","C1QA","C1QB",                           #resMAC-HLAIIhi
                 "HIF1A","NAMPT","PLAUR","ITGAX","HBEGF","OSM",      #moMAC-HBEGF+
                 "PSTPIP2","CXCL10","CXCL9","CCL2","CCL3","IL1B",    #moMAC-CXCL10+
                 "GPNMB","SPP1","APOC1","PLA2G7","CD68","CAPG",      #moFAM
                 "HMOX1","TREM2",                                    #moFAM
                 "C3","KCNQ3","ADGRB3","VASH1","CX3CR1",             #moMAC-C3+
                 "CLEC10A","FCER1A","CD1C",                          #cDC2
                 "TCF7L2","COTL1","FCGR3A",                          #ncMON
                 "FCN1",                                             #MON/ncMON
                 "VCAN","LYZ","CD36",                                #MON
                 "LAMP3","SLCO5A1","CCR7",#"EBI3","CCL19",           #mDC
                 "WDFY4","CADM1","CLEC9A","BATF3",                   #cDC1
                 "BCL11A","CLEC4C","IL3RA","PLD4",#"LILRA4",         #pDC
                 "S100A9","FCGR3B","S100A8","IFITM2",                #N
                 "TOP2A","MKI67",                                    #cycling
                 "NRXN1","GRIK2","CDH19","NCAM2"                     #SC/NEU
                 
)


pdf(file='Plots/Immune-NEU_Subclassl3_Marker_Dotplot.pdf',width=26,height=5.6)
DotPlot(KB.IMM, features = unique(imm.markers)) + RotatedAxis()
dev.off()
