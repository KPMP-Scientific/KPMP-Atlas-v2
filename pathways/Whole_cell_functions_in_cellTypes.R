rm(list = ls());
library(ggplot2)
library(reshape2)
library(gplots)
library(cowplot)
library(patchwork)
library(readxl)
library(gridExtra)
library(tidyr)
library(ggplotify)
library(ggpubr)

################# User-input needed ##############################
overall_directory = "D:/KPMP_disease_atlas/KPMP_data/"
working_directory = "D:/AACode/R_code/KPMP_atlas_v2_code/"
setwd(working_directory)
################# User-input needed ##############################

source('KPMP_atlas_v2_shared_code.R')

overall_data_directory = paste(overall_directory,"KPMP_atlas_v2_max5000cells_500DEGs/",sep='')
hfFunctions_directory = paste(overall_directory,"KPMP_atlas_v2_reference_pFs/",sep='')
dataset_cellTypes = c("SNboth_PODPEC","SNboth_PT","SNboth_TAL","SNboth_DCT","SNboth_CDCNT","SNboth_FIB","SNboth_EC","SNboth_Myeloid","SNboth_Lymphoid")

considered_disease_pf_groups = c("Metabolism","TM movement of ions, water and solutes",
                                 "Redox, iron and mitochondrial dynamics",
                                 "Cell structure dynamics (Cytoskeleton, PM and nucleus)","Cell adhesion","ECM remodeling",
                                 "Immune response and signaling","Gene expression")

stopifnot(length(which(!considered_disease_pf_groups %in% names(pfGroup_color_list)))==0)

dataset_abbreviation_dataset_list = list()
dataset_abbreviation_dataset_list[["SN"]] = "SNRNAseq"

max_minLog10Pvalue = 15
minColCount_subtype_pfGroups = 20

sn_diseaseOrSubtype_in_correct_order = c("HRT","Ref","CKD_high_eGFR","CKD_low_eGFR","AKI","CKD","border_row",
                                         
                                         "dPOD","POD","PEC",
                                         
                                         "PT-S1","PT-S2","PT-S3",
                                         "dPT-S1","dPT-S2","dPT-S3","dPT",
                                         "aPT2","aPT1","aPT-S1/S2","aPT-S3",
                                         "frPT-S1/S2","frPT-S3","cycPT",
                                         
                                         "C-TAL-A","C-TAL-B","C/M-TAL-A","C/M-TAL-B","M-TAL","MD",
                                         "dC-TAL","dM-TAL","aTAL1","aTAL2","frTAL","cycTAL",
                                         
                                         "DCT1","DCT2","dDCT","aDCT","frDCT",

                                         "MON","moMAC-CXCL10+","moMAC-HBEGF+","moFAM","moMAC-C3+","resMAC-LYVE1+","resMAC-HLAIIhi","MAST","ncMON","cycMAC","cDC1","cDC2","pDC","mDC","N","ERY",
                                         
                                         "dFIB","dIM-FIB","dOM-FIB","IM-FIB","C-FIB","C-MYOF","C-FIB-PATH","C-FIB-OSMRlo","C-FIB-OSMRhi","C/M-FIB","OM-FIB","pvFIB-RSPO3+","pvFIB-PI16+","pvFIB","pvMYOF","IM-pvMYOF"
)
datasetAbbreviation_diseaseOrSubtype_in_correct_order = list()
datasetAbbreviation_diseaseOrSubtype_in_correct_order[["SN"]] = sn_diseaseOrSubtype_in_correct_order

datasetCellType_vLineXintercepts_list = list()
datasetCellType_vLineXintercepts_list[["SNboth_PT"]] = c(0,3,4+c(0,3,7,11,13,14))
datasetCellType_vLineXintercepts_list[["SNboth_TAL"]] = c(0,3,4+c(0,6,8,10,11,12))
datasetCellType_vLineXintercepts_list[["SNboth_DCT"]] = c(0,3,4+c(0,2,3,4,5))
datasetCellType_vLineXintercepts_list[["SNboth_Myeloid"]] = c(0,3,4+c(0,5,7,8,9,10,15,15,16))
datasetCellType_vLineXintercepts_list[["SNboth_FIB"]] = c(0,3,4+c(0,3, 4, 9, 11, 15,16))
datasetCellType_vLineXintercepts_list[["SNboth_PODPEC"]] = c(0,3,4+c(0,2,3))

plot_only_disease_scps=FALSE
diseaseLabel = "SimD"
subClassLabel = "sl3"

Heatmaps = list()
hfGroup_disease_in_correct_order = c()
all_disease_hfGroups = c()
hfGroup_colors = c()
hfGroup_names = c()

pfGroup_heatmaps = list()
subtype_pfGroup_heatmaps = list()
pfGroup_heatmapsII = list()

indexDCT = length(dataset_cellTypes)
for (indexDCT in 1:length(dataset_cellTypes))
{#Begin
   dataset_cellType = dataset_cellTypes[indexDCT]

   splitStrings = strsplit(dataset_cellType,"_")[[1]]
   dataset_abbreviation = splitStrings[1]
   cellType = splitStrings[2]
   dataset = dataset_abbreviation_dataset_list[[dataset_abbreviation]]

   diseaseOrSubtype_in_correct_order = sn_diseaseOrSubtype_in_correct_order
   subtype_pathway_directory = paste(overall_data_directory,dataset_cellType,subClassLabel,"/",dataset_cellType,subClassLabel,"_pathways/",sep='')
   disease_pathway_directory = paste(overall_data_directory,dataset_cellType,diseaseLabel,"/",dataset_cellType,diseaseLabel,"_pathways/",sep='')
   
   hfFunctions_directory = subtype_pathway_directory
   hf_fileName = paste("WcfDefs_",dataset_cellType,subClassLabel,".txt",sep='')
   complete_hf_fileName = paste(hfFunctions_directory,hf_fileName,sep='')
   if (file.exists(complete_hf_fileName))
   {#Begin - Read hf definitions
       hf_definitions = read.csv(file=complete_hf_fileName,header=TRUE,stringsAsFactors = FALSE,sep='\t')
       
       hfFunctions_directory = disease_pathway_directory
       disease_hf_fileName = paste("WcfDefs_",dataset_cellType,diseaseLabel,".txt",sep='')
       complete_disease_hf_fileName = paste(hfFunctions_directory,disease_hf_fileName,sep='')
       if (file.exists(complete_disease_hf_fileName))
       {#Begin
          disease_hf_definitions = read.csv(file=complete_disease_hf_fileName,header=TRUE,stringsAsFactors = FALSE,sep='\t')
          indexColKeep = which(colnames(disease_hf_definitions) %in% colnames(hf_definitions))
          disease_hf_definitions = disease_hf_definitions[,indexColKeep]
          indexColKeep = which(colnames(hf_definitions) %in% colnames(disease_hf_definitions))
          hf_definitions = hf_definitions[,indexColKeep]
          hf_definitions = rbind(hf_definitions,disease_hf_definitions)
       }#End
   }#End - Read hf definitions
   
   {#Begin - Read and plot disease whole cell functions
      dataset_colors        = c("aquamarine4","red4" ,"orange"   ,"skyblue"  ,"dodgerblue4","midnightblue"       ,"tan4"  ,"magenta"     ,"forestgreen", "limegreen","orange")
      dataset_colors = replicate(length(dataset_colors),"midnightblue")
      names(dataset_colors) = c("SNboth_PODPEC"  ,"SNboth_PT","SNboth_DTLATL","SNboth_TAL" ,"SNboth_DCT"   ,"SNboth_CDCNT"   ,"SNboth_FIB","SNboth_EC","SNboth_Myeloid","SNboth_Lymphoid", "Bulk_SOX4KD")

      diseaseOrSubtype_in_correct_order = datasetAbbreviation_diseaseOrSubtype_in_correct_order[["SN"]]
      scps_in_correct_order = rev(c(names(pfGroup_color_list)))
     
      disease_directory = paste(overall_data_directory,dataset_cellType,diseaseLabel,"/",sep='')
      disease_pathway_directory = paste(disease_directory,dataset_cellType,diseaseLabel,"_pathways/",sep='')   
      diseasePathwayGroup_fileName = paste("WholeCellFunctions_",dataset_cellType,diseaseLabel,".txt",sep='')
      complete_disease_pathwayGroup_fileName = paste(disease_pathway_directory,diseasePathwayGroup_fileName,sep='')
      if (file.exists(complete_disease_pathwayGroup_fileName))
      {#Begin
         disease_pathwayGroups = read.csv(complete_disease_pathwayGroup_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
         dataset_names = unique(disease_pathwayGroups$Dataset.name)
         add = c()
         indexDN=3
         for (indexDN in 1:length(dataset_names))
         {#Begin
             dataset_name = dataset_names[indexDN]
             indexCurrentDN = which(disease_pathwayGroups$Dataset.name==dataset_name)
             dn_disease_pathwayGroups = disease_pathwayGroups[indexCurrentDN,]
             indexMissingScps = which(!scps_in_correct_order %in% dn_disease_pathwayGroups$SCP)
             missing_scps = scps_in_correct_order[indexMissingScps]
             indexM=1
             if (length(missing_scps)>0)
             {#Begin
                for (indexM in 1:length(missing_scps))
                {#Begin
                   missing_scp = missing_scps[indexM]
                   add_lines = as.data.frame(array(NA,c(length(dataset_name),dim(dn_disease_pathwayGroups)[2]),dimnames = list(dataset_name,colnames(dn_disease_pathwayGroups))))
                   add_lines$Dataset.name = dataset_name
                   add_lines$SCP = missing_scp
                   if (length(add[,1])>0) { add = rbind(add,add_lines) }
                   else { add = add_lines }
                }#End
             }#End
         }#End
         if (length(add[,1])>0) { disease_pathwayGroups = rbind(disease_pathwayGroups,add) }
         indexKeepScp = which(disease_pathwayGroups$SCP %in% considered_disease_pf_groups)
         disease_pathwayGroups = disease_pathwayGroups[indexKeepScp,]
         
         disease_pathwayGroups$Dataset.name_factor = factor(disease_pathwayGroups$Dataset.name,levels=diseaseOrSubtype_in_correct_order)
         disease_pathwayGroups$SCP_factor = factor(disease_pathwayGroups$SCP,levels=scps_in_correct_order)
         indexNA = which(is.na(disease_pathwayGroups$SCP_factor))
         if (length(indexNA)>0) { stop(paste("Not matchining pf groups for ",dataset_cellType,": ",unique(disease_pathwayGroups$SCP[indexNA]),sep='')) }
         disease_pathwayGroups$Label = round(disease_pathwayGroups$Minus.log10_pvalue*10)/10
         indexZero = which(disease_pathwayGroups$Minus.log10_pvalue==0)
         disease_pathwayGroups$Label[indexZero] = ""
         disease_pathwayGroups$Minus.log10_pvalue_fill = disease_pathwayGroups$Minus.log10_pvalue
         indexAboveCutoff = which(disease_pathwayGroups$Minus.log10_pvalue>=max_minLog10Pvalue)
         disease_pathwayGroups$Minus.log10_pvalue_fill[indexAboveCutoff] = max_minLog10Pvalue
         
         scp_colors = rev(unlist(pfGroup_color_list))
         indexKeep = which(names(scp_colors) %in% disease_pathwayGroups$SCP)
         scp_colors = scp_colors[indexKeep]
         
         indexColor = which(names(dataset_colors)==dataset_cellType)
         stopifnot(length(indexColor)==1)
         current_color = dataset_colors[indexColor]

         Heatmap = ggplot(disease_pathwayGroups,aes(y=SCP_factor,x=Dataset.name_factor,fill=Minus.log10_pvalue_fill,label=Label))
         Heatmap = Heatmap + geom_tile(color=current_color)# + geom_text(size=3)
         Heatmap = Heatmap + ggtitle(gsub("_"," ",dataset_cellType))
         Heatmap = Heatmap + theme_classic() + xlab("") + ylab("")
         Heatmap = Heatmap + theme(axis.line = element_blank())
         Heatmap = Heatmap + scale_x_discrete(position = "top")
         Heatmap = Heatmap + theme(axis.text.x = element_text(color="black",angle=90,vjust=0.5,hjust=0))
         Heatmap = Heatmap + theme(axis.text.y = element_text(color=scp_colors,face=2,size=9))
         Heatmap = Heatmap + theme(plot.title = element_text(color="black",vjust=0.5,hjust=0.5))
         Heatmap = Heatmap + scale_fill_gradient2(low="white",mid="white",na.value="white",high=current_color,midpoint=0,limits=c(0,max_minLog10Pvalue))
         Legend = get_legend(Heatmap)
         Heatmap = Heatmap + theme(legend.position = "none")
         pfGroup_heatmaps[[length(pfGroup_heatmaps)+1]] = Heatmap
         pfGroup_heatmaps[[length(pfGroup_heatmaps)+1]] = as_ggplot(Legend)
         
         disease_pathwayGroups$SCP_factor = factor(disease_pathwayGroups$SCP,levels=rev(scps_in_correct_order))
         disease_pathwayGroups$Dataset.name_factor = factor(disease_pathwayGroups$Dataset.name,levels=rev(diseaseOrSubtype_in_correct_order))
         rev_scp_colors = rev(scp_colors)

         Heatmap = ggplot(disease_pathwayGroups,aes(x=SCP_factor,y=Dataset.name_factor,fill=Minus.log10_pvalue_fill,label=Label))
         Heatmap = Heatmap + geom_tile(color=current_color)
         #Heatmap = Heatmap + ggtitle(gsub("_"," ",dataset_cellType))
         Heatmap = Heatmap + theme_classic() + xlab("") + ylab("")
         Heatmap = Heatmap + theme(axis.line = element_blank())
         Heatmap = Heatmap + scale_x_discrete(position = "top")
         Heatmap = Heatmap + guides(fill = guide_colorbar(barheight = unit(3, "mm"),ticks=FALSE, draw.llim = FALSE, draw.ulim = FALSE))
         Heatmap = Heatmap + theme(legend.position = "right", legend.title = element_blank())
         Heatmap = Heatmap + theme(legend.text = element_text(size=3))
         Heatmap = Heatmap + theme(legend.key = element_rect(fill="white",color="black"))
         Heatmap = Heatmap + theme(axis.text.y = element_text(color="black",angle=0,vjust=0.5,hjust=0,size=5,face=2))
         if (indexDCT==1) { Heatmap = Heatmap + theme(axis.text.x = element_text(color=rev_scp_colors,vjust=0.5,hjust=0,angle=90,face=2,size=8)) }
         if (indexDCT>1) { Heatmap = Heatmap + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) }
         #Heatmap = Heatmap + theme(plot.title = element_text(color="black",vjust=0.5,hjust=0.5))
         Heatmap = Heatmap + scale_fill_gradient2(low="white",mid="white",na.value="white",high=current_color,midpoint=0)#,limits=c(0,max_minLog10Pvalue))
         #Heatmap = Heatmap + theme(legend.position = "none")
         #Heatmap = Heatmap + theme(legend.position = "right", legend.title = element_blank())
         pfGroup_heatmapsII[[length(pfGroup_heatmapsII)+1]] = Heatmap
         
         #pfGroup_heatmapsII[[length(pfGroup_heatmapsII)+1]] = as.ggplot(legend)
      }#End
     
     
   }#End - Read and plot disease whole cell functions

   {#Begin - Read and plot disease whole cell functions II
     scps_in_correct_order = rev(scps_in_correct_order)
     dataset_colors        = c("aquamarine4","red4" ,"orange"   ,"skyblue"  ,"dodgerblue4","midnightblue"       ,"tan4"  ,"magenta"     ,"forestgreen", "limegreen","orange")
     names(dataset_colors) = c("SNboth_PODPEC"  ,"SNboth_PT","SNboth_DTLATL","SNboth_TAL" ,"SNboth_DCT"   ,"SNboth_CDCNT"   ,"SNboth_FIB","SNboth_EC","SNboth_Myeloid","SNboth_Lymphoid", "Bulk_SOX4KD")
     
     diseaseOrSubtype_in_correct_order = datasetAbbreviation_diseaseOrSubtype_in_correct_order[["SN"]]
     scps_in_correct_order = rev(c(names(pfGroup_color_list)))
     
     disease_directory = paste(overall_data_directory,dataset_cellType,diseaseLabel,"/",sep='')
     disease_pathway_directory = paste(disease_directory,dataset_cellType,diseaseLabel,"_pathways/",sep='')   
     diseasePathwayGroup_fileName = paste("WholeCellFunctions_",dataset_cellType,diseaseLabel,".txt",sep='')
     complete_disease_pathwayGroup_fileName = paste(disease_pathway_directory,diseasePathwayGroup_fileName,sep='')
     if (file.exists(complete_disease_pathwayGroup_fileName))
     {#Begin
       disease_pathwayGroups = read.csv(complete_disease_pathwayGroup_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
       indexMissingScps = which(!scps_in_correct_order %in% disease_pathwayGroups$SCP)
       missing_scps = scps_in_correct_order[indexMissingScps]
       dataset_names = unique(disease_pathwayGroups$Dataset.name)
       indexM=1
       if (length(missing_scps)>0)
       {#Begin
         for (indexM in 1:length(missing_scps))
         {#Begin
           missing_scp = missing_scps[indexM]
           add_lines = as.data.frame(array(NA,c(length(dataset_names),dim(disease_pathwayGroups)[2]),dimnames = list(dataset_names,colnames(disease_pathwayGroups))))
           add_lines$Dataset.name = dataset_names
           add_lines$SCP = missing_scp
           disease_pathwayGroups = rbind(disease_pathwayGroups,add_lines)
         }#End
       }#End
       indexKeepScp = which(disease_pathwayGroups$SCP %in% considered_disease_pf_groups)
       disease_pathwayGroups = disease_pathwayGroups[indexKeepScp,]

       disease_pathwayGroups$Dataset.name_factor = factor(disease_pathwayGroups$Dataset.name,levels=diseaseOrSubtype_in_correct_order)
       disease_pathwayGroups$SCP_factor = factor(disease_pathwayGroups$SCP,levels=scps_in_correct_order)
       indexNA = which(is.na(disease_pathwayGroups$SCP_factor))
       if (length(indexNA)>0) { stop(paste("Not matchining pf groups for ",dataset_cellType,": ",unique(disease_pathwayGroups$SCP[indexNA]),sep='')) }
       disease_pathwayGroups$Label = round(disease_pathwayGroups$Minus.log10_pvalue*10)/10
       indexZero = which(disease_pathwayGroups$Minus.log10_pvalue==0)
       disease_pathwayGroups$Label[indexZero] = ""
       disease_pathwayGroups$Minus.log10_pvalue_fill = disease_pathwayGroups$Minus.log10_pvalue
       indexAboveCutoff = which(disease_pathwayGroups$Minus.log10_pvalue>=max_minLog10Pvalue)
       disease_pathwayGroups$Minus.log10_pvalue_fill[indexAboveCutoff] = max_minLog10Pvalue
       
       scp_colors = rev(unlist(pfGroup_color_list))
       indexColor = which(names(dataset_colors)==dataset_cellType)
       stopifnot(length(indexColor)==1)
       current_color = dataset_colors[indexColor]
     }#End
   }#End - Read and plot disease whole cell functions II
   
   {#Begin - Read disease and subtype whole cell functions and plot with disease groups
     dataset_colors        = c("aquamarine4","red4" ,"orange"   ,"skyblue"  ,"dodgerblue4","midnightblue"       ,"tan4"  ,"magenta"     ,"forestgreen", "limegreen","orange")
     dataset_colors = replicate(length(dataset_colors),"midnightblue")
     names(dataset_colors) = c("SNboth_PODPEC"  ,"SNboth_PT","SNboth_DTLATL","SNboth_TAL" ,"SNboth_DCT"   ,"SNboth_CDCNT"   ,"SNboth_FIB","SNboth_EC","SNboth_Myeloid","SNboth_Lymphoid", "Bulk_SOX4KD")
     
     diseaseOrSubtype_in_correct_order = datasetAbbreviation_diseaseOrSubtype_in_correct_order[["SN"]]
     scps_in_correct_order = rev(c(names(pfGroup_color_list)))
     
     subtypePathwayGroup_fileName = paste("WholeCellFunctions_",dataset_cellType,subClassLabel,".txt",sep='')
     complete_subtype_pathwayGroup_fileName = paste(subtype_pathway_directory,subtypePathwayGroup_fileName,sep='')
     diseasePathway_fileName =  paste("WholeCellFunctions_",dataset_cellType,diseaseLabel,".txt",sep='')
     complete_disease_pathwayGroup_fileName = paste(disease_pathway_directory,diseasePathway_fileName,sep='')
     if (file.exists(complete_subtype_pathwayGroup_fileName))
     {#Begin
       subtype_pathwayGroups = read.csv(complete_subtype_pathwayGroup_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
       disease_pathwayGroups = read.csv(complete_disease_pathwayGroup_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
       indexKeepCol = which(colnames(disease_pathwayGroups) %in% colnames(subtype_pathwayGroups))
       disease_pathwayGroups = disease_pathwayGroups[,indexKeepCol]
       combined_pathwayGroups = rbind(subtype_pathwayGroups,disease_pathwayGroups)
       
       Col_names = colnames(combined_pathwayGroups)
       Col_length = length(Col_names)
       Row_names = ""
       Row_length = length(Row_names)
       border_row_line = as.data.frame(array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names)))
       border_row_line$Dataset.name = "border_row"
       border_row_line$SCP = combined_pathwayGroups$SCP[1]
       combined_pathwayGroups = rbind(combined_pathwayGroups,border_row_line)
       
       
       dataset_names = unique(combined_pathwayGroups$Dataset.name)
       add = c()
       indexDN=3
       for (indexDN in 1:length(dataset_names))
       {#Begin
         dataset_name = dataset_names[indexDN]
         indexCurrentDN = which(combined_pathwayGroups$Dataset.name==dataset_name)
         dn_combined_pathwayGroups = combined_pathwayGroups[indexCurrentDN,]
         indexMissingScps = which(!scps_in_correct_order %in% dn_combined_pathwayGroups$SCP)
         missing_scps = scps_in_correct_order[indexMissingScps]
         indexM=1
         if (length(missing_scps)>0)
         {#Begin
           for (indexM in 1:length(missing_scps))
           {#Begin
             missing_scp = missing_scps[indexM]
             add_lines = as.data.frame(array(NA,c(length(dataset_name),dim(dn_combined_pathwayGroups)[2]),dimnames = list(dataset_name,colnames(dn_combined_pathwayGroups))))
             add_lines$Dataset.name = dataset_name
             add_lines$SCP = missing_scp
             if (length(add[,1])>0) { add = rbind(add,add_lines) }
             else { add = add_lines }
           }#End
         }#End
       }#End
       if (length(add[,1])>0) { combined_pathwayGroups = rbind(combined_pathwayGroups,add) }

       cellSubTypes_count = length(unique(combined_pathwayGroups$Dataset.name))
       missing_columns_count = minColCount_subtype_pfGroups - cellSubTypes_count
       current_diseaseOrSubtype_in_correct_order = diseaseOrSubtype_in_correct_order
       if (missing_columns_count)
       {#Begin
         new_datasets = paste("Empty",1:missing_columns_count,sep='')
         current_diseaseOrSubtype_in_correct_order = c(current_diseaseOrSubtype_in_correct_order,new_datasets)
         Col_names = colnames(combined_pathwayGroups)
         Col_length = length(Col_names)
         Row_names = 1:missing_columns_count
         Row_length = length(Row_names)
         add = as.data.frame(array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names)))
         add$SCP = combined_pathwayGroups$SCP[1]
         add$Dataset.name = new_datasets
         add$Minus.log10_pvalue = 0
         combined_pathwayGroups = rbind(combined_pathwayGroups,add)
       }#End

       combined_pathwayGroups$Dataset.name_factor = factor(combined_pathwayGroups$Dataset.name,levels=current_diseaseOrSubtype_in_correct_order)
       combined_pathwayGroups$SCP_factor = factor(combined_pathwayGroups$SCP,levels=scps_in_correct_order)
       indexNA = which(is.na(combined_pathwayGroups$SCP_factor))
       if (length(indexNA)>0) { stop(paste("Not matchining pf groups for ",dataset_cellType,": ",unique(combined_pathwayGroups$SCP[indexNA]),sep='')) }
       combined_pathwayGroups$Label = round(combined_pathwayGroups$Minus.log10_pvalue*10)/10
       indexZero = which(combined_pathwayGroups$Minus.log10_pvalue==0)
       combined_pathwayGroups$Label[indexZero] = ""
       combined_pathwayGroups$Minus.log10_pvalue_fill = combined_pathwayGroups$Minus.log10_pvalue
       indexAboveCutoff = which(combined_pathwayGroups$Minus.log10_pvalue>=max_minLog10Pvalue)
       combined_pathwayGroups$Minus.log10_pvalue_fill[indexAboveCutoff] = max_minLog10Pvalue
       
       scp_colors = rev(unlist(pfGroup_color_list))
       indexKeep = which(names(scp_colors) %in% combined_pathwayGroups$SCP)
       scp_colors = scp_colors[indexKeep]
       
       hv_line_size = 0.5

       
       indexColor = which(names(dataset_colors)==dataset_cellType)
       stopifnot(length(indexColor)==1)
       current_color = dataset_colors[indexColor]
       
       Heatmap = ggplot(combined_pathwayGroups,aes(y=SCP_factor,x=Dataset.name_factor,fill=Minus.log10_pvalue_fill,label=Label))
       Heatmap = Heatmap + geom_tile(color=current_color)# + geom_text(size=3)
       Heatmap = Heatmap + ggtitle(gsub("_"," ",dataset_cellType))
       Heatmap = Heatmap + theme_classic() + xlab("") + ylab("")
       Heatmap = Heatmap + theme(axis.line = element_blank())
       Heatmap = Heatmap + scale_x_discrete(position = "top")
       Heatmap = Heatmap + theme(axis.text.x = element_text(color="black",angle=90,vjust=0.5,hjust=0))
       Heatmap = Heatmap + theme(axis.text.y = element_text(color=scp_colors,face=2,size=9))
       Heatmap = Heatmap + theme(plot.title = element_text(color="black",vjust=0.5,hjust=0.5))
       Heatmap = Heatmap + scale_fill_gradient2(low="white",mid="white",na.value="white",high=current_color,midpoint=0,limits=c(0,max_minLog10Pvalue))
       
       vLineXintercepts = c(datasetCellType_vLineXintercepts_list[[dataset_cellType]])
       segment_df = c()
       if (length(vLineXintercepts)>0)
       {#Begin
         vsegment_df = data.frame(x = 0.5 + vLineXintercepts,
                                  xend = 0.5 + vLineXintercepts,
                                  y = 0.5,
                                  yend = length(unique(combined_pathwayGroups$SCP)) + 0.5)           
         vsegment_df$Minus.log10_pvalue_fill = 0
         vsegment_df$Label = ""
         vsegment_df$SCP_factor = ""
         vsegment_df$Dataset.name_factor = ""
         
         hsegment_df = data.frame(y = c(0.5,length(unique(combined_pathwayGroups$SCP))+0.5),
                                  yend = c(0.5,length(unique(combined_pathwayGroups$SCP))+0.5),
                                  x = 0.5,
                                  xend = cellSubTypes_count + 0.5)
         hsegment_df$Minus.log10_pvalue_fill = 0
         hsegment_df$Label = ""
         hsegment_df$SCP_factor = ""
         hsegment_df$Dataset.name_factor = ""
         segment_df = rbind(hsegment_df,vsegment_df)
         
         Segments = geom_segment(data = segment_df, aes(x = x, xend = xend, y = y, yend = yend),linewidth = hv_line_size)
         Heatmap = Heatmap + Segments
       }#End
       indexColor=1
       Legend = get_legend(Heatmap)
       Heatmap = Heatmap + theme(legend.position = "none")
       subtype_pfGroup_heatmaps[[length(subtype_pfGroup_heatmaps)+1]] = Heatmap
       subtype_pfGroup_heatmaps[[length(subtype_pfGroup_heatmaps)+1]] = as_ggplot(Legend)
       
       datasets = unique(combined_pathwayGroups$Dataset.name)
       indexRemove = c(grep("Empty",datasets),grep("border_row",datasets),which(datasets==""))
       indexKeep = 1:length(datasets)
       indexKeep = indexKeep[!indexKeep %in% indexRemove]
       datasets = datasets[indexKeep]
       indexDataset = 1
       print(paste0("Calculate total min log10(p) for each subtype disease for ",cellType))
       for (indexDataset in 1:length(datasets))
       {#Begin
          dataset = datasets[indexDataset]
          indexCurrentDataset = which(combined_pathwayGroups$Dataset.name==dataset)
          all_minLog10Pvalues = sum(combined_pathwayGroups$Minus.log10_pvalue[indexCurrentDataset],na.rm = TRUE)
          print(paste0(dataset,": ",all_minLog10Pvalues))
       }#End
       print("")
       

       #pfGroup_heatmapsII[[length(pfGroup_heatmapsII)+1]] = as.ggplot(legend)
     }#End - Read and combined subtype pf groups
     
     
   }#End - Read disease and subtype whole cell functions and plot with disease groups
}#End

   
complete_pdf_fileName = paste(overall_data_directory,"WholeCellFunctions_accross_diseases_subtypes",".pdf",sep='') 
Generate_plots(subtype_pfGroup_heatmaps,complete_pdf_fileName,3,1)

complete_pdf_fileName = paste(overall_data_directory,"WholeCellFunctions_inFunctionalUnits_accross_diseases",".pdf",sep='') 
Generate_plots(pfGroup_heatmaps,complete_pdf_fileName,5,2)


complete_pdf_fileName = paste(overall_data_directory,"WholeCellFunctions_inFunctionalUnits_accross_diseases2",".pdf",sep='') 
#pdf(complete_pdf_fileName,width=8.27,height=11.69)

pdf(complete_pdf_fileName,width=3,height=7)

wrapped_rows <- lapply(pfGroup_heatmapsII, function(p) {
  p + theme(legend.position = "right", legend.title = element_blank()) +
    plot_layout(widths = c(5, 1))
})

# Stack those horizontally combined plots vertically
Final_plot <- wrap_plots(wrapped_rows, ncol = 1)
Final_plot = Final_plot + plot_layout(heights = rep(1, length(pfGroup_heatmapsII)))
Final_plot = Final_plot + plot_layout(ncol = 1) & theme(plot.margin = margin(0, 0, 0, 0))

print(Final_plot)
dev.off()





