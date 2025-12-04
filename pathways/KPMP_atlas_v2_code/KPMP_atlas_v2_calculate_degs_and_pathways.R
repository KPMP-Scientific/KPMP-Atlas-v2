#This script was written by Jens Hansen working for the Ravi Iyengarlab (Icahn School of Medicine at Mount Sinai, New York) and the Kidney Precision Medicine Project.

#Please acknoledge our work by citing Lake et. al bioRxiv 2025.09.26.678707.

#This script was tested in a Linux environment.
#It calls our application MBC PathNet that can be downloaded from mbc-ontology.org (github.com/SBCNY/Molecular-Biology-of-the-Cell). The directory of MBC PathNet can be specified in 'User input is needed here'.
#It also calls an executable file that was written in C#. The executable file is saved in the directory complete_calculateAvgPostHocPowerAnalysis_script_fileName that is specified in 'Information about necessary directories and files'.
#All other files that this script reads are also specified in 'Information about necessary directories and files'.
#If this script is run in PostHocPower mode (by setting calculateFinalDegs1_doPostHocPower2 in 'User input is needed here' to 2), the user should also adjust the parameter 'maxNoDegListsInMem' in 'Set parameter for PostHocPower analysis' to a value matching his memory availability.

#Progress updates are printed into the console or written into txt files into the directory 'AA_Documentations' that can be found in the related results directory (within the base_directory that can be defined in 'User input is needed here').

#To enable running  the script complete_calculateAvgPostHocPowerAnalysis_script_fileName:
#Install Microsoft package repository and keys
#wget https://packages.microsoft.com/config/ubuntu/22.04/packages-microsoft-prod.deb -O packages-microsoft-prod.deb
#sudo dpkg -i packages-microsoft-prod.deb
#rm packages-microsoft-prod.deb

# Update package lists
#sudo apt update

# Install .NET 8 runtime
#sudo apt install -y dotnet-runtime-8.0

unregister <- function() {
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)
  }    
  gc()
  unlink(".RData")
  rm(list = ls(all.names=TRUE));
  gc()

  ################# See 'User input is needed here' below, after 'Open libraries and document versions' #################
  
  Get_sessionInfo_summary_table = function()
  {#Begin - Get_sessionInfo_summary_table
    Col_names = c("Group","Library","Field","Entry")
    Col_length = length(Col_names)
    Row_names = 1
    Row_length = length(Row_names)
    sessionInfo_summary_baseLine = as.data.frame(array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names)))
    sessionInfo_summary = c();
    
    sessionInfos_fields_list = list( "R.version" = "combined"
                                     ,"platform" = "combined"
                                     ,"locale" = "combined"
                                     ,"running" = "combined"
                                     ,"RNGkind" = "combined"
                                     ,"basePkgs" = "combined"
                                     ,"otherPkgs" = c("Version","URL","Date/Publication","Built")
                                     ,"loadedOnly" = c("Version","URL","Date/Publication","Built")
                                     ,"matprod" = "combined"
                                     ,"BLAS" = "combined"
                                     ,"LAPACK" = "combined"
                                     ,"system.codepage" = "combined"
                                     ,"codepage" = "combined")
    
    r_sessionInfo = sessionInfo()
    r_sessionInfo_length = length(r_sessionInfo)
    indexSession=1
    for (indexSession in 1:r_sessionInfo_length)
    {#Begin
      current_infoName = names(r_sessionInfo)[indexSession]
      current_infos = r_sessionInfo[[current_infoName]]
      indexInfo=1
      fields_of_interest = sessionInfos_fields_list[[current_infoName]]
      if ((!is.null(fields_of_interest))&(length(fields_of_interest)>0))
      {#Begin
        for (indexInfo in 1:length(current_infos))
        {#Begin
          current_info = current_infos[indexInfo]
          if (fields_of_interest[1]=="combined")
          {#Begin
            sessionInfo_summary_line = sessionInfo_summary_baseLine
            sessionInfo_summary_line$Group = current_infoName
            if (!is.null(names(current_info)[1])){ sessionInfo_summary_line$Field = names(current_info)[1] }
            else { sessionInfo_summary_line$Field = "No given name" }
            sessionInfo_summary_line$Library = current_info[[1]]
            if (length(sessionInfo_summary_line$Library)>1) { rm(sessionInfos_fields_list) }
            if (length(sessionInfo_summary)>0) { sessionInfo_summary = rbind(sessionInfo_summary,sessionInfo_summary_line)}
            else { sessionInfo_summary = sessionInfo_summary_line }
          }#End
          else
          {#Begin
            for (indexField in 1:length(fields_of_interest))
            {#Begin
              field_of_interest = fields_of_interest[indexField]
              sessionInfo_summary_line = sessionInfo_summary_baseLine
              sessionInfo_summary_line$Group = current_infoName
              sessionInfo_summary_line$Field = field_of_interest
              sessionInfo_summary_line$Library = current_info[[1]][["Package"]]
              if (length(sessionInfo_summary_line$Library)>1) { rm(sessionInfos_fields_list) }
              if (!is.null(current_info[[1]][[field_of_interest]]))
              { sessionInfo_summary_line$Entry = current_info[[1]][[field_of_interest]] }
              if (length(sessionInfo_summary)>0) { sessionInfo_summary = rbind(sessionInfo_summary,sessionInfo_summary_line)}
              else { sessionInfo_summary = sessionInfo_summary_line }
            }#End
          }#End
        }#End
      }#End
    }#End
    sessionInfo_summary$Group = gsub("\t"," ",sessionInfo_summary$Group)
    sessionInfo_summary$Library = gsub("\t"," ",sessionInfo_summary$Library)
    sessionInfo_summary$Field = gsub("\t"," ",sessionInfo_summary$Field)
    sessionInfo_summary$Entry = gsub("\t"," ",sessionInfo_summary$Entry)
    sessionInfo_summary$Group = gsub("\r\n"," ",sessionInfo_summary$Group)
    sessionInfo_summary$Library = gsub("\r\n"," ",sessionInfo_summary$Library)
    sessionInfo_summary$Field = gsub("\r\n"," ",sessionInfo_summary$Field)
    sessionInfo_summary$Entry = gsub("\r\n"," ",sessionInfo_summary$Entry)
    sessionInfo_summary$Group = gsub("\n"," ",sessionInfo_summary$Group)
    sessionInfo_summary$Library = gsub("\n"," ",sessionInfo_summary$Library)
    sessionInfo_summary$Field = gsub("\n"," ",sessionInfo_summary$Field)
    sessionInfo_summary$Entry = gsub("\n"," ",sessionInfo_summary$Entry)
    return (sessionInfo_summary)
  }#End - Get_sessionInfo_summary_table

  {#Begin - Open libraries and document versions - BEGIN
    Col_names = c("Library","Version")
    Col_length = length(Col_names)
    Row_names = 1
    Row_length= length(Row_names)
    version_documentation_line = array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
    version_documentation_line = as.data.frame(version_documentation_line)
    version_documentations = c()
    
    libraries = c("future","matrixStats","progress","cowplot","scales","SingleCellExperiment","RColorBrewer","Seurat","SeuratDisk","Signac","SeuratObject","BPCells",
                  "harmony","gridExtra","ggplot2","ggridges","sctransform","mclust","dplyr","beeswarm","dendextend","fields","reshape2","gplots","ggnewscale","genefilter",
                  "clustree","Matrix","gridExtra","reticulate","umap","stringr","rhdf5","progress","doParallel","scales","ggbeeswarm","progress","readxl","openxlsx"
    )
    for (indexL in 1:length(libraries))
    {#Begin
      current_library = libraries[indexL]
      library(current_library,character.only=TRUE)
      new_version_documentation_line = version_documentation_line
      new_version_documentation_line$Library = current_library
      new_version_documentation_line$Version = packageVersion(current_library)
      if (length(version_documentations)==0)
      {#Begin
        version_documentations = new_version_documentation_line
      }#End
      else
      {#Begin
        version_documentations = rbind(version_documentations,new_version_documentation_line)
      }#End
    }#End
    sessionInfo_summary = Get_sessionInfo_summary_table()
  }#End - Open libraries and document versions
  
  {#Begin- User input is needed here
     calculateFinalDegs1_doPostHocPower2 = 2 #if run in PostHocPower mode the results will be saved in a directory ending with PHP (as specified with the variable 'add_to_results_directory')
     is_windows = FALSE #Script was only tested and DEGs were only calculated in a LINUX environment
     is_linux = !is_windows
     if (calculateFinalDegs1_doPostHocPower2==1)
     { max_cores_count = 38 }  #(for regular analysis, we suggest less cores, since there are only 25 tasks and uploading data into the cores takes some time as well)
     if (calculateFinalDegs1_doPostHocPower2==2)
     { max_cores_count = 38; }
     tmp_directory = "/data/tmp/"
     base_directory = paste("/data/KPMP_v2_atlas/KPMP_v2_data/") #ensure that directories end with "/"
     mbc_pathNet_directory = "/data/MBCO_windows_application/" #Download MBCPathNet from 'mbc-ontology.org'/'https://github.com/SBCNY/Molecular-Biology-of-the-Cell' and copy the content of the folder 'MBCO_windows_application' (windows is just a name) it into the specified directory,
                                                               #for use in LINUX, download mono as described in the ReadMe file found at the link above:
                                                                   #sudo apt update
                                                                   #sudo apt install mono-xsp4
                                                               #see iyengarlab.org/mbcpathnet for summary information about MBC PathNet
     #mbc_pathNet_directory = "D:/MBCO_windows_application/"
     directory_of_this_script = "/data/KPMP_disease_atlas/KPMP_atlas_v2_code/";
     name_of_this_script = "KPMP_atlas_v2_calculate_degs_and_pathways.R"
     working_directory = "/data/KPMP_v2_atlas/KPMP_atlas_v2_code/"; #has to contain the script 'KPMP_atlas_v2_shared_code.R'
     ##download the KPMP v2 atlas Seurat object and copy it into the 'SN_RNAseq_atlas_v2_2024November01' subdirectory as specified in 'Information about necessary directories and files' below
     ##adjust file name of and read function for the downloadable h5ad file in 'read_sn_rnaSeq'
     
     ##ensure that the user has read and write rights in directories
     ##Please also see Set parameter for PostHocPower analysis (if calculateFinalDegs1_doPostHocPower2==2) to adjust parameters to available memory on used computer
  }#End- User input is needed here
  
  {#Begin - Information about necessary directories and files
     complete_calculateAvgPostHocPowerAnalysis_script_fileName = paste(directory_of_this_script,"Average_DEGs_and_do_postHocPowerAnalysis_linux_x64/Average_DEGs_and_do_postHocPowerAnalysis",sep='') #no file extension
     #directories and files in base_directory
     #--KPMP_atlas_v2_reference_WcfDefs - these files are all assigned in assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list (see below)
     #----WcfDefs_SNboth_CDCNTSimD.xlsx    #SimD = simplified disease
     #----WcfDefs_SNboth_DCTSimD.xlsx    
     #----WcfDefs_SNboth_DCTsl3.xlsx       #sl3 = cell type/suptype level 3 (v2.subclass.l3)
     #----WcfDefs_SNboth_DTLATLSimD.xlsx
     #----WcfDefs_SNboth_ECSimD.xlsx
     #----WcfDefs_SNboth_FIBSimD.xlsx
     #----WcfDefs_SNboth_FIBsl3.xlsx
     #----WcfDefs_SNboth_LymphoidSimD.xlsx
     #----WcfDefs_SNboth_MyeloidSimD.xlsx
     #----WcfDefs_SNboth_Myeloidsl3.xlsx
     #----WcfDefs_SNboth_PODPECSimD.xlsx
     #----WcfDefs_SNboth_PODPECsl3.xlsx
     #----WcfDefs_SNboth_PTSimD.xlsx
     #----WcfDefs_SNboth_PTsl3.xlsx
     #----WcfDefs_SNboth_TALSimD.xlsx
     #----WcfDefs_SNboth_TALsl3.xlsx
     #--SN_RNAseq_atlas_v2_2024November01  --------------------------------- this directory needs to be filled by the user
     #----full_kidney_count_set_0424
     #------col_names
     #------idxptr
     #------index_data
     #------index_idx
     #------index_idx_offsets
     #------index_starts
     #------row_names
     #------shape
     #------storage_order
     #------val_data
     #------val_idx
     #------val_idx_offsets
     #------version
     #--Kidney_AtlasV2_Seurat_11012024.rds
     #--Metadata
     #----20250606_OpenAccessClinicalData.csv
  }#End - Information about necessary directories and files
 
  stopifnot(dir.exists(mbc_pathNet_directory))
  stopifnot(dir.exists(working_directory))
  setwd(working_directory)  
  source('KPMP_atlas_v2_shared_code.R')

  {#Begin - Define colors for whole cell functions
      wcfGroup_color_list = list()
      wcfGroup_color_list[["Metabolism"]] = "skyblue"
      wcfGroup_color_list[["TM movement of ions, water and solutes"]] = "seagreen3"
      wcfGroup_color_list[["Signaling networks reg. metabolism/TM movement"]]= "blue"
      wcfGroup_color_list[["Redox, iron and mitochondrial dynamics"]]= "purple4"
      wcfGroup_color_list[["Cell structure dynamics (Cytoskeleton, PM and nucleus)"]] = "red1"
      wcfGroup_color_list[["Cell adhesion"]] = "lightsalmon1"
      wcfGroup_color_list[["ECM remodeling"]] = "tan4"
      wcfGroup_color_list[["Signaling networks reg. cell/tissue structure"]] = "darkorange2"
      wcfGroup_color_list[["Immune response and signaling"]] = "gold2"
      wcfGroup_color_list[["Gene expression"]] = "magenta1"
      wcfGroup_color_list[["Intracellular degradation"]] = "magenta4"
      wcfGroup_color_list[["Cell cycle and DNA dynamics"]] = "gray50"
      wcfGroup_color_list[["Cell death"]] = "black"
      wcfGroup_color_list[["Vesicle traffic"]] = "plum"
  }#End - Define colors for whole cell functions
  
  findMarkers_logfc_threshold = 0.1
  findMarkers_adj_pvalue_cutoff = 0.05
  findMarkers_topDEGs = 500
  findMarkers_max_cells_per_ident = 5000
  pathway_pvalue_cutoff = 0.05
  simplified_diseases_of_interest = c("HRT","AKI","CKD")
  split_simplfied_disease_CKD_based_on_eGFR = TRUE
  if (split_simplfied_disease_CKD_based_on_eGFR)
  {#Begin
      simplified_diseases_of_interest = c("AKI","CKD_low_eGFR","CKD_high_eGFR")
  }#End

  if (!dir.exists(tmp_directory)) { dir.create(tmp_directory,recursive = TRUE) }
  Sys.setenv(TMPDIR = tmp_directory)
  

  if (calculateFinalDegs1_doPostHocPower2==1)
  { add_to_results_directory = paste0(findMarkers_topDEGs,"DEGs") }
  if (calculateFinalDegs1_doPostHocPower2==2)
  { add_to_results_directory = paste0(findMarkers_topDEGs,"DEGs_PHP") }
  #add_to_results_directory = "_soumya_DKD"
  #add_to_results_directory = "_cijiang"
  
  level1_subdirectory = paste("KPMP_atlas_v2_max",findMarkers_max_cells_per_ident,"cells_",add_to_results_directory,"/",sep='')
  
  #stopifnot(file.exists(complete_calculateAvgPostHocPowerAnalysis_script_fileName))
  
  mbc_pathNet_exe_fileName = "MBC_PathNet.exe"
  reference_wholeCellFunctionDef_directory = paste(base_directory,"KPMP_atlas_v2_reference_WcfDefs/",sep='')
  complete_level1_directory = paste(base_directory,level1_subdirectory,sep='')
  if (!dir.exists(complete_level1_directory)) { dir.create(complete_level1_directory) }
  documentation_directory = paste(base_directory,level1_subdirectory,"AA_Documentations/",sep='')
  if (!dir.exists(documentation_directory)) { dir.create(documentation_directory) }
  subtype_level = "sl3" # "sl2", "sl1", "sl3" (For most comparisons, sl1 will only have one subtype/state, containing all cells of the related cell type)
  delete_existing_degs = TRUE
  calculate_all_markers_for_subtypes = TRUE
  calculate_pathways_for_subtypes = TRUE
  calculate_all_markers_for_simplified_diseases = TRUE & (subtype_level=="sl3") #to avoid duplicated calculations
  calculate_pathways_for_simplified_diseases = TRUE & (subtype_level=="sl3") 
  read_sn_rnaSeq_altlas_v2_2024November = TRUE

  {#Begin - Write session info into documentation folder and copy this script into folder
    complete_this_scriptName = paste(directory_of_this_script,name_of_this_script,sep='')
    complete_copy_this_scriptName = paste(documentation_directory,name_of_this_script,sep='')
    file.copy(complete_this_scriptName,complete_copy_this_scriptName)
    complete_fileName = paste(documentation_directory,"SessionInfo_summary.txt",sep='')
    write.table(sessionInfo_summary,file=complete_fileName,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
  }#End - Write session info into documentation folder and copy this script into folder

  if (calculateFinalDegs1_doPostHocPower2==1)
  {#Begin - Set parameter for final analysis
     finalAnalysis_numberOfDatasets = 25   
     postHocPower_max_degListsCombinations = 0#100 #if set to zero postHocPower will be skipped
     postHocPower_minPercent_randomSeedNos = 1
     postHocPower_max_number_of_datasets = 0#100 #if set to zero postHocPower will be skipped
     random_seed_nos = 1:finalAnalysis_numberOfDatasets
     randomSeedNodeNoPostHocPower = -1
     maxNoDegListsInMem = 5000; #See explanation below
  }#End - Set parameter for final analysis
  if (calculateFinalDegs1_doPostHocPower2==2)
  {#Begin - Set parameter for PostHocPower analysis
     finalAnalysis_numberOfDatasets = 25  
     postHocPower_max_degListsCombinations = 100 #if set to zero postHocPower will be skipped
     postHocPower_minPercent_randomSeedNos = 1
     postHocPower_max_number_of_datasets = 100 #if set to zero postHocPower will be skipped
     random_seed_nos = 1:10000   # must be continuous and start at one
     randomSeedNodeNoPostHocPower = 4335 #-1 = no pseudo randomSeedNode
     maxNoDegListsInMem = 5000; #Determines how many generated DEGs lists (linked to the random_seed_nos) will be simultaneously saved in memory. If lower than the length of random_seed_nos, the script (complete_calculateAvgPostHocPowerAnalysis_script_fileName)
                                #will only keep this number of DEG lists in memory and replace existing DEG lists by DEG lists of interest, if necessary
                                #If this number is too high (for a given cell type), the computer might freeze and might need a cold restart. In this case the directory 'AA_Documentations' within the specified results directory will contain a progress report file that gives information
                                #about how many DEG lists were read before the available memory was exceeded. Though the number of storable DEG lists depends on the number of subtypes of a given cell type, we only offer one maxNoDegListsInMem for all cell types for simplicity.
                                #The number 5000 was evaluated on our computer with 400GB memory. Our compuer frooze after uploading > 7500 DEG lists when analyzing the PT DEGs.
  }#End - Set parameter for PostHocPower analysis
  
  dataset_assays = c("Multiome","SNonly")
  dataset_assays = c("SNboth")
  #dataset_assays = c("SNboth","Multiome","SNonly")
  
  datasetAssay_snRNAseqAssay_list = list()
  datasetAssay_snRNAseqAssay_list[["Multiome"]] = "10X Multiome"
  datasetAssay_snRNAseqAssay_list[["SNonly"]] = "10X snRNA-seq"
  datasetAssay_snRNAseqAssay_list[["SNboth"]] = "allSn"

  combined_results_directory = paste(base_directory,level1_subdirectory,"AA_combined_results/",sep='')
  if (!dir.exists(combined_results_directory)) { dir.create(combined_results_directory) }

  diseaseFocusLabel = "SimD"
  considered_pathway_fileName = "Pathways_sig_forWholeCellFunct.txt"
  
  ontologies = c("Mbco")#,"Go_bp","Reactome") #"Mbco","Go_bp","Go_mf","Go_cc","Reactome","Custom_1","Custom_2"

  cellTypeGroup_cellTypes_list = list()   
  cellTypeGroup_cellTypes_list[["PT"]] = "PT"
  cellTypeGroup_cellTypes_list[["TAL"]] = c("TAL")
  cellTypeGroup_cellTypes_list[["PODPEC"]] = c("PEC","POD")
  cellTypeGroup_cellTypes_list[["DTLATL"]] = c("DTL","ATL")
  cellTypeGroup_cellTypes_list[["DCT"]] = "DCT"
  cellTypeGroup_cellTypes_list[["CDCNT"]] = c("CNT","PC","IC")
  cellTypeGroup_cellTypes_list[["FIB"]] = "FIB"
  cellTypeGroup_cellTypes_list[["Myeloid"]] = "Myeloid"
  cellTypeGroup_cellTypes_list[["EC"]] = "EC"
  cellTypeGroup_cellTypes_list[["Lymphoid"]] = "Lymphoid"

  {#Begin - Define identity_cellCounts_line
     Col_names = c("Identity","Subtype_level","Dataset_assay","Cell_count")
     Col_length = length(Col_names)
     Row_names = ""
     Row_length = length(Row_names)
     identity_cellCounts_base_line = as.data.frame(array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names)))
     identity_cellCounts = c()
  }#End - Define identity_cellCounts_line

  assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list = list()
  for (indexDA in 1:length(dataset_assays))
  {#Begin - Fill assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list
    potential_assay = dataset_assays[indexDA]
    
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_PODPEC","sl3",sep='')]] = c("WcfDefs_SNboth_PODPECsl3.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_PT","sl3",sep='')]] = c("WcfDefs_SNboth_PTsl3.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_TAL","sl3",sep='')]] = c("WcfDefs_SNboth_TALsl3.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_DCT","sl3",sep='')]] = c("WcfDefs_SNboth_DCTsl3.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_Myeloid","sl3",sep='')]] = c("WcfDefs_SNboth_Myeloidsl3.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_FIB","sl3",sep='')]] = c("WcfDefs_SNboth_FIBsl3.xlsx")
    
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_PODPEC",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_PODPECSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_DTLATL",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_DTLATLSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_PT",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_PTSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_TAL",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_TALSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_DCT",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_DCTSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_CDCNT",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_CDCNTSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_EC",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_ECSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_Myeloid",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_MyeloidSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_Lymphoid",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_LymphoidSimD.xlsx")
    assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[paste(potential_assay,"_FIB",diseaseFocusLabel,sep='')]] = c("WcfDefs_SNboth_FIBSimD.xlsx")
  }#End - Fill assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list
  
  {#Begin - Check reference file column names
    assayCellTypes = names(assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list)
    indexACT = 1
    for (indexACT in 1:length(assayCellTypes))
    {#Begin
       assayCellType = assayCellTypes[indexACT]
       reference_pf_fileNames = assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[assayCellType]]
       indexRef=1
       for (indexRef in 1:length(reference_pf_fileNames))
       {#Begin
          reference_pf_fileName = reference_pf_fileNames[indexRef]
          complete_ref_pf_fileName = paste(reference_wholeCellFunctionDef_directory,reference_pf_fileName,sep='')
          reference_pf = read_excel(complete_ref_pf_fileName)
          if (!"SCP" %in% names(reference_pf)) { print(paste("'SCP' column is missing in ",reference_pf_fileName,sep='')) }
          if (!"Symbol" %in% names(reference_pf)) { print(paste("'Symbol' column is missing in ",reference_pf_fileName,sep='')) }
       }#End
    }#End
    rm(assayCellType);rm(assayCellTypes);rm(indexACT)
  }#End - Check reference file column names

  options(future.globals.maxSize= 20000000*1024^2)
  if (is_windows) { memory.limit(size=1000000) }

  dist.function = function(x) as.dist((1-cor((x),method="pearson")))
  hclust.function = function(x) hclust(x,method="average")

if (read_sn_rnaSeq_altlas_v2_2024November)
{#Begin - read_sn_rnaSeq
  #Adjust read functions to convert/read the downloadable f337b525-c8f7-4c96-8cfe-f258a9f5ca48.h5ad file
  print(paste("Read SNRNAseq/Multiome seurat object",sep=''))
  directory = paste(base_directory,"/SN_RNAseq_atlas_v2_2024November01/",sep='');
  destDirectory = paste(directory,"full_kidney_count_set_0424/",sep='')
  study_baseName = "SN RNAseq atlas v2 2024November01"
  seurat_object_fileName = "Kidney_AtlasV2_Seurat_11012024.rds" #f337b525-c8f7-4c96-8cfe-f258a9f5ca48.h5ad
  complete_seurat_object_fileName = paste(directory,seurat_object_fileName,sep='')
  seurat_object_all_assays = LoadSeuratRds(complete_seurat_object_fileName) 
  
  counts_matrix = open_matrix_dir(dir = destDirectory)
  #counts_matrix <- Azimuth:::ConvertEnsembleToSymbol(mat = counts_matrix, species = "human")
  
  seurat_object_all_assays[["RNA"]]$counts = counts_matrix
  seurat_object_all_assays = NormalizeData(seurat_object_all_assays,
                                           normalization.method = "LogNormalize",
                                           scale.factor = 10000,
                                           margin = 1)
  
  {#Begin - Split patient 3535 into two patients with NHT and HRT disease condition
    index3535 = which(seurat_object_all_assays$patient=="3535")
    indexNHT = which(seurat_object_all_assays$condition_level3=="NHT")
    indexHRT = which(seurat_object_all_assays$condition_level3=="HRT")
    index3535_nht = intersect(index3535,indexNHT)
    index3535_hrt = intersect(index3535,indexHRT)
    stopifnot(length(index3535_nht)>0)
    stopifnot(length(index3535_hrt)>0)
    stopifnot(length(index3535_nht) + length(index3535_hrt) == length(index3535))
    seurat_object_all_assays$condition_level3[index3535_hrt] = "NHT"
  }#End -  Split patient 3535 into two patients with NHT and HRT disease condition
  
  total_readCounts = unique(round(colSums(expm1(seurat_object_all_assays[["RNA"]]$data))))
  stopifnot(length(total_readCounts)==1)
  
  if (subtype_level=="cluster")
  { seurat_object_all_assays$Cell_subtype = paste(seurat_object_all_assays$v2.clusters,sep='') }
  if (subtype_level=="sl3")
  { seurat_object_all_assays$Cell_subtype = paste(seurat_object_all_assays$v2.subclass.l3,sep='') }
  if (subtype_level=="sl2")
  { seurat_object_all_assays$Cell_subtype = paste(seurat_object_all_assays$v2.subclass.l2,sep='') }
  if (subtype_level=="sl1")
  { seurat_object_all_assays$Cell_subtype = paste(seurat_object_all_assays$v2.subclass.l1,sep='') }
  if (is.null(seurat_object_all_assays$Cell_subtype)) { stop(paste(subtype_level," is not considered for SN RNAseq")) }
  seurat_object_all_assays$Cell_type = seurat_object_all_assays$v2.subclass.l1
  
  seurat_object_all_assays$Disease = seurat_object_all_assays$condition_level3
  
  seurat_object_all_assays$Cell_type = gsub("_","-",seurat_object_all_assays$Cell_type)
  seurat_object_all_assays$Cell_type = gsub("/","",seurat_object_all_assays$Cell_type)
  
  {#Begin - Get cellType_subtypes_list mapping v2.subclass.l1 to cell type
    cellType_subtypes_list = list()
    cellTypes_l1 = unique(seurat_object_all_assays$v2.subclass.l1)
    for (indexCT in 1:length(cellTypes_l1))
    {#Begin
      cellType = cellTypes_l1[indexCT]
      indexCurrentCellType = which(seurat_object_all_assays$v2.subclass.l1==cellType)
      cellType_subtypes_list[[cellType]] = unique(seurat_object_all_assays$Cell_type[indexCurrentCellType])
    }#End
  }#End - Get cellType_subtypes_list mapping v2.subclass.l1 to cell type

  seurat_object_all_assays$Cell_type = as.character(seurat_object_all_assays$Cell_type)
  seurat_object_all_assays$Cell_subtype = as.character(seurat_object_all_assays$Cell_subtype)
}#End - read_sn_rnaSeq

{#Begin - Add exact eGFR from clinical metadata spreadsheet
    complete_metadata_fileName = paste(base_directory,"Metadata/","20250606_OpenAccessClinicalData.csv",sep='')
    metadata = read.csv(complete_metadata_fileName,header=TRUE,stringsAsFactors = FALSE)
    seurat_object_all_assays$baseline_eGFR = -1
    patients = unique(seurat_object_all_assays$patient)
    indexP=1
    print("Add exact eGFR from clinical metadata spreadsheet")
    #pb = progress_bar(total=length(patients))
    for (indexP in 1:length(patients))
    {#Begin
       patient = patients[indexP]
       indexCurrentPatient_seurat = which(seurat_object_all_assays$patient==patient)
       indexCurrentPatient_metadata = which(metadata$Participant.ID==patient)
       stopifnot(length(indexCurrentPatient_metadata)<=1)
       if (length(indexCurrentPatient_metadata)==1)
       {#Begin
          seurat_object_all_assays$baseline_eGFR[indexCurrentPatient_seurat] = metadata$Baseline.eGFR..ml.min.1.73m2.[indexCurrentPatient_metadata]
       }#End
       #pb$tick()
    }#End
    #pb$terminate()
    indexUnknown = which(seurat_object_all_assays$baseline_eGFR==-1)
    patients_with_unknown_eGFR = unique(seurat_object_all_assays$patient[indexUnknown])
}#End - Add exact eGFR from clinical metadata spreadsheet

{#Begin - Get cellType_subtypes_list
  cellType_subtypes_list = list()
  cellTypes = unique(seurat_object_all_assays$Cell_type)
  for (indexCT in 1:length(cellTypes))
  {#Begin
    cellType = cellTypes[indexCT]
    indexCurrentCellType = which(seurat_object_all_assays$Cell_type==cellType)
    cellType_subtypes_list[[cellType]] = unique(seurat_object_all_assays$Cell_subtype[indexCurrentCellType])
  }#End
}#End - Get cellType_subtypes_list

{#Begin - Add simplified disease categories
  diseases = unique(seurat_object_all_assays$Disease)
  disease_simplifiedDisease_list = list()
  disease_simplifiedDisease_list[["NHT"]] = "ignore"
  disease_simplifiedDisease_list[["NHT-ON"]] = "ignore"
  disease_simplifiedDisease_list[["HRT-S"]] = "HRT"
  disease_simplifiedDisease_list[["HRT"]] = "HRT"
  disease_simplifiedDisease_list[["CKD"]] = "CKD"
  disease_simplifiedDisease_list[["RT-UCS"]] = "ignore"
  disease_simplifiedDisease_list[["DKD"]] = "CKD"
  disease_simplifiedDisease_list[["H-CKD"]] = "CKD"
  disease_simplifiedDisease_list[["AKI"]] = "AKI"
  disease_simplifiedDisease_list[["COV-AKI"]] = "AKI"
  disease_simplifiedDisease_list[["LD-HRT"]] = "HRT"
  disease_simplifiedDisease_list[["DM-R"]] = "ignore"
  seurat_object_all_assays$Simplified_disease = "error"
  indexD=1
  for (indexD in 1:length(diseases))
  {#Begin
    disease = diseases[indexD]
    indexCurrentDisease = which(seurat_object_all_assays$Disease==disease)
    stopifnot(disease %in% names(disease_simplifiedDisease_list))
    seurat_object_all_assays$Simplified_disease[indexCurrentDisease] = disease_simplifiedDisease_list[[disease]]
  }#End
  stopifnot(length(which(seurat_object_all_assays$Simplified_disease=="error"))==0)
}#End - Add simplified disease categories

if (split_simplfied_disease_CKD_based_on_eGFR)
{#Begin - Split simplified disease CKD based on eGFT
  indexCKD = which(seurat_object_all_assays$Simplified_disease=="CKD")
  indexHighEgfr = which(seurat_object_all_assays$baseline_eGFR >= 45)
  indexLowEgfr = which(seurat_object_all_assays$baseline_eGFR <= 44)
  indexCKDLowEgfr = intersect(indexCKD, indexLowEgfr)
  indexCKDHighEgfr = intersect(indexCKD, indexHighEgfr)
  seurat_object_all_assays$Simplified_disease[indexCKD] = "ignore"
  seurat_object_all_assays$Simplified_disease[indexCKDLowEgfr] = "CKD_low_eGFR"
  seurat_object_all_assays$Simplified_disease[indexCKDHighEgfr] = "CKD_high_eGFR"
  indexAKI = which(seurat_object_all_assays$Simplified_disease=="AKI")
}#End - Split simplified disease CKD based on eGFT

{#Begin - Print summary of diseases and simlified diseases
    print(paste0("Total participants: ",length(unique(seurat_object_all_assays$patient))))
    print("")
    conditions = unique(seurat_object_all_assays$condition)
    print("Conditions/diseases:")
    for (indexCon in 1:length(conditions))
    {#Begin
       condition = conditions[indexCon]
       indexCurrentCondition = which(seurat_object_all_assays$condition==condition)
       print(paste(condition,": ",length(unique(seurat_object_all_assays$patient[indexCurrentCondition]))," patients",sep=''))
    }#End

    simplified_diseases = unique(seurat_object_all_assays$Simplified_disease)
    indexSD=1
    print("")
    print("Simplified diseases:")
    for (indexSD in 1:length(simplified_diseases))
    {#Begin
       simplified_disease = simplified_diseases[indexSD]
       indexCurrentSimplifiedDisease = which(seurat_object_all_assays$Simplified_disease == simplified_disease)
       print(paste(simplified_disease,": ",length(unique(seurat_object_all_assays$patient[indexCurrentSimplifiedDisease]))," patients",sep=''))
    }#End
}#End - Print summary of diseases and simplified diseases
  
{#Begin - Define file name and directory additions  
    seuratDEGs_directoryAddition = "_seuratDEGs"
    pathways_directoryAddition = "_pathways"
    analysis_directoryAddition = "_analysis"
    randomSeedNodes_subdirectory = "RandomSeedNos/"
    wcfDef_fileName_start = "WcfDefs_"
    wcfDefGenerated_fileName_end = "_generated"
    noCutoff_subdirectory = "NoCutoffs/"
    postHocPower_analysis_subdirectory = "PostHocPowerAnalysis/";
    postHocPower_coefOfVar_fileName_beginning = "PHP_coefOfVar_noOfDrawnDatasets_";
    postHocPower_correlation_fileName_beginning = "PHP_correlation_noOfDrawnDatasets_vs_maxDrawnDatasets_";
    postHocPower_eachClusterGene_fileName_beginning = "PHP_clusterGene_inDependenceOfDatasetNos_";
    documentation_subdirectory = "AA_Documentations/"
}#End - Define file name and directory additions  
  
  
indexAssay=1  
for (indexAssay in 1:length(dataset_assays))
{#Begin - Analyze for given dataset assay
dataset_assay = dataset_assays[indexAssay]
snRNAseq_assay = datasetAssay_snRNAseqAssay_list[[dataset_assay]]

fixedDiseaseGroup_cellCounts_plots = list()


if (snRNAseq_assay!="allSn" )
{ assay_seurat_object = subset(seurat_object_all_assays, subset = assay == snRNAseq_assay) }
if (snRNAseq_assay=="allSn" )
{ assay_seurat_object = seurat_object_all_assays }

{#Begin - Get bgGenes
  rna_data_matrix = GetAssayData(assay_seurat_object, assay="RNA",layer = "data")
  indexNonZero = which(rowSums(rna_data_matrix)!=0)
  bg_Genes = rownames(rna_data_matrix)[indexNonZero]
  bg_Genes = toupper(bg_Genes)
}#End - Get bgGenes

identity_counts = c()
sig_documentations = c()
indexCellTypeGroup = 1

for (indexCellTypeGroup in 1:length(cellTypeGroup_cellTypes_list))    
{#Begin - Analyze for given cellType
  cellTypeGroup_of_interest = names(cellTypeGroup_cellTypes_list)[indexCellTypeGroup]
  report_term = paste(cellTypeGroup_of_interest," - ",dataset_assay," (",snRNAseq_assay,"): ",sep='')
  print(paste(report_term, "Start analysis",sep=''))
  assayCellTypeGroup = paste(dataset_assay,"_",cellTypeGroup_of_interest,sep='')
  assayCellTypeGroupSubtypeL = paste(assayCellTypeGroup,subtype_level,sep='')
  assayCellTypeGroupDisease = paste(assayCellTypeGroup,diseaseFocusLabel,sep='')
  current_cellTypes = cellTypeGroup_cellTypes_list[[cellTypeGroup_of_interest]]

  {#Begin - Generate all results directories and define fileNames
    shared_directory_label = paste(assayCellTypeGroupSubtypeL,sep='')
    results_directory = base_directory
    results_directory = paste(results_directory,level1_subdirectory,sep='')
    if (!dir.exists(results_directory)) { dir.create(results_directory) }
    report_directory = paste(results_directory,documentation_subdirectory,sep='')
    if (!dir.exists(report_directory)) { dir.create(report_directory) }
    sl_results_directory = paste(results_directory,assayCellTypeGroupSubtypeL,"/",sep='')
    if (!dir.exists(sl_results_directory)) { dir.create(sl_results_directory) }

    if (!dir.exists(sl_results_directory)) { dir.create(sl_results_directory) }
    subtype_degs_directory = paste(sl_results_directory,shared_directory_label,seuratDEGs_directoryAddition,"/",sep='')
    if (!dir.exists(subtype_degs_directory)) { dir.create(subtype_degs_directory) }
    subtype_degs_randomSeedNo_directory = paste(subtype_degs_directory,randomSeedNodes_subdirectory,sep='')
    if (!dir.exists(subtype_degs_randomSeedNo_directory)) { dir.create(subtype_degs_randomSeedNo_directory) }
    subtype_pathways_directory = paste(sl_results_directory,shared_directory_label,pathways_directoryAddition,"/",sep='')
    if (!dir.exists(subtype_pathways_directory)) { dir.create(subtype_pathways_directory) }
    subtype_analysis_directory = paste(sl_results_directory,shared_directory_label,analysis_directoryAddition,"/",sep='')
    if (!dir.exists(subtype_analysis_directory)) { dir.create(subtype_analysis_directory) }
    
    disease_results_directory = paste(results_directory,assayCellTypeGroupDisease,"/",sep='')
    
    shared_directory_label = paste(assayCellTypeGroupDisease,sep='')
    if (!dir.exists(disease_results_directory)) { dir.create(disease_results_directory) }
    disease_degs_directory = paste(disease_results_directory,shared_directory_label,seuratDEGs_directoryAddition,"/",sep='')
    if (!dir.exists(disease_degs_directory)) { dir.create(disease_degs_directory) }
    disease_degs_randomSeedNo_directory = paste(disease_degs_directory,randomSeedNodes_subdirectory,sep='')
    if (!dir.exists(disease_degs_randomSeedNo_directory)) { dir.create(disease_degs_randomSeedNo_directory) }
    disease_pathways_directory = paste(disease_results_directory,shared_directory_label,pathways_directoryAddition,"/",sep='')
    if (!dir.exists(disease_pathways_directory)) { dir.create(disease_pathways_directory) }
    disease_analysis_directory = paste(disease_results_directory,shared_directory_label,analysis_directoryAddition,"/",sep='')
    if (!dir.exists(disease_analysis_directory)) { dir.create(disease_analysis_directory) }
  }#End - Generate all results directories and define fileNames

  if (cellTypeGroup_of_interest=="All")
  { cellTypeGroup_seurat_object = assay_seurat_object }
  if (cellTypeGroup_of_interest!="All")
  { cellTypeGroup_seurat_object = subset(assay_seurat_object,subset = Cell_type %in% current_cellTypes) }
  
  stopifnot(length(which(!current_cellTypes %in% names(cellType_subtypes_list)))==0)
  stopifnot(length(unique(cellTypeGroup_seurat_object$Cell_type)) == length(current_cellTypes))
  stopifnot(length(unique(round(colSums(expm1(cellTypeGroup_seurat_object@assays$RNA@layers$data)))))==1)
  
  Generate_sig_documentation_lines_and_add = function(current_cellType_sig_documentations, allMarkers, i_assayCellTypeGroupSubtypeL)
  {#Begin - Generate_sig_documentation_lines
    Col_names = c("AssayCellTypeGroupSubtypeLDisease","Cluster","Sig_degs_count","Sig_cutoff")
    Col_length = length(Col_names)
    Row_names = ""
    Row_length = length(Row_names)
    sig_documentation_base_line = as.data.frame(array(NA,c(Row_length,Col_length),dimnames=list(Row_names,Col_names)))
    clusters = unique(allMarkers$cluster)
    indexCluster=1
    for (indexCluster in 1:length(clusters))
    {#Begin
      cluster = clusters[indexCluster]
      indexCurrentCluster = which(allMarkers$cluster==cluster)
      currentCluster_allMarkers = allMarkers[indexCurrentCluster,]
      sig_genes_count = dim(currentCluster_allMarkers)[[1]]
      sig_documentation_line = sig_documentation_base_line
      sig_documentation_line$AssayCellTypeGroupSubtypeLDisease = i_assayCellTypeGroupSubtypeL
      sig_documentation_line$Cluster = cluster
      sig_documentation_line$Sig_degs_count = sig_genes_count
      if (length(current_cellType_sig_documentations[,1])>0) { current_cellType_sig_documentations = rbind(current_cellType_sig_documentations,sig_documentation_line) }
      else { current_cellType_sig_documentations = sig_documentation_line }
    }#End
    return (current_cellType_sig_documentations)
  }#End - Generate_sig_documentation_lines
  
  current_cellType_sig_documentations = c()
  cellType_identity_lines = c()
  
  Calculate_pathways_using_mbc_pathNet = function(degs_directory, mbc_pathNet_directory, mbc_pathNet_exe_fileName, pathway_directory, ontologies, is_windows, is_linux)
  {#Begin - Calculate pathways using MBC PathNet
    custom_columnNames_fileName = "Custom_1_columnNames_userData.txt"
    complete_custom_columnNames_fileName = paste(degs_directory,custom_columnNames_fileName,sep='')
    input_parameterSettings_fileName = "MBC_pathNet_parameter_settings.txt"
    complete_input_parameterSettings_fileName = paste(mbc_pathNet_directory,input_parameterSettings_fileName,sep='')
    output_parameterSettings_fileName = "MBC_pathNet_parameter_settings.txt"
    complete_output_parameterSettings_fileName = paste(degs_directory,output_parameterSettings_fileName,sep='')

    allFiles = list.files(degs_directory, full.names=FALSE)
    allFiles = allFiles[!allFiles %in% c(output_parameterSettings_fileName, custom_columnNames_fileName)]
    if (length(allFiles)<2)
    {#Begin
        print(paste("Skip analysis of the data in ",degs_directory," with MBC PathNet using ",ontology,", random seed no: ",random_seed_no,": non MBCO files <2",sep=''))
    }#End
    if (length(allFiles)>=2)
    {#Begin - Calculate pathways using MBC PathNet
      current_working_directory = getwd()
      setwd(mbc_pathNet_directory)
      custom_1_columnNames <- data.frame(
        #Do not change, attributes have to match the attributes within MBC_pathNet
        `Attribute in the application` = c(
          "Dataset_name", #Index 1
          "Ncbi_official_gene_symbol", #Index 2
          "Time_point", #Index 3
          "Time_unit", #Index 4
          "Integration_group", #Index 5
          "Value_1st", #Index 6
          "Value_2nd", #Index 7
          "Dataset_color" ), #Index 8
        #Adapt to your data column names, each column name will be matched to the attribute at the same index,
        #indices with empty strings will be ignored
        `Column name in user data` = c(
          "Cluster", #Index 1
          "Gene", #Index 2
          "", #Index 3
          "", #Index 4
          "", #Index 5
          "Avg_log2FC_arithMean", #Index 6
          "P_val_adj_geomMean", #Index 7
          "" )) #Index 8
      complete_custom_columnNames_fileName = paste(degs_directory,"Custom_1_columnNames_userData.txt",sep='')
      colnames(custom_1_columnNames) = gsub("[.]"," ",colnames(custom_1_columnNames))
      write.table(custom_1_columnNames,file=complete_custom_columnNames_fileName,quote=FALSE,row.names=FALSE,col.names = TRUE,sep='\t')
      
      complete_analysis_finished_fileName = paste(degs_directory,"Analysis_finished.txt",sep='')
      indexOntology=1
      for (indexOntology in 1:length(ontologies))
      {#Begin
        ontology = ontologies[indexOntology]
        print(paste("Analyzing the data in '",degs_directory,"' with MBC PathNet using ",ontology," and generating selected charts and networks.",sep=''))
        complete_input_parameterSettings_fileName = paste(mbc_pathNet_directory,"MBC_pathNet_parameter_settings.txt",sep='')
        complete_output_parameterSettings_fileName = paste(degs_directory,"MBC_pathNet_parameter_settings.txt",sep='')
        mbcPathNet_parameter_lines <- readLines(complete_input_parameterSettings_fileName)

        indexUseCustomizedColors = grep("Bardiagram_options_class\tCustomized_colors\t",mbcPathNet_parameter_lines)
        stopifnot(length(indexUseCustomizedColors)==1)
        mbcPathNet_parameter_lines[indexUseCustomizedColors] = "Bardiagram_options_class\tCustomized_colors\tTrue"
        
        indexGenerateFigures = grep("Bardiagram_options_class\tGenerate_bardiagrams\t",mbcPathNet_parameter_lines)
        stopifnot(length(indexGenerateFigures)==1)
        mbcPathNet_parameter_lines[indexGenerateFigures] = "Bardiagram_options_class\tGenerate_bardiagrams\tFalse"
        
        indexGenerateNetworks = grep("MBCO_network_based_integration_options_class\tGenerate_scp_networks\t",mbcPathNet_parameter_lines)
        stopifnot(length(indexGenerateNetworks)==1)
        mbcPathNet_parameter_lines[indexGenerateNetworks] = "MBCO_network_based_integration_options_class\tGenerate_scp_networks\tFalse"
        
        indexKeepTopStandard_mbco = grep("MBCO_enrichment_pipeline_options_class\tKeep_top_standard_SCPs\tMbco\t",mbcPathNet_parameter_lines)
        stopifnot(length(indexKeepTopStandard_mbco)==1)
        mbcPathNet_parameter_lines[indexKeepTopStandard_mbco] = "MBCO_enrichment_pipeline_options_class\tKeep_top_standard_SCPs\tMbco\t-1;0;0;999;0"
        
        indexKeepTopDynamic_mbco = grep("MBCO_enrichment_pipeline_options_class\tKeep_top_dynamic_SCPs\tMbco\t",mbcPathNet_parameter_lines)
        stopifnot(length(indexKeepTopDynamic_mbco)==1)
        mbcPathNet_parameter_lines[indexKeepTopDynamic_mbco] = "MBCO_enrichment_pipeline_options_class\tKeep_top_dynamic_SCPs\tMbco\t-1;-1;0;0;-1"
        
        indexStandardPvalue_mbco = grep("MBCO_enrichment_pipeline_options_class\tMax_pvalue_standard_enrichment\tMbco\t",mbcPathNet_parameter_lines)
        stopifnot(length(indexStandardPvalue_mbco)==1)
        mbcPathNet_parameter_lines[indexStandardPvalue_mbco] = "MBCO_enrichment_pipeline_options_class\tMax_pvalue_standard_enrichment\tMbco\t0.05"

        writeLines(mbcPathNet_parameter_lines,complete_output_parameterSettings_fileName)
        exe_path = file.path(mbc_pathNet_directory, "MBC_PathNet.exe")
        shared_arguments = paste0(' --input-dir ', '"', degs_directory, '"',' --results-dir ',pathway_directory, ' --custom-1-column-names', ' --ontology ',ontology,sep='')
        if (is_windows) { cmd = paste0('"', exe_path, shared_arguments,sep='') }
        if (is_linux) { cmd = paste0('mono ',exe_path, shared_arguments,sep='') }
        system(cmd, wait=TRUE)
        if (file.exists(complete_analysis_finished_fileName))
        { print(paste("Enrichment analysis with MBC PathNet using ",ontology," successful. Enrichment results saved in '",pathway_directory,"'",sep='')) }
        if (!file.exists(complete_analysis_finished_fileName))
        {#Begin
          print(paste("For loop interrupted, since MBC PathNet analysis was not finished successfully.",sep=''))
          break;
        }#End
      }#End
      setwd(current_working_directory)
    }#End - Calculate pathways using MBC PathNet    
  }#End - Calculate pathways using MBC PathNet    

  {#Begin - Define task data frame for DEG calculation
     Task_col_names = c("FindMarkers_logfc_threshold","FindMarkers_adj_pvalue_cutoff","Degs_randomSeedNo_directory","Identity_attribute","CellTypeGroup_of_interest","Subtype_level","RandomSeedNo","Max_cells_per_ident","Delete_existing_degs","Max_cores_count","Min_percent_randomSeedNos_for_postHocPowerAnalysis","Average_with_csharp_script_and_do_postHocPowerAnalysis")
     Task_col_length = length(Task_col_names)
     Task_row_names = random_seed_nos
     Task_row_length = length(Task_row_names)
  }#End - Define task data frame for DEG calculation
       
  Calculate_and_write_markers_for_givenIdentityAttribute = function(core_seurat_object, core_task_line)
  {#Begin - Calculate_and_write_markers_for_givenIdentityAttribute
    i_findMarkers_logfc_threshold = core_task_line$FindMarkers_logfc_threshold
    i_findMarkers_adj_pvalue_cutoff = core_task_line$FindMarkers_adj_pvalue_cutoff
    i_degs_directory = core_task_line$Degs_randomSeedNo_directory
    i_identity_attribute = core_task_line$Identity_attribute
    i_random_seedNo = core_task_line$RandomSeedNo
    i_max_cells_per_ident = core_task_line$Max_cells_per_ident
    i_assayCellTypeGroupSubtypeL = core_task_line$AssayCellTypeGroupSubtypeL

    allMarkers = FindAllMarkers(core_seurat_object,assay="RNA",layer="data",group.by = i_identity_attribute, logfc.threshold=i_findMarkers_logfc_threshold,only.pos = TRUE,parallel=FALSE, max.cells.per.ident = i_max_cells_per_ident, random.seed = i_random_seedNo)
    #indexKeep = which(allMarkers$p_val_adj<1)
    #allMarkers = allMarkers[indexKeep,]
    if (!dir.exists(i_degs_directory)) { dir.create(i_degs_directory) }
    core_seurat_object$Identity_value = core_seurat_object[[i_identity_attribute]]
    i_identity_length = length(unique(core_seurat_object$Identity_value))
    
    i_random_seedNo_string = as.character(i_random_seedNo)

    complete_fileName = paste(i_degs_directory,"AllMarkers_",i_assayCellTypeGroupSubtypeL,"_ranSeedNo",i_random_seedNo_string,".txt",sep='')
    #if (length(allMarkers[,1])==0)
    #{#begin
       #allMarkers = data.frame(p_val = numeric(),avg_log2FC=numeric(),pct.1=numeric(),pct.2=numeric(),p_val_adj=numeric(),cluster=character(),gene=character())
    #}#End
    write.table(allMarkers,file=complete_fileName,row.names=FALSE,col.names = TRUE,quote=FALSE,sep='\t')
  }#End - Calculate_and_write_markers_for_givenIdentityAttribute_and_add_to_documentations

  Generate_considered_pathways = function(pathway_directory, considered_pathway_fileName, pathway_pvalue_cutoff, cellTypeGroup_of_interest)
  {#Begin - Generate_considered_pathways
    mbcPathNet_enrichment_fileName = "Mbco_human_standard_allPredictions.txt"
    complete_enrichment_fileName = paste(pathway_directory,mbcPathNet_enrichment_fileName,sep='')
    pathway_results = read.csv(complete_enrichment_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
    indexKeep = which(pathway_results$SCP.level==3)
    stopifnot(length(indexKeep)>0)
    pathway_results = pathway_results[indexKeep,]
    indexKeep = which(pathway_results$Pvalue<=pathway_pvalue_cutoff)
    pathway_results = pathway_results[indexKeep,]
    pathway_results$Cell_type = cellTypeGroup_of_interest
    
    complete_considered_enrichment_fileName = paste(pathway_directory,considered_pathway_fileName,sep='')
    write.table(pathway_results,file=complete_considered_enrichment_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
  }#End - Generate_considered_pathways
  
  Generate_and_write_physiological_functions_template = function(pathway_directory, considered_pathway_fileName, wcFunct_definition_fileName, reference_pfs_list=list())
  {#Begin
    complete_considered_enrichment_fileName = paste(pathway_directory,considered_pathway_fileName,sep='')
    pathway_results = read.csv(file=complete_considered_enrichment_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
    internal_cellTypeGroup_of_interest = unique(pathway_results$Cell_type)
    stopifnot(length(internal_cellTypeGroup_of_interest)==1)
    
    Col_names = c("Cell_type","Ontology","PF_group","Color","Hexadecimal_color","Comment","Reference","SCP","Symbol","Cell_subtypes (p-value)")
    Col_length = length(Col_names)
    indexPR=1
    preliminary_physiologial_functions = c()
    for (indexPR in 1:length(pathway_results[,1]))
    {#Begin - Generate preliminary physiological functions by reorganizing pathway enrichment results
      pathway_line = pathway_results[indexPR,]
      gene_symbols = strsplit(pathway_line$Overlapping.gene.symbols,",")[[1]]
      Row_names = 1:length(gene_symbols)
      Row_length = length(Row_names)
      physiological_block = as.data.frame(array("",c(Row_length,Col_length),dimnames=list(Row_names,Col_names)))
      physiological_block$Cell_type = internal_cellTypeGroup_of_interest
      physiological_block$Ontology = "Mbco_level3"
      physiological_block$PF_group = ""
      physiological_block$PieChart_slice_order = ""
      physiological_block$Color = ""
      physiological_block$Comment = ""
      physiological_block$Reference = ""
      physiological_block$Information_mapped_from = ""
      physiological_block$Hexadecimal_color = ""
      physiological_block$SCP = pathway_line$SCP
      physiological_block$Symbol = gene_symbols
      physiological_block[["Cell_subtypes (p-value)"]] = paste(pathway_line$Dataset.name," (",sprintf("%.1e", pathway_line$Pvalue),")",sep='')
      if (length(preliminary_physiologial_functions[,1])>0) { preliminary_physiologial_functions = rbind(preliminary_physiologial_functions, physiological_block) }
      else { preliminary_physiologial_functions = physiological_block }
    }#End - Generate preliminary physiological functions by reorganizing pathway enrichment results
    
    preliminary_physiologial_functions$Unique_identifier = paste(preliminary_physiologial_functions$SCP,"-",preliminary_physiologial_functions$Symbol,sep='')
    unique_identifiers = unique(preliminary_physiologial_functions$Unique_identifier)
    physiological_functions = c()
    indexUID=1
    for (indexUID in 1:length(unique_identifiers))
    {#Begin - Combine lines with same SCP gene annotations to physiological functions
      unique_identifier = unique_identifiers[indexUID]
      indexCurrentUID = which(preliminary_physiologial_functions$Unique_identifier==unique_identifier)
      uid_physiological_functions = preliminary_physiologial_functions[indexCurrentUID,]
      if (length(indexCurrentUID)>0)
      {#Begin
        uid_physiological_functions$`Cell_subtypes (p-value)` = paste(uid_physiological_functions$`Cell_subtypes (p-value)`,collapse=", ")
      }#End
      if (length(physiological_functions[,1])>0) { physiological_functions = rbind(physiological_functions,uid_physiological_functions[1,]) }
      else { physiological_functions = uid_physiological_functions[1,] }
    }#End - Combine lines with same SCP gene annotations to physiological functions
    
    if (length(reference_pfs_list)>0)
    {#Begin - Add annotations from references, if existent
      physiological_functions = physiological_functions[order(physiological_functions$Unique_identifier),]
      referenceNames = names(reference_pfs_list)
      indexR=1
      for (indexR in 1:length(referenceNames))
      {#Begin
        referenceName = referenceNames[indexR]
        reference_pfs = reference_pfs_list[[referenceName]]
        reference_pfs$Unique_identifier = paste(reference_pfs$SCP,"-",reference_pfs$Symbol,sep='')
        reference_pfs = reference_pfs[order(reference_pfs$Unique_identifier),]
        indexNotDuplicated = which(!duplicated(reference_pfs$Unique_identifier))
        reference_pfs = reference_pfs[indexNotDuplicated,]
        
        mapped_colnames = c("PF_group","Comment","Reference")
        indexM=3
        for (indexM in 1:length(mapped_colnames))
        {#Begin
          mapped_colname = mapped_colnames[indexM]
          stopifnot(mapped_colname %in% colnames(reference_pfs))
          {#Begin
            indexEmtpy_pf = unique(c(which(physiological_functions[[mapped_colname]]==""),
                                     which(is.na(physiological_functions[[mapped_colname]])),
                                     which(is.null(physiological_functions[[mapped_colname]]))))
            indexNotEmpty_ref = which((reference_pfs[[mapped_colname]]!="")&(!is.na(reference_pfs[[mapped_colname]])))
            indexUID_ref_in_pf = indexNotEmpty_ref[which(reference_pfs$Unique_identifier[indexNotEmpty_ref] %in% physiological_functions$Unique_identifier[indexEmtpy_pf])]
            indexUID_pf_in_ref = indexEmtpy_pf[which(physiological_functions$Unique_identifier[indexEmtpy_pf] %in% reference_pfs$Unique_identifier[indexNotEmpty_ref])]
            if (length(indexUID_ref_in_pf)>0)
            {#Begin
                indexNotEqual = which(physiological_functions$Unique_identifier[indexUID_pf_in_ref] != reference_pfs$Unique_identifier[indexUID_ref_in_pf])
                stopifnot(length(indexNotEqual)==0)
                stopifnot(all.equal(physiological_functions$SCP[indexUID_pf_in_ref],reference_pfs$SCP[indexUID_ref_in_pf]))
                stopifnot(all.equal(physiological_functions$Symbol[indexUID_pf_in_ref],reference_pfs$Symbol[indexUID_ref_in_pf]))
                physiological_functions[[mapped_colname]][indexUID_pf_in_ref] = reference_pfs[[mapped_colname]][indexUID_ref_in_pf]
                add_delimiter = ""
                indexRefExists = grep(referenceName,physiological_functions$Information_mapped_from)
                indexAddRef = indexUID_pf_in_ref[!indexUID_pf_in_ref %in% indexRefExists]
                if (length(indexAddRef)>0)
                {#Begin
                  if (max(nchar(physiological_functions$Information_mapped_from[indexAddRef]))>0) { add_delimiter = ", " }
                  physiological_functions$Information_mapped_from[indexAddRef] = paste(physiological_functions$Information_mapped_from[indexAddRef],add_delimiter,referenceName,sep='')
                }#End
            }#End
          }#End
        }#end
      }#End
    }#End - Add annotations from references, if existent
    
    physiological_functions$Unique_identifier = NULL
    #indexColor = which(physiological_functions$Color!="")
    #physiological_functions$Hexadecimal_color[indexColor] = col2hex(physiological_functions$Color[indexColor])
    
    complete_fileName = paste(pathway_directory,wcFunct_definition_fileName,sep='')
    write.table(physiological_functions,file=complete_fileName,row.names=FALSE,col.names=TRUE,quote=FALSE,sep='\t')
  }#End - Generate_and_write_physiological_functions_template
  
  Generate_matrix_from_dataframe = function(dataframe_input, rowHeadlineName, colHeadlineName, valueHeadlineName, empty_value)
  {#Begin
    library(purrr)
    library(dplyr)
    indexRowHeadlineName = which(colnames(dataframe_input)==rowHeadlineName)
    if (length(indexRowHeadlineName)!=1) { stop(paste("Generate_matrix_from_dataframe: Row headline name '",rowHeadlineName,"' does not exist or is duplicated",sep='')) }
    indexColHeadlineName = which(colnames(dataframe_input)==colHeadlineName)
    if (length(indexColHeadlineName)!=1) { stop(paste("Generate_matrix_from_dataframe: Col headline name '",colHeadlineName,"' does not exist or is duplicated",sep='')) }
    indexValueHeadlineName = which(colnames(dataframe_input)==valueHeadlineName)
    if (length(indexValueHeadlineName)!=1) { stop(paste("Generate_matrix_from_dataframe: Value headline name'",valueHeadlineName,"' does not exist or is duplicated",sep='')) }
    
    dataframe <- dataframe_input
    indexRow = which(colnames(dataframe)==rowHeadlineName)
    indexCol = which(colnames(dataframe)==colHeadlineName)
    indexValue = which(colnames(dataframe)==valueHeadlineName)
    colnames(dataframe)[indexRow] = "Row"
    colnames(dataframe)[indexCol] = "Col"
    colnames(dataframe)[indexValue] = "Value"
    dataframe$Row_col <- paste0(dataframe$Row, "_", dataframe$Col)
    if (anyDuplicated(dataframe$Row_col) > 0) { stop("Duplicated row and col names combination") }
    dataframe$Row_col <- NULL
    
    Col_names = unique(dataframe$Col)
    Col_names = Col_names[order(Col_names)]
    Col_length = length(Col_names)
    Row_names = unique(dataframe$Row)
    Row_names = Row_names[order(Row_names)]
    Row_length = length(Row_names)
    matrix = array(empty_value,c(Row_length,Col_length),dimnames=list(Row_names,Col_names))
    
    col_indices <- match(Col_names, colnames(matrix))
    print("Pre-check for unmatched columns to reduce checks in the loop")
    if (any(is.na(col_indices))) { stop("Some column names in Col_names are missing in the matrix") }
    
    print("Get row indices and values all at once for faster access")
    row_indices <- match(dataframe$Row, rownames(matrix))
    col_indices_df <- match(dataframe$Col, Col_names)
    
    print("Check for any unmatched rows early")
    if (any(is.na(row_indices))) { stop("Some row names in dataframe$Row are missing in the matrix") }
    
    print("Assign values directly by row and column indices")
    matrix[cbind(row_indices, col_indices[col_indices_df])] <- dataframe$Value
    
    print("Final NA check")
    if (any(is.na(matrix))) { stop("Matrix filled with NAs") }  
    
    return (as.matrix(matrix))
  }#End - Generate_matrix_from_dataframe
  
  Cluster_considered_pathways = function(pathway_directory, considered_pathway_fileName, wcFunct_definition_fileName)
  {#Begin - Cluster_considered_pathways
    complete_considered_pathway_fileName = paste(pathway_directory, considered_pathway_fileName,sep='')
    considered_pathways = read.csv(complete_considered_pathway_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
    internal_cellTypeGroup_of_interest = unique(considered_pathways$Cell_type) 
    stopifnot(length(internal_cellTypeGroup_of_interest)==1)
    
    if (exists("physioDefs")) { rm(physioDefs) }
    complete_physioDef_fileName = paste(pathway_directory, wcFunct_definition_fileName,sep='')
    if (file.exists(complete_physioDef_fileName)) { physioDefs = read.csv(complete_physioDef_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE) }
    
    pathway_matrix = Generate_matrix_from_dataframe(dataframe_input = considered_pathways,
                                                    colHeadlineName = "Dataset.name",
                                                    rowHeadlineName = "SCP",
                                                    valueHeadlineName = "Minus.log10_pvalue",
                                                    empty_value = 0)
    
    distance_correlation_method = "pearson"
    dist.function = function(x) as.dist((1-cor((x),method=distance_correlation_method))/2)
    hclust.function = function(x) hclust(x,method="average")
    
    row_dist = dist.function(t(pathway_matrix)) 
    row_hc = hclust.function(row_dist)
    row_dend = as.dendrogram(row_hc)
    
    col_dist = dist.function(pathway_matrix) 
    col_hc = hclust.function(col_dist)
    col_dend = as.dendrogram(col_hc)
    
    indexOriginalData_row = order.dendrogram(row_dend);
    indexOriginalData_col = order.dendrogram(col_dend);
    pathway_matrix_reordered = pathway_matrix[indexOriginalData_row,indexOriginalData_col]
    
    
    matrix_melt = melt(pathway_matrix_reordered)
    indexKeep = which(matrix_melt$value>0)
    matrix_melt = matrix_melt[indexKeep,]
    matrix_melt$Var1_factor = factor(matrix_melt$Var1, levels=unique(matrix_melt$Var1))
    matrix_melt$Var2_factor = factor(matrix_melt$Var2, levels=unique(matrix_melt$Var2))
    
    matrix_melt = matrix_melt[order(matrix_melt$Var1_factor),]
    pathways = unique(matrix_melt$Var1)
    
    pathways_with_multiple_pf_groups = c()
    pathway_colors = replicate(length(pathways),"gray40")
    
    if (exists("physioDefs"))
    {#Begin - Add physiological group colors
      names(pathway_colors) = pathways
      indexPathway = 1
      for (indexPathway in 1:length(pathways))
      {#Begin - Set pathway colors using pf group colors
        pathway = pathways[indexPathway]
        indexCurrentPathway = which(physioDefs$SCP==pathway)
        if (length(indexCurrentPathway)==0) { stop("Pathway is missing in physioDefinitions") }
        pf_group = unique(physioDefs$PF_group[indexCurrentPathway])
        if (length(pf_group)!=1)
        {#Begin
          pathways_with_multiple_pf_groups = c(pathways_with_multiple_pf_groups,as.character(pathway))
          pf_group=""
        }#End
        if ( (!is.na(pf_group))
             &(pf_group!=""))
        {#Begin
          color = unique(physioDefs$Color[indexCurrentPathway])
          stopifnot(length(color)==1)
          if (color!="")
          {#Begin
            pathway_colors[indexPathway] = color
          }#End
        }#End
      }#End - Set pathway colors using pf group colors
      if (length(pathways_with_multiple_pf_groups)>0)
      {#Begin
        stop(paste(paste(pathways_with_multiple_pf_groups,collapse=", ")," is/are annotated to multiple pf groups.",sep=''))
      }#End
    }#End - Add physiological group colors
    
    matrix_melt$Label = as.character(round(matrix_melt$value*10)/10)
    indexLarger10 = which(matrix_melt$value>=10)
    matrix_melt$Label[indexLarger10] = as.character(round(matrix_melt$value))
    Plots = list()
    
    Heatmap = ggplot(matrix_melt,aes(y=Var1_factor,x=Var2_factor,fill=Var1,label=Label))
    Heatmap = Heatmap + geom_tile()
    Heatmap = Heatmap + geom_text(size = 2, color="black")
    #Heatmap = Heatmap + scale_fill_gradient(low="white",high="orange")
    Heatmap = Heatmap + scale_fill_manual(values=pathway_colors)
    #Heatmap = Heatmap + ggnewscale::new_scale_fill()
    Heatmap = Heatmap + theme(legend.position = "none")
    #Heatmap = Heatmap + geom_tile(data = filter(matrix_melt, value>max_other), mapping = aes(fill = value > max_other))
    #Heatmap = Heatmap + scale_fill_manual("", values = "darkred", labels = paste("> ",max_other,sep=''), guide = guide_legend(order = 1))
    Heatmap = Heatmap + scale_x_discrete(position="top")
    Heatmap = Heatmap + theme(axis.text.x=element_text(angle = 90,vjust=1.5,hjust=0,face=2,size=12))
    Heatmap = Heatmap + theme(axis.text.y=element_text(color=pathway_colors,size=7))
    Heatmap = Heatmap + xlab("") + ylab("")
    Plots[[length(Plots)+1]] = Heatmap
    
    complete_pdf_fileName = paste(pathway_directory,"Pathway_heatmap_",internal_cellTypeGroup_of_interest,subtype_level,".pdf",sep='')
    Generate_plots(Plots,complete_pdf_fileName,rows_count=1,cols_count=1)
  }#End - Cluster_considered_pathways

  Generate_and_return_plots_for_averaged_degs = function(averaged_degs)
  {#Begin - Generate_plots_for_averaged_degs
      Plots = list()
      Plot = ggplot(averaged_degs,aes(x=cluster,y=Minus_log10_adj_pvalue_sd))
      Plot = Plot + geom_violin()
      Plot = Plot + ylab(paste("SD -log10(p)",sep=''))
      Plots[[length(Plots)+1]] = Plot

      quantile_cutoff = 0.75
      quantile_value_cutoff = quantile(averaged_degs$Minus_log10_adj_pvalue_sd, probs = quantile_cutoff, na.rm=TRUE)

      Plot = ggplot(averaged_degs,aes(x=cluster,y=Minus_log10_adj_pvalue_sd))
      Plot = Plot + geom_violin()
      Plot = Plot + ylim(0,quantile_value_cutoff)
      Plot = Plot + ylab(paste("SD -log10(p)\n",quantile_cutoff," quantile",sep=''))
      Plots[[length(Plots)+1]] = Plot

      averaged_degs$Coef_variation = averaged_degs$Minus_log10_adj_pvalue_sd / averaged_degs$Minus_log10_adj_pvalue
      Plot = ggplot(averaged_degs,aes(x=cluster,y=Coef_variation))
      Plot = Plot + geom_violin()
      Plot = Plot + ylab(paste("Coef var -log10(p)",sep=''))
      Plots[[length(Plots)+1]] = Plot

      quantile_cutoff = 0.75
      quantile_value_cutoff = quantile(averaged_degs$Coef_variation, probs = quantile_cutoff, na.rm=TRUE)
      
      Plot = ggplot(averaged_degs,aes(x=cluster,y=Coef_variation))
      Plot = Plot + geom_violin()
      Plot = Plot + ylim(0,quantile_value_cutoff)
      Plot = Plot + ylab(paste("Coef var -log10(p)\n",quantile_cutoff," quantile",sep=''))
      Plots[[length(Plots)+1]] = Plot

      return (Plots)      
  }#End - Generate_plots_for_averaged_degs
  
  Calculate_degs_using_parallel_cores_if_necessary_and_eventually_add_plots_to_avg_degs_plots = function(core_seurat_object, bg_genes, deg_tasks, current_cellType_sig_documentations, documentation_directory)
  {#Begin - Calculate_degs_using_parallel_cores
    if (exists("core_tasks")) { rm(core_tasks) }
    delete_existing_degs = unique(deg_tasks$Delete_existing_degs)
    max_cores_count = unique(deg_tasks$Max_cores_count)
    min_percent_randomSeedNos = unique(deg_tasks$Min_percent_randomSeedNos)
    max_degListsCombinations = unique(deg_tasks$Max_degListsCombinations)
    average_with_csharp_script_and_do_postHocPowerAnalysis = unique(deg_tasks$Average_with_csharp_script_and_do_postHocPowerAnalysis)
    i_assayCellTypeGroupSubtypeL = unique(deg_tasks$AssayCellTypeGroupSubtypeL)
    degs_directory = unique(deg_tasks$Degs_directory)
    degs_noCutoff_directory = unique(deg_tasks$Degs_noCutoff_directory)
    findMarkers_adj_pvalue_cutoff = unique(deg_tasks$FindMarkers_adj_pvalue_cutoff)
    findMarkers_topDEGs = unique(deg_tasks$FindMarkers_topDEGs)
    findMarkers_logfc_threshold = unique(deg_tasks$FindMarkers_logfc_threshold)
    identity_attribute = unique(deg_tasks$Identity_attribute)
    do_postHocPower_analysis = unique(deg_tasks$Do_postHocPower_analysis)
    finalAnalysis_numberOfDatasets = unique(deg_tasks$FinalAnalysis_numberOfDatasets)
    
    
    postHocPower_coefOfVar_fileName_beginning = unique(deg_tasks$PostHocPower_coefOfVar_fileName_beginning)
    postHocPower_correlation_fileName_beginning = unique(deg_tasks$PostHocPower_correlation_fileName_beginning)
    postHocPower_eachClusterGene_fileName_beginning = unique(deg_tasks$PostHocPower_eachClusterGene_fileName_beginning)
    postHocPower_analysis_subdirectory = unique(deg_tasks$PostHocPower_analysis_subdirectory)
    postHocPower_max_number_of_datasets = unique(deg_tasks$PostHocPower_max_number_of_datasets)
    
    stopifnot(length(delete_existing_degs)==1)
    stopifnot(length(max_cores_count)==1)
    stopifnot(length(max_degListsCombinations)==1)
    stopifnot(length(average_with_csharp_script_and_do_postHocPowerAnalysis)==1)
    stopifnot(length(i_assayCellTypeGroupSubtypeL)==1)
    stopifnot(length(degs_directory)==1)
    stopifnot(length(degs_noCutoff_directory)==1)
    stopifnot(length(findMarkers_adj_pvalue_cutoff)==1)
    stopifnot(length(findMarkers_topDEGs)==1)
    stopifnot(length(findMarkers_logfc_threshold)==1)
    stopifnot(length(postHocPower_coefOfVar_fileName_beginning)==1)
    stopifnot(length(postHocPower_correlation_fileName_beginning)==1)
    stopifnot(length(postHocPower_eachClusterGene_fileName_beginning)==1)
    stopifnot(length(postHocPower_analysis_subdirectory)==1)
    stopifnot(length(postHocPower_max_number_of_datasets)==1)
    stopifnot(length(identity_attribute)==1)

    identity_length = length(unique(core_seurat_object[[identity_attribute]][,1]))

    calculate_degs_and_eventually_delete_old_DEGs = TRUE
    if (calculate_degs_and_eventually_delete_old_DEGs)
    {#Begin - calculate_degs_and_eventually_delete_old_DEGs
        cores_count = max_cores_count
        length_tasks = length(deg_tasks[,1])
        if (length_tasks<cores_count) { cores_count = length_tasks }
        tasks_per_core = length_tasks/cores_count
        
        parallel_clusters = makeCluster(cores_count)
        on.exit(stopCluster(parallel_clusters))    
        clusterEvalQ(parallel_clusters, {
          library("Seurat")
          library("Matrix")
          library("dplyr")
          library("BPCells")
        }
        );
    
        if (delete_existing_degs)
        {#Begin - Delete existing DEGs
          print(paste("Delete existing DEGs and random DEGs",sep=''))
          complete_delete_fileNames = list.files(degs_directory,recursive = TRUE, full.names = TRUE)
          for (complete_delete_fileName in complete_delete_fileNames)
          { unlink(paste(complete_delete_fileName,sep='')) }
        }#End - Delete existing DEGs
        
            
        print(paste("Add tasks, seurat object and marker calculation function to each core",sep=''))
        pb = progress_bar$new(total=cores_count)
        for (indexCore in 1:cores_count)
        {#Begin - Fill cores
          startIndex = min(floor((indexCore-1) * tasks_per_core+1),length_tasks);
          endIndex = min(floor(indexCore * tasks_per_core),length_tasks)
          if (indexCore == cores_count) { endIndex = length_tasks }
          core_tasks = deg_tasks[startIndex:endIndex,]
          stopifnot(dim(core_tasks)[1]>0)
          clusterExport(parallel_clusters[indexCore],
                        varlist = c("Calculate_and_write_markers_for_givenIdentityAttribute", 
                                    "core_tasks", "core_seurat_object"),
                        envir = environment())    
          pb$tick()
          rm(core_tasks)
        }#End - Fill cores
        rm(core_seurat_object)
        pb$terminate()
        
        combined_core_tasks <- do.call('rbind', clusterEvalQ(parallel_clusters, core_tasks))
        stopifnot(length(combined_core_tasks[,1]) == length_tasks)
        stopifnot(length(which(!random_seed_nos %in% combined_core_tasks$RandomSeedNo))==0)
    
        #core_tasks = deg_tasks; indexCoreTask=1
        
        print(paste("Calculate marker genes for ",identity_length," clusters using ",cores_count," cores and ",length(unique(deg_tasks$RandomSeedNo))," random seed numbers.",sep=''))
        clusterEvalQ(parallel_clusters,
                     {#Begin - Parallel clusters - Calculate DEGs
                       
                       for (indexCoreTask in 1:length(core_tasks[,1]))
                       {#Begin
                         core_task_line = core_tasks[indexCoreTask,]
                         Calculate_and_write_markers_for_givenIdentityAttribute(core_seurat_object, core_task_line)
                       }#End
                     }#End - Parallel clusters - Calculate DEGs
        )
        print("Done")
        
        {#Begin - Close parallel clusters
        #  invisible(gc())
        #  parallel::stopCluster(parallel_clusters)
        #  invisible(gc())
        #  if (exists("parallel_clusters")) { rm(parallel_clusters) }
        }#End - Close parallel clusters
    }#End - calculate_degs_and_eventually_delete_old_DEGs
 
     average_with_csharp_script_and_do_postHocPowerAnalysis=TRUE
     if (average_with_csharp_script_and_do_postHocPowerAnalysis)
     {#Begin - average_with_csharp_script_and_do_postHocPowerAnalysis
        files_length = length(list.files(degs_randomSeedNo_directory))
        print(paste("Use C# script to average DEGs and do PostHocPower analysis using ",max_degListsCombinations," differnt combinations of deg lists by considering at least ",min_percent_randomSeedNos,"% of the random seed nos.",sep=''))
        print(paste("See 'Progress_report.txt' in '",report_directory,"'."))
        exe_path = file.path(complete_calculateAvgPostHocPowerAnalysis_script_fileName)
        stopifnot(file.exists(exe_path))
        shared_arguments = c("--overall_directory",complete_level1_directory,"--max_numberOfDatasets_postHocPower",postHocPower_max_number_of_datasets,"--numberOfDatasets_finalAnalysis",finalAnalysis_numberOfDatasets,"--assay_cellTypeGroup_subtypeLevel",i_assayCellTypeGroupSubtypeL,"--max_degListsCombinations",max_degListsCombinations,"--min_percentageOfDrawnDegLists",min_percent_randomSeedNos,"--randomSeedNodeNoPostHocPower",randomSeedNodeNoPostHocPower,"--maxNoDegListsInMem",maxNoDegListsInMem)
        system2(exe_path, args=shared_arguments)
        
        php_directory = paste(degs_directory,postHocPower_analysis_subdirectory,sep='')

        if ((postHocPower_max_number_of_datasets>0)&(postHocPower_max_degListsCombinations>0))
        {#Begin - Visualize results of postHocPower analysis
             Avg_degs_plots = list()
             
             {#Begin - Generate example gene/cluster vs mean -log10(p) plots
                complete_php_eachClusterGene_fileName = paste(php_directory,postHocPower_eachClusterGene_fileName_beginning,i_assayCellTypeGroupSubtypeL,".txt",sep='')
                eachClusterGene = read.csv(complete_php_eachClusterGene_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
                indexKeep = which(eachClusterGene$Cluster==eachClusterGene$Cluster[1])
                eachClusterGene = eachClusterGene[indexKeep,]
                
                eachClusterGene = eachClusterGene[order(eachClusterGene$Minus_log10_adj_pvalue_arithMean,decreasing=TRUE),]
                eachClusterGene$Unique_identifier = paste(eachClusterGene$Cluster,": gene ",eachClusterGene$Gene,sep='')
                eachClusterGene = eachClusterGene[order(eachClusterGene$Minus_log10_adj_pvalue_coefVar),]
                max_minusLog10Pvalue = max(eachClusterGene$Minus_log10_adj_pvalue_arithMean)
                indexMax = which(eachClusterGene$Minus_log10_adj_pvalue_arithMean==max_minusLog10Pvalue)
                uids_with_max_minusLog10Pvalue = unique(eachClusterGene$Unique_identifier[indexMax])
                uids_no_max_minusLog10Pvalue = unique(eachClusterGene$Unique_identifier)
                uids_no_max_minusLog10Pvalue = uids_no_max_minusLog10Pvalue[!uids_no_max_minusLog10Pvalue %in% uids_with_max_minusLog10Pvalue]
                length_uids_no_max_minusLog10Pvalue = length(uids_no_max_minusLog10Pvalue)
                indexSelectedUIDs = round(c(0.1,0.33,0.66,0.9)*length_uids_no_max_minusLog10Pvalue)
                selected_uids = uids_no_max_minusLog10Pvalue[indexSelectedUIDs]
                
                valueTypeOfInterest_array = c("Minus_log10_adj_pvalue","")
                
                indexUID=2
                for (indexUID in 1:length(selected_uids))
                {#Begin
                   selected_uid = selected_uids[indexUID]
                   indexSelected_uid = which(eachClusterGene$Unique_identifier==selected_uid)
                   selected_uid_eachClusterGene = eachClusterGene[indexSelected_uid,]
                
                   selected_uid_eachClusterGene$VoI_mean = selected_uid_eachClusterGene$Minus_log10_adj_pvalue_arithMean
                   selected_uid_eachClusterGene$VoI_sd = selected_uid_eachClusterGene$Minus_log10_adj_pvalue_sampleSd
                
                   ylabel = paste("Mean -log10(p)\n over <= ",unique(selected_uid_eachClusterGene$Number_of_iterations_for_same_number_of_datasets)," sets of cluster DEGs",sep='')
                
                   Plot = ggplot(selected_uid_eachClusterGene,aes(x=Number_of_datasets,y=VoI_mean))
                   Plot = Plot + geom_point() + geom_errorbar(aes(ymin=VoI_mean-VoI_sd,ymax=VoI_mean+VoI_sd))
                   Plot = Plot + theme_classic()
                   Plot = Plot + coord_cartesian(ylim=c(0,max(selected_uid_eachClusterGene$VoI_mean*1.2)))
                   Plot = Plot + ylab(ylabel) + ggtitle(gsub("_"," ",selected_uid))
                   Plot = Plot + xlab("# DEG lists within each set")
                   Plot = Plot + theme(plot.title = element_text(hjust=0.5,size=11,face=2))
                   Plot = Plot + theme(axis.title = element_text(hjust=0.5,size=9,face=1))
                   Plot = Plot + theme(axis.text = element_text(hjust=0.5,size=9,face=1))
                   Avg_degs_plots[[length(Avg_degs_plots)+1]] = Plot
                }#End
            }#End - Generate example gene/cluster vs mean -log10(p) plots
      
            {#Begin - Plot coefficient of variations
              base_attributes_of_interest_array = c( "Minus_log10_adj_pvalue"
                                                     ,"IsSig")
              complete_php_coefOfVar_fileName = paste(php_directory,postHocPower_coefOfVar_fileName_beginning,i_assayCellTypeGroupSubtypeL,".txt",sep='')
              coefOfVar = read.csv(complete_php_coefOfVar_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
              clusters = unique(coefOfVar$Cluster)
              indexC=3
              for (indexC in 1:length(clusters))
              {#Begin
                  cluster = clusters[indexC]
                  indexCurrentCluster = which(coefOfVar$Cluster==cluster)
                  cluster_coefOfVar = coefOfVar[indexCurrentCluster,]
                  indexBA=2
                  for (indexBA in 1:length(base_attributes_of_interest_array))
                  {#Begin
                      base_attributes_of_interest = base_attributes_of_interest_array[indexBA]
                      cluster_coefOfVar$VoI_mean = cluster_coefOfVar[[paste(base_attributes_of_interest,"_arithMean",sep='')]]
                      cluster_coefOfVar$VoI_sd = cluster_coefOfVar[[paste(base_attributes_of_interest,"_sampleSd",sep='')]]
                
                      ylabel = paste("Mean coefficient of variation\n over <= ",unique(cluster_coefOfVar$Number_of_genes)," genes",sep='')
                
                      max_y = max(cluster_coefOfVar$VoI_mean + cluster_coefOfVar$VoI_sd) * 1.2
                      
                      Plot = ggplot(cluster_coefOfVar,aes(x=Number_of_datasets,y=VoI_mean))
                      Plot = Plot + geom_point() + geom_errorbar(aes(ymin=VoI_mean-VoI_sd,ymax=VoI_mean+VoI_sd))
                      Plot = Plot + theme_classic()
                      Plot = Plot + coord_cartesian(ylim=c(0,max_y))
                      Plot = Plot + ylab(ylabel) + ggtitle(paste(cluster,": ",gsub("_"," ",base_attributes_of_interest),sep=''))
                      Plot = Plot + xlab("# DEG lists within each set")
                      Plot = Plot + theme(plot.title = element_text(hjust=0.5,size=11,face=2))
                      Plot = Plot + theme(axis.title = element_text(hjust=0.5,size=9,face=1))
                      Plot = Plot + theme(axis.text = element_text(hjust=0.5,size=9,face=1))
                      Avg_degs_plots[[length(Avg_degs_plots)+1]] = Plot
                  }#End
              }#End
            }#End - Plot coefficient of variations 
            
            {#Begin - Plot correlations
              base_attributes_of_interest_array = c( "Minus_log10_adj_pvalue"
                                                     ,"Is_sig")
              complete_php_correlation_fileName = paste(php_directory,postHocPower_correlation_fileName_beginning,i_assayCellTypeGroupSubtypeL,".txt",sep='')
              php_correlation = read.csv(complete_php_correlation_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
              indexKeep = which(php_correlation$ValueType_of_interest %in% base_attributes_of_interest_array)
              php_correlation = php_correlation[indexKeep,]
              php_correlation$Unique_identifier = paste(php_correlation$Cluster,": ",php_correlation$ValueType_of_interest,sep='')
              unique_identifiers = unique(php_correlation$Unique_identifier)
              unique_identifiers = unique_identifiers[order(unique_identifiers)]
              indexUID=16
              for (indexUID in 1:length(unique_identifiers))
              {#Begin
                  unique_identifier = unique_identifiers[indexUID]
                  indexCurrentUID = which(php_correlation$Unique_identifier==unique_identifier)
                  uid_correlation = php_correlation[indexCurrentUID,]
                  cluster = unique(uid_correlation$Cluster)
                  ylable = unique(uid_correlation$Correlation_valueType)
                  uid_correlation$VoI_mean = uid_correlation$Correlation_value_mean
                  uid_correlation$VoI_sd = uid_correlation$Correlation_value_sampleSd
                  
                  max_y = 1.2
      
                  Plot = ggplot(uid_correlation,aes(x=Number_of_drawn_datasets,y=VoI_mean))
                  Plot = Plot + geom_point() + geom_errorbar(aes(ymin=VoI_mean-VoI_sd,ymax=VoI_mean+VoI_sd))
                  Plot = Plot + theme_classic()
                  Plot = Plot + coord_cartesian(ylim=c(0,max_y))
                  Plot = Plot + ylab("Pearson correlation") + ggtitle(gsub("_"," ",unique_identifier))
                  Plot = Plot + xlab("# DEG lists within each set")
                  Plot = Plot + theme(plot.title = element_text(hjust=0.5,size=11,face=2))
                  Plot = Plot + theme(axis.title = element_text(hjust=0.5,size=9,face=1))
                  Plot = Plot + theme(axis.text = element_text(hjust=0.5,size=9,face=1))
                  Avg_degs_plots[[length(Avg_degs_plots)+1]] = Plot
              }#End
              
            }#End - Plot correlations
  
             complete_pdf_fileName = paste(php_directory,paste("PostHocPower_",i_assayCellTypeGroupSubtypeL,".pdf",sep=''))
             Generate_plots(Avg_degs_plots,complete_pdf_fileName,4,2)
             complete_pdf_fileName = paste(documentation_directory,paste("PostHocPower_",i_assayCellTypeGroupSubtypeL,".pdf",sep=''))
             Generate_plots(Avg_degs_plots,complete_pdf_fileName,4,2)
        }#End - Visualize results of postHocPower analysis
        
        {#Begin - Write noCutoff_degs and background genes into no cutoff directory
            noCutoff_data_fileName = list.files(degs_noCutoff_directory)
            indexBgGenes = grep("bgGenes",noCutoff_data_fileName)
            indexNoBgGenes = 1:length(noCutoff_data_fileName)
            indexNoBgGenes = indexNoBgGenes[!indexNoBgGenes %in% indexBgGenes]
            noCutoff_data_fileName = noCutoff_data_fileName[indexNoBgGenes]
            stopifnot(length(noCutoff_data_fileName)==1)
            complete_noCutoff_data_fileName = paste(degs_noCutoff_directory,noCutoff_data_fileName,sep='')
            noCutoff_data = read.csv(file=complete_noCutoff_data_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
            complete_bg_genes_fileName = gsub(".txt","_bgGenes.txt",complete_noCutoff_data_fileName)
            write.table(bg_Genes, file=complete_bg_genes_fileName,row.names = FALSE,col.names = FALSE,sep='\t',quote=FALSE)
        }#End - Write noCutoff_degs and background genes into no cutoff directory

        print(paste("Calculate and write significant DEGs for pathway enrichment analysis (adj.p-value <= ",findMarkers_adj_pvalue_cutoff," & significance rank <= ",findMarkers_topDEGs,")",sep=''))
        
        sig_avg_degs = noCutoff_data
        stopifnot(length(unique(sig_avg_degs$P_val_adj_geomMean))>1)
        stopifnot(length(unique(sig_avg_degs$Rank))>1)
        indexSig1 = which(sig_avg_degs$P_val_adj_geomMean<=findMarkers_adj_pvalue_cutoff)
        indexSig2 = which(sig_avg_degs$Rank<=findMarkers_topDEGs)
        indexKeep = intersect(indexSig1,indexSig2)
        sig_avg_degs = sig_avg_degs[indexKeep,]

        current_cellType_sig_documentations = Generate_sig_documentation_lines_and_add(current_cellType_sig_documentations, sig_avg_degs, i_assayCellTypeGroupSubtypeL)        

        subtype_marker_fileName = paste("SigAvgAllMarkers_",cellTypeGroup_of_interest,subtype_level,"_",identity_length,"clusters.txt",sep='')
        complete_marker_sig_fileName = paste(degs_directory,subtype_marker_fileName,sep='')
        write.table(sig_avg_degs,file=complete_marker_sig_fileName,row.names = FALSE,col.names = TRUE,sep='\t',quote=FALSE)
        complete_bg_genes_fileName = gsub(".txt","_bgGenes.txt",complete_marker_sig_fileName)
        write.table(bg_Genes, file=complete_bg_genes_fileName,row.names = FALSE,col.names = FALSE,sep='\t',quote=FALSE)
     }#End - average_with_csharp_script_and_do_postHocPowerAnalysis
 
     return (current_cellType_sig_documentations) 
  }#End - Calculate_degs_using_parallel_cores
  
  Avg_degs_plots = list()
  
  Get_identity_cellCounts = function(core_seurat_object,identity_attribute,subtype_level,dataset_assay,identity_cellCounts_base_line)
  {#Begin - Get_identity_cellCounts
     table_identity = table(core_seurat_object[[identity_attribute]])
     identities = names(table_identity)
     current_cellCounts = c()
     indexI = 1
     for (indexI in 1:length(identities))
     {#Begin
        identity = identities[indexI]
        identity_cellCounts_line = identity_cellCounts_base_line
        identity_cellCounts_line$Identity = names(table_identity)[indexI]
        identity_cellCounts_line$Cell_count = table_identity[indexI]
        identity_cellCounts_line$Subtype_level = subtype_level
        identity_cellCounts_line$Dataset_assay = dataset_assay
        if (length(current_cellCounts[,1])>0) { current_cellCounts = rbind(current_cellCounts,identity_cellCounts_line) }
        else { current_cellCounts = identity_cellCounts_line }
     }#End
     return (current_cellCounts)
  }#End - Get_identity_cellCounts
  
  if (calculate_all_markers_for_subtypes)
  {#Begin - Calculate and write all markers for each subtype after deleting all files in degs directory
    identity_attribute = "Cell_subtype"
    degs_directory = subtype_degs_directory
    degs_randomSeedNo_directory = subtype_degs_randomSeedNo_directory
    assayCellTypeGroupReference = assayCellTypeGroupSubtypeL
    core_seurat_object = cellTypeGroup_seurat_object

    add_cellCounts = Get_identity_cellCounts(core_seurat_object, identity_attribute, subtype_level, dataset_assay, identity_cellCounts_base_line)
    if (length(cellType_identity_lines[,1])>0) { cellType_identity_lines = rbind(cellType_identity_lines,add_cellCounts)
    } else { cellType_identity_lines = cellType_identity_lines }
    
    deg_tasks = as.data.frame(array(NA,c(Task_row_length,Task_col_length),dimnames = list(Task_row_names,Task_col_names)))
    deg_tasks$FindMarkers_logfc_threshold = findMarkers_logfc_threshold
    deg_tasks$FindMarkers_adj_pvalue_cutoff = findMarkers_adj_pvalue_cutoff
    deg_tasks$Degs_directory = degs_directory
    deg_tasks$Degs_randomSeedNo_directory = degs_randomSeedNo_directory
    deg_tasks$Degs_noCutoff_directory = paste(degs_directory,noCutoff_subdirectory,sep='')
    deg_tasks$CellTypeGroup_of_interest = cellTypeGroup_of_interest
    deg_tasks$Identity_attribute = identity_attribute
    deg_tasks$Subtype_level = subtype_level
    deg_tasks$AssayCellTypeGroupSubtypeL = assayCellTypeGroupReference
    deg_tasks$Max_degListsCombinations = postHocPower_max_degListsCombinations
    deg_tasks$FindMarkers_topDEGs = findMarkers_topDEGs
    deg_tasks$RandomSeedNo = random_seed_nos
    deg_tasks$Max_cells_per_ident = findMarkers_max_cells_per_ident
    deg_tasks$Delete_existing_degs = delete_existing_degs
    deg_tasks$Max_cores_count = max_cores_count
    deg_tasks$FinalAnalysis_numberOfDatasets = finalAnalysis_numberOfDatasets
    deg_tasks$Min_percent_randomSeedNos_for_postHocPowerAnalysis = postHocPower_minPercent_randomSeedNos
    deg_tasks$Average_with_csharp_script_and_do_postHocPowerAnalysis = TRUE
    deg_tasks$PostHocPower_coefOfVar_fileName_beginning = postHocPower_coefOfVar_fileName_beginning
    deg_tasks$PostHocPower_correlation_fileName_beginning = postHocPower_correlation_fileName_beginning
    deg_tasks$PostHocPower_eachClusterGene_fileName_beginning = postHocPower_eachClusterGene_fileName_beginning
    deg_tasks$PostHocPower_analysis_subdirectory = postHocPower_analysis_subdirectory
    deg_tasks$PostHocPower_max_number_of_datasets = postHocPower_max_number_of_datasets

    current_cellType_sig_documentations  = Calculate_degs_using_parallel_cores_if_necessary_and_eventually_add_plots_to_avg_degs_plots(core_seurat_object,bg_Genes,deg_tasks,current_cellType_sig_documentations,documentation_directory)
    
    rm(core_seurat_object)
    rm(deg_tasks)
  }#End - Calculate and write all markers for each subtype after deleting all files in degs directory

  if (calculate_pathways_for_subtypes)
  {#Begin - calculate_pathways_for_subtypes
      degs_directory = subtype_degs_directory
      pathway_directory = subtype_pathways_directory
      Calculate_pathways_using_mbc_pathNet(degs_directory, mbc_pathNet_directory, mbc_pathNet_exe_fileName, pathway_directory, ontologies, is_windows, is_linux)
      Generate_considered_pathways(pathway_directory, considered_pathway_fileName, pathway_pvalue_cutoff, cellTypeGroup_of_interest)

      wcFunct_definition_fileName = paste(wcfDef_fileName_start,assayCellTypeGroupSubtypeL,wcfDefGenerated_fileName_end,".txt",sep='')

      reference_pfs_list = list()
      reference_name = paste(dataset_assay,"_",cellTypeGroup_of_interest,subtype_level,sep='')
      reference_fileNames = assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[reference_name]]
      if (length(reference_fileNames)>0)
      {#Begin
          indexRef = 1
          for (indexRef in 1:length(reference_fileNames))
          {#Begin
            reference_fileName = reference_fileNames[indexRef]
            complete_reference_fileName = paste(reference_wholeCellFunctionDef_directory,reference_fileName,sep='')
            refs = read_excel(complete_reference_fileName)
            reference_pfs_list[[reference_fileName]] = refs
          }#End
      }#End
      Generate_and_write_physiological_functions_template(pathway_directory, considered_pathway_fileName, wcFunct_definition_fileName, reference_pfs_list)
  }#End - calculate_pathways_for_subtypes

  if (calculate_all_markers_for_simplified_diseases)
  {#Begin - Calculate and write all markers for each simplified disease after deleting all files in degs directory
    identity_attribute = "Simplified_disease"
    degs_directory = disease_degs_directory
    degs_randomSeedNo_directory = disease_degs_randomSeedNo_directory
    assayCellTypeGroupReference = assayCellTypeGroupDisease
    core_seurat_object = subset(cellTypeGroup_seurat_object, subset = Simplified_disease %in% simplified_diseases_of_interest)
    core_seurat_object = subset(core_seurat_object, subset = source == "KPMP")
    core_seurat_object = subset(core_seurat_object, subset = tissue_type == "Biopsy")
    
    add_cellCounts = Get_identity_cellCounts(core_seurat_object, identity_attribute, subtype_level, dataset_assay, identity_cellCounts_base_line)
    if (length(cellType_identity_lines[,1])>0) { cellType_identity_lines = rbind(cellType_identity_lines,add_cellCounts) }
    else { cellType_identity_lines = cellType_identity_lines }

    deg_tasks = as.data.frame(array(NA,c(Task_row_length,Task_col_length),dimnames = list(Task_row_names,Task_col_names)))
    deg_tasks$FindMarkers_logfc_threshold = findMarkers_logfc_threshold
    deg_tasks$FindMarkers_adj_pvalue_cutoff = findMarkers_adj_pvalue_cutoff
    deg_tasks$Degs_directory = degs_directory
    deg_tasks$Degs_randomSeedNo_directory = degs_randomSeedNo_directory
    deg_tasks$Degs_noCutoff_directory = paste(degs_directory,noCutoff_subdirectory,sep='')
    deg_tasks$CellTypeGroup_of_interest = cellTypeGroup_of_interest
    deg_tasks$Identity_attribute = identity_attribute
    deg_tasks$Subtype_level = subtype_level
    deg_tasks$FinalAnalysis_numberOfDatasets = finalAnalysis_numberOfDatasets
    deg_tasks$AssayCellTypeGroupSubtypeL = assayCellTypeGroupReference
    deg_tasks$Max_degListsCombinations = postHocPower_max_degListsCombinations
    deg_tasks$FindMarkers_topDEGs = findMarkers_topDEGs
    deg_tasks$RandomSeedNo = random_seed_nos
    deg_tasks$Max_cells_per_ident = findMarkers_max_cells_per_ident
    deg_tasks$Delete_existing_degs = delete_existing_degs
    deg_tasks$Max_cores_count = max_cores_count
    deg_tasks$Min_percent_randomSeedNos_for_postHocPowerAnalysis = postHocPower_minPercent_randomSeedNos
    deg_tasks$Average_with_csharp_script_and_do_postHocPowerAnalysis = TRUE
    deg_tasks$PostHocPower_coefOfVar_fileName_beginning = postHocPower_coefOfVar_fileName_beginning
    deg_tasks$PostHocPower_correlation_fileName_beginning = postHocPower_correlation_fileName_beginning
    deg_tasks$PostHocPower_eachClusterGene_fileName_beginning = postHocPower_eachClusterGene_fileName_beginning
    deg_tasks$PostHocPower_analysis_subdirectory = postHocPower_analysis_subdirectory
    deg_tasks$PostHocPower_max_number_of_datasets = postHocPower_max_number_of_datasets
    
    current_cellType_sig_documentations  = Calculate_degs_using_parallel_cores_if_necessary_and_eventually_add_plots_to_avg_degs_plots(core_seurat_object,bg_Genes,deg_tasks,current_cellType_sig_documentations,documentation_directory)

    rm(core_seurat_object)
    rm(deg_tasks)
  }#End - Calculate and write all markers for each simplified disease after deleting all files in degs directory
  
  if (calculate_pathways_for_simplified_diseases)
  {#Begin - calculate_pathways_for_simplified_diseases
    print(paste(report_term, "Calculate pathways for simplified diseases",sep=''))
    degs_directory = disease_degs_directory
    pathway_directory = disease_pathways_directory
    Calculate_pathways_using_mbc_pathNet(degs_directory, mbc_pathNet_directory, mbc_pathNet_exe_fileName, pathway_directory, ontologies, is_windows, is_linux)
    Generate_considered_pathways(pathway_directory, considered_pathway_fileName, pathway_pvalue_cutoff, cellTypeGroup_of_interest)

    wcFunct_definition_fileName = paste(wcfDef_fileName_start,dataset_assay,"_",cellTypeGroup_of_interest,diseaseFocusLabel,wcfDefGenerated_fileName_end,".txt",sep='')
    
    reference_pfs_list = list()
    reference_name = paste(dataset_assay,"_",cellTypeGroup_of_interest,diseaseFocusLabel,sep='')
    reference_fileNames = assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list[[reference_name]]
    if (length(reference_fileNames)>0)
    {#Begin
        fileNames = unlist(assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list)
        indexRef = 1
        for (indexRef in 1:length(reference_fileNames))
        {#Begin
          reference_fileName = reference_fileNames[indexRef]
          complete_reference_fileName = paste(reference_wholeCellFunctionDef_directory,reference_fileName,sep='')
          refs = read_excel(complete_reference_fileName)
          reference_pfs_list[[reference_fileName]] = refs
        }#End
    }#End
    Generate_and_write_physiological_functions_template(pathway_directory, considered_pathway_fileName, wcFunct_definition_fileName, reference_pfs_list)
    
  }#End - calculate_pathways_for_simplified_diseases
  
  complete_fileName = paste(subtype_analysis_directory,cellTypeGroup_of_interest,"_sig_degs_count.txt",sep='')
  write.table(current_cellType_sig_documentations,file=complete_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
  complete_fileName = paste(disease_analysis_directory,cellTypeGroup_of_interest,"_sig_degs_count.txt",sep='')
  write.table(current_cellType_sig_documentations,file=complete_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

  {if (length(identity_cellCounts[,1])>0)  { identity_cellCounts = rbind(identity_cellCounts,cellType_identity_lines) }
  else { identity_cellCounts = cellType_identity_lines }}

  {if (length(sig_documentations[,1])>0) { sig_documentations = rbind(sig_documentations,current_cellType_sig_documentations) }
  else { sig_documentations = current_cellType_sig_documentations }}
  
  gc()
  rm(cellTypeGroup_seurat_object)
  rm(assayCellTypeGroup)
  rm(cellTypeGroup_of_interest)
  rm(assayCellTypeGroupSubtypeL)
  rm(assayCellTypeGroupDisease)
  rm(current_cellTypes)
  #unlink(".RData")
}#End

if (calculateFinalDegs1_doPostHocPower2==1)
{#Begin - Add curated whole-cell functions to generated WcfDef from WcfDef processes in Excel and cross-compare WCF annotations accross cell types
  print(paste0(dataset_assay,": Add curated whole-cell functions to generated WcfDef from WcfDef processes in Excel and cross-compare WCF annotations accross cell types"))
  cellTypeGroups = names(cellTypeGroup_cellTypes_list)
  stopifnot(length(cellTypeGroups)>0)
  indexACTS = 2
  assayCellTypes = paste0(dataset_assay,"_",cellTypeGroups)
  
  wcfs_in_correct_order = names(wcfGroup_color_list)
  
  {#Begin - Define scp_annotation_base_line
     Col_names = c("AssayCellType","SCP","Whole_cell_function","Is_considered_elsewhere","Is_unselective_genes","Is_ignore_scp")
     Col_length = length(Col_names)
     Row_names = ""
     Row_length = length(Row_names)
     scp_annotation_base_line = as.data.frame(array(NA,c(Row_length,Col_length),dimnames = list(Row_names,Col_names)))
     scp_annotations = c()
  }#End - Define scp_annotation_base_line

  pathways_with_multiple_wcf_assignments=c()
  indexACTS = 7
  for (indexACTS in 1:length(assayCellTypes))
  {#Begin - Prepare supplied Wcf and copy into results directory
    assayCellType = assayCellTypes[indexACTS]
    input_fileName = paste0(wcfDef_fileName_start,assayCellType,subtype_level,".xlsx")
    complete_input_fileName = paste0(reference_wholeCellFunctionDef_directory,input_fileName)
    combined_pdDefs = c()
    if (file.exists(complete_input_fileName))
    {#Begin - Read subtype pfdefs
      pdDefs = read_excel(complete_input_fileName)
      pdDefs$Type = "Subtype"
      indexMissing = which(!pdDefs$PF_group %in% wcfs_in_correct_order)
      indexIsNotNa = which(!is.na(pdDefs$PF_group))
      indexMissing = intersect(indexMissing,indexIsNotNa)
      if (length(indexMissing)>0) { stop(paste0(input_fileName,": ",paste0(unique(pdDefs$PF_group[indexMissing]),collapse=","))) }
      if (length(combined_pdDefs[,1])>0)
      {#Begin 
        indexKeepCol = which(colnames(combined_pdDefs) %in% colnames(pdDefs))
        combined_pdDefs = combined_pdDefs[,indexKeepCol]
        indexKeepCol = which(colnames(pdDefs) %in% colnames(combined_pdDefs))
        pdDefs = pdDefs[,indexKeepCol]
        combined_pdDefs = rbind(combined_pdDefs,pdDefs)
      }#End
      else { combined_pdDefs = pdDefs }
    }#End - Read subtype pfdefs
    
    input_fileName = paste0(wcfDef_fileName_start,assayCellType,diseaseFocusLabel,".xlsx")
    complete_input_fileName = paste0(reference_wholeCellFunctionDef_directory,input_fileName)
    if (file.exists(complete_input_fileName))
    {#Begin - Read disease pfdefs
      pdDefs = read_excel(complete_input_fileName)
      pdDefs$Type = "Disease"
      indexMissing = which(!pdDefs$PF_group %in% wcfs_in_correct_order)
      indexIsNotNa = which(!is.na(pdDefs$PF_group))
      indexMissing = intersect(indexMissing,indexIsNotNa)
      if (length(indexMissing)>0) { stop(paste0(input_fileName,": ",paste0(unique(pdDefs$PF_group[indexMissing]),collapse=","))) }
      if (length(combined_pdDefs[,1])>0)
      {#Begin 
        indexKeepCol = which(colnames(combined_pdDefs) %in% colnames(pdDefs))
        combined_pdDefs = combined_pdDefs[,indexKeepCol]
        indexKeepCol = which(colnames(pdDefs) %in% colnames(combined_pdDefs))
        pdDefs = pdDefs[,indexKeepCol]
        combined_pdDefs = rbind(combined_pdDefs,pdDefs)
      }#End
      else { combined_pdDefs = pdDefs }
    }#End - Read disease pfdefs

    indexPFG=7
    combined_pdDefs$PieChart_slice_order = -1
    combined_pdDefs$Color = "white"
    combined_pdDefs$Hexadecimal_color = col2hex(combined_pdDefs$Color)
    for (indexPFG in 1:length(wcfs_in_correct_order))
    {#Begin - Set pie slice orders and colors
      pfGroup = wcfs_in_correct_order[indexPFG]
      indexCurrentPFGroup = which(combined_pdDefs$PF_group==pfGroup)
      if (length(indexCurrentPFGroup)>0)
      {#Begin
        current_combined =  combined_pdDefs[indexCurrentPFGroup,]
        current_combined$Color = wcfGroup_color_list[[pfGroup]]
        current_combined = current_combined[order(current_combined$SCP),]
        scps = unique(current_combined$SCP)
        for (indexScp in 1:length(scps))
        {#Begin
          scp = scps[indexScp]
          indexCurrentScp = which(current_combined$SCP==scp)
          current_combined$PieChart_slice_order[indexCurrentScp] = indexPFG * 100 + indexScp
        }#End
        combined_pdDefs[indexCurrentPFGroup,] = current_combined
      }#End
    }#End - Set pie slice orders and colors
    
    indexNoPf = c(which(combined_pdDefs$PF_group==""),which(is.na(combined_pdDefs$PF_group)))
    combined_pdDefs$PF_group[indexNoPf] = ""
    combined_pdDefs$Color[indexNoPf] = ""
    combined_pdDefs$Hexadecimal_color[indexNoPf] = ""
    
    indexPF = which(combined_pdDefs$PF_group!="")
    combined_pdDefs$Hexadecimal_color[indexPF] = col2hex(combined_pdDefs$Color[indexPF])
    indexNa = which(is.na(combined_pdDefs$Hexadecimal_color))
    combined_pdDefs$Hexadecimal_color[indexNa] = ""
    
    pathways = unique(combined_pdDefs$SCP)
    indexP=1
    for (indexP in 1:length(pathways))
    {#Begin - Add to scp annotations
      pathway = pathways[indexP]
      indexCurrentP = which(combined_pdDefs$SCP==pathway)
      scp_pdDefs = combined_pdDefs[indexCurrentP,]
      
      is_considered_elsewhere = length(grep("considered",scp_pdDefs$Comment)>0)
      is_unselective = length(grep("unselective",scp_pdDefs$Comment)>0)
      is_ignore = length(grep("ignore",scp_pdDefs$Comment)>0)

      scp_annotation_line = scp_annotation_base_line
      scp_annotation_line$Is_considered_elsewhere = is_considered_elsewhere
      scp_annotation_line$Is_unselective = is_unselective
      scp_annotation_line$Is_ignore_scp = is_ignore
      
      scp_annotation_line$AssayCellType = assayCellType
      scp_annotation_line$SCP = pathway
      
      local_pfGroups = unique(scp_pdDefs$PF_group)
      indexNone = which(local_pfGroups=="")
      local_pfGroups[indexNone] = "none"
      if (length(local_pfGroups)==1) { scp_annotation_line$Whole_cell_function = local_pfGroups }
      else { scp_annotation_line$Whole_cell_function = paste(local_pfGroups,collapse = ';') }
      rm(local_pfGroups)
      
      if (length(scp_annotations[,1])>0) { scp_annotations = rbind(scp_annotations,scp_annotation_line) }
      else { scp_annotations = scp_annotation_line }
    }#End - Add to scp annotations
    
    complete_scpAnnotations_fileName = paste(documentation_directory,"Scp_annotations_accross_cellTypes.txt",sep='')
    write.table(scp_annotations,file=complete_scpAnnotations_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
    
    
    pathways = unique(combined_pdDefs$SCP)
    indexP=1
    for (indexP in 1:length(pathways))
    {#Begin
      pathway = pathways[indexP]
      indexCurrentP = which(combined_pdDefs$SCP==pathway)
      indexConsideredElsewhere = c(grep("considered",combined_pdDefs$Comment),grep("unselective",combined_pdDefs$Comment))
      pathways_consideredElsewhere = unique(combined_pdDefs$SCP[indexConsideredElsewhere])
      indexPathwaysConsideredElsewhere = which(combined_pdDefs$SCP %in% pathways_consideredElsewhere)
      #indexExists = which(combined_pdDefs$Physiological_function!="")
      indexCurrentP = indexCurrentP[!indexCurrentP %in% indexPathwaysConsideredElsewhere]
      if (length(unique(combined_pdDefs$PF_group[indexCurrentP]))>1) { pathways_with_multiple_wcf_assignments = c(pathways_with_multiple_wcf_assignments, paste0(assayCellType, ": ", pathway)) }
    }#End
    

    indexSubtype = which(combined_pdDefs$Type == "Subtype")
    if (length(indexSubtype)>0)
    {#Begin
      pdDefs = combined_pdDefs[indexSubtype,]
      output_directory = paste0(base_directory,level1_subdirectory,assayCellType,subtype_level,"/",assayCellType,subtype_level,"_pathways/")
      complete_output_fileName = paste(output_directory,wcfDef_fileName_start,assayCellType,subtype_level,".txt",sep='')
      write.table(pdDefs,complete_output_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
    }#End
    indexDisease = which(combined_pdDefs$Type == "Disease")
    if (length(indexDisease)>0)
    {#Begin
      pdDefs = combined_pdDefs[indexDisease,]
      output_directory = paste0(base_directory,level1_subdirectory,assayCellType,diseaseFocusLabel,"/",assayCellType,diseaseFocusLabel,"_pathways/")
      complete_output_fileName = paste(output_directory,wcfDef_fileName_start,assayCellType,diseaseFocusLabel,".txt",sep='')
      write.table(pdDefs,complete_output_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
    }#End
  }#End - Prepare supplied Wcf and copy into results directory
  
  if (length(pathways_with_multiple_wcf_assignments)>0)
  {#Begin
     pathways_with_multiple_wcf_assignments = pathways_with_multiple_wcf_assignments[order(pathways_with_multiple_wcf_assignments)]
     print(paste0("Pathways with multiple whole cell function assignments: ",paste0(pathways_with_multiple_wcf_assignments,collaps=', ')))
  }#End
}#End - Add curated whole-cell functions to generated WcfDef from WcfDef processes in Excel and cross-compare WCF annotations accross cell types

if (calculateFinalDegs1_doPostHocPower2==1)
{#Begin - Check, if WcfDef generated by script are identical to WcfDef processes in Excel
  print(paste0(dataset_assay,": Check, if WcfDef generated by script are identical to WcfDef processes in Excel"))
  indexCellTypeGroup = 8
  for (indexCellTypeGroup in 1:length(cellTypeGroup_cellTypes_list))
  {#Begin
    cellTypeGroup_of_interest = names(cellTypeGroup_cellTypes_list)[indexCellTypeGroup]
    assayCellTypeGroup = paste(dataset_assay,"_",cellTypeGroup_of_interest,sep='')
    
    results_directory = paste(base_directory,level1_subdirectory,sep='')
    
    pathway_results_directory = paste0(results_directory,assayCellTypeGroup,subtype_level,"/",assayCellTypeGroup,subtype_level,"_pathways/",sep='')
    wcFunct_definition_fileName = paste(wcfDef_fileName_start,assayCellTypeGroup,subtype_level,wcfDefGenerated_fileName_end,".txt",sep='')
    complete_wcFunct_definition_fileName = paste(pathway_results_directory,wcFunct_definition_fileName,sep='')
    subtype_wcf_generated = read.csv(complete_wcFunct_definition_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)

    wcFunct_definition_fileName = paste(wcfDef_fileName_start,assayCellTypeGroup,subtype_level,".txt",sep='')
    complete_wcFunct_definition_fileName = paste(pathway_results_directory,wcFunct_definition_fileName,sep='')
    if (exists("subtype_wcf")) { rm(subtype_wcf) }
    if (file.exists(complete_wcFunct_definition_fileName))
    {#Begin
       subtype_wcf = read.csv(complete_wcFunct_definition_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
    }#End
    
    pathway_results_directory = paste0(results_directory,assayCellTypeGroup,diseaseFocusLabel,"/",assayCellTypeGroup,diseaseFocusLabel,"_pathways/",sep='')
    wcFunct_definition_fileName = paste(wcfDef_fileName_start,assayCellTypeGroup,diseaseFocusLabel,wcfDefGenerated_fileName_end,".txt",sep='')
    complete_wcFunct_definition_fileName = paste(pathway_results_directory,wcFunct_definition_fileName,sep='')
    disease_wcf_generated = read.csv(complete_wcFunct_definition_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
    
    wcFunct_definition_fileName = paste(wcfDef_fileName_start,assayCellTypeGroup,diseaseFocusLabel,".txt",sep='')
    complete_wcFunct_definition_fileName = paste(pathway_results_directory,wcFunct_definition_fileName,sep='')
    if (exists("disease_wcf")) { rm(disease_wcf) }
    if (file.exists(complete_wcFunct_definition_fileName))
    {#Begin
      disease_wcf = read.csv(complete_wcFunct_definition_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
    }#End
    

    wcf_list = list()
    if (exists("subtype_wcf")) { wcf_list[["Subtype"]] = subtype_wcf }
    if (exists("disease_wcf")) { wcf_list[["Disease"]] = disease_wcf }
    
    wcf_generated_list = list()
    wcf_generated_list[["Subtype"]] = subtype_wcf_generated
    wcf_generated_list[["Disease"]] = disease_wcf_generated

    clusterTypes = names(wcf_list)
    indexCT = 2
    for (indexCT in 1:length(clusterTypes))
    {#Begin
       clusterType = clusterTypes[indexCT]
       wcf = wcf_list[[clusterType]]
       wcf_generated = wcf_generated_list[[clusterType]]

       wcf = wcf[order(wcf$Symbol),]
       wcf = wcf[order(wcf$SCP),]
       wcf = wcf[order(wcf$PF_group),]
       
       wcf_generated = wcf_generated[order(wcf_generated$Symbol),]
       wcf_generated = wcf_generated[order(wcf_generated$SCP),]
       wcf_generated = wcf_generated[order(wcf_generated$PF_group),]

       indexMismatch = which(wcf_generated$PF_group != wcf$PF_group)
       stopifnot(length(indexMismatch)==0)
       indexMismatch = which(wcf_generated$SCP != wcf$SCP)
       stopifnot(length(indexMismatch)==0)
       indexMismatch = which(wcf_generated$Symbol != wcf$Symbol)
       stopifnot(length(indexMismatch)==0)
       indexMismatch = which(wcf_generated$Cell_subtypes..p.value. != wcf$Cell_subtypes..p.value.)
       if (cellTypeGroup_of_interest=="Lymphoid")
       {#Begin - Ignore lines containing subtype Naive Th due to non-standard characters
         indexContainsNaiveTh = intersect(grep("ve Th",wcf_generated$Cell_subtypes..p.value.),grep("Na",wcf_generated$Cell_subtypes..p.value.))
         indexMismatch = indexMismatch[!indexMismatch %in% indexContainsNaiveTh]
       }#End - Ignore lines containing subtype Naive Th due to non-standard characters
       stopifnot(length(indexMismatch)==0)
       print(paste0(assayCellTypeGroup,": ",clusterType,": No mismatches found",sep=''))
    }#End
  }#End
}#End - Check, if WcfDef generated by script are identical to WcfDef processes in Excel

if (calculateFinalDegs1_doPostHocPower2==1)
{#Begin - Add Whole cell functions to enrichment result tables and copy into supplementary directory
  Add_wcf_to_enrich_results_and_restore_original_colnames = function(enrich, wcf, pfGroup_color_list)
  {#Begin - Add_wcf_to_enrich_results
     enrich$Dataset.color = NULL
     colnames_in_correct_order = colnames(enrich)
     enrich$WholeCellFunction = "error"
     enrich$Wcf.color = "error"
     enrich$References.supporting.wcf.annotation = ""
     colnames_in_correct_order = c("WholeCellFunction","Wcf.color",colnames_in_correct_order,"References.supporting.wcf.annotation")
     enrich = enrich[,colnames_in_correct_order]

     
     scps = unique(enrich$SCP)
     indexScp=1
     for (indexScp in 1:length(scps))
     {#Begin
        scp = scps[indexScp]
        indexCurrentScp_in_enrich = which(enrich$SCP==scp)
        indexCurrentScp_in_wcf = which(wcf$SCP==scp)
        current_wcfName = unique(wcf$PF_group[indexCurrentScp_in_wcf])
        current_color = ""
        if ((current_wcfName!="")&(!is.na(current_wcfName)))
        { current_color = pfGroup_color_list[[current_wcfName]] }
        current_references = unique(wcf$Reference[indexCurrentScp_in_wcf])
        indexNotNa = which(!is.na(current_references))
        current_references = current_references[indexNotNa]
        if (length(current_references)>1) { current_references = paste0(current_references,collapse = ";") }
        stopifnot(length(current_wcfName)==1)
        stopifnot(length(current_color)==1)
        enrich$WholeCellFunction[indexCurrentScp_in_enrich] = current_wcfName
        enrich$Wcf.color[indexCurrentScp_in_enrich] = current_color
        if (length(current_references)>0)
        { enrich$References.supporting.wcf.annotation[indexCurrentScp_in_enrich] = current_references }
     }#End
     colnames(enrich) = gsub("[.]"," ",colnames(enrich))
     enrich = enrich[order(enrich$SCP),]
     enrich = enrich[order(enrich$WholeCellFunction),]
     indexNoWcf = which(enrich$WholeCellFunction=="")
     if (length(indexNoWcf)>0)
     {#Begin
        indexWcf = 1:length(enrich[,1])
        indexWcf = indexWcf[!indexWcf %in% indexNoWcf]
        enrich = rbind(enrich[indexWcf,],enrich[indexNoWcf,])
     }#End
     return (enrich)
  }#End - Add_wcf_to_enrich_results

  print(paste0(dataset_assay,": Add Whole cell functions to enrichment result tables and copy into supplementary directory"))
  results_directory = complete_level1_directory
  supplementary_data_directory = paste(results_directory,"AA_Supplementary_data/",sep='')
  if (!dir.exists(supplementary_data_directory)) { dir.create(supplementary_data_directory) }
  wb = createWorkbook()

  combined_subtype_pathways = c()
  combined_disease_pathways = c()
  
  indexCellTypeGroup = 1
  for (indexCellTypeGroup in 1:length(cellTypeGroup_cellTypes_list))
  {#Begin
     cellTypeGroup_of_interest = names(cellTypeGroup_cellTypes_list)[indexCellTypeGroup]
     assayCellTypeGroup = paste(dataset_assay,"_",cellTypeGroup_of_interest,sep='')
     enrich_fileName = "Mbco_human_standard_significantPredictions.txt";

     assayCellTypeGroupSubtypeLevel = paste(assayCellTypeGroup,subtype_level,sep='')
     assayCellTypeGroupDisease = paste(assayCellTypeGroup,diseaseFocusLabel,sep='')

     subtype_wcfDef_fileName = paste(wcfDef_fileName_start,assayCellTypeGroupSubtypeLevel,wcfDefGenerated_fileName_end,".txt",sep='')
     complete_subtype_enrichFileName = paste(results_directory,assayCellTypeGroupSubtypeLevel,"/",assayCellTypeGroupSubtypeLevel,pathways_directoryAddition,"/",enrich_fileName,sep='')
     subtype_enrich = read.csv(complete_subtype_enrichFileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
     complete_subtype_wcfDef_fileName = paste(results_directory,assayCellTypeGroupSubtypeLevel,"/",assayCellTypeGroupSubtypeLevel,pathways_directoryAddition,"/",subtype_wcfDef_fileName,sep='')
     subtype_wcfDef = read.csv(complete_subtype_wcfDef_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
     
     subtype_enrich = Add_wcf_to_enrich_results_and_restore_original_colnames(enrich=subtype_enrich, wcf=subtype_wcfDef, pfGroup_color_list = pfGroup_color_list)
     complete_output_subtype_fileName = paste(supplementary_data_directory,assayCellTypeGroupSubtypeLevel,"_pathways.txt",sep='')
     write.table(subtype_enrich,file=complete_output_subtype_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

     if (assayCellTypeGroupSubtypeLevel %in% names(assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list))
     {#Begin
        addWorksheet(wb, assayCellTypeGroupSubtypeLevel)
        writeData(wb,assayCellTypeGroupSubtypeLevel,subtype_enrich)
        
        subtype_enrich$Assay = dataset_assay
        subtype_enrich$Cell_type_group = cellTypeGroup_of_interest
        subtype_enrich$Subclass_level = subtype_level
        if (length(combined_subtype_pathways[,1])>0) { combined_subtype_pathways = rbind(combined_subtype_pathways,subtype_enrich) }
        else { combined_subtype_pathways = subtype_enrich }
     }#End
          
     disease_wcfDef_fileName = paste(wcfDef_fileName_start,assayCellTypeGroupDisease,wcfDefGenerated_fileName_end,".txt",sep='')
     complete_disease_enrichFileName = paste(results_directory,assayCellTypeGroupDisease,"/",assayCellTypeGroupDisease,pathways_directoryAddition,"/",enrich_fileName,sep='')
     disease_enrich = read.csv(complete_disease_enrichFileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
     complete_disease_wcfDef_fileName = paste(results_directory,assayCellTypeGroupDisease,"/",assayCellTypeGroupDisease,pathways_directoryAddition,"/",disease_wcfDef_fileName,sep='')
     disease_wcfDef = read.csv(complete_disease_wcfDef_fileName,header=TRUE,sep='\t',stringsAsFactors = FALSE)
     
     disease_enrich = Add_wcf_to_enrich_results_and_restore_original_colnames(enrich=disease_enrich, wcf=disease_wcfDef, pfGroup_color_list = pfGroup_color_list)
     complete_output_disease_fileName = paste(supplementary_data_directory,assayCellTypeGroupDisease,"_pathways.txt",sep='')
     write.table(disease_enrich,file=complete_output_disease_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
     
     if (assayCellTypeGroupDisease %in% names(assayCellTypeGroup_orderedSubtypeRPFLsFileNames_list))
     {#Begin
       addWorksheet(wb, assayCellTypeGroupDisease)
       writeData(wb,assayCellTypeGroupDisease,disease_enrich)
       
       disease_enrich$Assay = dataset_assay
       disease_enrich$Cell_type_group = cellTypeGroup_of_interest
       if (length(combined_disease_pathways[,1])>0) { combined_disease_pathways = rbind(combined_disease_pathways,disease_enrich) }
       else { combined_disease_pathways = disease_enrich }
     }#End

     rm(assayCellTypeGroup)
     rm(cellTypeGroup_of_interest)
  }#End
  
  complete_excel_fileName = paste(supplementary_data_directory,"Supplementary_table_20_sheets.xlsx",sep='')
  saveWorkbook(wb, complete_excel_fileName, overwrite = TRUE)
  
  colnames_subtypes = colnames(combined_subtype_pathways)
  first_colnames = c("Assay","Cell_type_group","Subclass_level","Ontology")
  stopifnot(length(which(!first_colnames %in% colnames_subtypes))==0)
  rest_colnames = colnames_subtypes[!colnames_subtypes %in% first_colnames]
  colnames_in_correct_order = c(first_colnames, rest_colnames)
  combined_subtype_pathways = combined_subtype_pathways[,colnames_in_correct_order]
  
  complete_combined_subtypePathways_fileName = paste(supplementary_data_directory,"Source_data_file.xlsx",sep='')
  write.xlsx(combined_subtype_pathways,complete_combined_subtypePathways_fileName)

  colnames_disease = colnames(combined_disease_pathways)
  first_colnames = c("Assay","Cell_type_group","Ontology")
  stopifnot(length(which(!first_colnames %in% colnames_disease))==0)
  rest_colnames = colnames_disease[!colnames_disease %in% first_colnames]
  colnames_in_correct_order = c(first_colnames, rest_colnames)
  combined_disease_pathways = combined_disease_pathways[,colnames_in_correct_order]

  complete_combined_diseasePathways_fileName = paste(supplementary_data_directory,"Supplementary_table_20.xlsx",sep='')
  write.xlsx(combined_disease_pathways,complete_combined_diseasePathways_fileName)
  rm(wb)
}#End - Add Whole cell functions to enrichment result tables and copy into supplementary directory

gc()
rm(assay_seurat_object)
rm(dataset_assay)
rm(snRNAseq_assay)

}#End - Analyze for given dataset assay

complete_identityCounts_fileName = paste(combined_results_directory,"Identity_cellCounts.txt",sep='')
write.table(identity_cellCounts,file=complete_identityCounts_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)


complete_sig_documentation_fileName = paste(combined_results_directory,"Sig_documentations.txt",sep='')
write.table(sig_documentations,file=complete_sig_documentation_fileName,row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

gc()
unlink(".RData")
rm(list = ls(all.names=TRUE));
gc()

  
