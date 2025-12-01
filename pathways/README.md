Code for pathway analyses on single nucleus Atlas V2 cell subtypes, states and disease conditions

- To start the analysis open the 'KPMP_atlas_v2_calculate_degs_and_pathways.R' file in the folder 'KPMP_atlas_v2_code' subfolder. Follow the instructions in the section 'User input is needed here'.
- As instructed in that file, the KPMP atlas v2 Single Nucleus seruat object has to be downloaded and copied into the 'SN_RNAseq_atlas_v2_2024November01' subdirectory.
- As instructed in that file, the application MBC PathNet has to be downloaded from 'mbc-ontology.org'/'https://github.com/SBCNY/Molecular-Biology-of-the-Cell' and copied into the specified directory.
- As instructed in that file, 'mono' has to be installed for running MBC PathNet in a Linux environment and NET 8 runtime has to installed for running the 'Average_DEGs_and_do_postHocPowerAnalysis_linux_x64' file in a LINUX environment
- The subfolder 'Average_DEGs_and_do_postHocPowerAnalysis' contains the source C# code that was published into the exe-file in 'Average_DEGs_and_do_postHocPowerAnalysis_linux_x64'. While the source code is not used, the exe file will be called from within the R-script. 
- The script 'Generate_pathwaySubtype_networks' can be run on Windows with Microsoft Visual Studio.

The R-script 'KPMP_atlas_v2_calculate_degs_and_pathways.R' will generate the data that is used by the 'Generate_pathwaySubtype_networks' script.
