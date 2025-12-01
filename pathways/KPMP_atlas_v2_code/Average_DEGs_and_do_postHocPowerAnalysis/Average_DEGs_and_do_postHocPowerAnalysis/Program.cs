/*
This script was written by Jens Hansen working for the Ravi Iyengar Lab (Icahn School of Medicine at Mount Sinai, New York) and the Kidney Precision Medicine Project 
 
Please acknoledge our work by citing Lake et. al bioRxiv 2025.09.26.678707.

It can be opened and published using Visual Studio in a Windows environment.
This script is used to publish a file set including an executable file that will be called by an R script running in a Linux environment.

Use 'Build' - 'Publish' to generate executable file for a Linux environment
Target location: './Average_DEGs_and_do_postHocPowerAnalysis_linux_x64'

Profile setting:
Configuration: 'Realeas | Any CPU'
Target framework: 'net8.0'
Deployment mode: 'Framework-dependent'
Target runtime: 'linux-x64'
Target_location: My_directory
File publish options: unselect 'Produce single file'
*/



using DEGs_and_postHocPower;
using Common_functions;
using System.Net.Http.Headers;
using System.Runtime.CompilerServices;
using Average_DEGs_and_do_postHocPowerAnalysis.Common_functions;
using Average_DEGs_and_do_postHocPowerAnalysis;
using System.Text;
using ReadWrite;

class Main_class
{
    public static void Main(params string[] arguments)
    {
        string cellType = "PTsl3";

        //System.Diagnostics.Debugger.Launch();
        string overall_directory = "D:/KPMP_v2_atlas/KPMP_v2_data/KPMP_atlas_v2_max5000cells_500DEGs_PHP/";
        string report_directory = overall_directory;
        //string assay_cellTypeGroup_subtype_level = "SNboth_PTsl3";
        string assay_cellTypeGroup_subtype_level = "SNboth_" + cellType;
        int arguments_length = arguments.Length;
        string argument;
        bool overall_directory_defined = false;
        bool assay_cellTypeGroup_subtype_level_defined = false;
        int maxDegListCombinations = 100;
        int min_percentage_of_degLists = 0; //min_number_of_datasets will be the max of 1 and the calculated number
        int numberOfDatasets_for_final_analysis = 100;
        float deg_max_adjusted_pvalue = 0.05F;
        float deg_max_rank = 500;
        int max_numberOfDatasets_postHocPower = 100;
        int randomSeedNo_number = -1;
        int maxNoDegListsInMemory = 100;
        List<string> unknown_arguments = new List<string>();
        string overall_directory_argumentLabel = "--overall_directory";
        string assay_cellTypeGroup_subtypeLevel_argumentLabel = "--assay_cellTypeGroup_subtypeLevel";
        string max_degListsCombinations_argumentLabel = "--max_degListsCombinations";
        string min_percentageOfDrawnDegLists_argumentLabel = "--min_percentageOfDrawnDegLists";
        string numberOfDatasets_finalAnalysis_argumentLabel = "--numberOfDatasets_finalAnalysis";
        string max_numberOfDatasets_postHocPower_argumentLabel = "--max_numberOfDatasets_postHocPower";
        string max_adjPvalue_argumentLabel = "--max_adj_pvalue";
        string max_rank_argumentLabel = "--max_rank";
        string randomSeedNo_argumentLabel = "--randomSeedNodeNoPostHocPower";
        string maxNoDegListsInMemory_argumentLabel = "--maxNoDegListsInMem";
        List<string> available_arguments_list = new List<string>();
        available_arguments_list.Add(overall_directory_argumentLabel);
        available_arguments_list.Add(assay_cellTypeGroup_subtypeLevel_argumentLabel);
        available_arguments_list.Add(max_degListsCombinations_argumentLabel);
        available_arguments_list.Add(min_percentageOfDrawnDegLists_argumentLabel);
        available_arguments_list.Add(numberOfDatasets_finalAnalysis_argumentLabel);
        available_arguments_list.Add(max_numberOfDatasets_postHocPower_argumentLabel);
        available_arguments_list.Add(max_adjPvalue_argumentLabel);
        available_arguments_list.Add(max_rank_argumentLabel);
        available_arguments_list.Add(randomSeedNo_argumentLabel);
        available_arguments_list.Add(maxNoDegListsInMemory_argumentLabel);
        for (int indexA = 0; indexA < arguments_length; indexA++)
        {
            argument = arguments[indexA];
            if ((argument.Equals(assay_cellTypeGroup_subtypeLevel_argumentLabel))
                && (indexA < arguments_length - 1))
            {
                indexA++;
                assay_cellTypeGroup_subtype_level = arguments[indexA];
                assay_cellTypeGroup_subtype_level_defined = true;
            }
            else if ((argument.Equals(overall_directory_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                overall_directory = arguments[indexA];
                overall_directory_defined = true;
            }
            else if ((argument.Equals(max_degListsCombinations_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!int.TryParse(arguments[indexA], out maxDegListCombinations))
                {
                    throw new Exception(max_degListsCombinations_argumentLabel + " is not followed by integer");
                }
            }
            else if ((argument.Equals(min_percentageOfDrawnDegLists_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!int.TryParse(arguments[indexA], out min_percentage_of_degLists))
                {
                    throw new Exception(min_percentageOfDrawnDegLists_argumentLabel + " is not followed by integer");
                }
            }
            else if ((argument.Equals(maxNoDegListsInMemory_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!int.TryParse(arguments[indexA], out maxNoDegListsInMemory))
                {
                    throw new Exception(maxNoDegListsInMemory_argumentLabel + " is not followed by integer");
                }
            }
            else if ((argument.Equals(randomSeedNo_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!int.TryParse(arguments[indexA], out randomSeedNo_number))
                {
                    throw new Exception(randomSeedNo_argumentLabel + " is not followed by integer");
                }
            }
            else if ((argument.Equals(numberOfDatasets_finalAnalysis_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!int.TryParse(arguments[indexA], out numberOfDatasets_for_final_analysis))
                {
                    throw new Exception(numberOfDatasets_finalAnalysis_argumentLabel + " is not followed by integer");
                }
            }
            else if ((argument.Equals(max_numberOfDatasets_postHocPower_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!int.TryParse(arguments[indexA], out max_numberOfDatasets_postHocPower))
                {
                    throw new Exception(max_numberOfDatasets_postHocPower_argumentLabel + " is not followed by integer");
                }
            }
            else if ((argument.Equals(max_adjPvalue_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!float.TryParse(arguments[indexA], out deg_max_adjusted_pvalue))
                {
                    throw new Exception(max_adjPvalue_argumentLabel + " is not followed by integer");
                }
            }
            else if ((argument.Equals(max_rank_argumentLabel))
                        && (indexA < arguments_length - 1))
            {
                indexA++;
                if (!float.TryParse(arguments[indexA], out deg_max_rank))
                {
                    throw new Exception(max_rank_argumentLabel + " is not followed by integer");
                }
            }
            else
            {
                unknown_arguments.Add(argument);
            }
        }
        if ((arguments.Length > 0) && (!assay_cellTypeGroup_subtype_level_defined)) { throw new Exception(assay_cellTypeGroup_subtypeLevel_argumentLabel + " argument is missing"); }
        if ((arguments.Length > 0) && (!overall_directory_defined)) { throw new Exception(overall_directory_argumentLabel + " has not been defined"); }
        if (unknown_arguments.Count > 0)
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Arguments ");
            foreach (string unknown_argument in unknown_arguments)
            {
                if (sb.Length > 0) { sb.Append(", "); }
                sb.Append("'" + unknown_argument + "'");
            }
            sb.Append(" is/are unknown or not followed by parameter. Available arguments:");
            string available_argument;
            int available_arguments_length = available_arguments_list.Count;
            for (int indexA=0; indexA<available_arguments_length; indexA++)
            {
                available_argument = available_arguments_list[indexA];
                if (indexA!=0) { sb.Append(", "); }
                sb.AppendFormat("{0}", available_argument);
            }
            throw new Exception(sb.ToString());
        }

        Seurat_data_processing_class seurat_processing = new Seurat_data_processing_class(overall_directory, assay_cellTypeGroup_subtype_level, randomSeedNo_number);
        seurat_processing.Options.Overall_directory = (string)overall_directory.Clone();
        seurat_processing.Options.Assay_cellTypeGroup_subtype_level = (string)assay_cellTypeGroup_subtype_level.Clone();
        seurat_processing.Options.PostHocPower_max_randomCombinations = maxDegListCombinations;
        seurat_processing.Options.PostHocPower_min_percentage_of_datasets = min_percentage_of_degLists;
        seurat_processing.Options.PostHocPower_max_number_of_datasets = max_numberOfDatasets_postHocPower;
        seurat_processing.Options.Number_of_datasets_final_data = numberOfDatasets_for_final_analysis;
        seurat_processing.Options.Consider_only_genes_that_were_significant_at_least_once_for_coefOfVar_and_pearsonCor = true;
        seurat_processing.Options.Deg_maxRank = deg_max_rank;
        seurat_processing.Options.Max_number_deg_lists_in_memory = maxNoDegListsInMemory;
        seurat_processing.Options.Deg_maxAdjPvalue = deg_max_adjusted_pvalue;
        seurat_processing.Options.Consider_only_degs_significant_in_final_for_postHocPower = true;
        seurat_processing.Write_options(randomSeedNo_number);

        ReadWriteClass.Write_one_type_array_as_single_column(available_arguments_list.ToArray(), seurat_processing.Options.Report_directory, "Available_arguments.txt");

        seurat_processing.Generate_for_all_runs();
        seurat_processing.Calculate_and_write_avgDEGs_using_first_random_seedNo_results_of_given_number_of_datasets();
        if ((seurat_processing.Options.PostHocPower_max_randomCombinations > 0)
            && (seurat_processing.Options.PostHocPower_max_number_of_datasets > 0)
            && (seurat_processing.Options.PostHocPower_min_percentage_of_datasets <= 100)
            && (seurat_processing.Options.PostHocPower_min_percentage_of_datasets >= 0))
        {
            seurat_processing.Do_postHocPower_analysis();
        }
    }
}

