/*
This script was written by Jens Hansen working for the Ravi Iyengar Lab (Icahn School of Medicine at Mount Sinai, New York) and the Kidney Precision Medicine Project 
 
Please acknoledge our work by citing Lake et. al bioRxiv 2025.09.26.678707.

It can be opened and debugged with Visual Studio in a Windows environment.

*/

using Common_functions;
using Enrichment_results;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using Trajectories;

class Main_class
{

    public static Enrichment_results_class Combine_scps_annotated_to_same_mergedScp_by_summing_up_minusLog10Pvalues_for_same_level(Enrichment_results_class sameLevel_enrichment_results, Dictionary<string, string> scp_mergedScp_dict)
    {
        #region Generate enich_pieChart by replacing scp by hfGroup and sum up duplicates
        Timeunit_enum timeunit = sameLevel_enrichment_results.Enrich[0].Timeunit;
        int scp_level = sameLevel_enrichment_results.Enrich[0].SCP_level;
        Enrichment_results_class enrich_mergedScp = sameLevel_enrichment_results.Deep_copy();
        int enrich_length = enrich_mergedScp.Enrich.Length;
        List<Enrichment_results_line_class> keep_enrich_lines = new List<Enrichment_results_line_class>();
        Enrichment_results_line_class collapsed_enrich_line;
        Enrichment_results_line_class enrich_line;
        for (int indexE = 0; indexE < enrich_length; indexE++)
        {
            enrich_line = enrich_mergedScp.Enrich[indexE];
            if (!enrich_line.Timeunit.Equals(timeunit)) { throw new Exception(); }
            if (!enrich_line.SCP_level.Equals(scp_level)) { throw new Exception(); }
            if (scp_mergedScp_dict.ContainsKey(enrich_line.SCP))
            {
                enrich_line.SCP = (string)scp_mergedScp_dict[enrich_line.SCP].Clone();
            }
            else
            {
                enrich_line.SCP = "Remove";
            }
            if (String.IsNullOrEmpty(enrich_line.SCP)) { throw new Exception(); }
        }
        enrich_mergedScp.Order_by_ontology_scp_datasetName_timepoint_upDownStatus();
        for (int indexE = 0; indexE < enrich_length; indexE++)
        {
            enrich_line = enrich_mergedScp.Enrich[indexE];
            if (!enrich_line.Timeunit.Equals(timeunit)) { throw new Exception(); }
            if ((indexE == enrich_length - 1)
                || (!enrich_line.SCP.Equals(enrich_mergedScp.Enrich[indexE + 1].SCP))
                || (!enrich_line.Dataset_name.Equals(enrich_mergedScp.Enrich[indexE + 1].Dataset_name))
                || (!enrich_line.Timepoint.Equals(enrich_mergedScp.Enrich[indexE + 1].Timepoint))
                || (!enrich_line.UpDown_status.Equals(enrich_mergedScp.Enrich[indexE + 1].UpDown_status)))
            {
                if (!enrich_line.SCP.Equals("Remove"))
                {
                    collapsed_enrich_line = new Enrichment_results_line_class();
                    switch (enrich_line.Ontology)
                    {
                        case Ontology_type_enum.Mbco:
                            collapsed_enrich_line.Ontology = Ontology_type_enum.Whole_cell_subfunction;
                            break;
                        case Ontology_type_enum.Whole_cell_subfunction:
                            collapsed_enrich_line.Ontology = Ontology_type_enum.Whole_cell_function;
                            break;
                        default:
                            throw new Exception();
                    }
                    collapsed_enrich_line.SCP_level = enrich_line.SCP_level;
                    collapsed_enrich_line.SCP = (string)enrich_line.SCP.Clone();
                    collapsed_enrich_line.Dataset_name = (string)enrich_line.Dataset_name.Clone();
                    collapsed_enrich_line.UpDown_status = enrich_line.UpDown_status;
                    collapsed_enrich_line.Timepoint = enrich_line.Timepoint;
                    collapsed_enrich_line.Timeunit = enrich_line.Timeunit;
                    collapsed_enrich_line.Minus_log10_pvalue = enrich_line.Minus_log10_pvalue;
                    keep_enrich_lines.Add(collapsed_enrich_line);
                }
            }
            else
            {
                enrich_mergedScp.Enrich[indexE + 1].Minus_log10_pvalue += enrich_line.Minus_log10_pvalue;
            }
        }
        enrich_mergedScp.Enrich = keep_enrich_lines.ToArray();
        #endregion

        return enrich_mergedScp;
    }

    public static void Generate_physiologicalFunctions_pfGroups_and_write_subtype_networks_for_given_cellType(int scp_level, string assay_cellType, string baseDirectory)
    {
        string input_enrichment_directory = baseDirectory + Global_directory_class.Get_enrichment_results_subdirectory(assay_cellType);
        string input_enrichment_fileName = Global_directory_class.Enrichment_results_input_fileName;
        Enrichment_results_class enrich = new Enrichment_results_class();
        enrich.Read(input_enrichment_directory, input_enrichment_fileName);
        enrich.Keep_only_lines_of_given_scpLevel(scp_level);

        string input_physiologicalFunctionsDefs_directory = baseDirectory + Global_directory_class.Get_physiologicalFunctionsDefinitions_subdirectory(assay_cellType);
        string input_physiologicalFunctionsDefs_fileName = Global_directory_class.Get_physiologicalFunctionDefinitions_fileName(assay_cellType);
        PhysiologicalFunctionsDefinition_class phyFunDefs = new PhysiologicalFunctionsDefinition_class();
        phyFunDefs.Generate_by_reading(input_physiologicalFunctionsDefs_directory, input_physiologicalFunctionsDefs_fileName);

        Dictionary<string, string> scp_wholeCellFunction_dict = phyFunDefs.Get_scp_wholeCellFunction_dict();
        Dictionary<string, int> pfGroup_pieChartNo_dict = phyFunDefs.Get_pfGroup_pieChartSliceNumber_dict();
        Dictionary<string, System.Drawing.Color> pfGroup_color_dict = phyFunDefs.Get_pfGroup_color_dict();

        Enrichment_results_class pfGroups = Combine_scps_annotated_to_same_mergedScp_by_summing_up_minusLog10Pvalues_for_same_level(enrich, scp_wholeCellFunction_dict);
        string output_pfGroups_directory = baseDirectory + Global_directory_class.Get_physiologicalFunctions_subdirectory(assay_cellType);
        string output_pfGroups_fileName = Global_directory_class.Get_pFGroup_fileName(assay_cellType);
        pfGroups.Write(output_pfGroups_directory, output_pfGroups_fileName);

        Trajectory_class trajectories = new Trajectory_class();
        trajectories.Generate_for_assay_cellType(assay_cellType);

        string output_nw_directory = baseDirectory + Global_directory_class.Get_physiologicalFunctions_subdirectory(assay_cellType);
        string output_nw_baseFileName = Global_directory_class.Get_pfGroup_nw_basefileName();
        KPMP_cluster_hierarchy_class cluster_hierarchy = new KPMP_cluster_hierarchy_class(assay_cellType);
        cluster_hierarchy.Options.Forced_reference_category_value = 100;
        cluster_hierarchy.Generate_from_trajectories_and_enrichment_results(trajectories, pfGroups, assay_cellType);
        cluster_hierarchy.Write_as_network_for_each_cellType(KPMP_category_for_cluster_enum.Scp_minusLog10Pvalues, pfGroup_pieChartNo_dict, pfGroup_color_dict, output_nw_directory, output_nw_baseFileName);

    }
    public static void Main(params string[] arguments)
    {
        #region User input needed
        string overall_directory = "D:/KPMP_v2_atlas/KPMP_v2_data/KPMP_atlas_v2_max5000cells_500DEGs/";
        #endregion

        string assay = "SNboth";
        string[] cellTypes_subclassl3 = new string[] { "PODPEC", "PT", "TAL", "DCT", "FIB", "Myeloid" };
        string[] cellTypes_simD = new string[] { "PODPEC", "PT", "DTLATL", "TAL", "DCT", "CDCNT", "EC", "FIB", "Lymphoid", "Myeloid" };
        string subclassLevel = "sl3";
        List<string> assayCellTypeSubclassLevels_list = new List<string>();
        foreach (string cellType_subclassl3 in cellTypes_subclassl3)
        {
            assayCellTypeSubclassLevels_list.Add(assay + "_" + cellType_subclassl3 + subclassLevel);
        }
        foreach (string cellType_simD in cellTypes_simD)
        {
            assayCellTypeSubclassLevels_list.Add(assay + "_" + cellType_simD + "SimD");
        }
        string[] assayCellTypeSubclassLevels = assayCellTypeSubclassLevels_list.ToArray();
        int scpLevel = 3;
        foreach (string assayCellTypeSubclassLevel in assayCellTypeSubclassLevels)
        {
            Generate_physiologicalFunctions_pfGroups_and_write_subtype_networks_for_given_cellType(scpLevel, assayCellTypeSubclassLevel, overall_directory + assayCellTypeSubclassLevel + "/");
        }
    }

}

