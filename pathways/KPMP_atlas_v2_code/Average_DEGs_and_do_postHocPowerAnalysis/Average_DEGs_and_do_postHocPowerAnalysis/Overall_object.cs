using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Average_DEGs_and_do_postHocPowerAnalysis.Common_functions;
using Common_functions;
using DEGs_and_postHocPower;
using ReadWrite;

namespace Average_DEGs_and_do_postHocPowerAnalysis
{
    class Seurat_data_processing_options_class
    {
        public string Overall_directory { get; set; }
        public string Assay_cellTypeGroup_subtype_level { get; set; }
        public int PostHocPower_max_randomCombinations { get; set; }
        public int PostHocPower_min_percentage_of_datasets { get; set; }
        public int PostHocPower_max_number_of_datasets { get; set; }
        public int Number_of_datasets_final_data { get; set; }
        public float Deg_maxAdjPvalue { get; set; }
        public bool Consider_only_degs_significant_in_final_for_postHocPower { get; set; }
        private float Deg_maxAdjPvalue_factor_realToPostHocPower { get; set; }
        public float Deg_maxAdjPvalue_for_consideration_in_postHocPower
        { get { return Deg_maxAdjPvalue * Deg_maxAdjPvalue_factor_realToPostHocPower; } }
        public float Deg_maxRank { get; set; }
        private float Deg_maxRank_factor_realToPostHocPower { get; set; }
        public float Deg_maxRank_for_consideration_in_postHocPower
        { get { return Deg_maxRank * Deg_maxRank_factor_realToPostHocPower; } }

        public string Report_directory {get;set;}
        public string Warning_directory { get; set; }
        public string RandomSeedNo_degs_subdirectory { get; set; }
        public string DEGs_directory { get { return Overall_directory + "/" + Assay_cellTypeGroup_subtype_level + "/" + Assay_cellTypeGroup_subtype_level + DEGs_subdirectory_name_ending +  "/"; } }
        public string DEGs_subdirectory_name_ending { get; set; }
        public string DEGsNoCutoff_subdirectory { get; set; }
        public string DEGsNoCutoff_directory { get { return DEGs_directory + DEGsNoCutoff_subdirectory; } }
        public string DEGsNoCutoff_fileName {  get { return DEGsAllMarker_fileName_beginning + Assay_cellTypeGroup_subtype_level + ".txt"; } }
        public string DEGsAllMarker_fileName_beginning { get; set; }
        public string RandomSeedNo_degs_directory { get { return DEGs_directory + RandomSeedNo_degs_subdirectory; } }
        public string PostHocPower_subdirectory { get; set; }
        public string PostHocPower_directory {  get { return DEGs_directory + PostHocPower_subdirectory; } }
        public string PostHocPower_coefOfVar_fileName { get { return PostHocPower_coefOfVar_fileName_beginning + Assay_cellTypeGroup_subtype_level + ".txt"; } }
        public string PostHocPower_correlation_fileName { get { return PostHocPower_correlation_fileName_beginning + Assay_cellTypeGroup_subtype_level + ".txt"; } }
        public string PostHocPower_eachClusterGene_fileName { get { return PostHocPower_eachClusterGene_fileName_beginning + Assay_cellTypeGroup_subtype_level + ".txt"; } }
        private string PostHocPower_correlation_fileName_beginning { get; set; }
        private string PostHocPower_coefOfVar_fileName_beginning { get; set; }
        private string PostHocPower_eachClusterGene_fileName_beginning { get; set; }
        public bool Use_data_normalized_towards_values_of_highest_datasetNos_for_postHocPower { get; set; }
        public bool Consider_only_genes_that_were_significant_at_least_once_for_coefOfVar_and_pearsonCor { get; set; }
        public int Max_number_deg_lists_in_memory { get; set; }
        public string RanSeedNo_label { get; set; }
        public Seurat_data_processing_options_class(string overall_directory, string assayCellTypeGrouplevel)
        {
            Overall_directory = (string)overall_directory.Clone();
            Assay_cellTypeGroup_subtype_level = (string)assayCellTypeGrouplevel.Clone();
            Report_directory = Overall_directory + "/AA_Documentations/";
            Warning_directory = Overall_directory + "/";
            PostHocPower_max_randomCombinations = -1;
            PostHocPower_min_percentage_of_datasets = -1;
            PostHocPower_max_number_of_datasets = -1;
            DEGs_subdirectory_name_ending = "_seuratDEGs";
            RandomSeedNo_degs_subdirectory = "RandomSeedNos/";
            PostHocPower_subdirectory = "PostHocPowerAnalysis/";
            DEGsNoCutoff_subdirectory = "NoCutoffs/";
            DEGsAllMarker_fileName_beginning = "AllMarkers_";
            PostHocPower_coefOfVar_fileName_beginning = "PHP_coefOfVar_noOfDrawnDatasets_";
            PostHocPower_correlation_fileName_beginning = "PHP_correlation_noOfDrawnDatasets_vs_maxDrawnDatasets_";
            PostHocPower_eachClusterGene_fileName_beginning = "PHP_clusterGene_inDependenceOfDatasetNos_";
            RanSeedNo_label = "_ranSeedNo";
            Deg_maxAdjPvalue = 0.05F;
            Deg_maxRank = 300;
            Deg_maxAdjPvalue_factor_realToPostHocPower = 4;
            Deg_maxRank_factor_realToPostHocPower = 3;
            Use_data_normalized_towards_values_of_highest_datasetNos_for_postHocPower = true;
            Max_number_deg_lists_in_memory = 99999999;
        }

        public string Get_fileName_for_given_oneBasedRRandomSeedNo(int oneBasedRRandomSeed_no)
        {
            return "AllMarkers_" + Assay_cellTypeGroup_subtype_level + RanSeedNo_label + oneBasedRRandomSeed_no + ".txt";
        }
        public void Write_options(string directory, string fileName, int randomSeedNode_no)
        {
            string complete_fileName = directory + fileName;
            char delimiter = Global_class.Tab;
            StreamWriter writer = new StreamWriter(complete_fileName);


            writer.WriteLine("Overall_directory" + delimiter + Overall_directory);
            writer.WriteLine("Assay_cellTypeGroup_subtype_level" + delimiter + Assay_cellTypeGroup_subtype_level);
            writer.WriteLine("PostHocPower_max_randomCombinations" + delimiter + PostHocPower_max_randomCombinations);
            writer.WriteLine("PostHocPower_min_percentage_of_datasets" + delimiter + PostHocPower_min_percentage_of_datasets);
            writer.WriteLine("PostHocPower_max_number_of_datasets" + delimiter + PostHocPower_max_number_of_datasets);
            writer.WriteLine("Number_of_datasets_final_data" + delimiter + Number_of_datasets_final_data);
            writer.WriteLine("Deg_maxAdjPvalue" + delimiter + Deg_maxAdjPvalue);
            writer.WriteLine("Consider_only_degs_significant_in_final_for_postHocPower" + delimiter + Consider_only_degs_significant_in_final_for_postHocPower);
            writer.WriteLine("Deg_maxAdjPvalue_for_consideration_in_postHocPower" + delimiter + Deg_maxAdjPvalue_for_consideration_in_postHocPower);
            writer.WriteLine("Deg_maxRank" + delimiter + Deg_maxRank);
            writer.WriteLine("Deg_maxRank_for_consideration_in_postHocPower" + delimiter + Deg_maxRank_for_consideration_in_postHocPower);
            writer.WriteLine("Report_directory" + delimiter + Report_directory);
            writer.WriteLine("Warning_directory" + delimiter + Warning_directory);
            writer.WriteLine("RandomSeedNo_degs_subdirectory" + delimiter + RandomSeedNo_degs_subdirectory);
            writer.WriteLine("DEGs_subdirectory_name_ending" + delimiter + DEGs_subdirectory_name_ending);
            writer.WriteLine("DEGs_directory" + delimiter + DEGs_directory);
            writer.WriteLine("DEGsNoCutoff_subdirectory" + delimiter + DEGsNoCutoff_subdirectory);
            writer.WriteLine("DEGsNoCutoff_directory" + delimiter + DEGsNoCutoff_directory);
            writer.WriteLine("DEGsNoCutoff_fileName" + delimiter + DEGsNoCutoff_fileName);
            writer.WriteLine("DEGsAllMarker_fileName_beginning" + delimiter + DEGsAllMarker_fileName_beginning);
            writer.WriteLine("RandomSeedNo_degs_directory" + delimiter + RandomSeedNo_degs_directory);
            writer.WriteLine("PostHocPower_subdirectory" + delimiter + PostHocPower_subdirectory);
            writer.WriteLine("PostHocPower_directory" + delimiter + PostHocPower_directory);
            writer.WriteLine("PostHocPower_coefOfVar_fileName" + delimiter + PostHocPower_coefOfVar_fileName);
            writer.WriteLine("PostHocPower_correlation_fileName" + delimiter + PostHocPower_correlation_fileName);
            writer.WriteLine("PostHocPower_eachClusterGene_fileName" + delimiter + PostHocPower_eachClusterGene_fileName);
            writer.WriteLine("Use_data_normalized_towards_values_of_highest_datasetNos_for_postHocPower" + delimiter + Use_data_normalized_towards_values_of_highest_datasetNos_for_postHocPower);
            writer.WriteLine("Consider_only_genes_that_were_significant_at_least_once_for_coefOfVar_and_pearsonCor" + delimiter + Consider_only_genes_that_were_significant_at_least_once_for_coefOfVar_and_pearsonCor);
            writer.WriteLine("Max_number_deg_lists_in_memory" + delimiter + Max_number_deg_lists_in_memory);
            if (randomSeedNode_no>0)
            {
                writer.WriteLine("Random_seed_node_number" + delimiter + randomSeedNode_no);
            }
            else
            {
                writer.WriteLine("Random_seed_node_number" + delimiter + "None (" + randomSeedNode_no + ")");
            }
            writer.Close();
        }
    }
    class Seurat_data_processing_class
    {
        public SCSN_deg_class[] SCSN_deg_instances { get; set; }
        public int[] Filled_scsn_deg_instances { get; set; }
        public Seurat_data_processing_options_class Options { get; set; }
        public SCSN_avgAnalysesWithinEachIteration_deg_class Final_scsn_avg_degs { get; set; }
        public Random_class Random_instance_draw_deg_lists { get; set; }
        public Random_class Random_instance_delete_deg_lists { get; set; }

        public Seurat_data_processing_class(string overall_directory, string assayCellTypeGrouplevel, int randomSeedNodes_minus1EqualsNone)
        {
            SCSN_deg_instances = new SCSN_deg_class[0];
            Options = new Seurat_data_processing_options_class(overall_directory, assayCellTypeGrouplevel);
            Final_scsn_avg_degs = new SCSN_avgAnalysesWithinEachIteration_deg_class(-1, -1);
            Random_instance_draw_deg_lists = new Random_class(randomSeedNodes_minus1EqualsNone);
            Random_instance_delete_deg_lists = new Random_class(-1);
        }
        private void Move_empty_files_to_subfolder()
        {
            string complete_deg_seedNo_directory = Options.RandomSeedNo_degs_directory;
            string empty_fileNames_directory = Options.RandomSeedNo_degs_directory + "Files_with_no_DEGs/";
            string[] complete_fileNames = Directory.GetFiles(complete_deg_seedNo_directory);
            string complete_fileName;
            string fileName;
            string complete_move_fileName;
            int fileNames_length = complete_fileNames.Length;
            SCSN_deg_instances = new SCSN_deg_class[fileNames_length];
            int empty_files = 0;
            for (int indexFN = 0; indexFN < fileNames_length; indexFN++)
            {
                complete_fileName = complete_fileNames[indexFN];
                StreamReader reader = new StreamReader(complete_fileName);
                string firstLine = reader.ReadLine();
                reader.Close();
                if (String.IsNullOrEmpty(firstLine))
                {
                    fileName = System.IO.Path.GetFileName(complete_fileName);
                    ReadWriteClass.Create_directory_if_it_does_not_exist(empty_fileNames_directory);
                    complete_move_fileName = empty_fileNames_directory + fileName;
                    System.IO.File.Move(complete_fileName, complete_move_fileName);
                    empty_files++;
                }
            }
            if (empty_files > 0)
            {
                string text = empty_files + " of " + fileNames_length + " were empty and have been moved to " + empty_fileNames_directory;
                Add_waring_to_report_and_warning_file(text);
            }
        }
        private void Add_headline_to_empty_files()
        {
            string complete_deg_seedNo_directory = Options.RandomSeedNo_degs_directory;
            string empty_fileNames_directory = Options.RandomSeedNo_degs_directory + "Files_with_no_DEGs/";
            string[] complete_fileNames = Directory.GetFiles(complete_deg_seedNo_directory);
            string complete_fileName;
            string fileName;
            string complete_move_fileName;
            int fileNames_length = complete_fileNames.Length;
            SCSN_deg_instances = new SCSN_deg_class[fileNames_length];
            int empty_files = 0;
            SCSN_deg_readWrite_options_class readWriteOptions = new SCSN_deg_readWrite_options_class("","");

            Add_to_report_file("Counting how many of " + fileNames_length + " files contain no data, and, if headers are missing, inserting the required headers into each empty file.");
            for (int indexFN = 0; indexFN < fileNames_length; indexFN++)
            {
                complete_fileName = complete_fileNames[indexFN];
                StreamReader reader = new StreamReader(complete_fileName);
                string firstLine = reader.ReadLine();
                if (String.IsNullOrEmpty(firstLine))
                {
                    reader.Close();
                    StreamWriter writer = new StreamWriter(complete_fileName);
                    StringBuilder sb = new StringBuilder();
                    foreach (string columnName in readWriteOptions.Key_columnNames)
                    {
                        if (sb.Length > 0) { sb.AppendFormat(Global_class.Tab.ToString()); }
                        sb.Append(columnName);
                    }
                    writer.WriteLine(sb.ToString());
                    writer.Close();
                    empty_files++;
                }
                else
                {
                    firstLine = reader.ReadLine();
                    reader.Close();
                    if (String.IsNullOrEmpty (firstLine)) { empty_files++; }
                }
            }
            if (empty_files > 0)
            {
                string text = empty_files + " of " + fileNames_length + " were empty. Headlines added to those files.";
                Add_waring_to_report_and_warning_file(text);
            }
            else
            {
                Add_to_report_file("All files contain data.");
            }
        }
        public void Write_options(int randomSeedNodeNo)
        {
            Options.Write_options(Options.Report_directory, "Selected_parameter_for_average_degs_and_do_postHocPowerAnalysis.txt", randomSeedNodeNo);
        }

        private SCSN_deg_class Generate_and_prepare_scsn_deg_instance(string directory, string fileName)
        {
            if (!System.IO.Directory.Exists(directory)) { throw new Exception(directory + " does not exist."); }
            if (!System.IO.File.Exists(directory + fileName)) { throw new Exception(fileName + " does not exist. Are the seurat random seed nos continuous and start at 1?"); }
            SCSN_deg_class add_scsn = new SCSN_deg_class();
            add_scsn.Generate_by_reading(directory, fileName);
            add_scsn.Remove_genes_above_max_adj_pvalue(Options.Deg_maxAdjPvalue_for_consideration_in_postHocPower);
            add_scsn.Remove_genes_above_max_significance_rank(Options.Deg_maxRank_for_consideration_in_postHocPower);
            add_scsn.Label_significant_genes(Options.Deg_maxAdjPvalue, Options.Deg_maxRank);
            return add_scsn;
        }
        public void Generate_for_all_runs()
        {
            Report_class.Delete_report_file_if_exists(Options.Report_directory, Options.Assay_cellTypeGroup_subtype_level);
            Report_class.Delete_warning_file_if_exists(Options.Report_directory, Options.Assay_cellTypeGroup_subtype_level);
            Report_class.Delete_warning_file_if_exists(Options.Warning_directory, Options.Assay_cellTypeGroup_subtype_level);

            //Move_empty_files_to_subfolder();
            Add_headline_to_empty_files();
            Filled_scsn_deg_instances = new int[0];

            string complete_deg_seedNo_directory = Options.RandomSeedNo_degs_directory;
            string[] complete_fileNames = Directory.GetFiles(complete_deg_seedNo_directory);
            //string fileName;
            int fileNames_length = complete_fileNames.Length;

            if (fileNames_length>Options.Max_number_deg_lists_in_memory)
            {
                Add_to_report_file("Up to " + Options.Max_number_deg_lists_in_memory + " of " + fileNames_length + " DEG files in " + Options.RandomSeedNo_degs_subdirectory + " will be uploaded during analysis. Afterwards the script will replace uploaded data by data that needs to be uploaded for next analysis steps.");
            }
            else
            {
                Add_to_report_file("All " + fileNames_length + " in " + Options.RandomSeedNo_degs_subdirectory + " will be uploaded during analysis.");
            }


            //string fileName;
            //SCSN_deg_instances = new SCSN_deg_class[fileNames_length];
            //SCSN_deg_class add_scsn;
            //int report_index = 500;
            //int files_with_no_DEGs_count = 0;
            //int store_files_length = Math.Min(Options.Max_number_deg_lists_in_memory, fileNames_length);
            //Add_to_report_file("Reading " + store_files_length + " of " + fileNames_length + " deg files.");
            //int oneBasedRIndex;
            //if (store_files_length<fileNames_length)
            //{
            //    Filled_scsn_deg_instances = new int[store_files_length];
            //}
            //for (int indexFN=0;indexFN < store_files_length; indexFN++)
            //{
            //    oneBasedRIndex = indexFN + 1;
            //    fileName = Options.Get_fileName_for_given_oneBasedRRandomSeedNo(oneBasedRIndex);
            //    add_scsn = Generate_and_prepare_scsn_deg_instance(complete_deg_seedNo_directory, fileName);
            //    if (add_scsn.DEGs.Length == 0) { files_with_no_DEGs_count++; }
            //    try
            //    {
            //        this.SCSN_deg_instances[indexFN] = add_scsn;
            //        if (store_files_length < fileNames_length)
            //        {
            //            Filled_scsn_deg_instances[indexFN] = indexFN;
            //        }
            //    }
            //    catch
            //    {
            //        throw new Exception("Next single cell/nucleus dataset could not be added to array any more. Decrease Max_number_deg_lists_in_memory");
            //    }
            //    if ((indexFN+1) % report_index == 0)
            //    {
            //        Add_to_report_file((indexFN + 1) + " files read.");
            //    }
            //}
            //if (files_with_no_DEGs_count > 0) { Add_waring_to_report_and_warning_file(files_with_no_DEGs_count + " read files of " + fileNames_length + " had no DEGs."); }
            //if (Options.Remove_empty_deg_lists)
            //{
            //    if (store_files_length<fileNames_length) { throw new Exception("Since Remove_empty_deg_lists is set to true, Max_number_deg_lists_in_memory has to be set to at least " + fileNames_length + ", i.e. the number of files in " + Options.RandomSeedNo_degs_directory + "."); }
            //    List<SCSN_deg_class> keep = new List<SCSN_deg_class>();
            //    SCSN_deg_class scsn_deg_instance;
            //    int scsn_length = this.SCSN_deg_instances.Length;
            //    for (int indexSCSN=0; indexSCSN<scsn_length; indexSCSN++)
            //    {
            //        scsn_deg_instance = this.SCSN_deg_instances [indexSCSN];
            //        if (scsn_deg_instance.DEGs.Length > 0)
            //        {
            //            keep.Add(scsn_deg_instance);
            //        }
            //    }
            //    this.SCSN_deg_instances = keep.ToArray();
            //    if (this.SCSN_deg_instances.Length!=scsn_length)
            //    {
            //        Add_waring_to_report_and_warning_file("DEG instances with no lines removed. " + this.SCSN_deg_instances.Length + " of " + scsn_length + " instances kept.");
            //    }
            //    else
            //    {
            //        Add_to_report_file("No DEG instances had zero lines and were removed.");
            //    }
            //}
        }
        private void Add_waring_to_report_and_warning_file(string warning_text)
        {
            Add_to_warning_file(warning_text);
            Add_to_report_file("Warning");
            Add_to_report_file(warning_text);
        }
        private void Add_to_report_file(string report_file_text)
        {
            Report_class.Add_to_report_file(Options.Report_directory, Options.Assay_cellTypeGroup_subtype_level, report_file_text);
        }
        private void Add_to_warning_file(string report_file_text)
        {
            Report_class.Add_to_warning_file(Options.Report_directory, Options.Assay_cellTypeGroup_subtype_level, report_file_text);
            Report_class.Add_to_warning_file(Options.Warning_directory, Options.Assay_cellTypeGroup_subtype_level, report_file_text);
        }
        private SCSN_avgAnalysesWithinEachIteration_deg_class Calculate_avg_DEGs_using_given_random_seedNo_results(int[] keep_indices, bool calculate_sds)
        {
            int keep_indices_length = keep_indices.Length;
            if (keep_indices_length > 0)
            {
            }
            else
            {
                keep_indices_length = SCSN_deg_instances.Length;
                keep_indices = new int[keep_indices_length];
                for (int indexKeep=0; indexKeep<keep_indices_length; indexKeep++)
                {
                    keep_indices[indexKeep] = indexKeep;
                }
                Add_to_report_file("Calculating average degs from " + keep_indices_length + " DEG files with different randomSeedNos");
            }
            int max_index = -1;
            foreach (int keep_index in keep_indices)
            {
                if (keep_index > max_index) {  max_index = keep_index; }
            }
            if (max_index >= this.SCSN_deg_instances.Length) { throw new Exception("Calcuate avg DEGs: Only " + SCSN_deg_instances.Length + " DEG lists are available. Max DEG list index (" + max_index + ") is too high. Adopt --max_numberOfDatasets or generate more DEG lists with additional random seed numbers."); }


            
            SCSN_deg_class scsn_degs;
            SCSN_avgAnalysesWithinEachIteration_deg_class scsn_avg_degs = new SCSN_avgAnalysesWithinEachIteration_deg_class(keep_indices_length, 0);
            int indexKeepIndex;
            int delete_index;
            int deleteFilledScSn_index;
            int oneBasedRIndex;
            string directory = Options.RandomSeedNo_degs_directory;
            string fileName;
            for (int indexIndex = 0; indexIndex < keep_indices_length; indexIndex++)
            {
                indexKeepIndex = keep_indices[indexIndex];
                scsn_degs = SCSN_deg_instances[indexKeepIndex];
                if (scsn_degs == null)
                {
                    if (Filled_scsn_deg_instances.Length == Options.Max_number_deg_lists_in_memory)
                    {
                        deleteFilledScSn_index = Random_instance_delete_deg_lists.Randomly_draw_non_overlapping_indices_from_number_of_available_indices(1, Filled_scsn_deg_instances.Length)[0];
                        delete_index = Filled_scsn_deg_instances[deleteFilledScSn_index];
                        if (SCSN_deg_instances[delete_index] == null) { throw new Exception(); }
                        SCSN_deg_instances[delete_index] = null;
                        Filled_scsn_deg_instances[deleteFilledScSn_index] = indexKeepIndex;
                    }
                    else
                    {
                        int filled_scsn_deg_instances_length = Filled_scsn_deg_instances.Length;
                        int[] new_filled_scsn_deg_instances = new int[filled_scsn_deg_instances_length + 1];
                        for (int indexFilled=0;indexFilled<filled_scsn_deg_instances_length;indexFilled++)
                        {
                            new_filled_scsn_deg_instances[indexFilled] = Filled_scsn_deg_instances[indexFilled];
                        }
                        new_filled_scsn_deg_instances[filled_scsn_deg_instances_length] = indexKeepIndex;
                        Filled_scsn_deg_instances = new_filled_scsn_deg_instances;
                        if (Filled_scsn_deg_instances.Length== Options.Max_number_deg_lists_in_memory)
                        {
                            Add_to_report_file("Maximum amount of files read and stored in memory (" + Filled_scsn_deg_instances.Length + " of " + SCSN_deg_instances.Length + " files). From now on, newly uploaded data will replace existing data in memory.");
                        }
                    }
                    if (!Filled_scsn_deg_instances.Distinct().ToArray().Length.Equals(Filled_scsn_deg_instances.Length)) { throw new Exception(); }
                    if (SCSN_deg_instances[indexKeepIndex] != null) { throw new Exception(); }
                    oneBasedRIndex = indexKeepIndex + 1;
                    fileName = Options.Get_fileName_for_given_oneBasedRRandomSeedNo(oneBasedRIndex);
                    SCSN_deg_instances[indexKeepIndex] = Generate_and_prepare_scsn_deg_instance(directory, fileName);
                    scsn_degs = SCSN_deg_instances[indexKeepIndex];
                }
                scsn_avg_degs.Add_other_dataset_by_summing_up_running_mean(scsn_degs);
            }
            scsn_avg_degs.Remove_empty_array_lines();
            if (calculate_sds)
            {
                for (int indexIndex = 0; indexIndex < keep_indices_length; indexIndex++)
                {
                    indexKeepIndex = keep_indices[indexIndex];
                    scsn_degs = SCSN_deg_instances[indexKeepIndex];
                    scsn_avg_degs.Add_other_dataset_by_summing_up_running_variances_after_mean_is_final(scsn_degs);
                }
                scsn_avg_degs.Calculate_sampleSD_and_finish_averaging_procedure();
            }
            else
            {
                scsn_avg_degs.Set_averaging_procedure_finished_to_true_if_no_sampleSD_calculated();
            }
            return scsn_avg_degs;
        }
        public void Calculate_and_write_avgDEGs_using_first_random_seedNo_results_of_given_number_of_datasets()
        {
            int[] dataset_indices = Enumerable.Range(0, Options.Number_of_datasets_final_data).ToArray();
            Final_scsn_avg_degs = Calculate_avg_DEGs_using_given_random_seedNo_results(dataset_indices, true);
            Final_scsn_avg_degs.Calculate_pValAdjGeomMean_from_minusLog10AdjPvalueArithMean();
            Final_scsn_avg_degs.Calculate_fractional_ranks_using_minusLog10AdjPvalueArithMean_and_avgLog2FcArithMean();
            Final_scsn_avg_degs.Write(Options.DEGsNoCutoff_directory, Options.DEGsNoCutoff_fileName);
        }

        public void Do_postHocPower_analysis()
        {
            ReadWriteClass.Delete_directory_recursive_if_it_exists(Options.PostHocPower_directory);

            int percentage_of_minimum_degLists = Options.PostHocPower_min_percentage_of_datasets;
            int number_of_available_datasets = this.SCSN_deg_instances.Length;
            int max_number_of_datasets = Options.PostHocPower_max_number_of_datasets;
            int max_randomizations_count = Options.PostHocPower_max_randomCombinations;
            int ideal_minimum_number_of_available_datasets = max_number_of_datasets * max_randomizations_count;
            int min_drawn_datasets = (int)Math.Max(1, Math.Round((float)percentage_of_minimum_degLists / 100F * (float)max_number_of_datasets));
            int total_drawn_datasets = 0;
            int number_of_datasets_for_real_data = Options.Number_of_datasets_final_data;
            for (int number_of_drawn_datasets = min_drawn_datasets; number_of_drawn_datasets<=max_number_of_datasets;number_of_drawn_datasets++)
            {
                total_drawn_datasets += number_of_drawn_datasets * max_randomizations_count;
            }


            if (max_number_of_datasets > number_of_available_datasets) { throw new Exception("--max_numberOfDatasets " + max_number_of_datasets + " was selected, but only " + number_of_available_datasets + " DEG dataset lists are available. To reduce amount of overlapping datasets between the analyses, estimated at least " + ideal_minimum_number_of_available_datasets + " DEG dataset lists should be available."); }

            Add_to_report_file("PostHoc Power Analysis: Drawing " + percentage_of_minimum_degLists + " to 100% of " + max_number_of_datasets + " from " + number_of_available_datasets + " existing DEG lists generated using different random seed numbers.");
            Add_to_report_file("Drawing for each iteration at max " + max_randomizations_count + " different combinations.");
            if (number_of_available_datasets < ideal_minimum_number_of_available_datasets)
            {
                string[] warning_texts = new string[] { number_of_available_datasets + " DEG dataset lists are available for Post hoc power analysis.",
                    "The pipeline draws " + max_randomizations_count + " different sets of a given number of DEG lists (i.e. from " + min_drawn_datasets,
                    "to " + max_number_of_datasets + ").",
                    "Each time, it first averages the DEG lists within each set, generating " + max_randomizations_count + " averaged DEG lists for each",
                    "given number of DEG lists.",
                    "One analysis then subjects the " + max_randomizations_count + " averaged DEG sets, derived from an equal number of initial DEG lists,",
                    "to a second averaging step. From the results of the second averaging step, it will calculate the",
                    "coefficient of variation to quantify if the initial number of drawn DEG lists is sufficient for",
                    "consistent results.",
                    "A second analysis will correlate the results of each first averaging with the final DEG list that",
                    "will be used for all downstream analyses. The final DEG list is generated by drawing and averaging",
                    number_of_datasets_for_real_data + " DEG lists. This saturation analysis will investigate if drawing less than " + number_of_datasets_for_real_data,
                    "initial DEG lists is already sufficient for consistent results, so that it can be concluded that",
                    "drawing more than " + number_of_datasets_for_real_data + " will not lead to more accuracy.",
                    "",
                    "Please consider that there are only " + number_of_available_datasets + " DEG datasets available, so that the same datasets might be",
                    "repetitively drawn, potentially biasing the results in favor of a too good outcome of the analysis.",
                    "As an intuitive/arbitrary estimate, there should be at least " + ideal_minimum_number_of_available_datasets + " DEG lists available to reduce such effects." };
                Add_to_warning_file(Options.Assay_cellTypeGroup_subtype_level + ":");
                Add_to_report_file("");
                Add_to_report_file("Warning - Begin");
                foreach (string warning_text in warning_texts)
                {
                    Add_to_report_file(warning_text);
                    Add_to_warning_file(warning_text);
                }
                Add_to_report_file("Warning - End");
                Add_to_report_file("");
                Add_to_warning_file("");
            }
            SCSN_avgWithinEachNumberOfDatasets_deg_class combined_avg_degs_eachNumberOfDatasets = new SCSN_avgWithinEachNumberOfDatasets_deg_class(Options.Consider_only_degs_significant_in_final_for_postHocPower, Options.Deg_maxAdjPvalue, Options.Deg_maxRank);
            SCSN_avgAnalysesWithinEachIteration_deg_class combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets;
            SCSN_avgAnalysesWithinEachIteration_deg_class current_avg_degs_eachIteration;
            SCSN_postHocPower_avg_class postHocPower = new SCSN_postHocPower_avg_class(this.Options.Consider_only_genes_that_were_significant_at_least_once_for_coefOfVar_and_pearsonCor);
            SCSN_postHocPower_correlation_class postHocPower_correlation = new SCSN_postHocPower_correlation_class(this.Options.Consider_only_genes_that_were_significant_at_least_once_for_coefOfVar_and_pearsonCor);
            int report_every_randoNo = 50;
            int current_max_randomizations_count;
            Math_class mfi = new Math_class();
            long binominal_coefficient;
            int[][] keep_indices_array;
            int[] keep_indices;

            for (int number_of_drawn_datasets = min_drawn_datasets; number_of_drawn_datasets <= max_number_of_datasets; number_of_drawn_datasets++)
            {
                combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets = new SCSN_avgAnalysesWithinEachIteration_deg_class(number_of_drawn_datasets,(int)1E7);
                current_max_randomizations_count = max_randomizations_count;
                binominal_coefficient = mfi.Get_binomial_coefficient(number_of_available_datasets, number_of_drawn_datasets);
                if ((binominal_coefficient>0)&&(current_max_randomizations_count > binominal_coefficient)) { current_max_randomizations_count = (int)binominal_coefficient; }

                HashSet<string> seen_sets = new HashSet<string>();
                keep_indices_array = new int[current_max_randomizations_count][];
                int randoNo = 0;
                while (randoNo < current_max_randomizations_count)
                {
                    keep_indices = Random_instance_draw_deg_lists.Randomly_draw_non_overlapping_indices_from_number_of_available_indices(number_of_drawn_datasets, number_of_available_datasets);
                    keep_indices = keep_indices.OrderBy(l => l).ToArray();
                    string key = string.Join(",", keep_indices);
                    if (!seen_sets.Contains(key))
                    {
                        seen_sets.Add(key);
                        keep_indices_array[randoNo] = keep_indices;
                        randoNo++;
                    }
                }

                if (number_of_drawn_datasets == number_of_available_datasets) { current_max_randomizations_count = 1; }
                Add_to_report_file("Drawing " + max_randomizations_count + " different combinations of " + number_of_drawn_datasets + " files from " + number_of_available_datasets + " potential datasets and averaging the results (first average).");
                for (randoNo = 0; randoNo < current_max_randomizations_count; randoNo++)
                {
                    keep_indices = keep_indices_array[randoNo];
                    bool calculate_sampleSDs = false;
                    if (number_of_drawn_datasets == min_drawn_datasets) { calculate_sampleSDs = true; }
                    current_avg_degs_eachIteration = Calculate_avg_DEGs_using_given_random_seedNo_results(keep_indices, calculate_sampleSDs);
                    current_avg_degs_eachIteration.Add_iteration_to_all_degs(randoNo);
                    combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets.Add_other_using_indexNextAvailableThis(current_avg_degs_eachIteration);
                    if ((randoNo+1) % report_every_randoNo == 0)
                    {
                        Add_to_report_file("Finished calculations for " + (randoNo + 1) + " drawn combinations.");
                    }
                }
                combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets.Remove_empty_array_lines();
                if (number_of_drawn_datasets == min_drawn_datasets)
                {
                    combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets.Write(Options.PostHocPower_directory, "DEGs_of_" + current_max_randomizations_count + "_iterations_with_" + number_of_drawn_datasets + "_drawn_datasets.txt");
                }
                Add_to_report_file("Averaging the avaraged results (second average).");
                combined_avg_degs_eachNumberOfDatasets.Calculate_mean_and_sampleSD_of_the_artihMeans_for_allValueTypesOfInterest_across_all_iterations_assuming_same_number_of_drawn_datasets(combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets, Final_scsn_avg_degs);
                //one average value for each cluster/gene/number_of_drawn_datasets
                Add_to_report_file("Calculating pearson correlations using the results after the first averaging step.");
                postHocPower_correlation.Calculate_pearson_correlation_for_valueTypeOfInterest_and_add_to_array(combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets, Final_scsn_avg_degs, AvgValueType_of_interest_enum.Minus_log10_adj_pvalue);
                postHocPower_correlation.Calculate_pearson_correlation_for_valueTypeOfInterest_and_add_to_array(combined_avg_degs_eachIteration_forSameNoOfDrawnDatasets, Final_scsn_avg_degs, AvgValueType_of_interest_enum.Is_sig);
            }
            //combined_avg_degs_eachNumberOfDatasets.Calculate_mean_and_sampleSD_of_the_artihMeans_for_allValueTypesOfInterest_across_all_iterations_assuming_same_number_of_drawn_datasets(Final_scsn_avg_degs);
            //postHocPower_correlation.Calculate_pearson_correlation_for_valueTypeOfInterest_and_add_to_array(Final_scsn_avg_degs, Final_scsn_avg_degs, AvgValueType_of_interest_enum.Minus_log10_adj_pvalue);
            //postHocPower_correlation.Calculate_pearson_correlation_for_valueTypeOfInterest_and_add_to_array(Final_scsn_avg_degs, Final_scsn_avg_degs, AvgValueType_of_interest_enum.Is_sig);

            //combined_avg_degs_eachNumberOfDatasets.Normalize_values_with_reference_to_values_of_highest_number_of_datasets();
            //one average value for each cluster/gene/number_of_drawn_datasets, normalized towards that value with the highest number of added datasets

            combined_avg_degs_eachNumberOfDatasets.Write(Options.PostHocPower_directory, Options.PostHocPower_eachClusterGene_fileName);
            //combined_avg_degs_eachNumberOfDatasets.Read(Options.PostHocPower_directory, Options.PostHocPower_eachClusterGene_fileName);

            postHocPower.Generate_from_avgWithinEachNumberOfDatasets_instance_and_add_to_array(combined_avg_degs_eachNumberOfDatasets);
            //one average value for each cluster/number_of_drawn_datasets
            postHocPower.Write_data(Options.PostHocPower_directory, Options.PostHocPower_coefOfVar_fileName);

            postHocPower_correlation.Write(Options.PostHocPower_directory, Options.PostHocPower_correlation_fileName);

            Add_to_report_file("PostHocPower Analysis finished.");
        }




    }
}
