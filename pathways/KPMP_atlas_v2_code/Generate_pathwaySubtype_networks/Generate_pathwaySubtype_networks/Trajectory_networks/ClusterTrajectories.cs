using Common_functions;
using ReadWrite;
using Network;
using Network_visualization;
using Enrichment_results;

namespace Trajectories
{
    interface IKPMP_cluster_name_conversion_line
    {
        string Cluster_name { get; set; }
        int Clustering_level { get; set; }
    }
    class KPMP_shared_conversion_class
    {
        private IKPMP_cluster_name_conversion_line[] Replace_hyphen_by_underscore(IKPMP_cluster_name_conversion_line[] clusterName_lines)
        {
            int clusterName_lines_length = clusterName_lines.Length;
            IKPMP_cluster_name_conversion_line clusterName_line;
            for (int indexCN = 0; indexCN < clusterName_lines_length; indexCN++)
            {
                clusterName_line = clusterName_lines[indexCN];
                clusterName_line.Cluster_name = clusterName_line.Cluster_name.Replace("-", "_");
            }
            return clusterName_lines;
        }
        private IKPMP_cluster_name_conversion_line[] Add_clustering_level_in_front_of_clusterName(IKPMP_cluster_name_conversion_line[] clusterName_lines)
        {
            int clusterName_lines_length = clusterName_lines.Length;
            IKPMP_cluster_name_conversion_line clusterName_line;
            for (int indexCN = 0; indexCN < clusterName_lines_length; indexCN++)
            {
                clusterName_line = clusterName_lines[indexCN];
                if (clusterName_line.Clustering_level <= 0) { throw new Exception(); }
                clusterName_line.Cluster_name = "L" + clusterName_line.Clustering_level + ":" + clusterName_line.Cluster_name;
            }
            return clusterName_lines;
        }
        private IKPMP_cluster_name_conversion_line[] Add_cellType_to_clusterName_if_only_number(IKPMP_cluster_name_conversion_line[] clusterName_lines, string cell_type)
        {
            IKPMP_cluster_name_conversion_line clusterName_line;
            int barcodes_length = clusterName_lines.Length;
            int cluster_no;
            for (int indexBarcode = 0; indexBarcode < barcodes_length; indexBarcode++)
            {
                clusterName_line = clusterName_lines[indexBarcode];
                if (int.TryParse(clusterName_line.Cluster_name, out cluster_no))
                {
                    clusterName_line.Cluster_name = cell_type + clusterName_line.Cluster_name;
                }
                //else if (numbers_added) { throw new Exception(); }
            }
            return clusterName_lines;
        }

        public IKPMP_cluster_name_conversion_line[] Adjust_clusterNames_including_level_addition(IKPMP_cluster_name_conversion_line[] clusterName_lines, string cell_type)
        {
            //clusterName_lines = Add_cellType_to_clusterName_if_only_number(clusterName_lines, cell_type);
            clusterName_lines = Add_clustering_level_in_front_of_clusterName(clusterName_lines);
            clusterName_lines = Replace_hyphen_by_underscore(clusterName_lines);
            return clusterName_lines;
        }

        public IKPMP_cluster_name_conversion_line[] Adjust_clusterNames_without_level_addition(IKPMP_cluster_name_conversion_line[] clusterName_lines, string cell_type)
        {
            clusterName_lines = Add_cellType_to_clusterName_if_only_number(clusterName_lines, cell_type);
            clusterName_lines = Replace_hyphen_by_underscore(clusterName_lines);
            return clusterName_lines;
        }
    }


    class Trajectory_line_class
    {
        public string Subtype_name { get; set; }
        public string Assay_cellType { get; set; }
        public int Lineage_no { get; set; }
        public int Position { get; set; }

        public Trajectory_line_class()
        {
            Subtype_name = "";
            Assay_cellType = "";
        }

        public Trajectory_line_class Deep_copy()
        {
            Trajectory_line_class copy = (Trajectory_line_class)this.MemberwiseClone();
            copy.Subtype_name = (string)this.Subtype_name.Clone();
            copy.Assay_cellType = (string)this.Assay_cellType.Clone();
            return copy;
        }
    }

    class Trajectory_class
    {
        public string Cell_type { get; set; }
        public Trajectory_line_class[] Trajectories { get; set; }

        public Trajectory_class()
        {
            Cell_type = "";
            Trajectories = new Trajectory_line_class[0];
        }

        public void Generate_for_assay_cellType(string assay_cellType)
        {
            Dictionary<int, string[]> lineageNo_subtypes_dict = new Dictionary<int, string[]>();
            switch (assay_cellType)
            {
                case "SNRNAseq_PTsl3":
                case "SNboth_PTsl3":
                case "SNonly_PTsl3":
                case "Multiome_PTsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { "aPT2", "aPT1", "aPT-S1/S2", "frPT-S1/S2" });
                    lineageNo_subtypes_dict.Add(2, new string[] { "aPT2", "aPT1", "aPT-S1/S2", "PT-S1", "PT-S2" });
                    lineageNo_subtypes_dict.Add(3, new string[] { "aPT2", "aPT1", "aPT-S3", "PT-S3" });
                    lineageNo_subtypes_dict.Add(4, new string[] { "aPT2", "aPT1", "frPT-S3" });
                    lineageNo_subtypes_dict.Add(5, new string[] { "aPT2", "aPT1", "frPT-S3", "aPT-S3" });
                    break;
                case "SNboth_TALsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { "aTAL1", "aTAL2", "C-TAL-A" });
                    lineageNo_subtypes_dict.Add(2, new string[] { "aTAL1", "aTAL2", "frTAL" });
                    break;
                case "SNboth_FIBsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { "C-FIB", "C-MYOF" });
                    lineageNo_subtypes_dict.Add(2, new string[] { "C-FIB", "C-FIB-OSMRlo", "C-FIB-OSMRhi" });
                    lineageNo_subtypes_dict.Add(3, new string[] { "C-FIB", "C-FIB-PATH" });
                    lineageNo_subtypes_dict.Add(4, new string[] { "pvFIB-RSPO3+", "pvFIB-PI16+", "pvFIB", "pvMYOF" });
                    break;
                case "SNboth_Myeloidsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { "MON", "moMAC-CXCL10+" });
                    lineageNo_subtypes_dict.Add(2, new string[] { "MON", "moMAC-C3+" });
                    lineageNo_subtypes_dict.Add(3, new string[] { "MON", "moMAC-HBEGF+", "moFAM" });
                    lineageNo_subtypes_dict.Add(4, new string[] { "resMAC-LYVE1+", "resMAC-HLAIIhi" });
                    break;
                case "SNboth_ECsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { });
                    break;
                case "SNboth_CDCNTsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { });
                    break;
                case "SNboth_DTLATLsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { });
                    break;
                case "SNboth_PODPECsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { });
                    break;
                case "SNBoth_DCTsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { });
                    break;
                case "SNBoth_Lymphoidsl3":
                    lineageNo_subtypes_dict.Add(1, new string[] { });
                    break;
                default:
                    lineageNo_subtypes_dict.Add(1, new string[] { });
                    break;
            }
            int[] lineageNos = lineageNo_subtypes_dict.Keys.ToArray();
            int lineageNo;
            int lineageNos_length = lineageNos.Length;
            string[] subtypes;
            string subtype;
            int subtypes_length;
            List<Trajectory_line_class> trajectory_lines = new List<Trajectory_line_class>();
            Trajectory_line_class trajectory_line;
            for (int indexLineageNo =0; indexLineageNo< lineageNos_length; indexLineageNo++)
            {
                lineageNo = lineageNos[indexLineageNo];
                subtypes = lineageNo_subtypes_dict[lineageNo];
                subtypes_length = subtypes.Length;
                for (int indexS=0; indexS<subtypes_length; indexS++)
                {
                    subtype = subtypes[indexS];
                    trajectory_line = new Trajectory_line_class();
                    trajectory_line.Lineage_no = lineageNo;
                    trajectory_line.Assay_cellType = (string)assay_cellType.Clone();
                    trajectory_line.Position = indexS;
                    trajectory_line.Subtype_name = (string)subtype.Clone();
                    trajectory_lines.Add(trajectory_line);
                }
            }
            this.Trajectories = trajectory_lines.ToArray();
        }
    }
    class KPMP_cluster_characterization_line_class
    {
        public string Dataset { get; set; }
        public string Cell_type { get; set; }
        public string Cell_subtype_name { get; set; }
        public KPMP_category_for_cluster_enum Category { get; set; }
        public string Category_entity { get; set; }
        public float Category_entity_value { get; set; }

        public KPMP_cluster_characterization_line_class()
        {
            this.Dataset = "";
            this.Cell_type = "";
            this.Category_entity = "";
            this.Cell_subtype_name = "";
        }

        public KPMP_cluster_characterization_line_class Deep_copy()
        {
            KPMP_cluster_characterization_line_class copy = (KPMP_cluster_characterization_line_class)this.MemberwiseClone();
            copy.Dataset = (string)this.Dataset.Clone();
            copy.Cell_type = (string)this.Cell_type.Clone();
            copy.Cell_subtype_name = (string)this.Cell_subtype_name.Clone();
            copy.Category_entity = (string)this.Category_entity.Clone();
            return copy;
        }
    }
    class KPMP_cluster_characterization_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ';'; } }
        public KPMP_cluster_characterization_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Dataset", "Cell_type", "Cluster_level", "Cluster_name", "Original_cluster_name", "Annotated_cluster_name", "Hierarchical_cluster_name", "Category", "Category_entity", "Category_entity_group", "Category_entity_value" };
            this.Key_columnNames = this.Key_propertyNames;
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
        }
    }
    class KPMP_cluster_hierarchy_line_class
    {
        public string Assay_cellType { get; set; }
        public string Source_cellSubtypeName { get; set; }
        public string Target_cellSubtypeName { get; set; }

        public KPMP_cluster_hierarchy_line_class()
        {
            Assay_cellType = "";
            Source_cellSubtypeName = "";
            Target_cellSubtypeName = "";
        }

        public KPMP_cluster_hierarchy_line_class Deep_copy()
        {
            KPMP_cluster_hierarchy_line_class copy = (KPMP_cluster_hierarchy_line_class)this.MemberwiseClone();
            copy.Assay_cellType = (string)this.Assay_cellType.Clone();
            copy.Source_cellSubtypeName = (string)this.Source_cellSubtypeName.Clone();
            copy.Target_cellSubtypeName = (string)this.Target_cellSubtypeName.Clone();
            return copy;
        }
    }
    class KPMP_cluster_hierarchy_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ';'; } }
        public KPMP_cluster_hierarchy_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Cell_type", "Parent_level", "Parent_clusterName", "Child_level", "Child_clusterName", "Parent_cluster_annotatedName", "Child_cluster_annotatedName", "Parent_cluster_hierarchicalName", "Child_cluster_hierarchicalName", "Shared_barcodes_count", "Parent_barcodes_count", "Child_barcodes_count", "Percent_child_barcodes_shared_with_parent", "Jaccard_index", "Is_parent_child_interaction" };
            this.Key_columnNames = this.Key_propertyNames;
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
        }
    }
    class KPMP_cluster_hierarchy_options_class
    {
        public float Minimum_jaccard_index { get; set; }
        public float Minimum_percent_of_child_barcodes_in_parent { get; set; }
        public System.Drawing.Color[] Category_colors { get; set; }
        public System.Drawing.Color MaxMinSize_color { get; set; }
        public Dictionary<string, System.Drawing.Color> CategoryEntity_color_dict { get; set; }
        public Dictionary<string, int> CategoryEntity_pieChartSliceOrder_dict { get; set; }
        public int Forced_reference_category_value { get; set; }


        public KPMP_cluster_hierarchy_options_class()
        {
            this.Minimum_jaccard_index = 0.0F;
            this.Minimum_percent_of_child_barcodes_in_parent = 80;
            this.Category_colors = new System.Drawing.Color[] { System.Drawing.Color.Orange,
                                                            System.Drawing.Color.RoyalBlue,
                                                            System.Drawing.Color.Red,
                                                            System.Drawing.Color.Navy,
                                                            System.Drawing.Color.Orchid,
                                                            System.Drawing.Color.Cyan,
                                                            System.Drawing.Color.Aquamarine,
                                                            System.Drawing.Color.Blue,
                                                            System.Drawing.Color.Plum,
                                                            System.Drawing.Color.DarkGoldenrod,
                                                            System.Drawing.Color.DarkOrchid,
                                                            System.Drawing.Color.Brown,
                                                            System.Drawing.Color.Snow,
                                                            System.Drawing.Color.Yellow,
                                                            System.Drawing.Color.MintCream,
                                                            System.Drawing.Color.Goldenrod,
                                                            System.Drawing.Color.DarkGoldenrod,
                                                            System.Drawing.Color.Lime,
                                                            System.Drawing.Color.Brown,
                                                            System.Drawing.Color.RosyBrown
            };
            Forced_reference_category_value = -1;
            CategoryEntity_color_dict = new Dictionary<string, System.Drawing.Color>();
            MaxMinSize_color = System.Drawing.Color.White;
            CategoryEntity_pieChartSliceOrder_dict = new Dictionary<string, int>();
        }

        public void Add_to_categoryEntity_color_dict(Dictionary<string, System.Drawing.Color> add_categoryEntity_color_dict)
        {
            string[] categories = add_categoryEntity_color_dict.Keys.ToArray();
            foreach (string category in categories)
            {
                if (!CategoryEntity_color_dict.ContainsKey(category)) { CategoryEntity_color_dict.Add(category, add_categoryEntity_color_dict[category]); }
                else if (!CategoryEntity_color_dict[category].Equals(add_categoryEntity_color_dict[category])) { throw new Exception(); }
            }
        }

        public void Add_to_categoryEntity_pieChartSliceOrder_dict(Dictionary<string, int> add_categoryEntity_pieChartSliceOrder_dict)
        {
            string[] categories = add_categoryEntity_pieChartSliceOrder_dict.Keys.ToArray();
            foreach (string category in categories)
            {
                if (!CategoryEntity_pieChartSliceOrder_dict.ContainsKey(category)) { CategoryEntity_pieChartSliceOrder_dict.Add(category, add_categoryEntity_pieChartSliceOrder_dict[category]); }
                else if (!CategoryEntity_pieChartSliceOrder_dict[category].Equals(add_categoryEntity_pieChartSliceOrder_dict[category])) { throw new Exception(); }
            }
        }

    }
    class KPMP_cluster_hierarchy_class
    {
        public KPMP_cluster_hierarchy_line_class[] Hierarchy { get; set; }
        public KPMP_cluster_characterization_line_class[] Cluster_characterizations { get; set; }
        public KPMP_cluster_hierarchy_options_class Options { get; set; }
        public bool Is_generationOfHierarchyFromBarcodes { get; set; }
        private string Cell_type { get; set; }

        public KPMP_cluster_hierarchy_class(string cell_type)
        {
            Cell_type = (string)cell_type.Clone();
            Is_generationOfHierarchyFromBarcodes = false;
            this.Options = new KPMP_cluster_hierarchy_options_class();
            this.Cluster_characterizations = new KPMP_cluster_characterization_line_class[0];
            this.Hierarchy = new KPMP_cluster_hierarchy_line_class[0];
        }

        private void Add_to_hierarchy(KPMP_cluster_hierarchy_line_class[] add_hierarchy)
        {
            int this_length = this.Hierarchy.Length;
            int add_length = add_hierarchy.Length;
            int new_length = this_length + add_length;
            KPMP_cluster_hierarchy_line_class[] new_hierarchy = new KPMP_cluster_hierarchy_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_hierarchy[indexNew] = this.Hierarchy[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_hierarchy[indexNew] = add_hierarchy[indexAdd];
            }
            this.Hierarchy = new_hierarchy;
        }

        private void Add_to_cluster_characterizations(KPMP_cluster_characterization_line_class[] add_cluster_characterizations)
        {
            int this_length = this.Cluster_characterizations.Length;
            int add_length = add_cluster_characterizations.Length;
            int new_length = this_length + add_length;
            KPMP_cluster_characterization_line_class[] new_cluster_characterizations = new KPMP_cluster_characterization_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_cluster_characterizations[indexNew] = this.Cluster_characterizations[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_cluster_characterizations[indexNew] = add_cluster_characterizations[indexAdd];
            }
            this.Cluster_characterizations = new_cluster_characterizations;
        }

        #region Check
        private void Check_that_category_does_not_exist(KPMP_category_for_cluster_enum category)
        {
            foreach (KPMP_cluster_characterization_line_class cluster_characterization_line in this.Cluster_characterizations)
            {
                if (cluster_characterization_line.Category.Equals(category))
                {
                    throw new Exception();
                }
            }
        }
        public void Check_for_correctness()
        {
            Dictionary<string, bool> existingCellType_dict = new Dictionary<string, bool>();
            foreach (KPMP_cluster_characterization_line_class char_line in this.Cluster_characterizations)
            {
                if (!existingCellType_dict.ContainsKey(char_line.Cell_subtype_name))
                {
                    existingCellType_dict.Add(char_line.Cell_subtype_name, false);
                }
            }
            foreach (KPMP_cluster_hierarchy_line_class hierarchy_line in this.Hierarchy)
            {
                if (!existingCellType_dict.ContainsKey(hierarchy_line.Source_cellSubtypeName))
                {
                    throw new Exception();
                }
                if (!existingCellType_dict.ContainsKey(hierarchy_line.Target_cellSubtypeName))
                {
                    throw new Exception();
                }
            }
        }
        #endregion

        #region Add to arrays
        private void Add_to_cluster_characterizations_after_checking_that_cluster_category_does_not_exist(KPMP_cluster_characterization_line_class[] add_cluster_characterizations)
        {
            int this_length = this.Cluster_characterizations.Length;
            int add_length = add_cluster_characterizations.Length;
            int new_length = this_length + add_length;

            KPMP_category_for_cluster_enum new_category;
            Dictionary<KPMP_category_for_cluster_enum, bool> category_dict = new Dictionary<KPMP_category_for_cluster_enum, bool>();
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                new_category = add_cluster_characterizations[indexAdd].Category;
                if (!category_dict.ContainsKey(new_category))
                {
                    category_dict.Add(new_category, true);
                    Check_that_category_does_not_exist(new_category);
                }
            }
            KPMP_cluster_characterization_line_class[] new_cluster_characterizations = new KPMP_cluster_characterization_line_class[new_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_cluster_characterizations[indexNew] = this.Cluster_characterizations[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_cluster_characterizations[indexNew] = add_cluster_characterizations[indexAdd];
            }
            this.Cluster_characterizations = new_cluster_characterizations;
        }
        #endregion

        #region Generate from trajectories
        private void Generate_unique_hierarchy_lines_from_trajectories(Trajectory_class trajectories)
        {
            trajectories.Trajectories = trajectories.Trajectories.OrderBy(l => l.Assay_cellType).ThenBy(l=>l.Lineage_no).ThenBy(l => l.Position).ToArray();
            int trajectories_length = trajectories.Trajectories.Length;
            Trajectory_line_class trajectory_line;
            string last_source = "";
            KPMP_cluster_hierarchy_line_class hierarchy_line;
            List<KPMP_cluster_hierarchy_line_class> new_hierarchy_lines = new List<KPMP_cluster_hierarchy_line_class>();
            KPMP_cluster_characterization_line_class new_cluster_characterization_line;
            List<KPMP_cluster_characterization_line_class> cluster_characterizations = new List<KPMP_cluster_characterization_line_class>();
            Dictionary<string, Dictionary<string, bool>> source_target_dict = new Dictionary<string, Dictionary<string, bool>>();
            for (int indexT = 0; indexT < trajectories_length; indexT++)
            {
                trajectory_line = trajectories.Trajectories[indexT];
                if (  (indexT==0)
                    || (!trajectory_line.Lineage_no.Equals(trajectories.Trajectories[indexT-1].Lineage_no)))
                {
                    last_source = "";
                }
                if (!String.IsNullOrEmpty(last_source))
                {
                    if (!source_target_dict.ContainsKey(last_source)) { source_target_dict.Add(last_source, new Dictionary<string, bool>()); }
                    if (!source_target_dict[last_source].ContainsKey(trajectory_line.Subtype_name))
                    {
                        source_target_dict[last_source].Add(trajectory_line.Subtype_name, true);
                        hierarchy_line = new KPMP_cluster_hierarchy_line_class();
                        hierarchy_line.Source_cellSubtypeName = (string)last_source.Clone();
                        hierarchy_line.Target_cellSubtypeName = (string)trajectory_line.Subtype_name.Clone();
                        hierarchy_line.Assay_cellType = (string)trajectory_line.Assay_cellType.Clone();
                        new_hierarchy_lines.Add(hierarchy_line);
                    }
                }
                last_source = trajectory_line.Subtype_name;
                new_cluster_characterization_line = new KPMP_cluster_characterization_line_class();
            }
            this.Hierarchy = new_hierarchy_lines.ToArray();
        }
        private void Generate_characterization_lines_from_enrichment_results(Enrichment_results_class enrich, string assay_cellType)
        {
            int enrich_length = enrich.Enrich.Length;
            KPMP_cluster_characterization_line_class[] cluster_char_lines = new KPMP_cluster_characterization_line_class[enrich_length];
            KPMP_cluster_characterization_line_class cluster_char_line;
            Enrichment_results_line_class enrich_results_line;
            for (int indexE=0; indexE<enrich_length; indexE++)
            {
                enrich_results_line = enrich.Enrich[indexE];
                cluster_char_line = new KPMP_cluster_characterization_line_class();
                cluster_char_line.Category_entity = (string)enrich_results_line.SCP.Clone();
                cluster_char_line.Category_entity_value = enrich_results_line.Minus_log10_pvalue;
                cluster_char_line.Category = KPMP_category_for_cluster_enum.Scp_minusLog10Pvalues;
                cluster_char_line.Cell_type = (string)assay_cellType.Clone();
                cluster_char_line.Cell_subtype_name = (string)enrich_results_line.Dataset_name.Clone();
                cluster_char_lines[indexE] = cluster_char_line;
            }
            this.Cluster_characterizations = cluster_char_lines;
        }
        public void Generate_from_trajectories_and_enrichment_results(Trajectory_class trajectories, Enrichment_results_class enrich, string assay_cellType)
        {
            Generate_unique_hierarchy_lines_from_trajectories(trajectories);
            Generate_characterization_lines_from_enrichment_results(enrich,assay_cellType);
            Check_for_correctness();
        }
        #endregion

        #region Write, copy
        public void Write(string directory, string fileName)
        {
            string hierarchy_fileName = Path.GetFileNameWithoutExtension(fileName) + "_hierarchy" + Path.GetExtension(fileName);
            KPMP_cluster_hierarchy_readWriteOptions_class hierarchy_readWriteOptions = new KPMP_cluster_hierarchy_readWriteOptions_class(directory, hierarchy_fileName);
            ReadWriteClass.WriteData(this.Hierarchy, hierarchy_readWriteOptions);
            string characterization_fileName = Path.GetFileNameWithoutExtension(fileName) + "_characterization" + Path.GetExtension(fileName);
            KPMP_cluster_characterization_readWriteOptions_class characterization_readWriteOptions = new KPMP_cluster_characterization_readWriteOptions_class(directory, characterization_fileName);
            ReadWriteClass.WriteData(this.Cluster_characterizations, characterization_readWriteOptions);
        }
        public void Read(string directory, string fileName)
        {
            string hierarchy_fileName = Path.GetFileNameWithoutExtension(fileName) + "_hierarchy" + Path.GetExtension(fileName);
            KPMP_cluster_hierarchy_readWriteOptions_class hierarchy_readWriteOptions = new KPMP_cluster_hierarchy_readWriteOptions_class(directory, hierarchy_fileName);
            this.Hierarchy = ReadWriteClass.ReadRawData_and_FillArray<KPMP_cluster_hierarchy_line_class>(hierarchy_readWriteOptions);
            string characterization_fileName = Path.GetFileNameWithoutExtension(fileName) + "_characterization" + Path.GetExtension(fileName);
            KPMP_cluster_characterization_readWriteOptions_class characterization_readWriteOptions = new KPMP_cluster_characterization_readWriteOptions_class(directory, characterization_fileName);
            this.Cluster_characterizations = ReadWriteClass.ReadRawData_and_FillArray<KPMP_cluster_characterization_line_class>(characterization_readWriteOptions);
        }
        public void Read_and_add_to_array(string directory, string fileName)
        {
            string hierarchy_fileName = Path.GetFileNameWithoutExtension(fileName) + "_hierarchy" + Path.GetExtension(fileName);
            KPMP_cluster_hierarchy_readWriteOptions_class hierarchy_readWriteOptions = new KPMP_cluster_hierarchy_readWriteOptions_class(directory, hierarchy_fileName);
            KPMP_cluster_hierarchy_line_class[] add_hierarchy = ReadWriteClass.ReadRawData_and_FillArray<KPMP_cluster_hierarchy_line_class>(hierarchy_readWriteOptions);
            Add_to_hierarchy(add_hierarchy);
            string characterization_fileName = Path.GetFileNameWithoutExtension(fileName) + "_characterization" + Path.GetExtension(fileName);
            KPMP_cluster_characterization_readWriteOptions_class characterization_readWriteOptions = new KPMP_cluster_characterization_readWriteOptions_class(directory, characterization_fileName);
            KPMP_cluster_characterization_line_class[] add_characterizations = ReadWriteClass.ReadRawData_and_FillArray<KPMP_cluster_characterization_line_class>(characterization_readWriteOptions);
            Add_to_cluster_characterizations(add_characterizations);
        }
        public void Add_other(KPMP_cluster_hierarchy_class other)
        {
            this.Add_to_cluster_characterizations(other.Cluster_characterizations);
            this.Add_to_hierarchy(other.Hierarchy);
        }
        private KPMP_cluster_characterization_line_class[] Get_cluster_characterization_lines_of_given_cellType_and_category_if_value_larger_zero(string cellType, KPMP_category_for_cluster_enum category)
        {
            List<KPMP_cluster_characterization_line_class> sameCellType_sameCategory_lines = new List<KPMP_cluster_characterization_line_class>();
            foreach (KPMP_cluster_characterization_line_class cluster_characterization_line in this.Cluster_characterizations)
            {
                if (  (cluster_characterization_line.Category.Equals(category))
                    &&(cluster_characterization_line.Cell_type.Equals(cellType))
                    &&(cluster_characterization_line.Category_entity_value>0))
                {
                    sameCellType_sameCategory_lines.Add(cluster_characterization_line);
                }
            }
            return sameCellType_sameCategory_lines.ToArray();
        }
        public void Write_as_network_for_each_cellType(KPMP_category_for_cluster_enum category, Dictionary<string, int> categoryEntity_pieChartNo_dict, Dictionary<string, System.Drawing.Color> categoryEntity_color_dict, string directory, string fileBaseName, params string[] addNodes)
        {
            this.Hierarchy = this.Hierarchy.OrderBy(l => l.Assay_cellType).ToArray();
            int hierarchy_length = this.Hierarchy.Length;
            KPMP_cluster_hierarchy_line_class hierarchy_line;
            List<KPMP_cluster_hierarchy_line_class> sameCellType_hierarchy_lines = new List<KPMP_cluster_hierarchy_line_class>();
            KPMP_cluster_characterization_line_class[] sameCellType_characterization_lines;
            for (int indexH = 0; indexH < hierarchy_length; indexH++)
            {
                hierarchy_line = this.Hierarchy[indexH];
                if ((indexH == 0)
                    || (!hierarchy_line.Assay_cellType.Equals(this.Hierarchy[indexH - 1].Assay_cellType)))
                {
                    sameCellType_hierarchy_lines.Clear();
                }
                sameCellType_hierarchy_lines.Add(hierarchy_line);
                if ((indexH == hierarchy_length - 1)
                    || (!hierarchy_line.Assay_cellType.Equals(this.Hierarchy[indexH + 1].Assay_cellType)))
                {
                    sameCellType_characterization_lines = this.Get_cluster_characterization_lines_of_given_cellType_and_category_if_value_larger_zero(hierarchy_line.Assay_cellType, category);
                    Write_as_network_for_same_cellType(sameCellType_hierarchy_lines.ToArray(), sameCellType_characterization_lines, categoryEntity_pieChartNo_dict, categoryEntity_color_dict, directory, fileBaseName + "_" + hierarchy_line.Assay_cellType);
                }
            }
        }
        private void Write_as_network_for_same_cellType(KPMP_cluster_hierarchy_line_class[] sameCellType_hierachy_lines, KPMP_cluster_characterization_line_class[] sameCellType_characterization_lines, Dictionary<string,int> categoryEntity_pieChartNo_dict, Dictionary<string,System.Drawing.Color> categoryEntity_color_dict, string directory, string fileBaseName)
        {
            if (sameCellType_characterization_lines.Length >= 1)
            {
                string cellType = (string)sameCellType_characterization_lines[0].Cell_type.Clone();
                KPMP_category_for_cluster_enum category = sameCellType_characterization_lines[0].Category;
                int categorization_lines_length = sameCellType_characterization_lines.Length;
                KPMP_cluster_characterization_line_class categorization_line;

                #region Check for correctness
                for (int indexCat=0; indexCat < categorization_lines_length; indexCat++)
                {
                    categorization_line = sameCellType_characterization_lines[indexCat];
                    if (!categorization_line.Cell_type.Equals(cellType)) { throw new Exception(); }
                    if (!categorization_line.Category.Equals(category)) { throw new Exception(); }
                }
                int hierarchy_lines_length = sameCellType_hierachy_lines.Length;
                KPMP_cluster_hierarchy_line_class hierarchy_line;
                for (int indexH=0; indexH<hierarchy_lines_length;indexH++)
                {
                    hierarchy_line = sameCellType_hierachy_lines[indexH];
                    if (!hierarchy_line.Assay_cellType.Equals(cellType)) { throw new Exception(); }
                }
                #endregion

                #region Generate nw
                NetworkBasis_class nw = new NetworkBasis_class();
                nw.Generate_and_overwrite_network(sameCellType_hierachy_lines, sameCellType_characterization_lines, categoryEntity_pieChartNo_dict, categoryEntity_color_dict);
                #endregion

                #region Generate and add legend nodes
                float reference_category_value = nw.UniqueNodes.Get_max_nodeColor_value();
                if (Options.Forced_reference_category_value != -1)
                {
                    reference_category_value = Options.Forced_reference_category_value;
                }
                float[] legend_sizes_fraction_of_reference = new float[] { 0.25F, 0.5F, 0.75F, 1F };
                float fraction;
                int fractions_length = legend_sizes_fraction_of_reference.Length;
                List<KPMP_cluster_characterization_line_class> legend_characterization_lines = new List<KPMP_cluster_characterization_line_class>();
                KPMP_cluster_characterization_line_class legend_characterization_line;
                float entity_value;
                for (int indexF=0; indexF<fractions_length; indexF++)
                {
                    fraction = legend_sizes_fraction_of_reference[indexF];
                    entity_value = reference_category_value * fraction;
                    legend_characterization_line = new KPMP_cluster_characterization_line_class();
                    legend_characterization_line.Cell_type = cellType;
                    legend_characterization_line.Cell_subtype_name = "Pie chart size ~ sum(-log10(p))\nthis area represents " + entity_value + "\n(circles with only one color are not size adjusted)";
                    legend_characterization_line.Category_entity = "Pie chart size slice a";
                    legend_characterization_line.Category = category;
                    legend_characterization_line.Category_entity_value = 0.5F * entity_value;
                    legend_characterization_lines.Add(legend_characterization_line);
                    if (indexF==0)
                    {
                        categoryEntity_pieChartNo_dict.Add(legend_characterization_line.Category_entity, 1);
                        categoryEntity_color_dict.Add(legend_characterization_line.Category_entity, System.Drawing.Color.Bisque);
                    }
                    legend_characterization_line = new KPMP_cluster_characterization_line_class();
                    legend_characterization_line.Cell_type = cellType;
                    legend_characterization_line.Cell_subtype_name = "Pie chart size ~ sum(-log10(p))\nthis area represents " + entity_value + "\n(circles with only one color are not size adjusted)";
                    legend_characterization_line.Category_entity = "Pie chart size slice b";
                    legend_characterization_line.Category = category;
                    legend_characterization_line.Category_entity_value = 0.5F * entity_value;
                    legend_characterization_lines.Add(legend_characterization_line);
                    if (indexF == 0)
                    {
                        categoryEntity_pieChartNo_dict.Add(legend_characterization_line.Category_entity, 1);
                        categoryEntity_color_dict.Add(legend_characterization_line.Category_entity, System.Drawing.Color.Aqua);
                    }
                }
                Dictionary<string, bool> categoryEntity_considered_dict = new Dictionary<string, bool>();
                for (int indexCat=0; indexCat<categorization_lines_length; indexCat++)
                {
                    categorization_line = sameCellType_characterization_lines[indexCat];
                    if (!categoryEntity_considered_dict.ContainsKey(categorization_line.Category_entity))
                    {
                        categoryEntity_considered_dict.Add(categorization_line.Category_entity, true);
                        legend_characterization_line = new KPMP_cluster_characterization_line_class();
                        legend_characterization_line.Cell_type = cellType;
                        legend_characterization_line.Cell_subtype_name = (string)categorization_line.Category_entity.Clone();
                        legend_characterization_line.Category_entity = (string)categorization_line.Category_entity.Clone();
                        legend_characterization_line.Category = category;
                        legend_characterization_line.Category_entity_value = reference_category_value * 0.2F;
                        legend_characterization_lines.Add(legend_characterization_line);
                    }
                }
                nw.Add_singleNodes_and_reindex_nw(legend_characterization_lines.ToArray(), categoryEntity_pieChartNo_dict, categoryEntity_color_dict);
                #endregion

                #region Generate visu_make
                RegularNW_generate_visualization_class visu_make = new RegularNW_generate_visualization_class();
                visu_make.Options.Headline = (string)cellType.Clone();
                visu_make.Options.Node_label = Regular_node_label_enum.Name;
                visu_make.Options.Node_shape = Regular_node_shape_enum.Ellipse;
                visu_make.Options.Node_color = Regular_node_color_enum.Selected_color_size_by_value;
                visu_make.Options.Add_border_lines = false;
                visu_make.Options.Reference_minLog10pvalue = reference_category_value;
                visu_make.Options.Reference_node_diameter = 200;
                visu_make.Options.Node_label_font_size = 25;
                #endregion

                Visualization_of_nw_basis visu = visu_make.Generate_visualization_instance(nw);
                yED_class yed = new yED_class();
                yed.Write_yED_file(visu, directory, fileBaseName + "_" + category);
            }
        }

        public KPMP_cluster_characterization_line_class[] Deep_copy_cluster_characterization_lines()
        {
            int characterization_lines_length = this.Cluster_characterizations.Length;
            KPMP_cluster_characterization_line_class[] copy_cluster_characterization_lines = new KPMP_cluster_characterization_line_class[characterization_lines_length];
            for (int indexClusterCharacterization = 0; indexClusterCharacterization < characterization_lines_length; indexClusterCharacterization++)
            {
                copy_cluster_characterization_lines[indexClusterCharacterization] = this.Cluster_characterizations[indexClusterCharacterization].Deep_copy();
            }
            return copy_cluster_characterization_lines;
        }

        public KPMP_cluster_hierarchy_class Deep_copy()
        {
            KPMP_cluster_hierarchy_class copy = (KPMP_cluster_hierarchy_class)this.MemberwiseClone();
            int hierarchy_length = this.Hierarchy.Length;
            copy.Hierarchy = new KPMP_cluster_hierarchy_line_class[hierarchy_length];
            for (int indexH = 0; indexH < hierarchy_length; indexH++)
            {
                copy.Hierarchy[indexH] = this.Hierarchy[indexH].Deep_copy();
            }
            copy.Cluster_characterizations = Deep_copy_cluster_characterization_lines();
            return copy;
        }
        #endregion
    }
}