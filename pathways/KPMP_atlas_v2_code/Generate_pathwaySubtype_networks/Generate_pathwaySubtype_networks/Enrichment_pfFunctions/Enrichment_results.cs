using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Reflection.Emit;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using Common_functions;
using ReadWrite;

namespace Enrichment_results
{
    class Enrichment_results_line_class
    {
        public Ontology_type_enum Ontology { get; set; }
        public int SCP_level { get; set; }
        public string Integration_group { get; set; }
        public string Dataset_name { get; set; }
        public float Timepoint { get; set; } 
        public Timeunit_enum Timeunit { get; set; }
        public UpDown_enum UpDown_status { get; set;  }
        public string Dataset_color_string { get; set; }
        public string SCP { get; set; }
        public int Bg_gene_symbols_count { get; set; }
        public int Experimental_gene_symbols_count {get;set;}
        public int SCP_gene_symbols_count { get; set; }
        public int Overlap_gene_symbols_count { get; set; }
        public float Pvalue { get; set; }
        public float Minus_log10_pvalue { get; set; }
        public float Fractional_rank { get; set; }
        public string[] Overlapping_gene_symbols { get; set; }
        public string ReadWrite_overlapping_gene_symbols
        {
          get { return ReadWriteClass.Get_writeLine_from_array<string>(Overlapping_gene_symbols,Enrichment_results_readWriteOptions_class.Array_delimiter); }
          set { Overlapping_gene_symbols = ReadWriteClass.Get_array_from_readLine<string>(value, Enrichment_results_readWriteOptions_class.Array_delimiter); }
        }
        public bool Significant { get; set; }

        public static Enrichment_results_line_class[] Order_by_ontology_scp_datasetName_timepoint_upDownStatus(Enrichment_results_line_class[] lines)
        {
            Timeunit_enum timeunit = lines[0].Timeunit;
            Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>>>> ontology_scp_datasetName_timepoint_upDownStatus_dict = new Dictionary<Ontology_type_enum, Dictionary<string, Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>>>>();
            Dictionary<string, Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>>> scp_datasetName_timepoint_upDownStatus_dict = new Dictionary<string, Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>>>();
            Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>> datasetName_timepoint_upDownStatus_dict = new Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>>();
            Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>> timepoint_upDownStatus_dict = new Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>();
            Dictionary<UpDown_enum, List<Enrichment_results_line_class>> upDownStatus_dict = new Dictionary<UpDown_enum, List<Enrichment_results_line_class>>();
            int lines_length = lines.Length;
            Enrichment_results_line_class line;
            for (int indexL = 0; indexL < lines_length; indexL++)
            {
                line = lines[indexL];
                if (!line.Timeunit.Equals(timeunit)) { throw new Exception(); }
                if (!ontology_scp_datasetName_timepoint_upDownStatus_dict.ContainsKey(line.Ontology))
                {
                    ontology_scp_datasetName_timepoint_upDownStatus_dict.Add(line.Ontology, new Dictionary<string, Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>>>());
                }
                if (!ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology].ContainsKey(line.SCP))
                {
                    ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology].Add(line.SCP, new Dictionary<string, Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>>());
                }
                if (!ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology][line.SCP].ContainsKey(line.Dataset_name))
                {
                    ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology][line.SCP].Add(line.Dataset_name, new Dictionary<float, Dictionary<UpDown_enum, List<Enrichment_results_line_class>>>());
                }
                if (!ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology][line.SCP][line.Dataset_name].ContainsKey(line.Timepoint))
                {
                    ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology][line.SCP][line.Dataset_name].Add(line.Timepoint, new Dictionary<UpDown_enum, List<Enrichment_results_line_class>>());
                }
                if (!ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology][line.SCP][line.Dataset_name][line.Timepoint].ContainsKey(line.UpDown_status))
                {
                    ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology][line.SCP][line.Dataset_name][line.Timepoint].Add(line.UpDown_status, new List<Enrichment_results_line_class>());
                }
                ontology_scp_datasetName_timepoint_upDownStatus_dict[line.Ontology][line.SCP][line.Dataset_name][line.Timepoint][line.UpDown_status].Add(line);
            }
            lines = null;
            Ontology_type_enum[] ontologies = ontology_scp_datasetName_timepoint_upDownStatus_dict.Keys.ToArray();
            Ontology_type_enum ontology;
            int ontologies_length = ontologies.Length;
            string[] scpNames;
            string scpName;
            int scpNames_length;
            string[] datasetNames;
            string datasetName;
            int datasetNames_length;
            float[] timepoints;
            float timepoint;
            int timepoints_length;
            UpDown_enum[] upDownStatuses;
            UpDown_enum upDownStatus;
            int upDownStatuses_length;

            ontologies = ontologies.OrderBy(l => l).ToArray();
            List<Enrichment_results_line_class> ordered_enrichment_results = new List<Enrichment_results_line_class>();
            for (int indexO = 0; indexO < ontologies_length; indexO++)
            {
                ontology = ontologies[indexO];
                scp_datasetName_timepoint_upDownStatus_dict = ontology_scp_datasetName_timepoint_upDownStatus_dict[ontology];
                scpNames = scp_datasetName_timepoint_upDownStatus_dict.Keys.ToArray();
                scpNames_length = scpNames.Length;
                scpNames = scpNames.OrderBy(l => l).ToArray();
                for (int indexSCP = 0; indexSCP < scpNames_length; indexSCP++)
                {
                    scpName = scpNames[indexSCP];
                    datasetName_timepoint_upDownStatus_dict = scp_datasetName_timepoint_upDownStatus_dict[scpName];
                    datasetNames = datasetName_timepoint_upDownStatus_dict.Keys.ToArray();
                    datasetNames_length = datasetNames.Length;
                    datasetNames = datasetNames.OrderBy(l => l).ToArray();
                    for (int indexDN = 0; indexDN < datasetNames_length; indexDN++)
                    {
                        datasetName = datasetNames[indexDN];
                        timepoint_upDownStatus_dict = datasetName_timepoint_upDownStatus_dict[datasetName];
                        timepoints = timepoint_upDownStatus_dict.Keys.ToArray();
                        timepoints_length = timepoints.Length;
                        timepoints = timepoints.OrderBy(l => l).ToArray();
                        for (int indexT = 0; indexT < timepoints_length; indexT++)
                        {
                            timepoint = timepoints[indexT];
                            upDownStatus_dict = timepoint_upDownStatus_dict[timepoint];
                            upDownStatuses = upDownStatus_dict.Keys.ToArray();
                            upDownStatuses_length = upDownStatuses.Length;
                            upDownStatuses = upDownStatuses.OrderBy(l => l).ToArray();
                            for (int indexUSD = 0; indexUSD < upDownStatuses_length; indexUSD++)
                            {
                                upDownStatus = upDownStatuses[indexUSD];
                                ordered_enrichment_results.AddRange(upDownStatus_dict[upDownStatus]);
                            }
                        }
                    }
                }
            }

            if (Global_class.Check_for_correct_ordering)
            {
                #region Check for correct ordering
                int ordered_length = ordered_enrichment_results.Count;
                if (ordered_length != lines_length) { throw new Exception(); }
                Enrichment_results_line_class previous_line;
                Enrichment_results_line_class current_line;
                for (int indexOrder = 1; indexOrder < ordered_length; indexOrder++)
                {
                    previous_line = ordered_enrichment_results[indexOrder - 1];
                    current_line = ordered_enrichment_results[indexOrder];
                    if (current_line.Ontology.CompareTo(previous_line.Ontology) < 0) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.SCP.CompareTo(previous_line.SCP) < 0)) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.SCP.Equals(previous_line.SCP))
                        && (current_line.Dataset_name.CompareTo(previous_line.Dataset_name) < 0)) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.SCP.Equals(previous_line.SCP))
                        && (current_line.Dataset_name.Equals(previous_line.Dataset_name))
                        && (current_line.Timepoint.CompareTo(previous_line.Timepoint) < 0)) { throw new Exception(); }
                    if ((current_line.Ontology.Equals(previous_line.Ontology))
                        && (current_line.SCP.Equals(previous_line.SCP))
                        && (current_line.Dataset_name.Equals(previous_line.Dataset_name))
                        && (current_line.Timepoint.Equals(previous_line.Timepoint))
                        && (current_line.UpDown_status.CompareTo(previous_line.UpDown_status) < 0)) { throw new Exception(); }
                }
            }
            #endregion

            return ordered_enrichment_results.ToArray();
        }

        public Enrichment_results_line_class()
        {
            SCP = "";
            Dataset_color_string = "";
            Dataset_name = "";
            Integration_group = "";
            Overlapping_gene_symbols = new string[0];
        }
        public Enrichment_results_line_class Deep_copy()
        {
            Enrichment_results_line_class copy = (Enrichment_results_line_class)this.MemberwiseClone();
            copy.SCP = (string)this.SCP.Clone();
            copy.Dataset_name = (string)this.Dataset_name.Clone();
            copy.Integration_group = (string)this.Integration_group.Clone();
            copy.Dataset_color_string = (string)this.Dataset_color_string.Clone();
            return copy;
        }
    }
    class Enrichment_results_readWriteOptions_class : ReadWriteOptions_base
    {
        public static char Array_delimiter { get { return ','; } }
        public Enrichment_results_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Ontology", "SCP_level", "Dataset_name", "Dataset_color_string", "SCP", "Bg_gene_symbols_count", "Experimental_gene_symbols_count", "SCP_gene_symbols_count", "Overlap_gene_symbols_count", "Pvalue", "Minus_log10_pvalue", "Fractional_rank", "ReadWrite_overlapping_gene_symbols", "Significant" };
            this.Key_columnNames   = new string[] { "Ontology", "SCP.level", "Dataset.name", "Dataset.color",        "SCP", "Bg.gene.symbols.count", "Experimental.gene.symbols.count", "SCP.gene.symbols.count", "Overlap.gene.symbols.count", "Pvalue", "Minus.log10_pvalue", "Fractional.rank", "Overlapping.gene.symbols", "Significant" };
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
            this.Add_to_existing_file = false;
        }
    }
    class Enrichment_results_class
    {
        public Enrichment_results_line_class[] Enrich { get; set; }
        public Enrichment_results_class()
        {
            this.Enrich = new Enrichment_results_line_class[0];
        }
        public void Order_by_ontology_scp_datasetName_timepoint_upDownStatus()
        {
            this.Enrich = Enrichment_results_line_class.Order_by_ontology_scp_datasetName_timepoint_upDownStatus(this.Enrich);
        }

        public void Keep_only_lines_of_given_scpLevel(int level)
        {
            List<Enrichment_results_line_class> keep = new List<Enrichment_results_line_class>();
            foreach (Enrichment_results_line_class enrich_line in this.Enrich)
            {
                if (enrich_line.SCP_level.Equals(level))
                {
                    keep.Add(enrich_line);
                }
            }
            this.Enrich = keep.ToArray();
        }

        public void Read(string directory, string fileName)
        {
            Enrichment_results_readWriteOptions_class readWriteOptions = new Enrichment_results_readWriteOptions_class(directory, fileName);
            this.Enrich = ReadWriteClass.ReadRawData_and_FillArray<Enrichment_results_line_class>(readWriteOptions);
        }
        public void Write(string directory, string fileName)
        {
            Enrichment_results_readWriteOptions_class readWriteOptions = new Enrichment_results_readWriteOptions_class(directory, fileName);
            ReadWriteClass.WriteData(this.Enrich, readWriteOptions);
        }
        public Enrichment_results_class Deep_copy()
        {
            Enrichment_results_class copy = (Enrichment_results_class)this.MemberwiseClone();
            int enrich_length = this.Enrich.Length;
            copy.Enrich = new Enrichment_results_line_class[enrich_length];
            for (int indexE=0; indexE<enrich_length; indexE++)
            {
                copy.Enrich[indexE] = this.Enrich[indexE].Deep_copy();
            }
            return copy;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class PhysiologicalFunctionsDefinition_line_class
    {
        public Ontology_type_enum HigherLevelFunction_ontology { get; set; }
        public Ontology_type_enum Ontology { get; set; }
        public string PF_group { get; set; }
        public string SCP { get; set; }
        public string SCP_abbreviation { get; set; }
        public int Index_in_SCP_distance_dendrogram { get; set; }
        public string Hexadecimal_color { get; set; }
        public int PieChart_slice_order { get; set; }

        public PhysiologicalFunctionsDefinition_line_class()
        {
            PF_group = "";
            SCP = "";
            SCP_abbreviation = "";
            Hexadecimal_color = "";
        }
        public PhysiologicalFunctionsDefinition_line_class Deep_copy()
        {
            PhysiologicalFunctionsDefinition_line_class copy = (PhysiologicalFunctionsDefinition_line_class)this.MemberwiseClone();
            copy.PF_group = (string)this.PF_group.Clone();
            copy.SCP = (string)this.SCP.Clone();
            copy.SCP_abbreviation = (string)this.SCP_abbreviation.Clone();
            copy.Hexadecimal_color = (string)this.Hexadecimal_color.Clone();
            return copy;
        }
    }

    class PhysiologicalFunctionsDefinition_readWriteOptions_class : ReadWriteOptions_base
    {
        public PhysiologicalFunctionsDefinition_readWriteOptions_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "PF_group", "SCP", "Hexadecimal_color", "PieChart_slice_order" };
            this.Key_columnNames = this.Key_propertyNames;
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
        }
    }

    class PhysiologicalFunctionsDefinition_class
    {
        public PhysiologicalFunctionsDefinition_line_class[] HigherFunctions { get; set; }

        public PhysiologicalFunctionsDefinition_class()
        {
            HigherFunctions = new PhysiologicalFunctionsDefinition_line_class[0];
        }

        private void Keep_only_selected_ontology(Ontology_type_enum ontology)
        {
            List<PhysiologicalFunctionsDefinition_line_class> keep = new List<PhysiologicalFunctionsDefinition_line_class>();
            foreach (PhysiologicalFunctionsDefinition_line_class higherFunction_line in HigherFunctions)
            {
                if (higherFunction_line.Ontology.Equals(ontology))
                {
                    keep.Add(higherFunction_line);
                }
            }
            this.HigherFunctions = keep.ToArray();
        }
        public void Generate_by_reading(string directory, string fileName)
        {
            Read(directory, fileName);
        }

        public Ontology_type_enum Get_ontology_and_check_if_only_one()
        {
            Ontology_type_enum ontology = this.HigherFunctions[0].Ontology;
            foreach (PhysiologicalFunctionsDefinition_line_class higherFunction_line in HigherFunctions)
            {
                if (!higherFunction_line.Ontology.Equals(ontology)) { throw new Exception(); }
            }
            return ontology;
        }
        public Dictionary<string, string> Get_scp_wholeCellFunction_dict()
        {
            Dictionary<string, string> scp_wholeCellFunction_dict = new Dictionary<string, string>();
            List<string> scps = new List<string>();
            int higherFunctions_length = this.HigherFunctions.Length;
            PhysiologicalFunctionsDefinition_line_class higherFunction_line;
            for (int indexHF = 0; indexHF < higherFunctions_length; indexHF++)
            {
                higherFunction_line = this.HigherFunctions[indexHF];
                if (!String.IsNullOrEmpty(higherFunction_line.PF_group))
                {
                    if ((indexHF != 0)
                        && (!higherFunction_line.Ontology.Equals(this.HigherFunctions[indexHF - 1].Ontology)))
                    {
                        throw new Exception();
                    }
                    if (!scp_wholeCellFunction_dict.ContainsKey(higherFunction_line.SCP))
                    {
                        scp_wholeCellFunction_dict.Add(higherFunction_line.SCP, higherFunction_line.PF_group);
                    }
                    else if (!scp_wholeCellFunction_dict[higherFunction_line.SCP].Equals(higherFunction_line.PF_group)) { throw new Exception(); }
                }
            }
            return scp_wholeCellFunction_dict;
        }
        public Dictionary<string, System.Drawing.Color> Get_pfGroup_color_dict()
        {
            Color_class color_conversion = new Color_class();
            Dictionary<string, System.Drawing.Color> higherFunction_color_dict = new Dictionary<string, System.Drawing.Color>();
            this.HigherFunctions = this.HigherFunctions.OrderBy(l => l.Ontology).ThenBy(l => l.PF_group).ThenByDescending(l => l.Hexadecimal_color).ToArray();
            int higherFunctions_length = this.HigherFunctions.Length;
            PhysiologicalFunctionsDefinition_line_class higherFunction_line;
            List<string> pf_groups_with_mismatching_colors = new List<string>();
            for (int indexHF = 0; indexHF < higherFunctions_length; indexHF++)
            {
                higherFunction_line = this.HigherFunctions[indexHF];
                if ((indexHF > 0)
                    && (higherFunction_line.PF_group.Equals(this.HigherFunctions[indexHF - 1].PF_group))
                    && (!higherFunction_line.Hexadecimal_color.Equals(this.HigherFunctions[indexHF - 1].Hexadecimal_color)))
                {
                    pf_groups_with_mismatching_colors.Add(higherFunction_line.PF_group);
                }
                if ((indexHF == higherFunctions_length - 1)
                    || (!higherFunction_line.PF_group.Equals(this.HigherFunctions[indexHF + 1].PF_group)))
                {
                    if (  (!String.IsNullOrEmpty(higherFunction_line.PF_group))
                        &&(!higherFunction_line.Hexadecimal_color.Equals("")))
                    {
                        higherFunction_color_dict.Add(higherFunction_line.PF_group, color_conversion.Get_closest_csharp_color_for_hexadecimal_color_if_exists(higherFunction_line.Hexadecimal_color));
                    }
                }
            }
            if (pf_groups_with_mismatching_colors.Count>0)
            {
                pf_groups_with_mismatching_colors = pf_groups_with_mismatching_colors.Distinct().OrderBy(l=>l).ToList();
                StringBuilder sb = new StringBuilder();
                sb.AppendFormat("Mismatching colors for the following whole cell functions: ");
                foreach (string pfGroup in pf_groups_with_mismatching_colors)
                {
                    if (sb.Length>0) { sb.Append(", "); }
                    sb.AppendFormat(pfGroup);
                }
                throw new Exception(sb.ToString());
            }
            return higherFunction_color_dict;
        }


        public Dictionary<string, int> Get_pfGroup_pieChartSliceNumber_dict()
        {
            Color_class color_conversion = new Color_class();
            Dictionary<string, int> hfGroup_pieChartSliceNumber_dict = new Dictionary<string, int>();
            this.HigherFunctions = this.HigherFunctions.OrderBy(l => l.Ontology).ThenBy(l=>l.PF_group).ThenBy(l => l.PieChart_slice_order).ToArray();
            int higherFunctions_length = this.HigherFunctions.Length;
            PhysiologicalFunctionsDefinition_line_class higherFunction_line;
            for (int indexHF = 0; indexHF < higherFunctions_length; indexHF++)
            {
                higherFunction_line = this.HigherFunctions[indexHF];
                if ((indexHF > 0)
                    && (!higherFunction_line.PieChart_slice_order.Equals(this.HigherFunctions[indexHF - 1].PieChart_slice_order)))
                {
                    //throw new Exception();
                }
                if ((indexHF == higherFunctions_length - 1)
                    || (!higherFunction_line.PF_group.Equals(this.HigherFunctions[indexHF + 1].PF_group)))
                {
                    if (!String.IsNullOrEmpty(higherFunction_line.PF_group))
                    {
                        hfGroup_pieChartSliceNumber_dict.Add(higherFunction_line.PF_group, higherFunction_line.PieChart_slice_order);
                    }
                }
            }
            return hfGroup_pieChartSliceNumber_dict;
        }

        private void Read(string directory, string fileName)
        {
            PhysiologicalFunctionsDefinition_readWriteOptions_class readWriteOptions = new PhysiologicalFunctionsDefinition_readWriteOptions_class(directory, fileName);
            this.HigherFunctions = ReadWriteClass.ReadRawData_and_FillArray<PhysiologicalFunctionsDefinition_line_class>(readWriteOptions);
        }

    }




}
