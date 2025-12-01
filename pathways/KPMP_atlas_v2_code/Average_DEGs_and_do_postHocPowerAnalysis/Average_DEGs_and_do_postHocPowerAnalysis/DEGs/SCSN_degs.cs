using System;
using System.Collections.Generic;
using System.Diagnostics.Metrics;
using System.Linq;
using System.Reflection;
using System.Reflection.Metadata.Ecma335;
using System.Runtime.Serialization.Formatters;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Schema;
using Common_functions;
using ReadWrite;

namespace DEGs_and_postHocPower
{

    enum AvgValueType_of_interest_enum {  E_m_p_t_y, Minus_log10_adj_pvalue, Avg_log2_fc, Rank, Is_sig }
    enum Statistical_measure_enum {  E_m_p_t_y, Coefficient_of_variation }
    enum Correlation_valueType_enum {  E_m_p_t_y, Pearson_coefficient, Jaccard_index }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class SCSN_deg_line_class
    {
        public string Cluster { get; set; }
        public string Gene { get; set; }
        public double Minus_log10_adj_pvalue { get; set; }
        public float Avg_log2FC { get; set; }
        public double P_val_adj { get; set; }
        public float Rank { get; set; }
        public int Is_significant { get; set; }
        public SCSN_deg_line_class()
        {
            Cluster = "";
            Gene = "";
            Rank = -1;
        }
        public SCSN_deg_line_class Deep_copy()
        {
            SCSN_deg_line_class copy = (SCSN_deg_line_class)MemberwiseClone();
            copy.Cluster = (string)this.Cluster.Clone();
            copy.Gene = (string)this.Gene.Clone();
            return copy;
        }
    }
    class SCSN_deg_readWrite_options_class : ReadWriteOptions_base
    {
        public SCSN_deg_readWrite_options_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Cluster", "Gene", "Avg_log2FC", "P_val_adj" };
            this.Key_columnNames = new string[] { "cluster", "gene", "avg_log2FC", "p_val_adj" };
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
            this.Report_unhandled_null_entries = true;
        }
    }
    class SCSN_deg_class
    {
        public SCSN_deg_line_class[] DEGs { get; set; }

        public SCSN_deg_class()
        {
            this.DEGs = new SCSN_deg_line_class[0];
        }

        private void Calculate_minusLog10Pvalue()
        {
            double max_minus_log10_adj_pvalue = Global_class.Max_considered_double;
            foreach (SCSN_deg_line_class deg_line in DEGs)
            {
                deg_line.Minus_log10_adj_pvalue = -Math.Log10(deg_line.P_val_adj);
                if (deg_line.Minus_log10_adj_pvalue> max_minus_log10_adj_pvalue)
                { deg_line.Minus_log10_adj_pvalue = max_minus_log10_adj_pvalue; }
            }
        }
        private void Calculate_fractional_ranks_based_on_isSignificant_minusLog10Pvalue_and_avgLog2FC_if_cutoffs_are_met(float rank_for_not_significant)
        {
            Math_class mfs = new Math_class();
            this.DEGs = this.DEGs.OrderBy(l => l.Cluster).ThenByDescending(l => l.Is_significant).ThenByDescending(l => l.Minus_log10_adj_pvalue).ThenByDescending(l => Math.Abs(l.Avg_log2FC)).ToArray();
            int degs_length = this.DEGs.Length;
            SCSN_deg_line_class deg_line;
            SCSN_deg_line_class inner_deg_line;
            List<float> current_ranks = new List<float>();
            float add_rank;
            float running_rank = -1;
            int firstIndex_sameValues = -1;
            for (int indexDeg = 0; indexDeg < degs_length; indexDeg++)
            {
                deg_line = this.DEGs[indexDeg];
                if ((indexDeg == 0)
                    || (!deg_line.Cluster.Equals(this.DEGs[indexDeg - 1].Cluster)))
                {
                    running_rank = 0F;
                }
                if ((indexDeg == 0)
                    || (!deg_line.Cluster.Equals(this.DEGs[indexDeg - 1].Cluster))
                    || (!deg_line.Is_significant.Equals(this.DEGs[indexDeg - 1].Is_significant))
                    || (!deg_line.Minus_log10_adj_pvalue.Equals(this.DEGs[indexDeg - 1].Minus_log10_adj_pvalue))
                    || (!Math.Abs(deg_line.Avg_log2FC).Equals(Math.Abs(this.DEGs[indexDeg - 1].Avg_log2FC))))
                {
                    current_ranks.Clear();
                    firstIndex_sameValues = indexDeg;
                }
                running_rank++;
                current_ranks.Add(running_rank);
                if ((indexDeg == degs_length - 1)
                    || (!deg_line.Cluster.Equals(this.DEGs[indexDeg + 1].Cluster))
                    || (!deg_line.Is_significant.Equals(this.DEGs[indexDeg + 1].Is_significant))
                    || (!deg_line.Minus_log10_adj_pvalue.Equals(this.DEGs[indexDeg + 1].Minus_log10_adj_pvalue))
                    || (!Math.Abs(deg_line.Avg_log2FC).Equals(Math.Abs(this.DEGs[indexDeg + 1].Avg_log2FC))))
                {
                    if (current_ranks.Count == 1) { deg_line.Rank = current_ranks[0]; }
                    else
                    {
                        add_rank = mfs.Get_average(current_ranks.ToArray());
                        for (int indexInner = firstIndex_sameValues; indexInner <= indexDeg; indexInner++)
                        {
                            inner_deg_line = this.DEGs[indexInner];
                            if (inner_deg_line.Is_significant == 1) { inner_deg_line.Rank = add_rank; }
                            else { inner_deg_line.Rank = rank_for_not_significant; }
                        }
                    }
                }
            }

        }
        private void Calculate_fractional_ranks_based_on_isSignificant_minusLog10Pvalue_and_avgLog2FC()
        {
            Math_class mfs = new Math_class();
            this.DEGs = this.DEGs.OrderBy(l => l.Cluster).ThenByDescending(l => l.Minus_log10_adj_pvalue).ThenByDescending(l => Math.Abs(l.Avg_log2FC)).ToArray();
            int degs_length = this.DEGs.Length;
            SCSN_deg_line_class deg_line;
            SCSN_deg_line_class inner_deg_line;
            List<float> current_ranks = new List<float>();
            float add_rank;
            float running_rank = -1;
            int firstIndex_sameValues = -1;
            for (int indexDeg = 0; indexDeg < degs_length; indexDeg++)
            {
                deg_line = this.DEGs[indexDeg];
                if ((indexDeg == 0)
                    || (!deg_line.Cluster.Equals(this.DEGs[indexDeg - 1].Cluster)))
                {
                    running_rank = 0F;
                }
                if ((indexDeg == 0)
                    || (!deg_line.Cluster.Equals(this.DEGs[indexDeg - 1].Cluster))
                    || (!deg_line.Minus_log10_adj_pvalue.Equals(this.DEGs[indexDeg - 1].Minus_log10_adj_pvalue))
                    || (!Math.Abs(deg_line.Avg_log2FC).Equals(Math.Abs(this.DEGs[indexDeg - 1].Avg_log2FC))))
                {
                    current_ranks.Clear();
                    firstIndex_sameValues = indexDeg;
                }
                running_rank++;
                current_ranks.Add(running_rank);
                if ((indexDeg == degs_length - 1)
                    || (!deg_line.Cluster.Equals(this.DEGs[indexDeg + 1].Cluster))
                    || (!deg_line.Minus_log10_adj_pvalue.Equals(this.DEGs[indexDeg + 1].Minus_log10_adj_pvalue))
                    || (!Math.Abs(deg_line.Avg_log2FC).Equals(Math.Abs(this.DEGs[indexDeg + 1].Avg_log2FC))))
                {
                    if (current_ranks.Count == 1) { deg_line.Rank = current_ranks[0]; }
                    else
                    {
                        add_rank = mfs.Get_average(current_ranks.ToArray());
                        for (int indexInner = firstIndex_sameValues; indexInner <= indexDeg; indexInner++)
                        {
                            inner_deg_line = this.DEGs[indexInner];
                            inner_deg_line.Rank = add_rank;
                        }
                    }
                }
            }

        }
        public void Generate_by_reading(string directory, string fileName)
        {
            Read(directory, fileName);
            Calculate_minusLog10Pvalue();
            Calculate_fractional_ranks_based_on_isSignificant_minusLog10Pvalue_and_avgLog2FC();
        }
        public void Remove_genes_below_min_avgLog2Fc(float min_avgLog2Fc)
        {
            List<SCSN_deg_line_class> keep = new List<SCSN_deg_line_class>();
            foreach (SCSN_deg_line_class deg_line in this.DEGs)
            {
                if (deg_line.Avg_log2FC >= min_avgLog2Fc)
                {
                    keep.Add(deg_line);
                }
            }
            this.DEGs = keep.ToArray();

        }
        public void Remove_genes_above_max_adj_pvalue(float max_adj_pvalue)
        {
            List<SCSN_deg_line_class> keep = new List<SCSN_deg_line_class>();
            foreach (SCSN_deg_line_class deg_line in this.DEGs)
            {
                if (deg_line.P_val_adj <= max_adj_pvalue)
                {
                    keep.Add(deg_line);
                }
            }
            this.DEGs = keep.ToArray();
        }
        public void Remove_genes_above_max_significance_rank(float significance_rank_cutoff)
        {
            List<SCSN_deg_line_class> keep = new List<SCSN_deg_line_class>();
            foreach (SCSN_deg_line_class deg_line in this.DEGs)
            {
                if (deg_line.Rank <= significance_rank_cutoff)
                {
                    keep.Add(deg_line);
                }
            }
            this.DEGs = keep.ToArray();
        }
        public void Label_significant_genes(float adj_pvalue_cutoff, float max_rank)
        {
            foreach (SCSN_deg_line_class deg_line in this.DEGs)
            {
                if (  (deg_line.P_val_adj<=adj_pvalue_cutoff)
                    &&(deg_line.Rank<=max_rank))
                {
                    deg_line.Is_significant = 1;
                }
                else
                {
                    deg_line.Is_significant = 0;
                }
            }
        }
        private void Read(string directory, string fileName)
        {
            SCSN_deg_readWrite_options_class readWriteOptions = new SCSN_deg_readWrite_options_class(directory, fileName);
            this.DEGs = ReadWriteClass.ReadRawData_and_FillArray<SCSN_deg_line_class>(readWriteOptions);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class SCSN_avgAnalysesWithinEachIteration_deg_line_class
    {
        #region Fields
        public string Cluster { get; set; }
        public string Gene { get; set; }
        public double Minus_log10_adj_pvalue_arithMean { get; set; }
        public double Minus_log10_adj_pvalue_popVariance { get; set; }
        public double Minus_log10_adj_pvalue_sampleSd { get; set; }
        public double Minus_log10_adj_pvalue_coefVar { get; set; }
        public double Avg_log2FC_arithMean { get; set; }
        public double Avg_log2FC_popVariance { get; set; }
        public double Avg_log2FC_sampleSd { get; set; }
        public double Avg_log2FC_coefVar { get; set; }
        public double Rank_arithMean { get; set; }
        public double Rank_popVariance { get; set; }
        public double Rank_sampleSd { get; set; }
        public double Rank_coefVar { get; set; }
        public double IsSig_arithMean { get; set; }
        public double IsSig_popVariance { get; set; }
        public double IsSig_sampleSd { get; set; }
        public double IsSig_coefVar { get; set; }
        public double P_val_adj_geomMean { get; set; }
        public int Number_of_datasets { get; set; }
        public int Iteration { get; set; }
        public float Rank { get; set; }
        #endregion

        public SCSN_avgAnalysesWithinEachIteration_deg_line_class()
        {
            Cluster = "";
            Gene = "";
            Iteration = -1;
            Minus_log10_adj_pvalue_arithMean = 0;
            Minus_log10_adj_pvalue_popVariance = 0;
            Minus_log10_adj_pvalue_sampleSd = -1;
            Avg_log2FC_arithMean = 0;
            Avg_log2FC_popVariance = 0;
            Avg_log2FC_sampleSd = -1;
            Rank_arithMean = 0;
            Rank_popVariance = 0;
            Rank_sampleSd = -1;
            IsSig_arithMean = 0;
            IsSig_popVariance = 0;
            IsSig_sampleSd = -1;
            Rank = -1;

        }
        public double Get_arithMean_for_value_of_interest(AvgValueType_of_interest_enum value_type_of_interest)
        {
            switch (value_type_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_arithMean;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_arithMean;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_arithMean;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_arithMean;
                default:
                    throw new Exception();
            }
        }
        public double Get_popVariance_for_value_of_interest(AvgValueType_of_interest_enum value_type_of_interest)
        {
            switch (value_type_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_popVariance;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_popVariance;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_popVariance;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_popVariance;
                default:
                    throw new Exception();
            }
        }
        public void Set_sampleSD_for_value_of_interest(double set_value, AvgValueType_of_interest_enum value_type_of_interest)
        {
            switch (value_type_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    Minus_log10_adj_pvalue_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    Avg_log2FC_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Rank:
                    Rank_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Is_sig:
                    IsSig_sampleSd = set_value;
                    break;
                default:
                    throw new Exception();
            }
        }
        public SCSN_avgAnalysesWithinEachIteration_deg_line_class Deep_copy()
        {
            SCSN_avgAnalysesWithinEachIteration_deg_line_class copy = (SCSN_avgAnalysesWithinEachIteration_deg_line_class)MemberwiseClone();
            copy.Cluster = (string)this.Cluster.Clone();
            copy.Gene = (string)this.Gene.Clone();
            return copy;
        }
    }
    class SCSN_avgAnalysesWithinEachIteration_deg_readWrite_options_class : ReadWriteOptions_base
    {
        public SCSN_avgAnalysesWithinEachIteration_deg_readWrite_options_class(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Number_of_datasets", "Iteration", "Cluster","Gene",
                                                    "Minus_log10_adj_pvalue_arithMean","Minus_log10_adj_pvalue_sampleSd","Minus_log10_adj_pvalue_coefVar",
                                                    "Avg_log2FC_arithMean","Avg_log2FC_sampleSd","Avg_log2FC_coefVar", "Minus_log10_adj_pvalue_arithMean",
                                                    "Rank_arithMean", "Rank_sampleSd", "Rank_coefVar",
                                                    "P_val_adj_geomMean", "Rank" };
            this.Key_columnNames = this.Key_propertyNames;
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
            this.Report_unhandled_null_entries = true;
        }
    }
    class SCSN_avgAnalysesWithinEachIteration_deg_class
    {
        public SCSN_avgAnalysesWithinEachIteration_deg_line_class[] DEGs { get; set; }
        private int IndexNextAvailableThis { get; set; }
        private int Array_block_length { get; set; }
        private int Number_of_datasets_added_for_mean { get; set; }
        private int Number_of_datasets_added_for_populationVariance{ get; set; }
        private int Anticipated_number_of_datasets_added { get; set; }
        public bool Averaging_procedure_finished { get; private set; }

        public SCSN_avgAnalysesWithinEachIteration_deg_class(int anticipated_number_of_datasets_added, int array_block_length)
        {
            Number_of_datasets_added_for_mean = 0;
            Number_of_datasets_added_for_populationVariance = 0;
            Averaging_procedure_finished = false;
            Anticipated_number_of_datasets_added = anticipated_number_of_datasets_added;
            IndexNextAvailableThis = 0;
            Array_block_length = array_block_length;
            this.DEGs = new SCSN_avgAnalysesWithinEachIteration_deg_line_class[0];
        }
        private void Check_for_duplicates()
        {
            Dictionary<int, Dictionary<string, Dictionary<string, bool>>> iteration_cluster_gene_dict = new Dictionary<int, Dictionary<string, Dictionary<string, bool>>>();
            foreach (SCSN_avgAnalysesWithinEachIteration_deg_line_class avg_line in this.DEGs)
            {
                if (avg_line != null)
                {
                    if (!iteration_cluster_gene_dict.ContainsKey(avg_line.Iteration))
                    {
                        iteration_cluster_gene_dict.Add(avg_line.Iteration, new Dictionary<string, Dictionary<string, bool>>());
                    }
                    if (!iteration_cluster_gene_dict[avg_line.Iteration].ContainsKey(avg_line.Cluster))
                    {
                        iteration_cluster_gene_dict[avg_line.Iteration].Add(avg_line.Cluster, new Dictionary<string, bool>());
                    }
                    iteration_cluster_gene_dict[avg_line.Iteration][avg_line.Cluster].Add(avg_line.Gene, true);
                }
            }
        }
        public void Check_for_correctness()
        {
            if (this.DEGs.Length!=Get_this_length_filled_with_data()) { throw new Exception(); }
            Check_for_duplicates();
        }

        private void Add_to_array(SCSN_avgAnalysesWithinEachIteration_deg_line_class[] add_degs)
        {
            int this_length = this.DEGs.Length;
            int add_length = add_degs.Length;
            int new_length = this_length + add_length;
            int indexNew = -1;
            SCSN_avgAnalysesWithinEachIteration_deg_line_class[] new_degs = new SCSN_avgAnalysesWithinEachIteration_deg_line_class[new_length];
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                new_degs[indexNew] = this.DEGs[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                new_degs[indexNew] = add_degs[indexAdd];
            }
            this.DEGs = new_degs;
        }
        private void Extend_array()
        {
            int this_length = this.DEGs.Length;
            int new_length = this_length + Array_block_length;
            SCSN_avgAnalysesWithinEachIteration_deg_line_class[] new_degs = new SCSN_avgAnalysesWithinEachIteration_deg_line_class[new_length];
            for (int indexThis=0; indexThis<this_length; indexThis++)
            {
                new_degs[indexThis] = this.DEGs[indexThis];
            }
            this.DEGs = new_degs;
        }
        private void Add_to_array_using_indexNextAvailableThis(SCSN_avgAnalysesWithinEachIteration_deg_line_class[] add_degs)
        {
            int this_length = this.DEGs.Length;
            int add_length = add_degs.Length;
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                if (IndexNextAvailableThis >= this_length)
                { 
                    if (add_length > Array_block_length) { Array_block_length = add_length; }
                    Extend_array();
                    this_length = this.DEGs.Length;
                }
                this.DEGs[IndexNextAvailableThis] = add_degs[indexAdd];
                IndexNextAvailableThis++;
            }
        }
        private int Get_this_length_filled_with_data()
        {
            return IndexNextAvailableThis;
        }
        public void Remove_empty_array_lines()
        {
            int this_length = this.DEGs.Length;
            int this_length_filled_with_data = Get_this_length_filled_with_data();
            if (this_length > this_length_filled_with_data)
            {
                int new_length = this_length_filled_with_data;
                SCSN_avgAnalysesWithinEachIteration_deg_line_class[] new_degs = new SCSN_avgAnalysesWithinEachIteration_deg_line_class[new_length];
                for (int indexThis = 0; indexThis < new_length; indexThis++)
                {
                    new_degs[indexThis] = this.DEGs[indexThis];
                }
                this.DEGs = new_degs;
            }
        }

        public void Add_other_using_indexNextAvailableThis(SCSN_avgAnalysesWithinEachIteration_deg_class other_degs)
        {
            other_degs.Check_for_correctness();
            if (DEGs.Length==0) { Averaging_procedure_finished = true; }
            if (!Averaging_procedure_finished) { throw new Exception("Add other: (Averaging_procedure_finsihed)"); }
            if (!other_degs.Averaging_procedure_finished) { throw new Exception("Add other: !other_degs.Averaging_procedure_finsihed"); }
            Add_to_array_using_indexNextAvailableThis(other_degs.DEGs);
        }
        public void Add_iteration_to_all_degs(int iteration_no)
        {
            foreach (SCSN_avgAnalysesWithinEachIteration_deg_line_class degs_line in this.DEGs)
            {
                if (degs_line.Iteration!=-1) { throw new Exception("Add_iteration_to_all_degs: iteration already set"); }
                degs_line.Iteration = iteration_no;
            }
        }
        //public void Add_other_dataset_by_summing_up_running_mean_and_variances(SCSN_deg_class other)
        //{
        //    Math_class mfs = new Math_class();
        //    if (Number_of_datasets_added == Anticipated_number_of_datasets_added) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }
        //    if (Averaging_procedure_finsihed) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }
        //    if ((DEGs.Length != 0) && (Number_of_datasets_added == 0)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }
        //    if ((DEGs.Length == 0) && (Number_of_datasets_added != 0)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }
        //    Number_of_datasets_added++;
        //    Dictionary<string, Dictionary<string, SCSN_deg_line_class>> other_cluster_gene_degLine_dict = new Dictionary<string, Dictionary<string, SCSN_deg_line_class>>();
        //    Dictionary<string, SCSN_deg_line_class> other_gene_degLine_dict = new Dictionary<string, SCSN_deg_line_class>();
        //    int other_length = other.DEGs.Length;
        //    SCSN_deg_line_class other_line;
        //    for (int indexO = 0; indexO < other_length; indexO++)
        //    {
        //        other_line = other.DEGs[indexO];
        //        if (!other_cluster_gene_degLine_dict.ContainsKey(other_line.Cluster))
        //        {
        //            other_cluster_gene_degLine_dict.Add(other_line.Cluster, new Dictionary<string, SCSN_deg_line_class>());
        //        }
        //        other_cluster_gene_degLine_dict[other_line.Cluster].Add(other_line.Gene, other_line);
        //    }
        //    int this_length = this.DEGs.Length;
        //    SCSN_avg_deg_line_class this_line;
        //    for (int indexThis = 0; indexThis < this_length; indexThis++)
        //    {
        //        this_line = this.DEGs[indexThis];
        //        if ((other_cluster_gene_degLine_dict.ContainsKey(this_line.Cluster))
        //            && (other_cluster_gene_degLine_dict[this_line.Cluster].ContainsKey(this_line.Gene)))
        //        {
        //            other_line = other_cluster_gene_degLine_dict[this_line.Cluster][this_line.Gene];
        //            if (!other_line.Cluster.Equals(this_line.Cluster)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }
        //            if (!other_line.Gene.Equals(this_line.Gene)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }

        //            this_line.Minus_log10_adj_pvalue_arithMean = mfs.Add_to_running_mean(this_line.Minus_log10_adj_pvalue_arithMean, other_line.Minus_log10_adj_pvalue, Anticipated_number_of_datasets_added);
        //            this_line.Minus_log10_adj_pvalue_squaredArithMean = mfs.Add_to_running_squared_mean(this_line.Minus_log10_adj_pvalue_squaredArithMean, other_line.Minus_log10_adj_pvalue, Anticipated_number_of_datasets_added);

        //            this_line.Avg_log2FC_arithMean = mfs.Add_to_running_mean(this_line.Avg_log2FC_arithMean, other_line.Avg_log2FC, Anticipated_number_of_datasets_added);
        //            this_line.Avg_log2FC_squaredArithMean = mfs.Add_to_running_squared_mean(this_line.Avg_log2FC_squaredArithMean, other_line.Avg_log2FC, Anticipated_number_of_datasets_added);

        //            other_cluster_gene_degLine_dict[this_line.Cluster][this_line.Gene] = new SCSN_deg_line_class();
        //        }
        //    }
        //    string[] other_clusters = other_cluster_gene_degLine_dict.Keys.ToArray();
        //    string other_cluster;
        //    int other_clusters_length = other_clusters.Length;
        //    string[] other_genes;
        //    string other_gene;
        //    int other_genes_length;
        //    SCSN_avg_deg_line_class add_this_line;
        //    List<SCSN_avg_deg_line_class> add = new List<SCSN_avg_deg_line_class>();
        //    for (int indexOC = 0; indexOC < other_clusters_length; indexOC++)
        //    {
        //        other_cluster = other_clusters[indexOC];
        //        other_gene_degLine_dict = other_cluster_gene_degLine_dict[other_cluster];
        //        other_genes = other_gene_degLine_dict.Keys.ToArray();
        //        other_genes_length = other_genes.Length;
        //        for (int indexOG = 0; indexOG < other_genes_length; indexOG++)
        //        {
        //            other_gene = other_genes[indexOG];
        //            other_line = other_gene_degLine_dict[other_gene];
        //            if (!String.IsNullOrEmpty(other_line.Cluster))
        //            {
        //                add_this_line = new SCSN_avg_deg_line_class();
        //                add_this_line.Cluster = (string)other_line.Cluster.Clone();
        //                add_this_line.Gene = (string)other_line.Gene.Clone();

        //                add_this_line.Minus_log10_adj_pvalue_arithMean = mfs.Add_to_running_mean(add_this_line.Minus_log10_adj_pvalue_arithMean, other_line.Minus_log10_adj_pvalue, Anticipated_number_of_datasets_added);
        //                add_this_line.Minus_log10_adj_pvalue_squaredArithMean = mfs.Add_to_running_squared_mean(add_this_line.Minus_log10_adj_pvalue_squaredArithMean, other_line.Minus_log10_adj_pvalue, Anticipated_number_of_datasets_added);

        //                add_this_line.Avg_log2FC_arithMean = mfs.Add_to_running_mean(add_this_line.Avg_log2FC_arithMean, other_line.Avg_log2FC, Anticipated_number_of_datasets_added);
        //                add_this_line.Avg_log2FC_squaredArithMean = mfs.Add_to_running_squared_mean(add_this_line.Avg_log2FC_squaredArithMean, other_line.Avg_log2FC, Anticipated_number_of_datasets_added);

        //                add.Add(add_this_line);
        //            }
        //        }
        //    }
        //    this.Add_to_array(add.ToArray());
        //}
        public void Add_other_dataset_by_summing_up_running_mean(SCSN_deg_class other)
        {
            Math_class mfs = new Math_class();
            if (Number_of_datasets_added_for_mean == Anticipated_number_of_datasets_added) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running"); }
            if (Averaging_procedure_finished) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean"); }
            if ((DEGs.Length != 0) && (Number_of_datasets_added_for_mean == 0)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean"); }
            //if ((DEGs.Length == 0) && (Number_of_datasets_added_for_mean != 0)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean"); }
            Number_of_datasets_added_for_mean++;
            Dictionary<string, Dictionary<string, SCSN_deg_line_class>> other_cluster_gene_degLine_dict = new Dictionary<string, Dictionary<string, SCSN_deg_line_class>>();
            Dictionary<string, SCSN_deg_line_class> other_gene_degLine_dict = new Dictionary<string, SCSN_deg_line_class>();
            int other_length = other.DEGs.Length;
            SCSN_deg_line_class other_line;
            for (int indexO = 0; indexO < other_length; indexO++)
            {
                other_line = other.DEGs[indexO];
                if (!other_cluster_gene_degLine_dict.ContainsKey(other_line.Cluster))
                {
                    other_cluster_gene_degLine_dict.Add(other_line.Cluster, new Dictionary<string, SCSN_deg_line_class>());
                }
                other_cluster_gene_degLine_dict[other_line.Cluster].Add(other_line.Gene, other_line);
            }
            int this_length = Get_this_length_filled_with_data();
            SCSN_avgAnalysesWithinEachIteration_deg_line_class this_line;
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                this_line = this.DEGs[indexThis];
                if ((other_cluster_gene_degLine_dict.ContainsKey(this_line.Cluster))
                    && (other_cluster_gene_degLine_dict[this_line.Cluster].ContainsKey(this_line.Gene)))
                {
                    other_line = other_cluster_gene_degLine_dict[this_line.Cluster][this_line.Gene];
                    if (!other_line.Cluster.Equals(this_line.Cluster)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }
                    if (!other_line.Gene.Equals(this_line.Gene)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_mean_and_variances"); }

                    this_line.Minus_log10_adj_pvalue_arithMean = mfs.Add_to_running_mean(this_line.Minus_log10_adj_pvalue_arithMean, other_line.Minus_log10_adj_pvalue, Anticipated_number_of_datasets_added);
                    this_line.Avg_log2FC_arithMean = mfs.Add_to_running_mean(this_line.Avg_log2FC_arithMean, other_line.Avg_log2FC, Anticipated_number_of_datasets_added);
                    this_line.Rank_arithMean = mfs.Add_to_running_mean(this_line.Rank_arithMean, other_line.Rank, Anticipated_number_of_datasets_added);
                    this_line.IsSig_arithMean = mfs.Add_to_running_mean(this_line.IsSig_arithMean, other_line.Is_significant, Anticipated_number_of_datasets_added);
                    other_cluster_gene_degLine_dict[this_line.Cluster][this_line.Gene] = new SCSN_deg_line_class();
                }
            }
            string[] other_clusters = other_cluster_gene_degLine_dict.Keys.ToArray();
            string other_cluster;
            int other_clusters_length = other_clusters.Length;
            string[] other_genes;
            string other_gene;
            int other_genes_length;
            List<SCSN_avgAnalysesWithinEachIteration_deg_line_class> add = new List<SCSN_avgAnalysesWithinEachIteration_deg_line_class>();
            for (int indexOC = 0; indexOC < other_clusters_length; indexOC++)
            {
                other_cluster = other_clusters[indexOC];
                other_gene_degLine_dict = other_cluster_gene_degLine_dict[other_cluster];
                other_genes = other_gene_degLine_dict.Keys.ToArray();
                other_genes_length = other_genes.Length;
                for (int indexOG = 0; indexOG < other_genes_length; indexOG++)
                {
                    other_gene = other_genes[indexOG];
                    other_line = other_gene_degLine_dict[other_gene];
                    if (!String.IsNullOrEmpty(other_line.Cluster))
                    {
                        this_line = new SCSN_avgAnalysesWithinEachIteration_deg_line_class();
                        this_line.Cluster = (string)other_line.Cluster.Clone();
                        this_line.Gene = (string)other_line.Gene.Clone();
                        
                        this_line.Minus_log10_adj_pvalue_arithMean = mfs.Add_to_running_mean(this_line.Minus_log10_adj_pvalue_arithMean, other_line.Minus_log10_adj_pvalue, Anticipated_number_of_datasets_added);
                        this_line.Avg_log2FC_arithMean = mfs.Add_to_running_mean(this_line.Avg_log2FC_arithMean, other_line.Avg_log2FC, Anticipated_number_of_datasets_added);
                        this_line.Rank_arithMean = mfs.Add_to_running_mean(this_line.Rank_arithMean, other_line.Rank, Anticipated_number_of_datasets_added);
                        this_line.IsSig_arithMean = mfs.Add_to_running_mean(this_line.IsSig_arithMean, other_line.Is_significant, Anticipated_number_of_datasets_added);
                        add.Add(this_line);
                    }
                }
            }
            //// All Genes from all clusters will be added before population variance is calcualted. If a gene is missing it will technically be added with values 0. This works for everything, except for the ranks
            this.Add_to_array_using_indexNextAvailableThis(add.ToArray());
            Check_for_duplicates();
        }
        public void Add_other_dataset_by_summing_up_running_variances_after_mean_is_final(SCSN_deg_class other)
        {
            if (Number_of_datasets_added_for_mean != Anticipated_number_of_datasets_added) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_variances_after_mean_is_final"); }
            if (Averaging_procedure_finished) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_variances_after_mean_is_final"); }
           // if (DEGs.Length == 0) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_variances_after_mean_is_final"); }
            Number_of_datasets_added_for_populationVariance++;
            Dictionary<string, Dictionary<string, SCSN_deg_line_class>> other_cluster_gene_degLine_dict = new Dictionary<string, Dictionary<string, SCSN_deg_line_class>>();
            Dictionary<string, SCSN_deg_line_class> other_gene_degLine_dict = new Dictionary<string, SCSN_deg_line_class>();
            int other_length = other.DEGs.Length;
            SCSN_deg_line_class other_line;
            for (int indexO = 0; indexO < other_length; indexO++)
            {
                other_line = other.DEGs[indexO];
                if (!other_cluster_gene_degLine_dict.ContainsKey(other_line.Cluster))
                {
                    other_cluster_gene_degLine_dict.Add(other_line.Cluster, new Dictionary<string, SCSN_deg_line_class>());
                }
                other_cluster_gene_degLine_dict[other_line.Cluster].Add(other_line.Gene, other_line);
            }
            int this_length = this.DEGs.Length;
            SCSN_avgAnalysesWithinEachIteration_deg_line_class this_line;
            Math_class mfs = new Math_class();
            other_line = new SCSN_deg_line_class();
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                this_line = this.DEGs[indexThis];
                if ((other_cluster_gene_degLine_dict.ContainsKey(this_line.Cluster))
                    && (other_cluster_gene_degLine_dict[this_line.Cluster].ContainsKey(this_line.Gene)))
                {
                    other_line = other_cluster_gene_degLine_dict[this_line.Cluster][this_line.Gene];
                    if (!other_line.Cluster.Equals(this_line.Cluster)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_variances_after_mean_is_final"); }
                    if (!other_line.Gene.Equals(this_line.Gene)) { throw new Exception("Exception in: Add_other_dataset_by_summing_up_running_variances_after_mean_is_final"); }
                }
                else
                {
                    other_line = new SCSN_deg_line_class();
                    other_line.Minus_log10_adj_pvalue = 0;
                    other_line.Avg_log2FC = 0;
                    other_line.Rank = 9999; //Results for ranks are not ok, since missing genes will be not be added to mean which means adding a rank 0
                    other_line.Is_significant = 0;
                }
                this_line.Minus_log10_adj_pvalue_popVariance = mfs.Add_to_running_population_variance(this_line.Minus_log10_adj_pvalue_popVariance, other_line.Minus_log10_adj_pvalue, this_line.Minus_log10_adj_pvalue_arithMean, Anticipated_number_of_datasets_added);
                this_line.Avg_log2FC_popVariance = mfs.Add_to_running_population_variance(this_line.Avg_log2FC_popVariance, other_line.Avg_log2FC, this_line.Avg_log2FC_arithMean, Anticipated_number_of_datasets_added);
                this_line.Rank_popVariance = mfs.Add_to_running_population_variance(this_line.Rank_popVariance, other_line.Rank, this_line.Rank_arithMean, Anticipated_number_of_datasets_added);
                this_line.IsSig_popVariance = mfs.Add_to_running_population_variance(this_line.IsSig_popVariance, other_line.Is_significant, this_line.IsSig_arithMean, Anticipated_number_of_datasets_added);
                if ((other_cluster_gene_degLine_dict.ContainsKey(this_line.Cluster))
                    && (other_cluster_gene_degLine_dict[this_line.Cluster].ContainsKey(this_line.Gene)))
                {
                    other_cluster_gene_degLine_dict[this_line.Cluster][this_line.Gene] = new SCSN_deg_line_class();
                }
            }
            string[] other_clusters = other_cluster_gene_degLine_dict.Keys.ToArray();
            string other_cluster;
            int other_clusters_length = other_clusters.Length;
            string[] other_genes;
            string other_gene;
            int other_genes_length;
            for (int indexOC = 0; indexOC < other_clusters_length; indexOC++)
            {
                other_cluster = other_clusters[indexOC];
                other_gene_degLine_dict = other_cluster_gene_degLine_dict[other_cluster];
                other_genes = other_gene_degLine_dict.Keys.ToArray();
                other_genes_length = other_genes.Length;
                for (int indexOG = 0; indexOG < other_genes_length; indexOG++)
                {
                    other_gene = other_genes[indexOG];
                    other_line = other_gene_degLine_dict[other_gene];
                    if (!String.IsNullOrEmpty(other_line.Cluster))
                    {
                        throw new Exception();
                    }
                }
            }
        }
        public void Calculate_pValAdjGeomMean_from_minusLog10AdjPvalueArithMean()
        {
            foreach (SCSN_avgAnalysesWithinEachIteration_deg_line_class avg_deg_line in this.DEGs)
            {
                avg_deg_line.P_val_adj_geomMean = Math.Pow(10, -avg_deg_line.Minus_log10_adj_pvalue_arithMean);
            }
        }
        public void Calculate_fractional_ranks_using_minusLog10AdjPvalueArithMean_and_avgLog2FcArithMean()
        {
            int degs_length = this.DEGs.Length;
            SCSN_avgAnalysesWithinEachIteration_deg_line_class avg_line;
            SCSN_avgAnalysesWithinEachIteration_deg_line_class inner_avg_line;
            this.DEGs = this.DEGs.OrderBy(l => l.Iteration).ThenBy(l => l.Cluster).ThenByDescending(l => l.Minus_log10_adj_pvalue_arithMean).ThenByDescending(l => Math.Abs(l.Avg_log2FC_arithMean)).ToArray();
            List<float> current_ranks = new List<float>();
            float running_rank = -1;
            float add_rank;
            int firstIndex_sameGroup = -1;
            Dictionary<string, bool> gene_dict = new Dictionary<string, bool>();
            Math_class mf = new Math_class();
            for (int indexDeg=0; indexDeg< degs_length; indexDeg++)
            {
                avg_line = this.DEGs[indexDeg];
                if (  (indexDeg==0)
                    || (!avg_line.Iteration.Equals(this.DEGs[indexDeg - 1].Iteration))
                    || (!avg_line.Cluster.Equals(this.DEGs[indexDeg - 1].Cluster)))
                {
                    running_rank = 0;
                    gene_dict.Clear();
                }
                gene_dict.Add(avg_line.Gene, true);//Check there are no duplicated genes
                if ((indexDeg == 0)
                    || (!avg_line.Iteration.Equals(this.DEGs[indexDeg - 1].Iteration))
                    || (!avg_line.Cluster.Equals(this.DEGs[indexDeg - 1].Cluster))
                    || (!avg_line.Minus_log10_adj_pvalue_arithMean.Equals(this.DEGs[indexDeg - 1].Minus_log10_adj_pvalue_arithMean))
                    || (!Math.Abs(avg_line.Avg_log2FC_arithMean).Equals(Math.Abs(this.DEGs[indexDeg - 1].Avg_log2FC_arithMean))))
                {
                    current_ranks.Clear();
                    firstIndex_sameGroup = indexDeg;
                }
                running_rank++;
                current_ranks.Add(running_rank);
                if ((indexDeg == degs_length-1)
                    || (!avg_line.Iteration.Equals(this.DEGs[indexDeg + 1].Iteration))
                    || (!avg_line.Cluster.Equals(this.DEGs[indexDeg + 1].Cluster))
                    || (!avg_line.Minus_log10_adj_pvalue_arithMean.Equals(this.DEGs[indexDeg + 1].Minus_log10_adj_pvalue_arithMean))
                    || (!Math.Abs(avg_line.Avg_log2FC_arithMean).Equals(Math.Abs(this.DEGs[indexDeg + 1].Avg_log2FC_arithMean))))
                {
                    add_rank = mf.Get_average(current_ranks.ToArray());
                    for (int indexInner = firstIndex_sameGroup; indexInner<=indexDeg; indexInner++)
                    {
                        inner_avg_line = this.DEGs[indexInner];
                        inner_avg_line.Rank = add_rank;
                    }
                }
            }
        }
        public void Calculate_sampleSD_and_finish_averaging_procedure()
        {
            AvgValueType_of_interest_enum[] avgValueOfInterest_array = new AvgValueType_of_interest_enum[] { AvgValueType_of_interest_enum.Minus_log10_adj_pvalue, AvgValueType_of_interest_enum.Rank, AvgValueType_of_interest_enum.Is_sig, AvgValueType_of_interest_enum.Avg_log2_fc };

            if (Averaging_procedure_finished) { throw new Exception("Error in Calculate_sampleSD"); }
            Averaging_procedure_finished = true;
            if (!Number_of_datasets_added_for_mean.Equals(Number_of_datasets_added_for_populationVariance)) { throw new Exception("Error in Calculate_sampleSD"); }
            if (!Number_of_datasets_added_for_mean.Equals(Anticipated_number_of_datasets_added)) { throw new Exception("Error in Calculate_sampleSD"); }
            double populationVariance;
            double sampleSD;
            int number_of_drawn_datasets_length;
            Math_class mf = new Math_class();
            foreach (SCSN_avgAnalysesWithinEachIteration_deg_line_class avg_line in this.DEGs)
            {
                avg_line.Number_of_datasets = Number_of_datasets_added_for_mean;
                foreach (AvgValueType_of_interest_enum avgValueOfInterst in avgValueOfInterest_array)
                {
                    populationVariance = avg_line.Get_popVariance_for_value_of_interest(avgValueOfInterst);
                    number_of_drawn_datasets_length = avg_line.Number_of_datasets;
                    sampleSD = mf.Calculate_sampleSD_from_population_variance(populationVariance, number_of_drawn_datasets_length);
                    avg_line.Set_sampleSD_for_value_of_interest(sampleSD, avgValueOfInterst);
                }
            }
        }
        public void Set_averaging_procedure_finished_to_true_if_no_sampleSD_calculated()
        {
            if (Averaging_procedure_finished) { throw new Exception("Error in Set_averaging_procedure_finished_to_true_if_no_sampleSD_calculated"); }
            Averaging_procedure_finished = true;
            if (!Number_of_datasets_added_for_populationVariance.Equals(0)) { throw new Exception("Error in Set_averaging_procedure_finished_to_true_if_no_sampleSD_calculated"); }
            if (!Number_of_datasets_added_for_mean.Equals(Anticipated_number_of_datasets_added)) { throw new Exception("Error in Set_averaging_procedure_finished_to_true_if_no_sampleSD_calculated"); }
            foreach (SCSN_avgAnalysesWithinEachIteration_deg_line_class avg_line in this.DEGs)
            {
                avg_line.Number_of_datasets = Number_of_datasets_added_for_mean;
            }
        }
        public void Keep_only_lines_with_isSigArtihMean_larger_zero()
        {
            List<SCSN_avgAnalysesWithinEachIteration_deg_line_class> keep = new List<SCSN_avgAnalysesWithinEachIteration_deg_line_class>();
            foreach (SCSN_avgAnalysesWithinEachIteration_deg_line_class avg_line in this.DEGs)
            {
                if (avg_line.IsSig_arithMean>0)
                {
                    keep.Add(avg_line);
                }
            }
            this.DEGs = keep.ToArray();
        }
        public SCSN_avgAnalysesWithinEachIteration_deg_class Deep_copy()
        {
            SCSN_avgAnalysesWithinEachIteration_deg_class copy = (SCSN_avgAnalysesWithinEachIteration_deg_class)this.MemberwiseClone();
            int degs_length = this.DEGs.Length;
            copy.DEGs = new SCSN_avgAnalysesWithinEachIteration_deg_line_class[degs_length];
            for (int indexDeg=0; indexDeg<degs_length;indexDeg++)
            {
                copy.DEGs[indexDeg] = this.DEGs[indexDeg].Deep_copy();
            }
            return copy;
        }
        public void Write(string directory, string fileName)
        {
            Check_for_correctness();
            if (!Averaging_procedure_finished) { throw new Exception("Averaging procedure not finsihed, no averaged data."); }
            SCSN_avgAnalysesWithinEachIteration_deg_readWrite_options_class readWriteOptions = new SCSN_avgAnalysesWithinEachIteration_deg_readWrite_options_class(directory, fileName);
            ReadWriteClass.WriteData(this.DEGs,readWriteOptions);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class SCSN_avgWithinEachNumberOfDatasets_deg_line_class
    {
        #region Fields
        public string Cluster { get; set; }
        public string Gene { get; set; }
        public double Minus_log10_adj_pvalue_arithMean { get; set; }
        public double Minus_log10_adj_pvalue_popVariance { get; set; }
        public double Minus_log10_adj_pvalue_sampleSd { get; set; }
        public double Minus_log10_adj_pvalue_coefVar { get; set; }
        public double Avg_log2FC_arithMean { get; set; }
        public double Avg_log2FC_sampleSd { get; set; }
        public double Avg_log2FC_coefVar { get; set; }
        public double Rank_arithMean { get; set; }
        public double Rank_sampleSd { get; set; }
        public double Rank_coefVar { get; set; }
        public double IsSig_arithMean { get; set; }
        public double IsSig_sampleSd { get; set; }
        public double IsSig_coefVar { get; set; }
        public double P_val_adj_geomMean { get; set; }
        public int Number_of_datasets { get; set; }
        public float Rank { get; set; }
        public int Number_of_iterations_for_same_number_of_datasets { get; set; }
        #endregion

        public SCSN_avgWithinEachNumberOfDatasets_deg_line_class()
        {
            Cluster = "";
            Gene = "";
            Number_of_iterations_for_same_number_of_datasets = 0;
            Minus_log10_adj_pvalue_arithMean = 0;
            Minus_log10_adj_pvalue_popVariance = 0;
            Avg_log2FC_arithMean = 0;
            Rank_arithMean = 0;
            IsSig_arithMean = 0;
        }

        #region Set values for valuesTypeOfInterest
        public void Set_arthMean_for_valueTypeOfInterest(double set_value, AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    Minus_log10_adj_pvalue_arithMean = set_value;
                    break;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    Avg_log2FC_arithMean = set_value;
                    break;
                case AvgValueType_of_interest_enum.Is_sig:
                    IsSig_arithMean = set_value;
                    break;
                case AvgValueType_of_interest_enum.Rank:
                    Rank_arithMean = set_value;
                    break;
                default:
                    throw new Exception();
            }
        }
        public void Set_sampleSD_for_valueTypeOfInterest(double set_value, AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    Minus_log10_adj_pvalue_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    Avg_log2FC_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Is_sig:
                    IsSig_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Rank:
                    Rank_sampleSd = set_value;
                    break;
                default:
                    throw new Exception();
            }
        }
        public void Set_coeffOfVar_for_valueTypeOfInterest(double set_value, AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    Minus_log10_adj_pvalue_coefVar = set_value;
                    break;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    Avg_log2FC_coefVar = set_value;
                    break;
                case AvgValueType_of_interest_enum.Is_sig:
                    IsSig_coefVar = set_value;
                    break;
                case AvgValueType_of_interest_enum.Rank:
                    Rank_coefVar = set_value;
                    break;
                default:
                    throw new Exception();
            }
        }
        public void Set_statisticalMeasure_for_valueTypeOfInterest(double set_value, AvgValueType_of_interest_enum valueType_of_interest, Statistical_measure_enum statistical_measure)
        {
            switch (statistical_measure)
            {
                case Statistical_measure_enum.Coefficient_of_variation:
                    Set_coeffOfVar_for_valueTypeOfInterest(set_value, valueType_of_interest);
                    break;
                default:
                    throw new Exception();
            }
        }

        #endregion

        #region Get values for valuesTypeOfInterest
        public double Get_statisticalMeasure_for_valueTypeOfInterest(AvgValueType_of_interest_enum valueType_of_interest, Statistical_measure_enum statistical_measure)
        {
            switch (statistical_measure)
            {
                case Statistical_measure_enum.Coefficient_of_variation:
                    return Get_coeffOfVar_for_valueTypeOfInterest(valueType_of_interest);
                default:
                    throw new Exception();
            }
        }
        public double Get_arthMean_for_valueTypeOfInterest(AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_arithMean;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_arithMean;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_arithMean;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_arithMean;
                default:
                    throw new Exception();
            }
        }
        public double Get_sampleSD_for_valueTypeOfInterest(AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_sampleSd;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_sampleSd;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_sampleSd;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_sampleSd;
                default:
                    throw new Exception();
            }
        }
        public double Get_coeffOfVar_for_valueTypeOfInterest(AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_coefVar;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_coefVar;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_coefVar;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_coefVar;
                default:
                    throw new Exception();
            }
        }
        #endregion

        public SCSN_avgWithinEachNumberOfDatasets_deg_line_class Deep_copy()
        {
            SCSN_avgWithinEachNumberOfDatasets_deg_line_class copy = (SCSN_avgWithinEachNumberOfDatasets_deg_line_class)MemberwiseClone();
            copy.Cluster = (string)this.Cluster.Clone();
            copy.Gene = (string)this.Gene.Clone();
            return copy;
        }

    }
    class SCSN_avgWithinEachNumberOfDatasets_deg_readWriteOptions_class : ReadWriteOptions_base
    {
        public SCSN_avgWithinEachNumberOfDatasets_deg_readWriteOptions_class(string directory, string fileName, bool data_is_normalized)
        {
            if (data_is_normalized) {  fileName = Path.GetFileNameWithoutExtension(fileName) + "_norm" + Path.GetExtension(fileName); }
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Number_of_datasets", "Number_of_iterations_for_same_number_of_datasets", "Cluster","Gene",
                                                    "Minus_log10_adj_pvalue_arithMean","Minus_log10_adj_pvalue_sampleSd","Minus_log10_adj_pvalue_coefVar",
                                                    "Avg_log2FC_arithMean","Avg_log2FC_sampleSd","Avg_log2FC_coefVar", 
                                                    "Rank_arithMean", "Rank_sampleSd", "Rank_coefVar",
                                                    "IsSig_arithMean", "IsSig_sampleSd", "IsSig_coefVar",
                                                    "P_val_adj_geomMean","Rank" };
            this.Key_columnNames = this.Key_propertyNames;
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
            this.Report_unhandled_null_entries = true;
        }
    }
    class SCSN_avgWithinEachNumberOfDatasets_deg_options_class
    {
        public bool Consider_only_genes_significant_in_final_degs { get; set; }
        public float Max_adj_pvalue { get; set; }
        public float Max_rank { get; set; }
        public SCSN_avgWithinEachNumberOfDatasets_deg_options_class(bool consider_only_genes_significant_in_final_degs, float max_adj_pavlue, float max_rank)
        {
            this.Consider_only_genes_significant_in_final_degs = consider_only_genes_significant_in_final_degs;
            this.Max_adj_pvalue = max_adj_pavlue;
            this.Max_rank = max_rank;
        }
    }
    class SCSN_avgWithinEachNumberOfDatasets_deg_class
    {
        public SCSN_avgWithinEachNumberOfDatasets_deg_line_class[] DEGs { get; set; }
        public bool Data_is_normalized_with_highest_datasetNumber_as_reference { get; set; }
        public SCSN_avgWithinEachNumberOfDatasets_deg_options_class Options { get; set; }

        public SCSN_avgWithinEachNumberOfDatasets_deg_class(bool consider_only_genes_significant_in_final_degs, float max_adj_pavlue, float max_rank)
        {
            Data_is_normalized_with_highest_datasetNumber_as_reference = false;
            this.DEGs = new SCSN_avgWithinEachNumberOfDatasets_deg_line_class[0];
            Options = new SCSN_avgWithinEachNumberOfDatasets_deg_options_class(consider_only_genes_significant_in_final_degs, max_adj_pavlue, max_rank);
        }
        private void Check_for_duplicates()
        {
            Dictionary<int, Dictionary<string, Dictionary<string, bool>>> numberOfDatasets_cluster_gene_dict = new Dictionary<int, Dictionary<string, Dictionary<string, bool>>>();
            foreach (SCSN_avgWithinEachNumberOfDatasets_deg_line_class avg_line in this.DEGs)
            {
                if (!numberOfDatasets_cluster_gene_dict.ContainsKey(avg_line.Number_of_datasets))
                {
                    numberOfDatasets_cluster_gene_dict.Add(avg_line.Number_of_datasets, new Dictionary<string, Dictionary<string, bool>>());
                }
                if (!numberOfDatasets_cluster_gene_dict[avg_line.Number_of_datasets].ContainsKey(avg_line.Cluster))
                {
                    numberOfDatasets_cluster_gene_dict[avg_line.Number_of_datasets].Add(avg_line.Cluster, new Dictionary<string, bool>());
                }
                numberOfDatasets_cluster_gene_dict[avg_line.Number_of_datasets][avg_line.Cluster].Add(avg_line.Gene, true);
            }
        }
        private void Add_to_array_and_check_if_number_of_datasets_does_not_already_exist(SCSN_avgWithinEachNumberOfDatasets_deg_line_class[] add_degs)
        {
            int this_length = this.DEGs.Length;
            int add_length = add_degs.Length;
            int new_length = this_length + add_length;
            int indexNew = -1;
            SCSN_avgWithinEachNumberOfDatasets_deg_line_class[] new_degs = new SCSN_avgWithinEachNumberOfDatasets_deg_line_class[new_length];
            SCSN_avgWithinEachNumberOfDatasets_deg_line_class thisAdd_line;
            Dictionary<int, bool> numberOfConsideredDatasets_dict = new Dictionary<int, bool>();
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                thisAdd_line = this.DEGs[indexThis];
                new_degs[indexNew] = thisAdd_line;
                if (!numberOfConsideredDatasets_dict.ContainsKey(thisAdd_line.Number_of_datasets))
                {
                    numberOfConsideredDatasets_dict.Add(thisAdd_line.Number_of_datasets, true);
                }
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                thisAdd_line = add_degs[indexAdd];
                new_degs[indexNew] = thisAdd_line;
                if (numberOfConsideredDatasets_dict.ContainsKey(thisAdd_line.Number_of_datasets)) { throw new Exception("Number of considered datasets already exists in avgWithinEachNumberOfDatasets"); }
            }
            this.DEGs = new_degs;
        }
        public void Add_other(SCSN_avgWithinEachNumberOfDatasets_deg_class other_degs)
        {
            if (Data_is_normalized_with_highest_datasetNumber_as_reference) { throw new Exception("Add other: data already normalized"); }
            Add_to_array_and_check_if_number_of_datasets_does_not_already_exist(other_degs.DEGs);
        }
        public void Calculate_mean_and_sampleSD_of_the_artihMeans_for_allValueTypesOfInterest_across_all_iterations_assuming_same_number_of_drawn_datasets(SCSN_avgAnalysesWithinEachIteration_deg_class scsn_avg, SCSN_avgAnalysesWithinEachIteration_deg_class final_scsn_avg)
        {
            scsn_avg.Check_for_correctness();
            if (Data_is_normalized_with_highest_datasetNumber_as_reference) { throw new Exception("Exception in: Add_average_values_from_dataset_by_summing_up_running_mean_assuming_same_number_of_drawn_datasets"); }
            if (!scsn_avg.Averaging_procedure_finished) { throw new Exception("Exception in: Add_average_values_from_dataset_by_summing_up_running_mean_assuming_same_number_of_drawn_datasets"); }
            int number_of_drawn_datasets = scsn_avg.DEGs[0].Number_of_datasets;
            Dictionary<int, bool> iterationNos_dict = new Dictionary<int, bool>();
            foreach (SCSN_avgAnalysesWithinEachIteration_deg_line_class deg_line in scsn_avg.DEGs)
            {
                if (deg_line.Number_of_datasets != number_of_drawn_datasets) { throw new Exception("Error in Add_average_values_from_dataset_by_summing_up_running_mean_assuming_same_number_of_drawn_datasets"); }
                if (!iterationNos_dict.ContainsKey(deg_line.Iteration)) { iterationNos_dict.Add(deg_line.Iteration, true); }
            }
            int iterations_length = iterationNos_dict.Keys.Count;
            double one_div_by_iterations_length = 1.0 / (double)iterations_length;

            Math_class mfs = new Math_class();

            int number_of_datasets = scsn_avg.DEGs[0].Number_of_datasets;

            Dictionary<string, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>> cluster_gene_valueType_runningMean_dict = new Dictionary<string, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>>();
            Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>> gene_valueType_runningMean_dict = new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>();
            Dictionary<AvgValueType_of_interest_enum, double> valueType_runningMean_dict = new Dictionary<AvgValueType_of_interest_enum, double>();
            double runningMean;

            Dictionary<string, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>> cluster_gene_valueType_runningPopVariance_dict = new Dictionary<string, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>>();
            Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>> gene_valueType_runningPopVariance_dict = new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>();
            Dictionary<AvgValueType_of_interest_enum, double> valueType_runningPopVariance_dict = new Dictionary<AvgValueType_of_interest_enum, double>();
            double runningPopVariance;


            int degs_length;
            SCSN_avgAnalysesWithinEachIteration_deg_line_class scsn_avg_line;

            AvgValueType_of_interest_enum[] avgValueType_of_interest_array = new AvgValueType_of_interest_enum[] { AvgValueType_of_interest_enum.Minus_log10_adj_pvalue, AvgValueType_of_interest_enum.Avg_log2_fc, AvgValueType_of_interest_enum.Rank, AvgValueType_of_interest_enum.Is_sig };
            AvgValueType_of_interest_enum avgValueType_of_interest;
            int avgValueType_of_interest_array_length = avgValueType_of_interest_array.Length;

            Dictionary<string, Dictionary<string, bool>> cluster_gene_consider_dict = new Dictionary<string, Dictionary<string, bool>>();

            #region Fill cluster gene consider dict
            if (Options.Consider_only_genes_significant_in_final_degs)
            {
                int final_degs_length = final_scsn_avg.DEGs.Length;
                for (int indexFinal = 0; indexFinal < final_degs_length; indexFinal++)
                {
                    scsn_avg_line = final_scsn_avg.DEGs[indexFinal];
                    if (  (scsn_avg_line.P_val_adj_geomMean <= Options.Max_adj_pvalue)
                        &&(scsn_avg_line.Rank <= Options.Max_rank))
                    {
                        if (!cluster_gene_consider_dict.ContainsKey(scsn_avg_line.Cluster))
                        {
                            cluster_gene_consider_dict.Add(scsn_avg_line.Cluster, new Dictionary<string, bool>());
                        }
                        if (!cluster_gene_consider_dict[scsn_avg_line.Cluster].ContainsKey(scsn_avg_line.Gene))
                        {
                            cluster_gene_consider_dict[scsn_avg_line.Cluster].Add(scsn_avg_line.Gene, true);
                        }
                    }
                }
            }
            else
            {
                degs_length = scsn_avg.DEGs.Length;
                for (int indexDEG = 0; indexDEG < degs_length; indexDEG++)
                {
                    scsn_avg_line = scsn_avg.DEGs[indexDEG];

                    if (!cluster_gene_consider_dict.ContainsKey(scsn_avg_line.Cluster))
                    {
                        cluster_gene_consider_dict.Add(scsn_avg_line.Cluster, new Dictionary<string, bool>());
                    }
                    if (!cluster_gene_consider_dict[scsn_avg_line.Cluster].ContainsKey(scsn_avg_line.Gene))
                    {
                        cluster_gene_consider_dict[scsn_avg_line.Cluster].Add(scsn_avg_line.Gene, true);
                    }
                }
            }
            #endregion

            Dictionary<string, Dictionary<string, List<SCSN_avgAnalysesWithinEachIteration_deg_line_class>>> cluster_gene_avgLines_dict = new Dictionary<string, Dictionary<string, List<SCSN_avgAnalysesWithinEachIteration_deg_line_class>>>();
            Dictionary<string, List<SCSN_avgAnalysesWithinEachIteration_deg_line_class>> gene_avgLines_dict;

            #region Calcuate means by adding to running means and fill cluster_gene_avgLine_dict and iteration_dict
            double add_value;
            degs_length = scsn_avg.DEGs.Length;
            for (int indexDEG = 0; indexDEG < degs_length; indexDEG++)
            {
                scsn_avg_line = scsn_avg.DEGs[indexDEG];
                if (cluster_gene_consider_dict[scsn_avg_line.Cluster].ContainsKey(scsn_avg_line.Gene))
                {
                    #region Fill cluster_gene_avgLines_dict and iteration_dict
                    if (!cluster_gene_avgLines_dict.ContainsKey(scsn_avg_line.Cluster))
                    {
                        cluster_gene_avgLines_dict.Add(scsn_avg_line.Cluster, new Dictionary<string, List<SCSN_avgAnalysesWithinEachIteration_deg_line_class>>());
                    }
                    if (!cluster_gene_avgLines_dict[scsn_avg_line.Cluster].ContainsKey(scsn_avg_line.Gene))
                    {
                        cluster_gene_avgLines_dict[scsn_avg_line.Cluster].Add(scsn_avg_line.Gene, new List<SCSN_avgAnalysesWithinEachIteration_deg_line_class>());
                    }
                    cluster_gene_avgLines_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene].Add(scsn_avg_line);
                    #endregion

                    if (!cluster_gene_valueType_runningMean_dict.ContainsKey(scsn_avg_line.Cluster))
                    {
                        cluster_gene_valueType_runningMean_dict.Add(scsn_avg_line.Cluster, new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>());
                    }
                    if (!cluster_gene_valueType_runningMean_dict[scsn_avg_line.Cluster].ContainsKey(scsn_avg_line.Gene))
                    {
                        cluster_gene_valueType_runningMean_dict[scsn_avg_line.Cluster].Add(scsn_avg_line.Gene, new Dictionary<AvgValueType_of_interest_enum, double>());
                    }
                    for (int indexAvgVT = 0; indexAvgVT < avgValueType_of_interest_array_length; indexAvgVT++)
                    {
                        avgValueType_of_interest = avgValueType_of_interest_array[indexAvgVT];
                        if (!cluster_gene_valueType_runningMean_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene].ContainsKey(avgValueType_of_interest))
                        {
                            cluster_gene_valueType_runningMean_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene].Add(avgValueType_of_interest, 0);
                        }
                        add_value = scsn_avg_line.Get_arithMean_for_value_of_interest(avgValueType_of_interest);
                        cluster_gene_valueType_runningMean_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene][avgValueType_of_interest]
                            = mfs.Add_to_running_mean(cluster_gene_valueType_runningMean_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene][avgValueType_of_interest],
                                                      add_value,
                                                      iterations_length);
                    }
                }
            }
            #endregion

            #region Calcuate population variances by adding to running population variances
            string[] clusters = cluster_gene_avgLines_dict.Keys.ToArray();
            string cluster;
            int clusters_length = clusters.Length;
            string[] genes;
            string gene;
            int genes_length;
            List<SCSN_avgAnalysesWithinEachIteration_deg_line_class> avgLines;
            SCSN_avgAnalysesWithinEachIteration_deg_line_class avgLine_for_missing_gene;
            int avgLines_count;

            for (int indexCluster=0; indexCluster<clusters_length;indexCluster++)
            {
                cluster = clusters[indexCluster];
                gene_avgLines_dict = cluster_gene_avgLines_dict[cluster];
                genes = gene_avgLines_dict.Keys.ToArray();
                genes_length = genes.Length;
                for (int indexGene=0; indexGene<genes_length;indexGene++)
                {
                    gene = genes[indexGene];
                    avgLines = gene_avgLines_dict[gene];
                    avgLines_count = avgLines.Count;
                    for (int indexAdd=avgLines_count;indexAdd<iterations_length;indexAdd++)
                    {
                        avgLine_for_missing_gene = new SCSN_avgAnalysesWithinEachIteration_deg_line_class();
                        avgLine_for_missing_gene.Cluster = (string)cluster.Clone();
                        avgLine_for_missing_gene.Gene = (string)gene.Clone();
                        avgLine_for_missing_gene.Minus_log10_adj_pvalue_arithMean = 0;
                        avgLine_for_missing_gene.Rank_arithMean = 9999;
                        avgLine_for_missing_gene.IsSig_arithMean = 0;
                        avgLine_for_missing_gene.Avg_log2FC_arithMean = 0;
                        avgLines.Add(avgLine_for_missing_gene);
                    }
                    for (int indexAvgLine=0; indexAvgLine<avgLines_count;indexAvgLine++)
                    {
                        scsn_avg_line = avgLines[indexAvgLine];
                        if (!cluster_gene_valueType_runningPopVariance_dict.ContainsKey(scsn_avg_line.Cluster))
                        {
                            cluster_gene_valueType_runningPopVariance_dict.Add(scsn_avg_line.Cluster, new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>());
                        }
                        if (!cluster_gene_valueType_runningPopVariance_dict[scsn_avg_line.Cluster].ContainsKey(scsn_avg_line.Gene))
                        {
                            cluster_gene_valueType_runningPopVariance_dict[scsn_avg_line.Cluster].Add(scsn_avg_line.Gene, new Dictionary<AvgValueType_of_interest_enum, double>());
                        }
                        for (int indexAvgVT = 0; indexAvgVT < avgValueType_of_interest_array_length; indexAvgVT++)
                        {
                            avgValueType_of_interest = avgValueType_of_interest_array[indexAvgVT];
                            if (!cluster_gene_valueType_runningPopVariance_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene].ContainsKey(avgValueType_of_interest))
                            {
                                cluster_gene_valueType_runningPopVariance_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene].Add(avgValueType_of_interest, 0);
                            }
                            add_value = scsn_avg_line.Get_arithMean_for_value_of_interest(avgValueType_of_interest);
                            cluster_gene_valueType_runningPopVariance_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene][avgValueType_of_interest]
                                = mfs.Add_to_running_population_variance(cluster_gene_valueType_runningPopVariance_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene][avgValueType_of_interest],
                                                                         add_value,
                                                                         cluster_gene_valueType_runningMean_dict[scsn_avg_line.Cluster][scsn_avg_line.Gene][avgValueType_of_interest],
                                                                         iterations_length);
                        }
                    }

                }
            }
            #endregion

            #region Generate new avg_lines
            SCSN_avgWithinEachNumberOfDatasets_deg_line_class avg_line;
            List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class> avg_lines = new List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>();

            clusters = cluster_gene_valueType_runningMean_dict.Keys.ToArray();
            clusters_length = clusters.Length;
            double sampleSD;
            double coeffVar;

            for (int indexC=0; indexC<clusters_length;indexC++)
            {
                cluster= clusters[indexC];
                gene_valueType_runningMean_dict = cluster_gene_valueType_runningMean_dict[cluster];
                gene_valueType_runningPopVariance_dict = cluster_gene_valueType_runningPopVariance_dict[cluster];
                genes = gene_valueType_runningMean_dict.Keys.ToArray();
                genes_length = genes.Length;
                for (int indexG=0; indexG<genes_length;indexG++)
                {
                    gene = genes[indexG];
                    valueType_runningMean_dict = gene_valueType_runningMean_dict[gene];
                    valueType_runningPopVariance_dict = gene_valueType_runningPopVariance_dict[gene];
                    avg_line = new SCSN_avgWithinEachNumberOfDatasets_deg_line_class();
                    avg_line.Cluster = (string)cluster.Clone();
                    avg_line.Gene = (string)gene.Clone();
                    avg_line.Number_of_datasets = number_of_datasets;
                    avg_line.Number_of_iterations_for_same_number_of_datasets = iterations_length;
                    for (int indexAVT=0; indexAVT<avgValueType_of_interest_array_length;indexAVT++)
                    {
                        avgValueType_of_interest = avgValueType_of_interest_array[indexAVT];
                        runningMean = valueType_runningMean_dict[avgValueType_of_interest];
                        runningPopVariance = valueType_runningPopVariance_dict[avgValueType_of_interest];
                        sampleSD = mfs.Calculate_sampleSD_from_population_variance(runningPopVariance, iterations_length);
                        coeffVar = mfs.Calculate_coefficient_of_variance(runningMean, sampleSD);
                        avg_line.Set_arthMean_for_valueTypeOfInterest(runningMean, avgValueType_of_interest);
                        avg_line.Set_sampleSD_for_valueTypeOfInterest(sampleSD, avgValueType_of_interest);
                        avg_line.Set_coeffOfVar_for_valueTypeOfInterest(coeffVar, avgValueType_of_interest);
                    }
                    avg_lines.Add(avg_line);
                }
            }
            #endregion

            this.Add_to_array_and_check_if_number_of_datasets_does_not_already_exist(avg_lines.ToArray());
            Check_for_duplicates();
        }
        public void Normalize_values_with_reference_to_values_of_highest_number_of_datasets()
        {
            if (Data_is_normalized_with_highest_datasetNumber_as_reference) { throw new Exception("Data is alredy normalized with highest datasetNumber as reference"); }
            Data_is_normalized_with_highest_datasetNumber_as_reference = true;
            Dictionary<string, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>> cluster_gene_avgVTOfInterst_refValue_dict = new Dictionary<string, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>>();

            AvgValueType_of_interest_enum[] avgValueType_of_interest_array = new AvgValueType_of_interest_enum[] { AvgValueType_of_interest_enum.Minus_log10_adj_pvalue, AvgValueType_of_interest_enum.Avg_log2_fc };
            AvgValueType_of_interest_enum avgValueType_of_interest;
            int avgValueType_of_interest_array_length = avgValueType_of_interest_array.Length;

            int highest_number_of_datasets = -1;
            int data_length = this.DEGs.Length;
            SCSN_avgWithinEachNumberOfDatasets_deg_line_class avg_deg_line;
            for (int indexAvg = 0; indexAvg < data_length; indexAvg++)
            {
                avg_deg_line = this.DEGs[indexAvg];
                if (highest_number_of_datasets < avg_deg_line.Number_of_datasets) { highest_number_of_datasets = avg_deg_line.Number_of_datasets; }
            }
            double one_div_by_highest_number_of_datasets = 1.0 / (double)highest_number_of_datasets;
            for (int indexAvg = 0; indexAvg < data_length; indexAvg++)
            {
                avg_deg_line = this.DEGs[indexAvg];
                if (highest_number_of_datasets == avg_deg_line.Number_of_datasets)
                {
                    if (!cluster_gene_avgVTOfInterst_refValue_dict.ContainsKey(avg_deg_line.Cluster))
                    {
                        cluster_gene_avgVTOfInterst_refValue_dict.Add(avg_deg_line.Cluster, new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>());
                    }
                    if (!cluster_gene_avgVTOfInterst_refValue_dict[avg_deg_line.Cluster].ContainsKey(avg_deg_line.Gene))
                    {
                        cluster_gene_avgVTOfInterst_refValue_dict[avg_deg_line.Cluster].Add(avg_deg_line.Gene, new Dictionary<AvgValueType_of_interest_enum, double>());
                    }
                    for (int indexAVT=0; indexAVT<avgValueType_of_interest_array_length;indexAVT++)
                    {
                        avgValueType_of_interest = avgValueType_of_interest_array[indexAVT];
                        cluster_gene_avgVTOfInterst_refValue_dict[avg_deg_line.Cluster][avg_deg_line.Gene].Add(avgValueType_of_interest, avg_deg_line.Get_arthMean_for_valueTypeOfInterest(avgValueType_of_interest));
                    }
                }
            }
            double reference_value;
            double current_mean;
            double normalized_mean;
            double current_sampleSD;
            double normalized_sampleSD;
            Math_class mf = new Math_class();
            for (int indexAvg = 0; indexAvg < data_length; indexAvg++)
            {
                avg_deg_line = this.DEGs[indexAvg];
                for (int indexAVT = 0; indexAVT < avgValueType_of_interest_array_length; indexAVT++)
                {
                    avgValueType_of_interest = avgValueType_of_interest_array[indexAVT];
                    reference_value = cluster_gene_avgVTOfInterst_refValue_dict[avg_deg_line.Cluster][avg_deg_line.Gene][avgValueType_of_interest];

                    current_mean = avg_deg_line.Get_arthMean_for_valueTypeOfInterest(avgValueType_of_interest);
                    normalized_mean = mf.Normalize_value(current_mean, reference_value);
                    avg_deg_line.Set_arthMean_for_valueTypeOfInterest(normalized_mean, avgValueType_of_interest);

                    current_sampleSD = avg_deg_line.Get_sampleSD_for_valueTypeOfInterest(avgValueType_of_interest);
                    normalized_sampleSD = mf.Normalize_value(current_sampleSD, reference_value);
                    avg_deg_line.Set_sampleSD_for_valueTypeOfInterest(normalized_sampleSD, avgValueType_of_interest);
                }
            }
        }
        public void Keep_only_lines_with_isSigArithMean_larger_zero()
        {
            List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class> keep = new List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>();
            foreach (SCSN_avgWithinEachNumberOfDatasets_deg_line_class line in this.DEGs)
            {
                if (line.IsSig_arithMean>0)
                {
                    keep.Add(line);
                }
            }
            this.DEGs = keep.ToArray();
        }
        public SCSN_avgWithinEachNumberOfDatasets_deg_class Deep_copy()
        {
            SCSN_avgWithinEachNumberOfDatasets_deg_class copy = (SCSN_avgWithinEachNumberOfDatasets_deg_class)this.MemberwiseClone();
            int degs_length = this.DEGs.Length;
            copy.DEGs = new SCSN_avgWithinEachNumberOfDatasets_deg_line_class[degs_length];
            for (int indexDeg=0; indexDeg<degs_length; indexDeg++)
            {
                copy.DEGs[indexDeg] = this.DEGs[indexDeg].Deep_copy();
            }
            return copy;
        }
        public void Write(string directory, string fileName)
        {
            SCSN_avgWithinEachNumberOfDatasets_deg_readWriteOptions_class readWriteOptions = new SCSN_avgWithinEachNumberOfDatasets_deg_readWriteOptions_class(directory, fileName, Data_is_normalized_with_highest_datasetNumber_as_reference);
            ReadWriteClass.WriteData(this.DEGs, readWriteOptions);
        }

        public void Read(string directory, string fileName)
        {
            SCSN_avgWithinEachNumberOfDatasets_deg_readWriteOptions_class readWriteOptions = new SCSN_avgWithinEachNumberOfDatasets_deg_readWriteOptions_class(directory, fileName, Data_is_normalized_with_highest_datasetNumber_as_reference);
            this.DEGs = ReadWriteClass.ReadRawData_and_FillArray<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>(readWriteOptions);
        }

    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class SCSN_postHocPower_avg_line_class
    {
        public string Cluster { get; set; }
        public double Minus_log10_adj_pvalue_arithMean { get; set; }
        public double Minus_log10_adj_pvalue_sampleSd { get; set; }
        public double Minus_log10_adj_pvalue_coefVar { get; set; }
        public double Avg_log2FC_arithMean { get; set; }
        public double Avg_log2FC_sampleSd { get; set; }
        public double Avg_log2FC_coefVar { get; set; }
        public double Rank_arithMean { get; set; }
        public double Rank_sampleSd { get; set; }
        public double Rank_coefVar { get; set; }
        public double IsSig_arithMean { get; set; }
        public double IsSig_sampleSd { get; set; }
        public double IsSig_coefVar { get; set; }
        public int Number_of_datasets { get; set; }
        public int Number_of_genes { get; set; }
        public Statistical_measure_enum Source_statistical_measure { get; set; }
        public float Rank { get; set; }
        public bool Avg_data_was_normalized { get; set; }
        public SCSN_postHocPower_avg_line_class()
        {
            Cluster = "";
            Source_statistical_measure = Statistical_measure_enum.E_m_p_t_y;
        }
        #region Set values for valuesTypeOfInterest
        public void Set_arthMean_for_valueTypeOfInterest(double set_value, AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    Minus_log10_adj_pvalue_arithMean = set_value;
                    break;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    Avg_log2FC_arithMean = set_value;
                    break;
                case AvgValueType_of_interest_enum.Is_sig:
                    IsSig_arithMean = set_value;
                    break;
                case AvgValueType_of_interest_enum.Rank:
                    Rank_arithMean = set_value;
                    break;
                default:
                    throw new Exception();
            }
        }
        public void Set_sampleSD_for_valueTypeOfInterest(double set_value, AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    Minus_log10_adj_pvalue_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    Avg_log2FC_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Is_sig:
                    IsSig_sampleSd = set_value;
                    break;
                case AvgValueType_of_interest_enum.Rank:
                    Rank_sampleSd = set_value;
                    break;
                default:
                    throw new Exception();
            }
        }
        public void Set_coeffOfVar_for_valueTypeOfInterest(double set_value, AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    Minus_log10_adj_pvalue_coefVar = set_value;
                    break;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    Avg_log2FC_coefVar = set_value;
                    break;
                case AvgValueType_of_interest_enum.Is_sig:
                    IsSig_coefVar = set_value;
                    break;
                case AvgValueType_of_interest_enum.Rank:
                    Rank_coefVar = set_value;
                    break;
                default:
                    throw new Exception();
            }
        }
        #endregion

        #region Get values for valuesTypeOfInterest
        public double Get_arthMean_for_valueTypeOfInterest(AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_arithMean;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_arithMean;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_arithMean;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_arithMean;
                default:
                    throw new Exception();
            }
        }
        public double Get_sampleSD_for_valueTypeOfInterest(AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_sampleSd;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_sampleSd;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_sampleSd;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_sampleSd;
                default:
                    throw new Exception();
            }
        }
        public double Get_coeffOfVar_for_valueTypeOfInterest(AvgValueType_of_interest_enum valueType_of_interest)
        {
            switch (valueType_of_interest)
            {
                case AvgValueType_of_interest_enum.Minus_log10_adj_pvalue:
                    return Minus_log10_adj_pvalue_coefVar;
                case AvgValueType_of_interest_enum.Avg_log2_fc:
                    return Avg_log2FC_coefVar;
                case AvgValueType_of_interest_enum.Is_sig:
                    return IsSig_coefVar;
                case AvgValueType_of_interest_enum.Rank:
                    return Rank_coefVar;
                default:
                    throw new Exception();
            }
        }
        #endregion

        public SCSN_postHocPower_avg_line_class Deep_copy()
        {
            SCSN_postHocPower_avg_line_class copy = (SCSN_postHocPower_avg_line_class)MemberwiseClone();
            copy.Cluster = (string)this.Cluster.Clone();
            return copy;
        }
    }
    class SCSN_postHocPower_avg_options_class
    {
        public bool Only_considere_genes_that_were_significant_at_least_once { get; set; }
        public SCSN_postHocPower_avg_options_class()
        {
            Only_considere_genes_that_were_significant_at_least_once = true;
        }
    }
    class SCSN_postHocPower_avg_readWriteOptions : ReadWriteOptions_base
    {
        public SCSN_postHocPower_avg_readWriteOptions(string directory, string fileName, bool data_is_normalized)
        {
            if (data_is_normalized) { fileName = Path.GetFileNameWithoutExtension(fileName) + "_norm" + Path.GetExtension(fileName); }
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Cluster",
                                                    "Minus_log10_adj_pvalue_arithMean", "Minus_log10_adj_pvalue_sampleSd", "Minus_log10_adj_pvalue_coefVar",
                                                    "Rank_arithMean", "Rank_sampleSd", "Rank_coefVar",
                                                    "IsSig_arithMean", "IsSig_sampleSd", "IsSig_coefVar",
                                                    "Avg_log2FC_arithMean", "Avg_log2FC_sampleSd", "Avg_log2FC_coefVar",
                                                    "Avg_data_was_normalized","Source_statistical_measure",
                                                    "Number_of_datasets", "Number_of_genes" };
            this.Key_columnNames = this.Key_propertyNames;
            this.Delimiter = Global_class.Tab;
            this.File_has_headline = true;
            this.Report_unhandled_null_entries = true;
        }
    }
    class SCSN_postHocPower_avg_class
    {
        public SCSN_postHocPower_avg_options_class Options { get; set; }
        public SCSN_postHocPower_avg_line_class[] PostHocPower_lines { get; set; }
        public bool Data_was_normalized_before_postHocPowerResults_generation { get; set; }

        public SCSN_postHocPower_avg_class(bool consider_only_genes_that_were_significant_at_least_once)
        {
            this.PostHocPower_lines = new SCSN_postHocPower_avg_line_class[0];
            Data_was_normalized_before_postHocPowerResults_generation = false;
            this.Options = new SCSN_postHocPower_avg_options_class();
            this.Options.Only_considere_genes_that_were_significant_at_least_once = consider_only_genes_that_were_significant_at_least_once;
        }
        private void Add_to_array_and_check_if_number_of_added_datasets_do_not_exist(SCSN_postHocPower_avg_line_class[] add_postHocPower_lines)
        {
            int this_length = this.PostHocPower_lines.Length;
            int add_length = add_postHocPower_lines.Length;
            int new_length = this_length + add_length;
            SCSN_postHocPower_avg_line_class[] new_postHocPower_lines = new SCSN_postHocPower_avg_line_class[new_length];
            int indexNew = -1;
            SCSN_postHocPower_avg_line_class thisAdd_line;
            Dictionary<int,bool> this_numberOfAddedDatasets_dist = new Dictionary<int,bool>();
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                indexNew++;
                thisAdd_line = this.PostHocPower_lines[indexThis];
                if (!this_numberOfAddedDatasets_dist.ContainsKey(thisAdd_line.Number_of_datasets))
                {
                    this_numberOfAddedDatasets_dist.Add(thisAdd_line.Number_of_datasets, true);
                }
                new_postHocPower_lines[indexNew] = thisAdd_line;
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                indexNew++;
                thisAdd_line = add_postHocPower_lines[indexAdd];
                if (this_numberOfAddedDatasets_dist.ContainsKey(thisAdd_line.Number_of_datasets)) { throw new Exception("Added number of added datasets already exists"); }
                new_postHocPower_lines[indexNew] = thisAdd_line;
            }
            this.PostHocPower_lines = new_postHocPower_lines;
        }
        public void Generate_from_avgWithinEachNumberOfDatasets_instance_and_add_to_array(SCSN_avgWithinEachNumberOfDatasets_deg_class scsn_avg_input)
        {
            SCSN_avgWithinEachNumberOfDatasets_deg_class scsn_avg = scsn_avg_input.Deep_copy();
            if (Options.Only_considere_genes_that_were_significant_at_least_once)
            {
                scsn_avg.Keep_only_lines_with_isSigArithMean_larger_zero();
            }


            if (PostHocPower_lines.Length==0) { Data_was_normalized_before_postHocPowerResults_generation = scsn_avg.Data_is_normalized_with_highest_datasetNumber_as_reference; }
            if (Data_was_normalized_before_postHocPowerResults_generation != scsn_avg.Data_is_normalized_with_highest_datasetNumber_as_reference) { throw new Exception("Exception in: Generate_from_avgWithinEachNumberOfDatasets_instance_and_add_to_array"); }
            int number_of_drawn_datasets = scsn_avg.DEGs[0].Number_of_datasets;

            Math_class mfs = new Math_class();

            int number_of_datasets = scsn_avg.DEGs[0].Number_of_datasets;

            Dictionary<int, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>> numberOfAddedDatasets_cluster_valueType_runningMean_dict = new Dictionary<int, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>>();
            Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>> cluster_valueType_runningMean_dict = new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>();
            Dictionary<AvgValueType_of_interest_enum, double> valueType_runningMean_dict = new Dictionary<AvgValueType_of_interest_enum, double>();
            double runningMean;

            Dictionary<int, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>> numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict = new Dictionary<int, Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>>();
            Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>> cluster_valueType_runningPopVariance_dict = new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>();
            Dictionary<AvgValueType_of_interest_enum, double> valueType_runningPopVariance_dict = new Dictionary<AvgValueType_of_interest_enum, double>();
            double runningPopVariance;

            Dictionary<string, Dictionary<string, bool>> cluster_gene_dict = new Dictionary<string, Dictionary<string, bool>>();
            Dictionary<int, Dictionary<string, Dictionary<string, bool>>> numberOfDatasets_cluster_gene_dict = new Dictionary<int, Dictionary<string, Dictionary<string, bool>>>();

            int degs_length = scsn_avg.DEGs.Length;
            SCSN_avgWithinEachNumberOfDatasets_deg_line_class scsn_avg_line;
            for (int indexDEG = 0; indexDEG < degs_length; indexDEG++)
            {
                scsn_avg_line = scsn_avg.DEGs[indexDEG];
                if (!cluster_gene_dict.ContainsKey(scsn_avg_line.Cluster)) { cluster_gene_dict.Add(scsn_avg_line.Cluster, new Dictionary<string, bool>()); }
                if (!cluster_gene_dict[scsn_avg_line.Cluster].ContainsKey(scsn_avg_line.Gene)) { cluster_gene_dict[scsn_avg_line.Cluster].Add(scsn_avg_line.Gene, true); }

                #region Check for right input data by checking for duplicates within the same number of datasets
                if (!numberOfDatasets_cluster_gene_dict.ContainsKey(scsn_avg_line.Number_of_datasets)) { numberOfDatasets_cluster_gene_dict.Add(scsn_avg_line.Number_of_datasets, new Dictionary<string, Dictionary<string, bool>>()); }
                if (!numberOfDatasets_cluster_gene_dict[scsn_avg_line.Number_of_datasets].ContainsKey(scsn_avg_line.Cluster)) { numberOfDatasets_cluster_gene_dict[scsn_avg_line.Number_of_datasets].Add(scsn_avg_line.Cluster, new Dictionary<string, bool>()); 
                numberOfDatasets_cluster_gene_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster].Add(scsn_avg_line.Gene, true); }
                #endregion
            }
            Statistical_measure_enum source_statistical_measure = Statistical_measure_enum.Coefficient_of_variation;
            AvgValueType_of_interest_enum[] avgValueType_of_interest_array = new AvgValueType_of_interest_enum[] { AvgValueType_of_interest_enum.Minus_log10_adj_pvalue, AvgValueType_of_interest_enum.Avg_log2_fc, AvgValueType_of_interest_enum.Rank, AvgValueType_of_interest_enum.Is_sig };
            AvgValueType_of_interest_enum avgValueType_of_interest;
            int avgValueType_of_interest_array_length = avgValueType_of_interest_array.Length;


            Dictionary<int, Dictionary<string, List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>>> numberOfDatasets_cluster_avgLines_dict = new Dictionary<int, Dictionary<string, List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>>>();
            Dictionary<string, List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>> cluster_avgLines_dict;

            #region Calcuate means over all genes expressed in same cluster by adding to running means and fill cluster_gene_avgLine_dict
            double add_value;
            int genes_count;
            double max_add_value = 0;
            for (int indexDEG = 0; indexDEG < degs_length; indexDEG++)
            {
                scsn_avg_line = scsn_avg.DEGs[indexDEG];

                #region Fill cluster_gene_avgLines_dict and iteration_dict
                if (!numberOfDatasets_cluster_avgLines_dict.ContainsKey(scsn_avg_line.Number_of_datasets))
                {
                    numberOfDatasets_cluster_avgLines_dict.Add(scsn_avg_line.Number_of_datasets, new Dictionary<string, List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>>());
                }
                if (!numberOfDatasets_cluster_avgLines_dict[scsn_avg_line.Number_of_datasets].ContainsKey(scsn_avg_line.Cluster))
                {
                    numberOfDatasets_cluster_avgLines_dict[scsn_avg_line.Number_of_datasets].Add(scsn_avg_line.Cluster, new List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class>());
                }
                numberOfDatasets_cluster_avgLines_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster].Add(scsn_avg_line);
                #endregion

                if (!numberOfAddedDatasets_cluster_valueType_runningMean_dict.ContainsKey(scsn_avg_line.Number_of_datasets))
                {
                    numberOfAddedDatasets_cluster_valueType_runningMean_dict.Add(scsn_avg_line.Number_of_datasets, new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>());
                }
                if (!numberOfAddedDatasets_cluster_valueType_runningMean_dict[scsn_avg_line.Number_of_datasets].ContainsKey(scsn_avg_line.Cluster))
                {
                    numberOfAddedDatasets_cluster_valueType_runningMean_dict[scsn_avg_line.Number_of_datasets].Add(scsn_avg_line.Cluster, new Dictionary<AvgValueType_of_interest_enum, double>());
                }
                genes_count = cluster_gene_dict[scsn_avg_line.Cluster].Keys.Count;
                for (int indexAvgVT = 0; indexAvgVT < avgValueType_of_interest_array_length; indexAvgVT++)
                {
                    avgValueType_of_interest = avgValueType_of_interest_array[indexAvgVT];
                    if (!numberOfAddedDatasets_cluster_valueType_runningMean_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster].ContainsKey(avgValueType_of_interest))
                    {
                        numberOfAddedDatasets_cluster_valueType_runningMean_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster].Add(avgValueType_of_interest, 0);
                    }
                    add_value = scsn_avg_line.Get_statisticalMeasure_for_valueTypeOfInterest(avgValueType_of_interest, source_statistical_measure);
                    if (add_value > max_add_value) { max_add_value = add_value; }
                    numberOfAddedDatasets_cluster_valueType_runningMean_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster][avgValueType_of_interest]
                        = mfs.Add_to_running_mean(numberOfAddedDatasets_cluster_valueType_runningMean_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster][avgValueType_of_interest],
                                                    add_value,
                                                    genes_count);
                }
            }
            #endregion

            #region Calculate popultion variances by adding to running popVariances
            int[] numberOfDatasets_array = numberOfDatasets_cluster_avgLines_dict.Keys.ToArray();
            int numberOfDatasets;
            int numberOfDatasets_array_length = numberOfDatasets_array.Length;
            string[] clusters;
            string cluster;
            int clusters_length;
            List<SCSN_avgWithinEachNumberOfDatasets_deg_line_class> avg_lines;
            SCSN_avgWithinEachNumberOfDatasets_deg_line_class add_missingGene_line;
            int avg_lines_count;
            for (int indexNumberOfDatasets = 0; indexNumberOfDatasets < numberOfDatasets_array_length; indexNumberOfDatasets++)
            {
                numberOfDatasets = numberOfDatasets_array[indexNumberOfDatasets];
                cluster_avgLines_dict = numberOfDatasets_cluster_avgLines_dict[numberOfDatasets];
                clusters = cluster_avgLines_dict.Keys.ToArray();
                clusters_length = clusters.Length;
                for (int indexCluster=0; indexCluster<clusters_length;indexCluster++)
                {
                    cluster = clusters[indexCluster];
                    avg_lines = cluster_avgLines_dict[cluster];
                    avg_lines_count = avg_lines.Count;
                    genes_count = cluster_gene_dict[cluster].Count;
                    for (int indexAdd=avg_lines_count; indexAdd<genes_count; indexAdd++)
                    {
                        add_missingGene_line = new SCSN_avgWithinEachNumberOfDatasets_deg_line_class();
                        for (int indexAvgVT = 0; indexAvgVT < avgValueType_of_interest_array_length; indexAvgVT++)
                        {
                            avgValueType_of_interest = avgValueType_of_interest_array[indexAvgVT];
                            add_missingGene_line.Set_statisticalMeasure_for_valueTypeOfInterest(0, avgValueType_of_interest, source_statistical_measure);
                        }
                        avg_lines.Add(add_missingGene_line);
                    }
                    for (int indexAvg = 0; indexAvg < avg_lines_count; indexAvg++)
                    {
                        scsn_avg_line = avg_lines[indexAvg];
                        if (!numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict.ContainsKey(scsn_avg_line.Number_of_datasets))
                        {
                            numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict.Add(scsn_avg_line.Number_of_datasets, new Dictionary<string, Dictionary<AvgValueType_of_interest_enum, double>>());
                        }
                        if (!numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict[scsn_avg_line.Number_of_datasets].ContainsKey(scsn_avg_line.Cluster))
                        {
                            numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict[scsn_avg_line.Number_of_datasets].Add(scsn_avg_line.Cluster, new Dictionary<AvgValueType_of_interest_enum, double>());
                        }
                        for (int indexAvgVT = 0; indexAvgVT < avgValueType_of_interest_array_length; indexAvgVT++)
                        {
                            avgValueType_of_interest = avgValueType_of_interest_array[indexAvgVT];
                            if (!numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster].ContainsKey(avgValueType_of_interest))
                            {
                                numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster].Add(avgValueType_of_interest, 0);
                            }
                            add_value = scsn_avg_line.Get_statisticalMeasure_for_valueTypeOfInterest(avgValueType_of_interest, source_statistical_measure);
                            numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster][avgValueType_of_interest]
                                = mfs.Add_to_running_population_variance(numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster][avgValueType_of_interest],
                                                                         add_value,
                                                                         numberOfAddedDatasets_cluster_valueType_runningMean_dict[scsn_avg_line.Number_of_datasets][scsn_avg_line.Cluster][avgValueType_of_interest],
                                                                         genes_count);
                        }
                    }
                }
            }
            #endregion

            #region Generate new avg_lines
            SCSN_postHocPower_avg_line_class php_avg_line;
            List<SCSN_postHocPower_avg_line_class> php_avg_lines = new List<SCSN_postHocPower_avg_line_class>();

            int[] numberOfAddedDatasets_array = numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict.Keys.ToArray();
            int numberOfAddedDatasets;
            int numberOfAddedDatasets_array_length = numberOfAddedDatasets_array.Length;
            double sampleSD;
            double coeffVar;

            for (int indexNAD=0; indexNAD<numberOfAddedDatasets_array_length;indexNAD++)
            {
                numberOfAddedDatasets = numberOfAddedDatasets_array[indexNAD];
                cluster_valueType_runningMean_dict = numberOfAddedDatasets_cluster_valueType_runningMean_dict[numberOfAddedDatasets];
                cluster_valueType_runningPopVariance_dict = numberOfAddedDatasets_cluster_valueType_runningPopVariance_dict[numberOfAddedDatasets];
                clusters = cluster_valueType_runningMean_dict.Keys.ToArray();
                clusters_length = clusters.Length;
                for (int indexC=0; indexC<clusters_length;indexC++)
                {
                    cluster = clusters[indexC];
                    genes_count = cluster_gene_dict[cluster].Keys.Count;
                    valueType_runningMean_dict = cluster_valueType_runningMean_dict[cluster];
                    valueType_runningPopVariance_dict = cluster_valueType_runningPopVariance_dict[cluster];
                    php_avg_line = new SCSN_postHocPower_avg_line_class();
                    php_avg_line.Cluster = (string)cluster.Clone();
                    php_avg_line.Number_of_datasets = numberOfAddedDatasets;
                    php_avg_line.Number_of_genes = genes_count;
                    php_avg_line.Avg_data_was_normalized = Data_was_normalized_before_postHocPowerResults_generation;
                    for (int indexAVT=0; indexAVT<avgValueType_of_interest_array_length;indexAVT++)
                    {
                        avgValueType_of_interest = avgValueType_of_interest_array[indexAVT];
                        runningMean = valueType_runningMean_dict[avgValueType_of_interest];
                        runningPopVariance = valueType_runningPopVariance_dict[avgValueType_of_interest];
                        sampleSD = mfs.Calculate_sampleSD_from_population_variance(runningPopVariance, genes_count);
                        if ((genes_count == 1)&&(sampleSD!=0)) { throw new Exception(); }
                        coeffVar = mfs.Calculate_coefficient_of_variance(runningMean, sampleSD);
                        php_avg_line.Source_statistical_measure = source_statistical_measure;
                        php_avg_line.Set_arthMean_for_valueTypeOfInterest(runningMean, avgValueType_of_interest);
                        php_avg_line.Set_sampleSD_for_valueTypeOfInterest(sampleSD, avgValueType_of_interest);
                        php_avg_line.Set_coeffOfVar_for_valueTypeOfInterest(coeffVar, avgValueType_of_interest);
                    }
                    php_avg_lines.Add(php_avg_line);
                }
            }
            #endregion

            Add_to_array_and_check_if_number_of_added_datasets_do_not_exist(php_avg_lines.ToArray());
        }
        public void Write_data(string directory, string fileName)
        {
            SCSN_postHocPower_avg_readWriteOptions readWriteOptions = new SCSN_postHocPower_avg_readWriteOptions(directory, fileName, Data_was_normalized_before_postHocPowerResults_generation);
            ReadWriteClass.WriteData(this.PostHocPower_lines, readWriteOptions);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class SCSN_postHocPower_correlation_line_class
    {
        #region Fields
        public string Cluster { get; set; }
        public double Correlation_value_mean { get; set; }
        public double Correlation_value_sampleSd { get; set; }
        public Correlation_valueType_enum Correlation_valueType { get; set; }
        public AvgValueType_of_interest_enum ValueType_of_interest { get; set; }
        public int Number_of_drawn_datasets { get; set; }
        #endregion

        public SCSN_postHocPower_correlation_line_class()
        {
            Cluster = "";
            Number_of_drawn_datasets = 0;
        }
    }

    class SCSN_readWriteOptions : ReadWriteOptions_base
    {
        public SCSN_readWriteOptions(string directory, string fileName)
        {
            this.File = directory + fileName;
            this.Key_propertyNames = new string[] { "Cluster","Correlation_value_mean", "Correlation_value_sampleSd", "Correlation_valueType","ValueType_of_interest","Number_of_drawn_datasets" };
            this.Key_columnNames = this.Key_propertyNames;
            this.File_has_headline = true;
            this.Delimiter = Global_class.Tab;
        }
    }
    class SCSN_postHoPower_correlation_options_class
    {
        public bool Consider_only_genes_that_were_predicted_with_significance_at_least_once { get; set; }

    }

    class SCSN_postHocPower_correlation_class
    {
        public SCSN_postHocPower_correlation_line_class[] Correlations { get; set; }
        public SCSN_postHoPower_correlation_options_class Options { get; set; }

        public SCSN_postHocPower_correlation_class(bool consider_only_genes_that_were_predicted_with_significance_at_least_once)
        {
            this.Correlations = new SCSN_postHocPower_correlation_line_class[0];
            Options = new SCSN_postHoPower_correlation_options_class();
            Options.Consider_only_genes_that_were_predicted_with_significance_at_least_once = consider_only_genes_that_were_predicted_with_significance_at_least_once;
        }

        private void Add_to_array_and_check_if_valueTypeOfInterest_and_numberOfDatasets_combination_does_not_exist(SCSN_postHocPower_correlation_line_class[] add_correlations)
        {
            int this_length = this.Correlations.Length;
            int add_length = add_correlations.Length;
            int new_length = this_length + add_length;
            SCSN_postHocPower_correlation_line_class thisAdd_line;
            SCSN_postHocPower_correlation_line_class[] new_correlations = new SCSN_postHocPower_correlation_line_class[new_length];
            int indexNew = -1;
            Dictionary<int, Dictionary<AvgValueType_of_interest_enum, bool>> existing_numberOfDDatasets_valueType_dict = new Dictionary<int, Dictionary<AvgValueType_of_interest_enum, bool>>();
            for (int indexThis = 0; indexThis < this_length; indexThis++)
            {
                thisAdd_line = this.Correlations[indexThis];
                indexNew++;
                new_correlations[indexNew] = thisAdd_line;
                if (!existing_numberOfDDatasets_valueType_dict.ContainsKey(thisAdd_line.Number_of_drawn_datasets))
                {
                    existing_numberOfDDatasets_valueType_dict.Add(thisAdd_line.Number_of_drawn_datasets, new Dictionary<AvgValueType_of_interest_enum, bool>());
                }
                if (!existing_numberOfDDatasets_valueType_dict[thisAdd_line.Number_of_drawn_datasets].ContainsKey(thisAdd_line.ValueType_of_interest))
                {
                    existing_numberOfDDatasets_valueType_dict[thisAdd_line.Number_of_drawn_datasets].Add(thisAdd_line.ValueType_of_interest, true);
                }
            }
            for (int indexAdd = 0; indexAdd < add_length; indexAdd++)
            {
                thisAdd_line = add_correlations[indexAdd];
                indexNew++;
                new_correlations[indexNew] = thisAdd_line;
                if (  (existing_numberOfDDatasets_valueType_dict.ContainsKey(thisAdd_line.Number_of_drawn_datasets))
                    &&(existing_numberOfDDatasets_valueType_dict[thisAdd_line.Number_of_drawn_datasets].ContainsKey(thisAdd_line.ValueType_of_interest)))
                {
                    throw new Exception("AvgValueType exists already.");
                }
            }
            this.Correlations = new_correlations;
        }

        public void Calculate_pearson_correlation_for_valueTypeOfInterest_and_add_to_array(SCSN_avgAnalysesWithinEachIteration_deg_class scsn_avg_input, SCSN_avgAnalysesWithinEachIteration_deg_class reference_scsn_avg, AvgValueType_of_interest_enum valueType_of_interest)
        {
            SCSN_avgAnalysesWithinEachIteration_deg_class scsn_avg = scsn_avg_input.Deep_copy();
            if (Options.Consider_only_genes_that_were_predicted_with_significance_at_least_once)
            {
                scsn_avg.Keep_only_lines_with_isSigArtihMean_larger_zero();
            }

            SCSN_avgAnalysesWithinEachIteration_deg_line_class avg_line;
            Dictionary<string, Dictionary<string, double>> maxDatasetNos_cluster_gene_valueArthMean_dict = new Dictionary<string, Dictionary<string, double>>();
            Dictionary<string, double> maxDatasetNos_gene_valueArthMean_dict = new Dictionary<string, double>();
            Dictionary<string, Dictionary<string, bool>> cluster_gene_dict = new Dictionary<string, Dictionary<string, bool>>();
            int reference_avg_length = reference_scsn_avg.DEGs.Length;
            int max_datasets = -1;
            bool multiple_genes_for_same_cluster = false;
            for (int indexAvg = 0; indexAvg < reference_avg_length; indexAvg++)
            {
                avg_line = reference_scsn_avg.DEGs[indexAvg];
                if ((max_datasets == -1)
                    && (avg_line.Number_of_datasets > max_datasets))
                {
                    maxDatasetNos_cluster_gene_valueArthMean_dict.Clear();
                    max_datasets = avg_line.Number_of_datasets;
                    multiple_genes_for_same_cluster = false;
                }
                if (avg_line.Number_of_datasets == max_datasets)
                {
                    if (!maxDatasetNos_cluster_gene_valueArthMean_dict.ContainsKey(avg_line.Cluster))
                    {
                        maxDatasetNos_cluster_gene_valueArthMean_dict.Add(avg_line.Cluster, new Dictionary<string, double>());
                    }
                    if (avg_line.Get_arithMean_for_value_of_interest(valueType_of_interest) != 0)
                    {
                        if (!maxDatasetNos_cluster_gene_valueArthMean_dict[avg_line.Cluster].ContainsKey(avg_line.Gene))
                        { maxDatasetNos_cluster_gene_valueArthMean_dict[avg_line.Cluster].Add(avg_line.Gene, avg_line.Get_arithMean_for_value_of_interest(valueType_of_interest)); }
                        else { multiple_genes_for_same_cluster = true; }
                    }
                }
                if (!cluster_gene_dict.ContainsKey(avg_line.Cluster))
                {
                    cluster_gene_dict.Add(avg_line.Cluster, new Dictionary<string, bool>());
                }
                if (!cluster_gene_dict[avg_line.Cluster].ContainsKey(avg_line.Gene))
                {
                    cluster_gene_dict[avg_line.Cluster].Add(avg_line.Gene, true);
                }
            }
            if (multiple_genes_for_same_cluster) { throw new Exception(); }

            scsn_avg.DEGs = scsn_avg.DEGs.OrderBy(l => l.Cluster).ThenByDescending(l => l.Number_of_datasets).ThenBy(l=>l.Iteration).ThenBy(l => l.Gene).ToArray();
            int scsn_avg_length = scsn_avg.DEGs.Length;
            Dictionary<string, Dictionary<string, double>> maxDataset_cluster_gene_arithMean_dict = new Dictionary<string, Dictionary<string, double>>();
            double[] maxDatasetNo_values = new double[0];
            double[] thisDatasetNo_values = new double[0];
            string[] allGenes = new string[0];
            int allGenes_length = -1;
            string allGene;
            int indexAllGenes = -1;
            int allGenes_compare = -2;
            double pearson_correlation_coefficient;
            double pearson_correlation_coefficient_instable;
            Correlation_coefficient_class correlation = new Correlation_coefficient_class();
            SCSN_postHocPower_correlation_line_class new_postHocPower_correlation_line;
            List<SCSN_postHocPower_correlation_line_class> new_postHocPower_correlation_lines = new List<SCSN_postHocPower_correlation_line_class>();
            List<double> current_cluster_numberOfDatasets_pearson_correlations = new List<double>();
            Math_class mfs = new Math_class();
            double pearson_correlation_mean;
            double pearson_correlation_sampleSD;
            for (int indexAvg = 0; indexAvg < scsn_avg_length; indexAvg++)
            {
                avg_line = scsn_avg.DEGs[indexAvg];
                if ((indexAvg == 0)
                    || (!avg_line.Cluster.Equals(scsn_avg.DEGs[indexAvg - 1].Cluster)))
                {
                    allGenes = cluster_gene_dict[avg_line.Cluster].Keys.ToArray();
                    allGenes = allGenes.OrderBy(l => l).ToArray();
                    allGenes_length = allGenes.Length;
                    maxDatasetNos_gene_valueArthMean_dict = maxDatasetNos_cluster_gene_valueArthMean_dict[avg_line.Cluster];
                    maxDatasetNo_values = new double[allGenes_length];
                    for (indexAllGenes = 0; indexAllGenes < allGenes_length; indexAllGenes++)
                    {
                        allGene = allGenes[indexAllGenes];
                        if (maxDatasetNos_gene_valueArthMean_dict.ContainsKey(allGene)) { maxDatasetNo_values[indexAllGenes] = maxDatasetNos_gene_valueArthMean_dict[allGene]; }
                        else { maxDatasetNo_values[indexAllGenes] = 0; }
                    }
                }
                if ((indexAvg == 0)
                    || (!avg_line.Cluster.Equals(scsn_avg.DEGs[indexAvg - 1].Cluster))
                    || (!avg_line.Number_of_datasets.Equals(scsn_avg.DEGs[indexAvg - 1].Number_of_datasets)))
                {
                    current_cluster_numberOfDatasets_pearson_correlations.Clear();
                }
                if ((indexAvg == 0)
                    || (!avg_line.Cluster.Equals(scsn_avg.DEGs[indexAvg - 1].Cluster))
                    || (!avg_line.Number_of_datasets.Equals(scsn_avg.DEGs[indexAvg - 1].Number_of_datasets))
                    || (!avg_line.Iteration.Equals(scsn_avg.DEGs[indexAvg - 1].Iteration)))
                {
                    thisDatasetNo_values = new double[allGenes_length];
                    indexAllGenes = 0;
                }
                allGenes_compare = -2;
                while ((indexAllGenes < allGenes_length) && (allGenes_compare < 0))
                {
                    allGene = allGenes[indexAllGenes];
                    allGenes_compare = allGene.CompareTo(avg_line.Gene);
                    if (allGenes_compare < 0)
                    {
                        thisDatasetNo_values[indexAllGenes] = 0;
                        indexAllGenes++;
                    }
                    else if (allGenes_compare == 0)
                    {
                        thisDatasetNo_values[indexAllGenes] = avg_line.Get_arithMean_for_value_of_interest(valueType_of_interest);
                        indexAllGenes++;
                    }
                }
                if ((indexAvg == scsn_avg_length - 1)
                    || (!avg_line.Cluster.Equals(scsn_avg.DEGs[indexAvg + 1].Cluster))
                    || (!avg_line.Number_of_datasets.Equals(scsn_avg.DEGs[indexAvg + 1].Number_of_datasets))
                    || (!avg_line.Iteration.Equals(scsn_avg.DEGs[indexAvg + 1].Iteration)))
                {
                    while (indexAllGenes < allGenes_length)
                    {
                        thisDatasetNo_values[indexAllGenes] = 0;
                        indexAllGenes++;
                    }
                    pearson_correlation_coefficient_instable = correlation.Get_pearson_correlation_coefficient_numerically_instable(thisDatasetNo_values, maxDatasetNo_values);
                    pearson_correlation_coefficient = correlation.Get_pearson_correlation_coefficient_stable(thisDatasetNo_values, maxDatasetNo_values);
                    if ((Global_class.Is_developer_mode) && (pearson_correlation_coefficient_instable!=0))
                    {
                        if (Math.Round(100 * pearson_correlation_coefficient) != Math.Round(100 * pearson_correlation_coefficient_instable)) { throw new Exception(); }
                    }
                    current_cluster_numberOfDatasets_pearson_correlations.Add(pearson_correlation_coefficient);
                }
                if ((indexAvg == scsn_avg_length - 1)
                    || (!avg_line.Cluster.Equals(scsn_avg.DEGs[indexAvg + 1].Cluster))
                    || (!avg_line.Number_of_datasets.Equals(scsn_avg.DEGs[indexAvg + 1].Number_of_datasets)))
                {

                    mfs.Get_mean_and_sample_sd(current_cluster_numberOfDatasets_pearson_correlations.ToArray(), out pearson_correlation_mean, out pearson_correlation_sampleSD);
                    new_postHocPower_correlation_line = new SCSN_postHocPower_correlation_line_class();
                    new_postHocPower_correlation_line.Cluster = (string)avg_line.Cluster.Clone();
                    new_postHocPower_correlation_line.Number_of_drawn_datasets = avg_line.Number_of_datasets;
                    new_postHocPower_correlation_line.Correlation_value_mean = pearson_correlation_mean;
                    new_postHocPower_correlation_line.Correlation_value_sampleSd = pearson_correlation_sampleSD;
                    new_postHocPower_correlation_line.Correlation_valueType = Correlation_valueType_enum.Pearson_coefficient;
                    new_postHocPower_correlation_line.ValueType_of_interest = valueType_of_interest;
                    new_postHocPower_correlation_lines.Add(new_postHocPower_correlation_line);
                }
            }
            this.Add_to_array_and_check_if_valueTypeOfInterest_and_numberOfDatasets_combination_does_not_exist(new_postHocPower_correlation_lines.ToArray());
        }

        public void Write(string directory, string fileName)
        {
            SCSN_readWriteOptions readWriteOptions = new SCSN_readWriteOptions(directory, fileName);
            ReadWriteClass.WriteData(this.Correlations, readWriteOptions);
        }
    }
}
