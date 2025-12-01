using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace Common_functions
{

    enum Ontology_type_enum { E_m_p_t_y, Mbco, Mbco_physiological_function, Mbco_pf_group }
    enum Timeunit_enum { E_m_p_t_y, hrs, min }
    enum UpDown_enum { E_m_p_t_y, Up, Down }
    enum Network_interaction_type_enum { E_m_p_t_y, Activation, Inhibition, Directed, Complex, Interaction }
    enum Network_databases_enum { E_m_p_t_y, Bind, Biocarta, BioGRID, Biogrid_filtered, DIP, Dip_filtered, figeys, Figeys, Hprd, HPRD, Humap, Innatedb_filtered, Iref, Intact_filtered, Kegg, Mint, Mips, Bioplex_ppi, InnateDB_mod, INtAct_mod, KEA_mod, KEGG, MINT, MIPS, murphy, Pdzbase, pdzbase, Ppid, ppid, Predicted_ppid, Snavi, SNAVI_mod, Stelzl, Vidal, Vidal_hi_ii_14, Bioplex_interactionlist_v4, Seth, Metabolites, Interactome }
    enum Network_direction_type_enum { E_m_p_t_y = 0, Directed_forward = 1, Directed_backwards = 2, Undirected_double = 3, Undirected_single = 4, Bipartite_forward = 5, Bipartite_backwards = 6 }
    enum KPMP_category_for_cluster_enum { E_m_p_t_y, Scp_minusLog10Pvalues }



    class Global_class
    {
        const bool is_developer_mode = true;

        #region Const
        //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
        private const string empty_entry = "E_m_p_t_y";  //check enums, Empty has to be the same!!
        private const string space_text = "S_p_a_c_e";
        #endregion

        public static bool Is_developer_mode { get { return is_developer_mode; } }
        public static string Empty_entry { get { return empty_entry; } }
        public static string Space_text { get { return space_text; } }
        public static char Tab { get { return '\t'; } }
        public static double Max_considered_double {  get { return 1E150; } }
    }
    class Global_directory_class
    {
    }
    class Math_class
    {
        public Math_class()
        {
        }

        public long Get_binomial_coefficient(int n, int k)
        {
            if (k < 0 || k > n)
                return 0;
            if (k == 0 || k == n)
                return 1;

            k = Math.Min(k, n - k); // Use symmetry: C(n, k) == C(n, n-k)

            long result = 1;
            for (int i = 1; i <= k; i++)
            {
                result *= (n - (k - i));
                result /= i;
            }

            return result;
        }

        public float Get_sum(float[] values)
        {
            int values_length = values.Length;
            float sum = values[0];
            for (int indexV = 1; indexV < values_length; indexV++)
            {
                sum += values[indexV];
            }
            return sum;
        }
        public float Get_average(float[] values)
        {
            float values_length = values.Length;
            return Get_sum(values) / values_length;
        }
        private double Round_value(double sample_value)
        {
            if (sample_value == 0) { return sample_value; }
            else
            {
                double rounding_factor = 1E10;
                double check_value = sample_value * rounding_factor;
                sample_value = Math.Round(rounding_factor * sample_value) / rounding_factor;
                return sample_value;
            }
        }
        public double Get_average(double[] values)
        {
            int values_length = values.Length;
            double sum = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                sum += values[indexV];
            }
            return sum / (double)values_length;
        }
        public void Get_mean_and_population_sd(double[] values, out double mean, out double sd)
        {
            mean = Get_average(values);
            int values_length = values.Length;
            double root_term = 0;
            for (int indexV = 0; indexV < values_length; indexV++)
            {
                root_term += Math.Pow(values[indexV] - mean,2);
            }
            sd = Math.Sqrt(root_term / values_length);
        }
        public void Get_mean_and_sample_sd(double[] values, out double mean, out double sd)
        {
            int values_length = values.Length;
            mean = Get_average(values);
            double sum_squared_diff = 0;
            for (int i = 0; i < values_length; i++)
            {
                sum_squared_diff += Math.Pow(values[i] - mean, 2);
            }
            if (sum_squared_diff > 0) { sd = Math.Sqrt(sum_squared_diff / (values_length - 1)); }
            else { sd = 0; }
        }

        public static double Round_value_to_significant_digits(double sample_value)
        {
            //if (sample_value== 0) { return 0; }
            //else
            //{
            //    int considered_digits_count = 15;
            //    double sign = Math.Sign(sample_value);
            //    double abs_sample_value = Math.Abs(sample_value);
            //    int order_of_magnitude = (int)Math.Floor(Math.Log10(abs_sample_value));
            //    double rounding_factor = Math.Pow(10, order_of_magnitude - considered_digits_count + 1);
            //    double rounded_value = sign * Math.Round(sample_value / rounding_factor) * rounding_factor;
            //    return rounded_value;
            //}
            return sample_value;
        }

        public double Add_to_running_population_variance(double running_population_variance, double sample_value, double final_mean, int samples_length)
        {
            sample_value = Round_value_to_significant_digits(sample_value);
            final_mean = Round_value_to_significant_digits(final_mean);
            double factor_term = (sample_value - final_mean);
            factor_term = Round_value_to_significant_digits(factor_term);
            double product_term = factor_term * factor_term;
            double add = product_term / samples_length;
            add = Round_value_to_significant_digits(add);
            running_population_variance += add;
            if (double.IsInfinity(running_population_variance)) { throw new Exception("double.IsInfinity(running_population_variance)"); }
            if (double.IsNaN(running_population_variance)) { throw new Exception("double.IsNaN(running_population_variance)"); }
            return running_population_variance;
        }
        public double Add_to_running_mean(double running_mean, double sample_value, int samples_length)
        {
            double add = sample_value / (double)samples_length;
            add = Round_value_to_significant_digits(add);
            running_mean += add;
            if (double.IsInfinity(running_mean)) { throw new Exception("(double.IsInfinity(running_mean))"); }
            if (double.IsNaN(running_mean)) { throw new Exception("(double.IsNaN(running_mean))"); }
            return running_mean;
        }
        public double Calculate_sampleSD_from_population_variance(double population_variance, int samples_length)
        {
            if (Round_value_to_significant_digits(population_variance) == 0) { return 0; }
            else
            {
                double sample_sd = Math.Sqrt((population_variance / (samples_length - 1)) * samples_length);
                if (double.IsNaN(sample_sd)) { throw new Exception("double.IsNaN(sample_sd)"); }
                if (double.IsInfinity(sample_sd)) { throw new Exception("double.IsNaN(sample_sd)"); }
                if ((samples_length == 1) && (sample_sd != 0)) { throw new Exception(); }
                sample_sd = Round_value_to_significant_digits(sample_sd);
                return sample_sd;
            }
        }
        public double Calculate_coefficient_of_variance(double mean_value, double sd_value)
        {
            double coeff = sd_value / mean_value;
            if ((mean_value==0) && (sd_value==0)) { coeff = 0; }
            if (double.IsNaN(coeff)) { throw new Exception("double.IsNaN(coeff)"); }
            coeff = Round_value_to_significant_digits(coeff);
            return coeff;
        }
        public double Normalize_value(double original_value, double reference_value)
        {
            double normalized_value = 100 * original_value / reference_value;
            normalized_value = Round_value_to_significant_digits(normalized_value);
            return normalized_value;
        }
    }
    class Random_class
    {
        private Random Random_instance { get; set; }


        public Random_class(int randomSeedNodes_minus1equalsNone)
        {
            if (randomSeedNodes_minus1equalsNone > 0)
            {
                Random_instance = new Random(randomSeedNodes_minus1equalsNone);
            }
            else
            {
                Random_instance = new Random();
            }
        }

        public int[] Randomly_draw_non_overlapping_indices_from_number_of_available_indices(int number_of_drawn_indices, int number_of_available_indices)
        {
            if (number_of_drawn_indices > number_of_available_indices)
            {
                throw new ArgumentException("Cannot draw more unique indices than available.");
            }

            int[] pool = Enumerable.Range(0, number_of_available_indices).ToArray();

            // Fisher-Yates shuffle
            for (int indexPool = number_of_available_indices - 1; indexPool > 0; indexPool--)
            {
                int indexPool_switch = Random_instance.Next(indexPool + 1);
                (pool[indexPool], pool[indexPool_switch]) = (pool[indexPool_switch], pool[indexPool]);
            }
            return pool.Take(number_of_drawn_indices).ToArray();
        }
    }



}
