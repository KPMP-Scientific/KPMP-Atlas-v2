using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Common_functions
{

    enum Ontology_type_enum {  E_m_p_t_y, Mbco, Whole_cell_subfunction, Whole_cell_function }
    enum Timeunit_enum {  E_m_p_t_y, hrs, min }
    enum UpDown_enum {  E_m_p_t_y, Up, Down }
    enum Network_interaction_type_enum { E_m_p_t_y, Activation, Inhibition, Directed, Complex, Interaction }
    enum Network_databases_enum { E_m_p_t_y, Bind, Biocarta, BioGRID, Biogrid_filtered, DIP, Dip_filtered, figeys, Figeys, Hprd, HPRD, Humap, Innatedb_filtered, Iref, Intact_filtered, Kegg, Mint, Mips, Bioplex_ppi, InnateDB_mod, INtAct_mod, KEA_mod, KEGG, MINT, MIPS, murphy, Pdzbase, pdzbase, Ppid, ppid, Predicted_ppid, Snavi, SNAVI_mod, Stelzl, Vidal, Vidal_hi_ii_14, Bioplex_interactionlist_v4, Seth, Metabolites, Interactome }
    enum Network_direction_type_enum { E_m_p_t_y = 0, Directed_forward = 1, Directed_backwards = 2, Undirected_double = 3, Undirected_single = 4, Bipartite_forward = 5, Bipartite_backwards = 6 }
    enum KPMP_category_for_cluster_enum { E_m_p_t_y, Scp_minusLog10Pvalues }



    class Global_class
    {
        const bool check_for_correct_ordering = true;

        #region Const
        //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
        private const string empty_entry = "E_m_p_t_y";  //check enums, Empty has to be the same!!
        private const string space_text = "S_p_a_c_e";
        #endregion

        public static bool Check_for_correct_ordering {  get { return check_for_correct_ordering; } }
        public static string Empty_entry {  get {  return empty_entry; } }
        public static string Space_text {  get { return space_text; } }
        public static char Tab {  get { return '\t'; } }
    }
    class Global_directory_class
    {
        public static string Enrichment_results_input_fileName {  get { return "Pathways_sig_forWholeCellFunct.txt"; } }
        public static string NW_subdirectory { get { return "Wcf_nws/"; } }
        public static string Get_enrichment_results_subdirectory(string assay_cellType)
        {
            string enrichment_subdirectory = assay_cellType + "_pathways/";
            return enrichment_subdirectory;
        }
        public static string Get_physiologicalFunctionsDefinitions_subdirectory(string assay_cellType)
        {
            string physiologicalFunctions_subdirectory = assay_cellType + "_pathways/";
            return physiologicalFunctions_subdirectory;
        }
        public static string Get_physiologicalFunctions_subdirectory(string assay_cellType)
        {
            string physiologicalFunctions_subdirectory = assay_cellType + "_pathways/";
            return physiologicalFunctions_subdirectory;
        }
        public static string Get_pFGroup_subdirectory(string assay_cellType)
        {
            string physiologicalFunctions_subdirectory = assay_cellType + "_pathways/";
            return physiologicalFunctions_subdirectory;
        }
        public static string Get_physiologicalFunctionDefinitions_fileName(string assay_cellType)
        {
            string physiologicalFunctions_fileName = "WcfDefs_" + assay_cellType + ".txt";
            return physiologicalFunctions_fileName;
        }
        public static string Get_physiologicalFunctions_fileName(string assay_cellType)
        {
            string physiologicalFunctions_fileName = "WholeCellSubFunctions_" + assay_cellType + ".txt";
            return physiologicalFunctions_fileName;
        }
        public static string Get_pFGroup_fileName(string assay_cellType)
        {
            string physiologicalFunctions_fileName = "WholeCellFunctions_" + assay_cellType + ".txt";
            return physiologicalFunctions_fileName;
        }
        public static string Get_pfGroup_nw_basefileName()
        {
            string physiologicalFunctions_subdirectory = "Trajectories";
            return physiologicalFunctions_subdirectory;
        }
    }

    public class Color_class
    {
        System.Drawing.Color[] All_csharp_colors { get; set; }
        public Color_class()
        {
            All_csharp_colors = Get_all_csharp_colors();
        }

        private static string Get_hexadecimal_sign(int number)
        {
            string sign = "no value";
            switch (number)
            {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                    sign = number.ToString();
                    break;
                case 10:
                    sign = "A";
                    break;
                case 11:
                    sign = "B";
                    break;
                case 12:
                    sign = "C";
                    break;
                case 13:
                    sign = "D";
                    break;
                case 14:
                    sign = "E";
                    break;
                case 15:
                    sign = "F";
                    break;
                default:
                    throw new Exception();
            }
            return sign;
        }
        public static string Convert_into_two_digit_hexadecimal(int number)
        {
            if ((number > 255) || (number < 0))
            {
                throw new Exception();
            }
            else
            {
                int multiples_of_16 = (int)Math.Floor((double)number / (double)16);
                int modulus = number % 16;
                return Get_hexadecimal_sign(multiples_of_16) + Get_hexadecimal_sign(modulus);
            }
        }
        public static string Get_hexadecimal_code(int red, int green, int blue)
        {
            return "#" + Convert_into_two_digit_hexadecimal(red) + Convert_into_two_digit_hexadecimal(green) + Convert_into_two_digit_hexadecimal(blue);
        }
        public static string Get_hexadecimal_black()
        {
            return Get_hexadecimal_code(0, 0, 0);
        }


        ///Inset was written by Chat GPT - Begin
        public (int R, int G, int B) ContinuousToRGB(double value)
        {
            if (value < 0 || value > 1)
                throw new ArgumentOutOfRangeException("Value must be in the range [0, 1]");

            double hue = value;
            double saturation = 1.0;
            double brightness = 1.0;

            return HSVToRGB(hue * 360, saturation, brightness);
        }

        public string ContinuousToHexadecimalColor(double value)
        {
            int red; int green; int blue;
            (red, green, blue) = ContinuousToRGB(value);
            return RGBToHex(red, green, blue);
        }

        public string RGBToHex(int r, int g, int b)
        {
            // Ensure the RGB values are within the valid range
            if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255)
            {
                throw new ArgumentException("RGB values must be in the range 0-255.");
            }

            // Convert each component to hexadecimal and concatenate
            return "#" + r.ToString("X2") + g.ToString("X2") + b.ToString("X2");
        }

        private (int R, int G, int B) HSVToRGB(double h, double s, double v)
        {
            double r = 0, g = 0, b = 0;

            int i = (int)(h / 60.0) % 6;
            double f = (h / 60.0) - Math.Floor(h / 60.0);
            double p = v * (1.0 - s);
            double q = v * (1.0 - f * s);
            double t = v * (1.0 - (1.0 - f) * s);

            switch (i)
            {
                case 0:
                    r = v; g = t; b = p; break;
                case 1:
                    r = q; g = v; b = p; break;
                case 2:
                    r = p; g = v; b = t; break;
                case 3:
                    r = p; g = q; b = v; break;
                case 4:
                    r = t; g = p; b = v; break;
                case 5:
                    r = v; g = p; b = q; break;
            }

            return ((int)(r * 255.0), (int)(g * 255.0), (int)(b * 255.0));
        }
        ///Inset was written by Chat GPT - END

        #region Static functions for call in property fields
        public static string Get_color_string(System.Drawing.Color color)
        {
            string color_string = color.ToString().Replace("Color ", "").Replace("[", "").Replace("]", "");
            return color_string;
        }

        public static System.Drawing.Color Set_color_from_string(string color_string)
        {
            System.Drawing.Color return_color;
            return_color = System.Drawing.Color.FromName(color_string);
            return return_color;
        }
        #endregion

        private System.Drawing.Color[] Get_all_csharp_colors()
        {
            Dictionary<System.Drawing.Color, bool> selectable_colors_dict = new Dictionary<System.Drawing.Color, bool>();
            System.Drawing.Color add_color;
            foreach (System.Reflection.PropertyInfo property in typeof(System.Drawing.Color).GetProperties())
            {
                if (property.PropertyType == typeof(System.Drawing.Color))
                {
                    add_color = (System.Drawing.Color)property.GetValue(null, null);
                    if ((!selectable_colors_dict.ContainsKey(add_color))
                        && (!add_color.Equals(System.Drawing.Color.Transparent))
                        && (!add_color.Equals(System.Drawing.Color.White)))
                    {
                        selectable_colors_dict.Add(add_color, false);
                    }
                }
            }
            return selectable_colors_dict.Keys.ToArray();
        }

        private System.Drawing.Color Get_closest_csharp_color(int input_red, int input_green, int input_blue)
        {
            int all_colors_length = All_csharp_colors.Length;
            System.Drawing.Color current_color;
            int csharp_red = -1;
            int csharp_green = -1;
            int csharp_blue = -1;
            float current_distance;
            float minimum_distance = 999999999;
            System.Drawing.Color selected_csharp_color = System.Drawing.Color.Gray;
            for (int indexColor = 0; indexColor < all_colors_length; indexColor++)
            {
                current_color = All_csharp_colors[indexColor];
                csharp_blue = int.Parse(current_color.B.ToString());
                csharp_red = int.Parse(current_color.R.ToString());
                csharp_green = int.Parse(current_color.G.ToString());
                current_distance = (float)Math.Sqrt(Math.Pow(input_red - csharp_red, 2) + Math.Pow(input_blue - csharp_blue, 2) + Math.Pow(input_green - csharp_green, 2));
                if (current_distance < minimum_distance)
                {
                    minimum_distance = current_distance;
                    selected_csharp_color = current_color;
                }
            }
            return selected_csharp_color;
        }

        public System.Drawing.Color Get_closest_csharp_color_for_hexadecimal_color_if_exists(string color_string)
        {
            System.Drawing.Color closest_color = System.Drawing.Color.FromName(color_string);
            if ((color_string.Substring(0, 1).Equals("#"))
                && (color_string.Length == 7))
            {
                try
                {
                    int red = int.Parse(color_string.Substring(1, 2), System.Globalization.NumberStyles.AllowHexSpecifier);
                    int green = int.Parse(color_string.Substring(3, 2), System.Globalization.NumberStyles.AllowHexSpecifier);
                    int blue = int.Parse(color_string.Substring(5, 2), System.Globalization.NumberStyles.AllowHexSpecifier);
                    closest_color = Get_closest_csharp_color(red, green, blue);

                }
                catch { }
            }
            return closest_color;
        }

        public static string Get_hexadecimal_code_for_color(System.Drawing.Color color)
        {
            return "#" + color.R.ToString("X2") + color.G.ToString("X2") + color.B.ToString("X2");
        }
    }

}
