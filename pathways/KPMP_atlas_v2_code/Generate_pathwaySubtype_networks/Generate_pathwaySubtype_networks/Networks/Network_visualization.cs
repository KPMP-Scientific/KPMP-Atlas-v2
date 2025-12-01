using Common_functions;
using System.Drawing;
using Network_visualization;


namespace Network
{
    enum Regular_node_label_enum { E_m_p_t_y, Name }
    enum Regular_node_color_enum { E_m_p_t_y, Selected_color_size_by_value }
    enum Regular_node_shape_enum { E_m_p_t_y, Ellipse }

    class Visualization_of_regularNW_options : Visualization_of_nw_options
    {
        #region Fields

        public Regular_node_color_enum Node_color { get; set; }
        public Regular_node_label_enum Node_label { get; set; }
        public Regular_node_shape_enum Node_shape { get; set; }
        public bool Add_border_lines { get; set; }
        public float Node_label_font_size { get; set; }
        public int Reference_node_diameter { get; set; }
        public float Reference_minLog10pvalue { get; set; }
        #endregion

        public Visualization_of_regularNW_options()
        {
            Node_color = Regular_node_color_enum.Selected_color_size_by_value;
            Node_label = Regular_node_label_enum.Name;
            Node_shape = Regular_node_shape_enum.Ellipse;
            Node_label_font_size = 20;
            Add_border_lines = false;
        }

    }

    //////////////////////////////////////////////////////

    class Visualization_of_ontologyNW_class : Visualization_of_nw_basis
    {
        #region Fields
        Visualization_of_regularNW_options RegularOptions { get; set; }
        #endregion

        public Visualization_of_ontologyNW_class(Visualization_of_regularNW_options options)
        {
            RegularOptions = options;
            Options = options;
        }

        #region Nodes
        private void Set_nodeId_label_and_compartment(ref Visualization_of_nw_node_line visNode_line, Node_line_class un_line)
        {
            visNode_line.Id = (string)un_line.Name.Clone();
            visNode_line.Label_alignement = "center";
            visNode_line.Model_name = "sides";
            visNode_line.Model_position = "s";
            visNode_line.FontSize = RegularOptions.Node_label_font_size;
            switch (RegularOptions.Node_label)
            {
                case Regular_node_label_enum.Name:
                    visNode_line.Label = (string)un_line.Name.Clone();
                    break;
                default:
                    throw new Exception();
            }
        }
        private void Set_nodeColor(ref Visualization_of_nw_node_line visNode_line, Node_line_class un_line)
        {
            switch (RegularOptions.Node_color)
            {
                case Regular_node_color_enum.Selected_color_size_by_value:
                    int node_color_values = un_line.Node_color_values.Length;
                    if (node_color_values == 0) { throw new Exception(); }
                    Color_specification_line_class node_color_specification;
                    List<Color_specification_line_class> value_color_specifications = new List<Color_specification_line_class>();
                    for (int indexNodeValueColors = 0; indexNodeValueColors < node_color_values; indexNodeValueColors++)
                    {
                        node_color_specification = new Color_specification_line_class();
                        node_color_specification.Fill_color = un_line.Node_color_values[indexNodeValueColors].Node_color;
                        node_color_specification.Size = un_line.Node_color_values[indexNodeValueColors].Node_color_value;
                        node_color_specification.PieChart_sliceOrderNo = un_line.Node_color_values[indexNodeValueColors].PieChart_sliceOrderNo;
                        value_color_specifications.Add(node_color_specification);
                    }
                    visNode_line.Color_specifications = value_color_specifications.ToArray();
                    if (Options.Add_border_lines)
                    {
                        visNode_line.Border_style_color = Color_class.Get_hexadecimal_code_for_color(Color.Black);
                    }
                    else
                    {
                        visNode_line.Border_style_color = Color_class.Get_hexadecimal_code_for_color(value_color_specifications[0].Fill_color);
                    }
                    break;
                default:
                    throw new Exception();
            }
        }

        private void Set_nodeShape(ref Visualization_of_nw_node_line visNode_line, Node_line_class un_line)
        {
            switch (RegularOptions.Node_shape)
            {
                case Regular_node_shape_enum.Ellipse:
                    visNode_line.Shape_type = "ellipse";
                    visNode_line.Geometry_heigth = RegularOptions.Reference_node_diameter;
                    visNode_line.Geometry_width = RegularOptions.Reference_node_diameter;
                    break;
                default:
                    throw new Exception();
            }
        }

        private void Generate_pre_nodes(Node_class uniqueNodes)
        {
            int un_length = uniqueNodes.UNs.Length;
            Node_line_class un_line;
            Visualization_of_nw_node_line visNode_line;
            List<Visualization_of_nw_node_line> visNodes_list = new List<Visualization_of_nw_node_line>();
            for (int indexUN = 0; indexUN < un_length; indexUN++)
            {
                un_line = uniqueNodes.UNs[indexUN];
                visNode_line = new Visualization_of_nw_node_line();
                Set_nodeId_label_and_compartment(ref visNode_line, un_line);
                Set_nodeColor(ref visNode_line, un_line);
                Set_nodeShape(ref visNode_line, un_line);
                visNodes_list.Add(visNode_line);
            }
            VisNodes = visNodes_list.ToArray();
        }

        private void Adjust_nodeSizes_if_specified()
        {
            if (RegularOptions.Node_color.Equals(Regular_node_color_enum.Selected_color_size_by_value))
            {
                float max_color_values_for_node = -1;
                float min_color_values_for_node = -1;
                int visNodes_length = this.VisNodes.Length;
                Visualization_of_nw_node_line node_line;
                float currentNode_color_values;
                for (int indexVisu = 0; indexVisu < visNodes_length; indexVisu++)
                {
                    node_line = this.VisNodes[indexVisu];
                    currentNode_color_values = 0;
                    foreach (Color_specification_line_class color_specification_line in node_line.Color_specifications)
                    {
                        currentNode_color_values += color_specification_line.Size;
                    }
                    if ((max_color_values_for_node == -1)
                        || (currentNode_color_values > max_color_values_for_node))
                    {
                        max_color_values_for_node = currentNode_color_values;
                    }
                    if ((min_color_values_for_node == -1)
                        || (currentNode_color_values < min_color_values_for_node))
                    {
                        min_color_values_for_node = currentNode_color_values;
                    }
                }
                double reference_node_area = Math.Pow(0.5 * RegularOptions.Reference_node_diameter, 2) * Math.PI;
                double reference_minLog10P = RegularOptions.Reference_minLog10pvalue;
                double this_node_area;
                double this_diameter;
                for (int indexVisu = 0; indexVisu < visNodes_length; indexVisu++)
                {
                    node_line = this.VisNodes[indexVisu];
                    if (!node_line.Shape_type.Equals("ellipse")) { throw new Exception(); }
                    currentNode_color_values = 0;
                    foreach (Color_specification_line_class color_specification_line in node_line.Color_specifications)
                    {
                        currentNode_color_values += color_specification_line.Size;
                    }
                    this_node_area = reference_node_area * currentNode_color_values / reference_minLog10P;
                    this_diameter = 2 * Math.Sqrt(this_node_area / Math.PI);
                    node_line.Geometry_heigth = (float)this_diameter;
                    node_line.Geometry_width = (float)this_diameter;
                }
            }
        }

        private void Generate_nodes(NetworkBasis_class onto_nw)
        {
            Generate_pre_nodes(onto_nw.UniqueNodes);
            Adjust_nodeSizes_if_specified();
        }
        #endregion

        #region Edges
        private void Generate_edges(NetworkBasis_class nw, params Dictionary<string, Dictionary<string, float>>[] source_target_weight_dicts)
        {
            if (source_target_weight_dicts.Length > 1) { throw new Exception(); }
            int nw_count = nw.NW.Length;
            NetworkBasis_line_class nw_line;
            nw.UniqueNodes.Order_by_index();
            Node_line_class source_node_line;
            Node_line_class target_node_line;
            int targets_count;
            int indexTarget;
            List<Visualization_of_nw_edge_line> visEdges_list = new List<Visualization_of_nw_edge_line>();
            for (int indexNW = 0; indexNW < nw_count; indexNW++)
            {
                source_node_line = nw.UniqueNodes.Get_individual_node_line(indexNW);
                nw_line = nw.NW[indexNW];
                targets_count = nw_line.Targets.Length;
                for (int indexTargetIndex = 0; indexTargetIndex < targets_count; indexTargetIndex++)
                {
                    indexTarget = nw_line.Targets[indexTargetIndex];
                    target_node_line = nw.UniqueNodes.Get_individual_node_line(indexTarget);
                    Visualization_of_nw_edge_line edge_line = new Visualization_of_nw_edge_line();
                    edge_line.Arrow_color = Color_class.Get_hexadecimal_black();
                    edge_line.Edge_id = source_node_line.Name + "_to_" + target_node_line.Name;
                    edge_line.Source_id = (string)source_node_line.Name.Clone();
                    edge_line.Target_id = (string)target_node_line.Name.Clone();
                    edge_line.Arrow_source_end = "none";
                    edge_line.Arrow_target_end = "standard";
                    edge_line.Arrow_type = "line";
                    if (source_target_weight_dicts.Length == 1)
                    {
                        edge_line.Arrow_width = source_target_weight_dicts[0][source_node_line.Name][target_node_line.Name].ToString();
                    }
                    visEdges_list.Add(edge_line);
                }
            }
            VisEdges = visEdges_list.ToArray();
        }
        #endregion

        public void Generate(NetworkBasis_class nw, params Dictionary<string, Dictionary<string, float>>[] source_target_weight_dicts)
        {
            Generate_nodes(nw);
            Generate_edges(nw, source_target_weight_dicts);
            Generate_resource_lines_and_connect_to_nodes();
        }

    }

    //////////////////////////////////////////////////////

    class RegularNW_generate_visualization_class : NetworkBasis_class
    {
        public Visualization_of_regularNW_options Options { get; set; }

        public RegularNW_generate_visualization_class()
        {
            Options = new Visualization_of_regularNW_options();
        }

        public Visualization_of_nw_basis Generate_visualization_instance(NetworkBasis_class nw, params Dictionary<string, Dictionary<string, float>>[] source_target_weight_dicts)
        {
            if (source_target_weight_dicts.Length > 1) { throw new Exception(); }

            base.Deep_copy_into_this(nw);
            Visualization_of_ontologyNW_class visu = new Visualization_of_ontologyNW_class(this.Options);
            visu.Generate(this, source_target_weight_dicts);
            return visu;
        }

    }
}
