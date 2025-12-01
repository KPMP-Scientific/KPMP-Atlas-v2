using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Reflection;
using System.Drawing;
using Common_functions;
using ReadWrite;
using Trajectories;
using Network_visualization;
using System.Xml.Serialization;

namespace Network
{
    class Node_color_value_line_class
    {
        public Color Node_color { get; set; }
        public float Node_color_value { get; set; }
        public int PieChart_sliceOrderNo { get; set; }

        public Node_color_value_line_class Deep_copy()
        {
            Node_color_value_line_class copy = (Node_color_value_line_class)this.MemberwiseClone();
            return copy;
        }
    }

    class Node_line_class
    {
        #region Fields
        public string Name { get; internal set; }
        public Node_color_value_line_class[] Node_color_values { get; set; }
        public int Index { get; internal set; }
        public int IndexOld { get; internal set; }
        #endregion

        public Node_line_class()
        {
            Name = Global_class.Empty_entry;
            Index = -1;
            IndexOld = -1;
            Node_color_values = new Node_color_value_line_class[0];
        }

        public Node_line_class Deep_copy()
        {
            Node_line_class newLine = (Node_line_class)this.MemberwiseClone();
            newLine.Name = (string)Name.Clone();
            int node_color_values_length = this.Node_color_values.Length;
            newLine.Node_color_values = new Node_color_value_line_class[node_color_values_length];
            for (int indexCV = 0; indexCV < node_color_values_length; indexCV++)
            {
                newLine.Node_color_values[indexCV] = this.Node_color_values[indexCV].Deep_copy();
            }
            return newLine;
        }
    }

    class Node_writeOptions_class : ReadWriteOptions_base
    {
        const char array_delimiter = ',';

        public static char Get_array_delimiter()
        {
            return array_delimiter;
        }

        public Node_writeOptions_class(string directory, string fileName, bool report)
        {
            string completeFileName = directory + fileName;
            File = completeFileName;
            Key_propertyNames = new string[] { "Index", "Name", "Name2", "Classification", "NW_classification", "DE", "ReadWrite_timepoints", "ReadWrite_log2_fc_changes", "Minus_log10_p_value", "Seed_node_neighbor_ratio", "Degree", "Node_colors" };
            Key_columnNames = Key_propertyNames;
            File_has_headline = true;
            Delimiter = Global_class.Tab;
        }
    }

    class Node_class
    {
        #region Fields
        public Node_line_class[] UNs { get; set;}
        public bool Index_changes_adopted { get; set; }
        #endregion

        public Node_class()
        {
            UNs = new Node_line_class[0];
            Index_changes_adopted = true;
        }

        public Dictionary<int, int> Get_oldIndex_index_dict()
        {
            Dictionary<int, int> oldIndex_newIndex_dict = new Dictionary<int, int>();
            int un_count = this.UNs.Length;
            Node_line_class un_line;
            for (int indexUN = 0; indexUN < un_count; indexUN++)
            {
                un_line = this.UNs[indexUN];
                if (un_line.IndexOld != -1)
                {
                    oldIndex_newIndex_dict.Add(un_line.IndexOld, un_line.Index);
                }
            }
            return oldIndex_newIndex_dict;
        }

        public float Get_max_nodeColor_value()
        {
            int un_count = this.UNs.Length;
            float max_nodeColor_value = -1;
            float current_nodeColor_value;
            Node_line_class un_line;
            for (int indexUN = 0; indexUN < un_count; indexUN++)
            {
                un_line = this.UNs[indexUN];
                current_nodeColor_value = 0;
                foreach (Node_color_value_line_class color_line in un_line.Node_color_values)
                {
                    current_nodeColor_value += color_line.Node_color_value;
                }
                if (max_nodeColor_value < current_nodeColor_value)
                { 
                    max_nodeColor_value = current_nodeColor_value;
                }
            }
            return max_nodeColor_value;
        }
        public int Get_max_index()
        {
            int un_count = this.UNs.Length;
            int max_index = -1;
            Node_line_class un_line;
            for (int indexUN = 0; indexUN < un_count; indexUN++)
            {
                un_line = this.UNs[indexUN];
                if (un_line.Index > max_index)
                {
                    max_index = un_line.Index;
                }
            }
            return max_index;
        }

        private void Add_to_array_and_reindex_alphabetically(Node_line_class[] add_UNs)
        {
            int this_UNs_length = this.UNs.Length;
            int add_UNs_length = add_UNs.Length;
            int new_UNs_length = this_UNs_length + add_UNs_length;
            Node_line_class[] new_UNs = new Node_line_class[new_UNs_length];
            int indexNew = -1;
            for (int indexThis = 0; indexThis < this_UNs_length; indexThis++)
            {
                indexNew++;
                new_UNs[indexNew] = this.UNs[indexThis];
            }
            for (int indexAdd = 0; indexAdd < add_UNs_length; indexAdd++)
            {
                indexNew++;
                new_UNs[indexNew] = add_UNs[indexAdd];
            }
            this.UNs = new_UNs;
            Index_nodes_alphabetically_and_set_index_old();
        }

        #region Order Nodes
        public void Order_by_name()
        {
            Dictionary<string, List<Node_line_class>> name_dict = new Dictionary<string, List<Node_line_class>>();
            int un_length = UNs.Length;
            Node_line_class node_line;
            for (int indexNode = 0; indexNode < un_length; indexNode++)
            {
                node_line = this.UNs[indexNode];
                if (!name_dict.ContainsKey(node_line.Name))
                {
                    name_dict.Add(node_line.Name, new List<Node_line_class>());
                }
                name_dict[node_line.Name].Add(node_line);
            }
            string[] names = name_dict.Keys.ToArray();
            string name;
            int names_length = names.Length;
            names = names.OrderBy(l => l).ToArray();
            List<Node_line_class> ordered_nodes = new List<Node_line_class>();
            for (int indexName = 0; indexName < names_length; indexName++)
            {
                name = names[indexName];
                ordered_nodes.AddRange(name_dict[name]);
            }
            UNs = ordered_nodes.ToArray();
        }
        public void Order_by_indexOld()
        {
            UNs = UNs.OrderBy(l => l.IndexOld).ToArray();
        }
        public void Order_by_index()
        {
            UNs = UNs.OrderBy(l => l.Index).ToArray();
        }
        public int Order_by_index_and_get_last_source_node_index_plus_one()
        {
            int pos;
            int lastSourceNodeIndex;
            Order_by_index_and_get_lastSourceIndex_and_firstTargetIndex(out lastSourceNodeIndex, out pos);
            return lastSourceNodeIndex + 1;
        }
        public void Order_by_index_and_get_lastSourceIndex_and_firstTargetIndex(out int lastSourceIndex, out int firstTargetIndex)
        {
            Order_by_index();
            lastSourceIndex = UNs.Length - 1;
            if (!lastSourceIndex.Equals(UNs[lastSourceIndex].Index)) { throw new Exception(); }
            firstTargetIndex = 0;
            if (!firstTargetIndex.Equals(UNs[firstTargetIndex].Index)) { throw new Exception(); }
        }
        #endregion

        #region Clear
        public void Clear_unique_nodes()
        {
            UNs = new Node_line_class[0];
            Index_changes_adopted = true;
        }
        #endregion

        public void Generate_and_overwrite_nodes_index_alphabetically(KPMP_cluster_characterization_line_class[] cluster_characterization_lines, Dictionary<string,int> categoryEntity_pieCharSliceOrderNo_dict, Dictionary<string,System.Drawing.Color> categoryEntity_color_dict)
        {
            this.UNs = new Node_line_class[0];
            Add_nodes_and_reindex_alphabetically(cluster_characterization_lines, categoryEntity_pieCharSliceOrderNo_dict, categoryEntity_color_dict);
        }

        public void Add_nodes_and_reindex_alphabetically(KPMP_cluster_characterization_line_class[] cluster_characterization_lines, Dictionary<string, int> categoryEntity_pieCharSliceOrderNo_dict, Dictionary<string, System.Drawing.Color> categoryEntity_color_dict)
        {
            KPMP_category_for_cluster_enum category = cluster_characterization_lines[0].Category;
            int cluster_characterization_lines_length = cluster_characterization_lines.Length;
            cluster_characterization_lines = cluster_characterization_lines.OrderBy(l => l.Cell_subtype_name).ThenBy(l => l.Category_entity).ToArray();
            KPMP_cluster_characterization_line_class cluster_characterization_line;
            Node_line_class un_line;
            List<Node_line_class> nodes = new List<Node_line_class>();
            Node_color_value_line_class node_color_line;
            List<Node_color_value_line_class> node_color_lines = new List<Node_color_value_line_class>();
            for (int indexChar = 0; indexChar < cluster_characterization_lines_length; indexChar++)
            {
                cluster_characterization_line = cluster_characterization_lines[indexChar];
                if (!cluster_characterization_line.Category.Equals(category)) { throw new Exception(); }
                if ((indexChar != 0)
                    && (cluster_characterization_line.Cell_subtype_name.Equals(cluster_characterization_lines[indexChar - 1].Cell_subtype_name))
                    && (cluster_characterization_line.Category_entity.Equals(cluster_characterization_lines[indexChar - 1].Category_entity)))
                {
                    throw new Exception();
                }
                if ((indexChar == 0)
                    || (!cluster_characterization_line.Cell_subtype_name.Equals(cluster_characterization_lines[indexChar - 1].Cell_subtype_name)))
                {
                    node_color_lines.Clear();
                }
                node_color_line = new Node_color_value_line_class();
                node_color_line.Node_color_value = cluster_characterization_line.Category_entity_value;
                if (categoryEntity_pieCharSliceOrderNo_dict.ContainsKey(cluster_characterization_line.Category_entity))
                { node_color_line.PieChart_sliceOrderNo = categoryEntity_pieCharSliceOrderNo_dict[cluster_characterization_line.Category_entity]; }
                if (categoryEntity_color_dict.ContainsKey(cluster_characterization_line.Category_entity))
                { node_color_line.Node_color = categoryEntity_color_dict[cluster_characterization_line.Category_entity]; }
                else { node_color_line.Node_color = System.Drawing.Color.Black; }
                node_color_lines.Add(node_color_line);
                if ((indexChar == cluster_characterization_lines_length - 1)
                    || (!cluster_characterization_line.Cell_subtype_name.Equals(cluster_characterization_lines[indexChar + 1].Cell_subtype_name)))
                {
                    un_line = new Node_line_class();
                    un_line.Node_color_values = node_color_lines.ToArray();
                    un_line.Name = (string)cluster_characterization_line.Cell_subtype_name.Clone();
                    nodes.Add(un_line);
                }
            }
            Add_to_array_and_reindex_alphabetically(nodes.ToArray());
        }
        public Dictionary<string, int> Get_sourceName_index_dict()
        {
            Dictionary<string, int> sourceName_index_dict = new Dictionary<string, int>();
            foreach (Node_line_class node_line in this.UNs)
            {
                sourceName_index_dict.Add(node_line.Name, node_line.Index);
            }
            return sourceName_index_dict;
        }


        #region Get individual nodes
        public Node_line_class Get_individual_node_line(int indexUN)
        {
            Node_line_class un_line = UNs[indexUN];
            if (un_line.Index != indexUN)
            {
                throw new Exception(typeof(Node_class).Name + " Get_individual_node_line: Indexes do not match ( " + un_line.Index + " <-> " + indexUN + "), UniqueNodes not ordered by index");
            }
            else
            {
                return un_line; ;
            }
        }
        #endregion

        #region Check
        public bool Small_correctness_check()
        {
            bool everything_correct = true;
            if (UNs.Length == 0)
            {
                everything_correct = false;
                throw new Exception();
            }
            if ((Index_changes_adopted == false))
            {
                everything_correct = false;
                throw new Exception();
            }
            //Test for duplicated nodes
            Order_by_name();
            int un_count = UNs.Length;
            Node_line_class un_line;
            for (int indexUN = 0; indexUN < un_count; indexUN++)
            {
                un_line = UNs[indexUN];
                if ((indexUN != 0) && (un_line.Name.Equals(UNs[indexUN - 1].Name)))
                {
                    everything_correct = false;
                    throw new Exception();
                }
            }
            return everything_correct;
        }
        #endregion

        public void Index_nodes_alphabetically_and_set_index_old()
        {
            if (UNs.Length > 0)
            {
                UNs = UNs.OrderBy(l => l.Name).ToArray();
                int newUN_count = UNs.Length;
                for (int i = 0; i < newUN_count; i++)
                {
                    UNs[i].IndexOld = UNs[i].Index;
                    UNs[i].Index = i;
                }
            }
            Index_changes_adopted = false;
        }

        #region Copy and output
        public void Deep_copy_into_this_without_UN(Node_class other)
        {
            Index_changes_adopted = other.Index_changes_adopted;
            UNs = new Node_line_class[0];
        }

        public void Deep_copy_into_this(Node_class other)
        {
            Deep_copy_into_this_without_UN(other);
            int un_length = other.UNs.Length;
            this.UNs = new Node_line_class[un_length];
            for (int indexUN = 0; indexUN < un_length; indexUN++)
            {
                this.UNs[indexUN] = other.UNs[indexUN].Deep_copy();
            }
        }
        public Node_class Deep_copy_without_UN()
        {
            Node_class newUN = (Node_class)this.MemberwiseClone();
            newUN.UNs = new Node_line_class[0];
            return newUN;
        }
        public Node_class Deep_copy()
        {
            Small_correctness_check();
            Node_class newUN = Deep_copy_without_UN();
            int uns_length = this.UNs.Length;
            Node_line_class[] newUNs = new Node_line_class[uns_length];
            for (int indexNode=0; indexNode<uns_length; indexNode++)
            {
                newUN.UNs[indexNode] = this.UNs[indexNode].Deep_copy();
            }
            return newUN;
        }
        #endregion
    }

    ////////////////////////////////////////////////////////////////
}
