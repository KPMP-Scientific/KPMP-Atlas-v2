// private NetworkBasis_class Generate_subnetwork(int[] subnetwork_indexes_old, List<Node_line_class> seed_nodes_old)
// does not consider duplicated seed_nodes_old


using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Trajectories;

namespace Network
{
    class NetworkBasis_line_class
    {
        #region Fields
        public int[] Targets { get; set; }
        #endregion

        #region Constructor
        public NetworkBasis_line_class()
        {
            Targets = new int[0];
        }
        #endregion

        public void Add_not_existing_targets_and_order(params int[] addTargets)
        {
            List<int> new_targets = new List<int>();
            new_targets.AddRange(Targets);
            new_targets.AddRange(addTargets);
            this.Targets = new_targets.Distinct().OrderBy(l => l).ToArray();
        }
        public NetworkBasis_line_class Deep_copy()
        {
            NetworkBasis_line_class newLine = (NetworkBasis_line_class)this.MemberwiseClone();
            int targets_length = this.Targets.Length;
            newLine.Targets = new int[targets_length];
            for (int indexT=0; indexT<targets_length; indexT++)
            {
                newLine.Targets[indexT] = this.Targets[indexT];
            }
            return newLine;
        }

    }

    class NetworkBasis_class
    {
        #region Fields
        public NetworkBasis_line_class[] NW { get; protected set; }
        public Node_class UniqueNodes { get; private set; }
        #endregion

        #region Constructors
        public NetworkBasis_class()
        {
            UniqueNodes = new Node_class();
            NW = new NetworkBasis_line_class[0];
        }
        #endregion

        #region Generate
        public void Generate_and_overwrite_network(KPMP_cluster_hierarchy_line_class[] cluster_hierarchy_lines, KPMP_cluster_characterization_line_class[] cluster_characterization_lines, Dictionary<string,int> categoryEntity_pieChartOrderNo_dict, Dictionary<string,System.Drawing.Color> categoryEntity_color_dict)
        {
            this.UniqueNodes.Generate_and_overwrite_nodes_index_alphabetically(cluster_characterization_lines, categoryEntity_pieChartOrderNo_dict, categoryEntity_color_dict);
            Dictionary<string,int> sourceName_index_dict = this.UniqueNodes.Get_sourceName_index_dict();
            int cluster_hierarchy_lines_length = cluster_hierarchy_lines.Length;
            cluster_hierarchy_lines = cluster_hierarchy_lines.OrderBy(l => l.Source_cellSubtypeName).ToArray();
            KPMP_cluster_hierarchy_line_class cluster_hierarchy_line;
            NetworkBasis_line_class nwLine;
            int indexClusterHierarchy = 0;
            Node_line_class un_line;
            int sourceNameCompare;
            int un_lastSourceNode_index = UniqueNodes.Order_by_index_and_get_last_source_node_index_plus_one();
            NetworkBasis_line_class[] nw_array = new NetworkBasis_line_class[un_lastSourceNode_index];
            int added_hierarchy_sourceNames_count = 0;
            for (int indexUN = 0; indexUN < un_lastSourceNode_index; indexUN++)
            {
                un_line = UniqueNodes.UNs[indexUN];
                nwLine = new NetworkBasis_line_class();
                sourceNameCompare = -2;
                while ((indexClusterHierarchy < cluster_hierarchy_lines_length) && (sourceNameCompare <= 0))
                {
                    cluster_hierarchy_line = cluster_hierarchy_lines[indexClusterHierarchy];
                    sourceNameCompare = cluster_hierarchy_line.Source_cellSubtypeName.CompareTo(un_line.Name);
                    if (sourceNameCompare < 0)
                    {
                        indexClusterHierarchy++;
                    }
                    else if (sourceNameCompare == 0)
                    {
                        nwLine.Add_not_existing_targets_and_order(sourceName_index_dict[cluster_hierarchy_line.Target_cellSubtypeName]);
                        added_hierarchy_sourceNames_count++;
                        indexClusterHierarchy++;
                    }
                }
                nw_array[indexUN] = nwLine;
            }
            NW = nw_array;
        }
        #endregion

        private void Reindex_network()
        {
            if (UniqueNodes.UNs.Length > 0)
            {
                int nw_count = NW.Length;
                int un_count = UniqueNodes.UNs.Length;
                List<int> transformIndex = new List<int>();
                Dictionary<int, int> oldIndex_newIndex_dict = this.UniqueNodes.Get_oldIndex_index_dict();
                int max_index = this.UniqueNodes.Get_max_index();
                NetworkBasis_line_class nw_line;
                NetworkBasis_line_class[] new_nw = new NetworkBasis_line_class[max_index+1];
                for (int indexUN = 0; indexUN <= max_index; indexUN++)
                {
                    new_nw[indexUN] = new NetworkBasis_line_class();
                }
                int targets_count;
                List<int> new_targets = new List<int>();
                int target;
                for (int indexOldNW = 0; indexOldNW < nw_count; indexOldNW++)
                {
                    nw_line = this.NW[indexOldNW];
                    if (oldIndex_newIndex_dict.ContainsKey(indexOldNW))
                    {
                        targets_count = nw_line.Targets.Length;
                        new_targets = new List<int>();
                        for (int indexT = 0; indexT < targets_count; indexT++)
                        {
                            target = nw_line.Targets[indexT];
                            if (oldIndex_newIndex_dict.ContainsKey(target))
                            {
                                new_targets.Add(oldIndex_newIndex_dict[target]);
                            }
                        }
                        nw_line.Targets = new_targets.OrderBy(l => l).ToArray();
                        if (new_nw[oldIndex_newIndex_dict[indexOldNW]].Targets.Length > 0) { throw new Exception(); }
                        new_nw[oldIndex_newIndex_dict[indexOldNW]] = nw_line;
                    }
                }
                NW = new_nw;
                UniqueNodes.Index_changes_adopted = true;
            }
            else
            {
                throw new Exception();
            }
        }

        public void Add_singleNodes_and_reindex_nw(KPMP_cluster_characterization_line_class[] cluster_characterization_lines, Dictionary<string, int> categoryEntity_pieChartOrderNo_dict, Dictionary<string, System.Drawing.Color> categoryEntity_color_dict)
        {
            this.UniqueNodes.Add_nodes_and_reindex_alphabetically(cluster_characterization_lines, categoryEntity_pieChartOrderNo_dict, categoryEntity_color_dict);
            Reindex_network();
        }

        #region Copy, read, write
        public void Deep_copy_into_this(NetworkBasis_class other)
        {
            UniqueNodes.Deep_copy_into_this(other.UniqueNodes);
            int nw_length = other.NW.Length;
            NW = new NetworkBasis_line_class[nw_length];
            for (int indexNW = 0; indexNW < nw_length; indexNW++)
            {
                this.NW[indexNW] = other.NW[indexNW].Deep_copy();
            }
        }
        public NetworkBasis_class Shallow_copy()
        {
            NetworkBasis_class copy = (NetworkBasis_class)this.MemberwiseClone();
            return copy;
        }
        #endregion
    }
}
