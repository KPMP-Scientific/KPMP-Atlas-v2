using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Average_DEGs_and_do_postHocPowerAnalysis.Common_functions
{
    class Report_class
    {
        public static string Get_complete_reportFileName(string report_directory, string assayCellTypeSubtypeLevel)
        {
            return report_directory + "Progress_report_" + assayCellTypeSubtypeLevel + ".txt";
        }
        public static string Get_complete_warining_fileName(string report_directory, string assayCellTypeSubtypeLevel)
        {
            return report_directory + "Warnings_" + assayCellTypeSubtypeLevel + ".txt";
        }
        public static void Delete_report_file_if_exists(string report_directory, string assayCellTypeSubtypeLevel)
        {
            string complete_report_fileName = Get_complete_reportFileName(report_directory, assayCellTypeSubtypeLevel);
            if (File.Exists(complete_report_fileName)) { File.Delete(complete_report_fileName); }
        }
        public static void Add_to_report_file(string report_directory, string assayCellTypeSubtypeLevel, string report_text)
        {
            string complete_report_fileName = Get_complete_reportFileName(report_directory, assayCellTypeSubtypeLevel);
            bool add_to_existing_file = false;
            if (File.Exists(complete_report_fileName)) { add_to_existing_file = true; }
            StreamWriter writer = new StreamWriter(complete_report_fileName, add_to_existing_file);
            string complete_report_text = assayCellTypeSubtypeLevel + ": " + DateTime.Now.ToString("HH:mm") + ": " + report_text;
            writer.WriteLine(complete_report_text);
            writer.Close();
        }
        public static void Delete_warning_file_if_exists(string report_directory, string assayCellTypeSubtypeLevel)
        {
            string complete_warning_fileName = Get_complete_warining_fileName(report_directory, assayCellTypeSubtypeLevel);
            if (File.Exists(complete_warning_fileName)) { File.Delete(complete_warning_fileName); }
        }
        public static void Add_to_warning_file(string report_directory, string assayCellTypeSubtypeLevel, string report_text)
        {
            string complete_warning_fileName = Get_complete_warining_fileName(report_directory, assayCellTypeSubtypeLevel);
            bool add_to_existing_file = false;
            if (File.Exists(complete_warning_fileName)) { add_to_existing_file = true; }
            StreamWriter writer = new StreamWriter(complete_warning_fileName, add_to_existing_file);
            writer.WriteLine(report_text);
            writer.Close();
        }
    }
}
