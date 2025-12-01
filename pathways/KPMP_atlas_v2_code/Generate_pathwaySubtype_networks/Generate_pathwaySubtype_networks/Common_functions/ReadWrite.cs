using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.IO.Compression;
using System.Reflection;
using System.ComponentModel;
using System.Diagnostics;
using System.Text;
using Common_functions;

namespace ReadWrite
{
    public abstract class ReadWriteOptions_base
    {
        private bool add_to_existing_file { get; set; }
        public string File { get; set; }
        public string[] Key_columnNames { get; set; }
        public string[] Key_propertyNames { get; set; }
        public int[] Key_columnIndexes { get; set; }
        public bool File_has_headline { get; set; }
        public char Delimiter { get; set; }
        public int Empty_integer_value { get; set; }
        public string Empty_string_value { get; set; }
        public bool Report_unhandled_null_entries { get; set; }
        public bool Add_to_existing_file
        {
            get { return add_to_existing_file; }
            set
            {
                this.add_to_existing_file = value;
                if ((this.add_to_existing_file == true) && (System.IO.File.Exists(this.File)))
                {
                    File_has_headline = false;
                }
            }
        }
        public bool Check_for_duplicated_columnNames_in_file { get; set; }

        public ReadWriteOptions_base()
        {
            Key_columnIndexes = new int[0];
            Key_columnNames = new string[0];
            Key_propertyNames = new string[0];
            File = "error";
            Report_unhandled_null_entries = true;
            Empty_string_value = Global_class.Empty_entry;
            Add_to_existing_file = false;
            Check_for_duplicated_columnNames_in_file = true;
        }
    }

    class ReadWriteClass
    {
        public static string Get_writeLine_from_array<T>(T[] array, char delimiter)
        {
            StringBuilder stringBuild = new StringBuilder();
            int array_length = array.Length;
            for (int i = 0; i < array_length; i++)
            {
                if (i == 0) { stringBuild.AppendFormat("{0}", array[i]); }
                else { stringBuild.AppendFormat("{0}{1}", delimiter, array[i]); }
            }
            return stringBuild.ToString();
        }
        public static T[] Get_array_from_readLine<T>(string readLine, char delimiter)
        {
            if (String.IsNullOrEmpty(readLine)) { return new T[0]; }
            else
            {
                string[] split = readLine.Split(delimiter);
                int split_length = split.Length;
                if (string.IsNullOrEmpty(split[split_length - 1])) { split_length--; }
                var TType = typeof(T);
                T[] array = new T[split_length];
                for (int i = 0; i < split_length; i++)
                {
                    array[i] = (T)Convert.ChangeType(split[i], TType);
                }
                return array;
            }
        }
        private static string[] Get_and_modify_columnNames(string headline, ReadWriteOptions_base Options)
        {
            List<string> columnNamesList = new List<string>();
            columnNamesList.AddRange(headline.Split(Options.Delimiter));
            string[] columnNames = columnNamesList.ToArray();
            if ((Options.Check_for_duplicated_columnNames_in_file) && (columnNames.Distinct().ToArray().Length != columnNames.Length))
            {
                columnNames = columnNames.OrderBy(l => l).ToArray();
                int columnNames_length = columnNames.Length;
                string this_columnName;
                string previous_columnName;
                List<string> duplicated_columnNames = new List<string>();
                for (int indexC = 1; indexC < columnNames_length; indexC++)
                {
                    this_columnName = columnNames[indexC];
                    previous_columnName = columnNames[indexC - 1];
                    if (this_columnName.Equals(previous_columnName))
                    {
                        duplicated_columnNames.Add(this_columnName);
                    }
                }
                throw new Exception();
            }
            return columnNames;
        }
        private static int[] Get_columnIndexes_of_given_columnNames<T>(string[] columnNames, params string[] given_columnNames)
        {
            int given_length = given_columnNames.Length;
            if (given_length == 0)
            {
                throw new Exception(typeof(T).Name + ": Get_columnIndexes_of_given_columnNames: no given_columNames");
            }
            int[] columnIndexes = new int[given_length];
            for (int i = 0; i < given_length; i++)
            {
                int index = Array.IndexOf(columnNames, given_columnNames[i]);
                if (index >= 0) { columnIndexes[i] = index; }
                else
                {
                    string missing = given_columnNames[i];
                    throw new Exception(typeof(T).Name + ": columnName \"" + missing + "\" does not exist");
                }
            }
            return columnIndexes;
        }

        private static int[] Get_propertyIndexes_of_corresponding_given_columnNames<T>(PropertyInfo[] propInfo, string[] propertyNames, string[] given_columnNames, string[] search_given_columnNames)
        {
            int search_length = search_given_columnNames.Length;
            int[] columnNames_indexes = new int[search_length];
            if (search_length == 0)
            {
                throw new Exception(typeof(T).Name + ": Get_propertyIndexes_of_corresponding_given_columnNames: no search_given_columnNames");
            }
            List<string> missing_columnNames = new List<string>();
            for (int i = 0; i < search_length; i++)
            {
                int index = Array.IndexOf(given_columnNames, search_given_columnNames[i]);
                if (index >= 0) { columnNames_indexes[i] = index; }
                else
                {
                    missing_columnNames.Add(given_columnNames[i]);
                }
            }
            if (missing_columnNames.Count > 0) { throw new Exception(typeof(T).Name + ": Get_propertyIndexes_of_corresponding_given_columnNames: missing column names"); }
            string[] corresponding_propertyNames = new string[search_length];
            for (int indexS = 0; indexS < search_length; indexS++)
            {
                corresponding_propertyNames[indexS] = propertyNames[columnNames_indexes[indexS]];
            }
            int[] propertyIndexes = Get_propertyIndexes<T>(propInfo, corresponding_propertyNames);
            return propertyIndexes;
        }
        private static int[] Get_propertyIndexes<T>(PropertyInfo[] propInfo, string[] key_propertyNames)
        {
            int key_length = key_propertyNames.Length;
            int[] propertyIndexes = new int[key_length];
            string[] propInfo_names = new string[propInfo.Length];

            List<string> missing_propertyNames = new List<string>();

            for (int i = 0; i < propInfo.Length; i++)
            {
                propInfo_names[i] = propInfo[i].Name;
            }

            for (int i = 0; i < key_length; i++)
            {
                int index = Array.IndexOf(propInfo_names, key_propertyNames[i]);
                if (index >= 0) { propertyIndexes[i] = index; }
                if (index < 0)
                {
                    missing_propertyNames.Add(key_propertyNames[i]);
                }
            }
            if (missing_propertyNames.Count > 0) { throw new Exception(typeof(T).Name + ": Get_propertyIndexes: missing property names"); }
            return propertyIndexes;
        }
        public static void Create_subdirectory_if_it_does_not_exist(string actDirectory, string sub_directory_name)
        {
            if (!Directory.Exists(actDirectory + sub_directory_name))
            {
                DirectoryInfo dir = new DirectoryInfo(actDirectory);
                dir.CreateSubdirectory(sub_directory_name);
            }
        }
        public static void Create_directory_if_it_does_not_exist(string input_directory)
        {
            string directory = Path.GetDirectoryName(input_directory);
            if (!String.IsNullOrEmpty(directory))
            {
                if (!Directory.Exists(directory))
                {
                    Directory.CreateDirectory(directory);
                }
            }
        }
        public static List<T> ReadRawData_and_FillList<T>(ReadWriteOptions_base options) where T : class
        {
            FileInfo file = new FileInfo(options.File);
            StreamReader stream = file.OpenText();
            List<T> Data = ReadRawData_and_FillList<T>(stream, options, options.File);
            stream.Close();
            return Data;
        }
        public static List<T> ReadRawData_and_FillList<T>(StreamReader stream, ReadWriteOptions_base options, string file_name) where T : class
        {
            PropertyInfo[] propInfo = typeof(T).GetProperties();

            #region Determine columns to be safed and invalidLine_defining columns and properties
            //Read headline, if it exists, determine indexes of columns to be safed in list
            //Begin
            string[] columnNames = { Global_class.Empty_entry };
            int[] columnIndexes;
            int[] invalidLine_defining_columnIndexes = new int[0];
            int[] invalidLine_defining_popertyIndexes = new int[0];
            int[] propertyIndexes;

            if (options.File_has_headline)
            {
                string headline = stream.ReadLine();
                columnNames = Get_and_modify_columnNames(headline, options);
                columnIndexes = Get_columnIndexes_of_given_columnNames<T>(columnNames, options.Key_columnNames);
                options.Key_columnIndexes = columnIndexes;
            }
            else { columnIndexes = options.Key_columnIndexes; }
            propertyIndexes = Get_propertyIndexes<T>(propInfo, options.Key_propertyNames);
            if (columnIndexes.Length != propertyIndexes.Length)
            {
                throw new Exception(typeof(T).Name + "Length columnIndexes (Key_columnNames/columnIndexes) != propertyIndexes (Key_propertyNames)\"");
            }
            //End
            #endregion

            #region Generate and fill list
            List<T> Data = new List<T>();
            var TType = typeof(T);

            int invalidLine_defining_columnIndexes_length = invalidLine_defining_columnIndexes.Length;
            string inputLine;
            int readLines = 0;
            int safedLines = 0;
            int colIndex;
            int propIndex;
            bool safeLine;
            bool safeLine_and_stop = false;
            bool report_check_lineDelimiters = false;
            bool valid;
            string invalidLineDefiningColumnEntry;
            int line_count = 0;
            bool safeCondition_encountered = false;

            while ((!safeLine_and_stop) && ((inputLine = stream.ReadLine()) != null))
            {
                if ((inputLine.Length > 0) && (!inputLine.Substring(0, 5).Equals("-----")))
                {
                    line_count++;
                    //line_count++;
                    string[] columnEntries = inputLine.Split(options.Delimiter);
                    if (columnEntries.Length == 1)
                    {
                        report_check_lineDelimiters = true;
                    }
                    safeLine = true;
                    valid = true;
                    for (int indexIndex = 0; indexIndex < invalidLine_defining_columnIndexes_length; indexIndex++)
                    {
                        invalidLineDefiningColumnEntry = columnEntries[invalidLine_defining_columnIndexes[indexIndex]];
                        try
                        {
                            var obj = Convert.ChangeType(invalidLineDefiningColumnEntry, propInfo[invalidLine_defining_popertyIndexes[indexIndex]].PropertyType);
                            valid = true;
                        }
                        catch (InvalidCastException)
                        {
                            valid = false;
                        }
                        catch (FormatException)
                        {
                            valid = false;
                        }
                        catch (OverflowException)
                        {
                            valid = false;
                        }
                        catch (ArgumentNullException)
                        {
                            valid = false;
                        }
                    }
                    if (((safeLine)
                           || (safeLine_and_stop))
                        && (valid))
                    {
                        T newLine = (T)Activator.CreateInstance(TType);
                        for (int i = 0; i < columnIndexes.Length; i++)
                        {
                            colIndex = columnIndexes[i];
                            propIndex = propertyIndexes[i];
                            if (columnEntries[colIndex].Contains(Global_class.Space_text.ToString()))
                            {
                                columnEntries[colIndex] = (string)columnEntries[colIndex].Replace(Global_class.Space_text.ToString(), "");
                            }
                            if (columnEntries[colIndex] == "#DIV/0!") { columnEntries[colIndex] = "NaN"; }
                            if (propInfo[propIndex].PropertyType.IsEnum)
                            {

                                try { propInfo[propIndex].SetValue(newLine, Enum.Parse(propInfo[propIndex].PropertyType, columnEntries[colIndex]), null); }
                                catch
                                {
                                    columnEntries[colIndex] = char.ToUpper(columnEntries[colIndex][0]) + columnEntries[colIndex].ToLower().Substring(1);
                                    propInfo[propIndex].SetValue(newLine, Enum.Parse(propInfo[propIndex].PropertyType, columnEntries[colIndex]), null);
                                }
                            }
                            else if (string.IsNullOrEmpty(columnEntries[colIndex]))
                            {
                                if (propInfo[propIndex].PropertyType == typeof(int))
                                {
                                    propInfo[propIndex].SetValue(newLine, options.Empty_integer_value, null);
                                }
                                else if (propInfo[propIndex].PropertyType == typeof(string))
                                {
                                    propInfo[propIndex].SetValue(newLine, "", null);
                                }
                                else if (options.Report_unhandled_null_entries)
                                {
                                    throw new Exception(typeof(ReadWriteClass).Name + "ReadRawData_and_FillList: " + options.File + " unhandled null entry");
                                }
                            }
                            else
                            {
                                if ((columnEntries[colIndex] != "") && ((columnEntries[colIndex] != "NA") || (propInfo[propIndex].PropertyType == typeof(string))))
                                {
                                    {
                                        if (columnEntries[colIndex].Equals("Inf"))
                                        {
                                            if (propInfo[propIndex].PropertyType.Equals(typeof(double)))
                                            {
                                                propInfo[propIndex].SetValue(newLine, double.PositiveInfinity, null);
                                            }
                                            else if (propInfo[propIndex].PropertyType.Equals(typeof(float)))
                                            {
                                                propInfo[propIndex].SetValue(newLine, float.PositiveInfinity, null);
                                            }
                                            else if (propInfo[propIndex].PropertyType.Equals(typeof(int)))
                                            {
                                                throw new Exception();
                                            }
                                        }
                                        else if (columnEntries[colIndex].Equals("-Inf"))
                                        {
                                            if (propInfo[propIndex].PropertyType.Equals(typeof(double)))
                                            {
                                                propInfo[propIndex].SetValue(newLine, double.NegativeInfinity, null);
                                            }
                                            else if (propInfo[propIndex].PropertyType.Equals(typeof(float)))
                                            {
                                                propInfo[propIndex].SetValue(newLine, float.NegativeInfinity, null);
                                            }
                                            else if (propInfo[propIndex].PropertyType.Equals(typeof(int)))
                                            {
                                                throw new Exception();
                                            }
                                        }
                                        else
                                        {
                                            propInfo[propIndex].SetValue(newLine, Convert.ChangeType(columnEntries[colIndex], propInfo[propIndex].PropertyType), null);
                                        }
                                    }
                                }
                            }
                        }
                        Data.Add(newLine);
                        safedLines = safedLines + 1;
                    }
                    readLines = readLines + 1;
                }
            }
            #endregion

            #region Final report
            if (report_check_lineDelimiters)
            {
                throw new Exception(typeof(T).Name + "{0}: only one column entry: Check lineDelimiters");
            }
            #endregion

            return Data;
        }
        public static T[] ReadRawData_and_FillArray<T>(ReadWriteOptions_base Options) where T : class
        {
            return ReadRawData_and_FillList<T>(Options).ToArray();
        }
        public static void WriteData<T>(List<T> Data, ReadWriteOptions_base Options) where T : class
        {
            ReadWriteClass.Create_directory_if_it_does_not_exist(Options.File);
            StreamWriter writer = new StreamWriter(Options.File, false);
            WriteData(Data, Options, writer);
            writer.Close();
        }
        public static void WriteData<T>(T[] Data, ReadWriteOptions_base Options) where T : class
        {
            string directory = Path.GetDirectoryName(Options.File);
            Create_directory_if_it_does_not_exist(directory + "\\");
            StreamWriter writer = new StreamWriter(Options.File, Options.Add_to_existing_file);
            WriteData(Data, Options, writer);
            writer.Close();
        }
        public static void WriteData<T>(T[] Data, ReadWriteOptions_base Options, StreamWriter writer) where T : class
        {
            WriteData(Data.ToList(), Options, writer);
        }
        public static void WriteData<T>(List<T> Data, ReadWriteOptions_base Options, StreamWriter writer) where T : class
        {
            PropertyInfo[] propInfo = typeof(T).GetProperties();
            PropertyInfo prop;

            int[] propertyIndexes = Get_propertyIndexes<T>(propInfo, Options.Key_propertyNames);

            //Generate and write Headline
            int propertyIndexes_length = propertyIndexes.Length;
            if ((Options.File_has_headline == true))
            {
                char headline_delimiter = Options.Delimiter;
                StringBuilder headline = new StringBuilder();
                for (int index = 0; index < propertyIndexes_length; index++)
                {
                    if (index < propertyIndexes_length - 1)
                    {
                        headline.AppendFormat("{0}{1}", Options.Key_columnNames[index], headline_delimiter);
                    }
                    else
                    {
                        headline.AppendFormat("{0}", Options.Key_columnNames[index]);
                    }
                }
                writer.WriteLine(headline);
            }

            //Generate and write lines
            char line_delimiter = Options.Delimiter;
            StringBuilder line = new StringBuilder();
            int data_count = Data.Count;
            for (int lineIndex = 0; lineIndex < data_count; lineIndex++)
            {
                line.Clear();
                for (int index = 0; index < propertyIndexes_length; index++)
                {
                    prop = propInfo[propertyIndexes[index]];
                    if (index < propertyIndexes_length - 1) { line.AppendFormat("{0}{1}", prop.GetValue(Data[lineIndex], null), line_delimiter); }
                    else { line.AppendFormat("{0}", prop.GetValue(Data[lineIndex], null)); }
                }
                writer.WriteLine(line);
            }
            writer.Close();
        }
    }
}