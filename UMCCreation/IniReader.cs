using System;
using System.Collections.Generic;
using System.IO;

namespace UMCCreation
{
    public class IniReader
    {
        // considered code from http://www.codeproject.com/Articles/1966/An-INI-file-handling-class-using-C
        // Discarded that choice since the kernel function call is apparently expensive, and reading the file is easy.

        public string _path;

        /// <summary>
        /// INIFile Constructor.
        /// </summary>
        /// <param name="path"></param>
        public IniReader(string path)
        {
            _path = path;
            ParseIniFile();
        }

        public int ReadInteger(string szSection, string szKey, int iDefaultValue)
        {
            string result = ReadString(szSection, szKey, iDefaultValue.ToString()).Trim();
            int iResult = int.Parse(result);
            return iResult;
        }

        public float ReadFloat(string szSection, string szKey, float fltDefaultValue)
        {
            string result = ReadString(szSection, szKey, fltDefaultValue.ToString()).Trim();
            float fltResult = float.Parse(result);
            return fltResult;
        }

        public double ReadDouble(string szSection, string szKey, double dblDefaultValue)
        {
            string result = ReadString(szSection, szKey, dblDefaultValue.ToString()).Trim();
            double dblResult = double.Parse(result);
            return dblResult;
        }

        public bool ReadBoolean(string szSection, string szKey, bool bolDefaultValue)
        {
            string result = ReadString(szSection, szKey, bolDefaultValue.ToString()).Trim();
            bool bolResult = bool.Parse(result);
            return bolResult;
        }

        public string ReadString(string szSection, string szKey, string szDefaultValue)
        {
            var result = szDefaultValue;
            if (_iniFileValues.ContainsKey(szSection) && _iniFileValues[szSection].ContainsKey(szKey))
            {
                result = _iniFileValues[szSection][szKey];
            }
            return result;
        }

        /// <summary>
        /// Storing the whole file in memory, since it is a small file.
        /// </summary>
        private Dictionary<string, Dictionary<string, string>> _iniFileValues = new Dictionary<string, Dictionary<string, string>>();

        /// <summary>
        /// Read the whole file into the <see cref="_iniFileValues"/> dictionary
        /// </summary>
        private void ParseIniFile()
        {
            using (var reader = new StreamReader(new FileStream(_path, FileMode.Open, FileAccess.Read, FileShare.Read)))
            {
                var section = "";
                _iniFileValues.Add(section, new Dictionary<string, string>());
                while (!reader.EndOfStream)
                {
                    var line = reader.ReadLine();
                    if (string.IsNullOrWhiteSpace(line))
                    {
                        continue;
                    }
                    // Strip out comments
                    if (line.Contains(";"))
                    {
                        var posc = line.IndexOf(';');
                        line = line.Substring(posc).Trim();
                        if (string.IsNullOrWhiteSpace(line))
                        {
                            continue;
                        }
                    }
                    if (line.StartsWith("["))
                    {
                        section = line.Trim('[', ']');
                        if (!_iniFileValues.ContainsKey(section))
                        {
                            _iniFileValues.Add(section, new Dictionary<string, string>());
                        }
                        continue;
                    }
                    var pos = line.IndexOf("=", StringComparison.InvariantCulture);
                    var key = line.Substring(0, pos);
                    var value = line.Substring(pos + 1);
                    _iniFileValues[section][key] = value;
                }
            }
        }
    }
}
