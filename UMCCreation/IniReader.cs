using System;
using System.Collections.Generic;
using System.Data;
using System.Globalization;
using System.IO;

namespace UMCCreation
{
    public class IniReader
    {
        // considered code from http://www.codeproject.com/Articles/1966/An-INI-file-handling-class-using-C
        // Discarded that choice since the kernel function call is apparently expensive, and reading the file is easy.

        public string IniFilePath { get; }

        /// <summary>
        /// INIFile Constructor.
        /// </summary>
        /// <param name="path"></param>
        public IniReader(string path)
        {
            IniFilePath = path;
            ParseIniFile();
        }

        public int ReadInteger(string szSection, string szKey, int iDefaultValue)
        {
            var result = ReadString(szSection, szKey, iDefaultValue.ToString()).Trim();
            var iResult = int.Parse(result);
            return iResult;
        }

        public float ReadFloat(string szSection, string szKey, float fltDefaultValue)
        {
            var result = ReadString(szSection, szKey, fltDefaultValue.ToString(CultureInfo.InvariantCulture)).Trim();
            var fltResult = float.Parse(result);
            return fltResult;
        }

        public double ReadDouble(string szSection, string szKey, double dblDefaultValue)
        {
            var result = ReadString(szSection, szKey, dblDefaultValue.ToString(CultureInfo.InvariantCulture)).Trim();
            var dblResult = double.Parse(result);
            return dblResult;
        }

        public bool ReadBoolean(string szSection, string szKey, bool bolDefaultValue)
        {
            var result = ReadString(szSection, szKey, bolDefaultValue.ToString()).Trim();
            var bolResult = bool.Parse(result);
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
        private readonly Dictionary<string, Dictionary<string, string>> _iniFileValues = new Dictionary<string, Dictionary<string, string>>();

        /// <summary>
        /// Read the whole file into the <see cref="_iniFileValues"/> dictionary
        /// </summary>
        private void ParseIniFile()
        {
            using (var reader = new StreamReader(new FileStream(IniFilePath, FileMode.Open, FileAccess.Read, FileShare.Read)))
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
