using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using UMCCreation;

namespace UMCCreationTests
{
    [TestFixture]
    public class IniReaderTests
    {
        [Test]
        public void TestReadIniFile()
        {
            var path = @"D:\Users\gibb166\Documents\GitTransfer\LC-MS-Feature-Finder\bin\Debug\ExampleSettings.ini";
            var iniReader = new IniReader(path);
            string isos_file = iniReader.ReadString("Files", "InputFileName", "");
            string output_dir = iniReader.ReadString("Files", "OutputDirectory", ".");
            float isotopicFit = iniReader.ReadFloat("DataFilters", "MaxIsotopicFit", 1);
            int intensityFilter = iniReader.ReadInteger("DataFilters", "MinimumIntensity", 500);
            float mflt_mono_mass_start = iniReader.ReadFloat("DataFilters", "MonoMassStart", 0);
            float mflt_mono_mass_end = iniReader.ReadFloat("DataFilters", "MonoMassEnd", float.MaxValue);
            bool mbln_process_chunks = iniReader.ReadBoolean("DataFilters", "ProcessDataInChunks", false);
            int maxPoints = iniReader.ReadInteger("DataFilters", "MaxDataPointsPerChunk", int.MaxValue);
            int chunkSize = iniReader.ReadInteger("DataFilters", "ChunkSize", 3000);
            int imsMinScan = iniReader.ReadInteger("DataFilters", "IMSMinScan", 0);
            int imsMaxScan = iniReader.ReadInteger("DataFilters", "IMSMaxScan", 50000);
            int lcMinScan = iniReader.ReadInteger("DataFilters", "LCMinScan", 0);
            int lcMaxScan = iniReader.ReadInteger("DataFilters", "LCMaxScan", 50000);
            Console.WriteLine();
        }
    }
}
