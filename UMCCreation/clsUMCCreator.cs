using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace UMCCreation
{
    /// <summary>
    /// IsotopePeak object
    /// </summary>
    public class clsIsotopePeak
    {
        public int mint_original_index;
        public int mint_umc_index;
        public int mint_lc_scan;
        public short mshort_charge;
        public double mdbl_abundance;
        public double mdbl_mz;
        public float mflt_fit;
        public double mdbl_average_mass;
        public double mdbl_mono_mass;
        public double mdbl_max_abundance_mass;
        public double mdbl_i2_abundance;
        public float mflt_ims_drift_time;
    }

    /// <summary>
    /// For status reporting
    /// </summary>
    public enum enmStatus
    {
        IDLE = 0,
        LOADING,
        CLUSTERING,
        SUMMARIZING,
        FAILED,
        COMPLETE,
        CHUNKING
    }

    /// <summary>
    /// Entry point for UMCCreation
    /// </summary>
    public class clsUMCCreator
    {
        private int mint_percent_done;
        private bool mbln_process_chunks;
        private float mflt_mono_mass_start;
        private float mflt_mono_mass_end;
        private int mint_mono_mass_overlap;

        private UMCCreator mobj_umc_creator;

        private StreamWriter mfile_logFile;
        private string mstr_baseFileName;

        /// <summary>
        /// Creates the baseFileName that will be used for properly naming the output files.
        ///	- Desired Output format: baseFileName_LCMSFeatures.txt and baseFileName_LCMSFeatureToPeakMap.txt
        /// </summary>
        /// <param name="directoryName"></param>
        /// <param name="inputFileName"></param>
        /// <returns></returns>
        private string CreateBaseFileName(string directoryName, string inputFileName)
        {
            // Remove any directory information from the filename
            var inputFileNameTok = Path.GetFileName(inputFileName);

            // Remove the "_isos.csv" file extension
            inputFileNameTok = inputFileNameTok.Replace("_isos.csv", "");

            // Append the Input File Name to the Output Directory to create the Base File Name
            return Path.Combine(directoryName, inputFileNameTok);
        }

        private void createLogFile()
        {
            string logFileName = mstr_baseFileName + "_FeatureFinder_Log.txt";
            mfile_logFile =
                new StreamWriter(new FileStream(logFileName, FileMode.Create, FileAccess.Write, FileShare.Read));
        }

        private void log(string textToLog)
        {
            mfile_logFile.WriteLine("{0:MM/dd/yyyy HH:mm:ss}\t{1}", DateTime.Now, textToLog);
            mfile_logFile.Flush();
            System.Console.WriteLine(textToLog);
        }

        private void log(string textToLog, int numToLog)
        {
            mfile_logFile.WriteLine("{0:MM/dd/yyyy HH:mm:ss}\t{1}{2}", DateTime.Now, textToLog, numToLog);
            mfile_logFile.Flush();
            System.Console.WriteLine(textToLog + numToLog);
        }

        private void log(string textToLog, float numToLog)
        {
            mfile_logFile.WriteLine("{0:MM/dd/yyyy HH:mm:ss}\t{1}{2,4:F4}", DateTime.Now, textToLog, numToLog);
            mfile_logFile.Flush();
            System.Console.WriteLine(textToLog + numToLog);
        }

        public short PercentComplete
        {
            get
            {
                switch (Status)
                {
                    case enmStatus.LOADING:
                        if (mobj_umc_creator != null)
                            return mobj_umc_creator.GetPercentComplete();
                        return 0;
                    case enmStatus.CLUSTERING:
                        if (mobj_umc_creator != null)
                            return mobj_umc_creator.GetPercentComplete();
                        return 0;
                    case enmStatus.SUMMARIZING:
                        if (mobj_umc_creator != null)
                            return mobj_umc_creator.GetPercentComplete();
                        return 0;
                    case enmStatus.COMPLETE:
                        return 100;
                }
                return 0;
            }
        }

        public enmStatus Status { get; private set; }

        public string Message { get; private set; }

        public string OutputFileName { get; set; }

        public string OptionsFileName { get; set; }

        public string FileName { get; set; }

        public int MinUMCLength { get; set; }

        public int MinScan
        {
            get { return mobj_umc_creator.mint_lc_min_scan; }
        }

        public int MaxScan
        {
            get { return mobj_umc_creator.mint_lc_max_scan; }
        }

        public clsUMCCreator()
        {
            mobj_umc_creator = new UMCCreator();
        }

        /// <summary>
        /// Anuj added this method to load all the necessary details for the LC ms feature finder
        /// Program options, data filters and output directories are all loaded here
        /// The file is read and data filters applied and finally UMCs are found.
        /// </summary>
        public void LoadFindUMCs()
        {
            Status = enmStatus.LOADING;
            LoadProgramOptions();

            if (mbln_process_chunks)
            {
                double mflt_mono_mass_chunk_start = mflt_mono_mass_start;
                float chunk_size = mobj_umc_creator.GetSegmentSize();

                log("Processing with Chunks ...");
                int iChunk = 0;
                int UMC_count = 0;

                while (Status != enmStatus.COMPLETE)
                {
                    //here's where you need to read a file process maxPoints and then continue from where you left off
                    Status = enmStatus.CHUNKING;

                    //here's where we have to figuure out how to laod the isos sqlite database file into the mvect_isotope_peaks
                    //vector using custom reading code. Once we have that information, we can process the rest of the data as is
                    //and write the UMC's out to the final output file. We can create the output file in APPEND mode if necessary

                    Status = enmStatus.LOADING;
                    if (mflt_mono_mass_chunk_start + chunk_size > mflt_mono_mass_end)
                    {
                        Status = enmStatus.COMPLETE;
                        continue;
                    }

                    //set the range of mass in which we'll operate
                    mobj_umc_creator.SetMassRange((float) mflt_mono_mass_chunk_start,
                        (float) mflt_mono_mass_chunk_start + chunk_size);

                    //instead of READCSV, here's where we need to read from the sqLite database and get a list of peaks
                    //also you need to do this read in such a manner that you get close to the maxPoints per segment range
                    //int numPeaks = mobj_umc_creator->LoadPeaksFromDatabase();

                    mflt_mono_mass_chunk_start = mflt_mono_mass_chunk_start + chunk_size;
                    log("Processing one Chunk");
                    int numPeaks = mobj_umc_creator.ReadCSVFile();

                    log("Total number of peaks we'll consider = ", numPeaks);
                    if (numPeaks == 0)
                    {
                        Status = enmStatus.COMPLETE;
                        log("Nothing read in");
                        continue;
                    }

                    //if ( false ){
                    //since mbln_process is true, the file reading will stop when it has loaded maxPoints of peaks that pass filters
                    //so the mono mass of the last peak loaded will be the end mass 
                    Status = enmStatus.CLUSTERING;
                    mobj_umc_creator.CreateUMCsSinglyLinkedWithAll();
                    Status = enmStatus.SUMMARIZING;

                    Message = "Filtering out short clusters";
                    mobj_umc_creator.RemoveShortUMCs(MinUMCLength);
                    Message = "Calculating UMC statistics";
                    mobj_umc_creator.CalculateUMCs();

                    //now here's where we'll have to perform some magic of writing to temporary files and reloading those files in 
                    //conjunction with loading data from the original isos files. The isos file data can be directly loaded using ReadCSV
                    //but it'll have to stack the peaks on top of existing UMC peaks. Maybe we could just serialize all objects
                    //to a temporary file and reload from there? Something to consider
                    // mobj_umc_creator->SerializeObjects();
                    // objectsSerialized = true;

                    PrintUMCsToFile(iChunk, UMC_count);
                    UMC_count += mobj_umc_creator.GetNumUmcs();
                    //}
                    iChunk++;
                }

                // Combine several partial UMC files into one file.
            }
            else
            {
                log("Processing without Chunks...");
                Status = enmStatus.LOADING;
                int numPeaks = mobj_umc_creator.ReadCSVFile();
                log("Total number of peaks we'll consider = ", numPeaks);

                Status = enmStatus.CLUSTERING;
                log("Creating UMCs...");
                mobj_umc_creator.CreateUMCsSinglyLinkedWithAll();

                Status = enmStatus.SUMMARIZING;
                log("Filtering out short UMCs...");
                mobj_umc_creator.RemoveShortUMCs(MinUMCLength);

                log("Calculating UMC statistics...");
                mobj_umc_creator.CalculateUMCs();
                Status = enmStatus.COMPLETE;

                log("Total number of UMCs = ", mobj_umc_creator.GetNumUmcs());

                log("Writing output files...");
                PrintUMCsToFile();
            }

            mfile_logFile.Close();
        }

        public void FindUMCs()
        {
            Status = enmStatus.CLUSTERING;
            Message = "Clustering Isotope Peaks";
            mobj_umc_creator.CreateUMCsSinglyLinkedWithAll();
            Status = enmStatus.SUMMARIZING;
            Message = "Filtering out short clusters";
            mobj_umc_creator.RemoveShortUMCs(MinUMCLength);
            Message = "Calculating UMC statistics";
            mobj_umc_creator.CalculateUMCs();
            Status = enmStatus.COMPLETE;
        }

        private void LoadFindUMCs(bool is_pek_file)
        {
            System.Console.WriteLine("Loading UMCs");
            Status = enmStatus.LOADING;

            if (is_pek_file)
            {
                Message = "Loading PEK file";
                mobj_umc_creator.ReadPekFileMemoryMapped(FileName);
            }
            else
            {
                Message = "Loading CSV file";
                System.Console.WriteLine(Message);
                mobj_umc_creator.ReadCSVFile(FileName);
            }

            Status = enmStatus.CLUSTERING;
            Message = "Clustering Isotope Peaks";
            mobj_umc_creator.CreateUMCsSinglyLinkedWithAll();
            Status = enmStatus.SUMMARIZING;
            Message = "Filtering out short clusters";
            mobj_umc_creator.RemoveShortUMCs(MinUMCLength);
            Message = "Calculating UMC statistics";
            mobj_umc_creator.CalculateUMCs();
        }

        public void LoadFindUMCsPEK()
        {
            LoadFindUMCs(true);
            Status = enmStatus.COMPLETE;
        }

        public void LoadFindUMCsCSV()
        {
            LoadFindUMCs(false);
            Status = enmStatus.COMPLETE;
        }

        public void ResetStatus()
        {
            Status = enmStatus.IDLE;
            mobj_umc_creator.Reset();
        }

        public void SetIsotopePeaks(ref clsIsotopePeak[] isotope_peaks)
        {
            int numPeaks = isotope_peaks.Length;
            List<IsotopePeak> vectPeaks = new List<IsotopePeak>(numPeaks);

            for (int pkNum = 0; pkNum < numPeaks; pkNum++)
            {
                clsIsotopePeak isoPk;
                IsotopePeak pk = new IsotopePeak();
                isoPk = isotope_peaks[pkNum];

                pk.mint_original_index = isoPk.mint_original_index;
                pk.mint_umc_index = isoPk.mint_umc_index;
                pk.mint_lc_scan = isoPk.mint_lc_scan;
                pk.mshort_charge = isoPk.mshort_charge;
                pk.mdbl_abundance = isoPk.mdbl_abundance;
                pk.mdbl_mz = isoPk.mdbl_mz;
                pk.mflt_fit = isoPk.mflt_fit;
                pk.mflt_ims_drift_time = isoPk.mflt_ims_drift_time;
                pk.mdbl_average_mass = isoPk.mdbl_average_mass;
                pk.mdbl_mono_mass = isoPk.mdbl_mono_mass;
                pk.mdbl_max_abundance_mass = isoPk.mdbl_max_abundance_mass;
                pk.mdbl_i2_abundance = isoPk.mdbl_i2_abundance;

                vectPeaks.Add(pk);
            }
            mobj_umc_creator.SetPeks(ref vectPeaks);
        }

        public bool LoadProgramOptions()
        {
            string settings_file;
            string logText;
            bool success = false;

            settings_file = OptionsFileName;
            IniReader iniReader = new IniReader(settings_file);

            //first load incoming and outgoing filenames and folder options
            string isos_file = iniReader.ReadString("Files", "InputFileName", "");
            string output_dir = iniReader.ReadString("Files", "OutputDirectory", ".");
            mobj_umc_creator.InputFileName = isos_file;
            mobj_umc_creator.OutputDirectory = output_dir;

            mstr_baseFileName = CreateBaseFileName(output_dir, isos_file);

            createLogFile();

            logText = "Loading settings from INI file: " + settings_file;
            log(logText);

            //next load data filters
            float isotopicFit = iniReader.ReadFloat("DataFilters", "MaxIsotopicFit", 1);
            if (isotopicFit == 0)
            {
                isotopicFit = 1;
            }

            int intensityFilter = iniReader.ReadInteger("DataFilters", "MinimumIntensity", 500);
            mflt_mono_mass_start = iniReader.ReadFloat("DataFilters", "MonoMassStart", 0);
            mflt_mono_mass_end = iniReader.ReadFloat("DataFilters", "MonoMassEnd", float.MaxValue);
            mbln_process_chunks = iniReader.ReadBoolean("DataFilters", "ProcessDataInChunks", false);

            //mint_mono_mass_overlap = iniReader.ReadInteger("DataFilters", "MonoMassSegmentOverlapDa", 2);
            if (mflt_mono_mass_end == 0)
            {
                if (mbln_process_chunks)
                {
                    mflt_mono_mass_end = mflt_mono_mass_start + 250;
                }
            }

            int maxPoints = iniReader.ReadInteger("DataFilters", "MaxDataPointsPerChunk", int.MaxValue);
            if (maxPoints == 0)
            {
                maxPoints = int.MaxValue;
            }

            int chunkSize = iniReader.ReadInteger("DataFilters", "ChunkSize", 3000);
            if (chunkSize == 0)
            {
                chunkSize = 3000;
            }

            int imsMinScan = iniReader.ReadInteger("DataFilters", "IMSMinScan", 0);

            int imsMaxScan = iniReader.ReadInteger("DataFilters", "IMSMaxScan", 50000);
            if (imsMaxScan == 0)
            {
                imsMaxScan = int.MaxValue;
            }

            int lcMinScan = iniReader.ReadInteger("DataFilters", "LCMinScan", 0);
            int lcMaxScan = iniReader.ReadInteger("DataFilters", "LCMaxScan", 50000);
            if (lcMaxScan == 0)
            {
                lcMaxScan = int.MaxValue;
            }

            mobj_umc_creator.SetFilterOptions(isotopicFit, intensityFilter, lcMinScan, lcMaxScan, imsMinScan, imsMaxScan,
                mflt_mono_mass_start, mflt_mono_mass_end, mbln_process_chunks, maxPoints, mint_mono_mass_overlap,
                chunkSize);

            //next load the UMC creation options
            float monoMassWeight = iniReader.ReadFloat("UMCCreationOptions", "MonoMassWeight", 0.01f);
            float monoMassConstraint = iniReader.ReadFloat("UMCCreationOptions", "MonoMassConstraint", 50);
            bool monoMassPPM = iniReader.ReadBoolean("UMCCreationOptions", "MonoMassConstraintIsPPM", true);
            float imsDriftWeight = iniReader.ReadFloat("UMCCreationOptions", "IMSDriftTimeWeight", 0.1f);
            float logAbundanceWeight = iniReader.ReadFloat("UMCCreationOptions", "LogAbundanceWeight", 0.1f);
            float netWeight = iniReader.ReadFloat("UMCCreationOptions", "NETWeight", 0.01f);
            float fitWeight = iniReader.ReadFloat("UMCCreationOptions", "FitWeight", 0.01f);
            float avgMassWeight = iniReader.ReadFloat("UMCCreationOptions", "AvgMassWeight", 0.01f);
            float avgMassConstr = iniReader.ReadFloat("UMCCreationOptions", "AvgMassConstraint", 10);
            bool avgMassPPM = iniReader.ReadBoolean("UMCCreationOptions", "AvgMassConstraintIsPPM", true);
            float scanWeight = iniReader.ReadFloat("UMCCreationOptions", "ScanWeight", 0);
            float maxDist = iniReader.ReadFloat("UMCCreationOptions", "MaxDistance", 0.1f);
            bool useGeneric = iniReader.ReadBoolean("UMCCreationOptions", "UseGenericNET", true);
            MinUMCLength = iniReader.ReadInteger("UMCCreationOptions", "MinFeatureLengthPoints", 2);
            bool useCharge = iniReader.ReadBoolean("UMCCreationOptions", "UseCharge", false);

            //this one is not sent over for now
            bool useWeightedEuclidean = iniReader.ReadBoolean("UMCCreationOptions", "UseWeightedEuclidean", false);

            log("Data Filters - ");
            log(" Minimum LC scan = ", lcMinScan);
            log(" Maximum LC scan = ", lcMaxScan);
            log(" Minimum IMS scan = ", imsMinScan);
            log(" Maximum IMS scan = ", imsMaxScan);
            log(" Maximum fit = ", isotopicFit);
            log(" Minimum intensity = ", intensityFilter);
            log(" Mono mass start = ", mflt_mono_mass_start);
            log(" Mono mass end = ", mflt_mono_mass_end);
            log(" Require matching charge state = " + useCharge.ToString());

            //load all the umc creation options
            mobj_umc_creator.SetOptionsEx(monoMassWeight, monoMassConstraint, monoMassPPM, avgMassWeight, avgMassConstr,
                avgMassPPM, logAbundanceWeight, scanWeight, netWeight, fitWeight, maxDist, useGeneric, imsDriftWeight,
                useCharge);

            return success;
        }

        public int GetUmcMapping(ref int[] isotope_peaks_index, ref int[] umc_index)
        {
            int numMappings = mobj_umc_creator.mmultimap_umc_2_peak_index.Sum(i => i.Value.Count);
            isotope_peaks_index = new int[numMappings];
            umc_index = new int[numMappings];

            int mappingNum = 0;
            foreach (var item in mobj_umc_creator.mmultimap_umc_2_peak_index)
            {
                int currentUmcNum = item.Key;
                foreach (var p in item.Value)
                {
                    int pkIndex = p;
                    isotope_peaks_index[mappingNum] = pkIndex;
                    umc_index[mappingNum] = currentUmcNum;
                    mappingNum++;
                }
            }
            return numMappings;
        }

        public void SetLCMinMaxScans(int min, int max)
        {
            mobj_umc_creator.SetLCMinMaxScan(min, max);

            //if (max <= min)
            //cout << " Max scan must be greater than min scan" << endl;
            //throw new exception("Max scan must be greater than min scan");
        }

        public void SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan,
            float wt_fit, float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net,
            float wt_ims_drift_time, bool use_cs)
        {
            mobj_umc_creator.SetOptions(wt_mono_mass, wt_avg_mass, wt_log_abundance, wt_scan, wt_fit, wt_net,
                mono_constraint, avg_constraint, max_dist, use_net, wt_ims_drift_time, use_cs);
        }

        /// <summary>
        /// Calls methods for creating UMCs and writing them to files
        /// </summary>
        /// <returns></returns>
        public bool PrintUMCsToFile()
        {
            return mobj_umc_creator.CreateFeatureFiles(mstr_baseFileName);
        }

        /// <summary>
        /// This overloaded version of the method is used to support chunking
        ///		- Chunking was only thought to be necessary because we ran out of memory on 32-bit machines.
        ///		- Running this program on a 64-bit machine with enough memory will allow it to process very large isos files
        /// </summary>
        /// <param name="chunkIndex"></param>
        /// <param name="featureStartIndex"></param>
        /// <returns></returns>
        public bool PrintUMCsToFile(int chunkIndex, int featureStartIndex)
        {
            string baseFileName = mstr_baseFileName + "_chunk" + chunkIndex;
            return mobj_umc_creator.CreateFeatureFiles(baseFileName, featureStartIndex);
        }

        //MaxIsotopicFit=0.15
        //MinimumIntensity=0
        //MonoMassStart=0
        //MonoMassEnd=0
        //ProcessDataInMonoMassSegments=False
        //MaxDataPointsPerMonoMassSegment=1000000
        //MonoMassSegmentOverlapDa=2
        public void SetFilterOptions(float isotopic_fit, int min_intensity, int min_lc_scan, int max_lc_scan,
            int min_ims_scan, int max_ims_scan, float mono_mass_start, float mono_mass_end, bool process_mass_seg,
            int maxDataPoints, int monoMassSegOverlap, float segmentSize)
        {
            mobj_umc_creator.SetFilterOptions(isotopic_fit, min_intensity, min_lc_scan, max_lc_scan, min_ims_scan,
                max_ims_scan, mono_mass_start, mono_mass_end, process_mass_seg, maxDataPoints, monoMassSegOverlap,
                segmentSize);
        }

        public void SetOptionsEx(float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
            float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm, float wt_log_abundance, float wt_scan,
            float wt_net, float wt_fit, double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs)
        {
            mobj_umc_creator.SetOptionsEx(wt_mono_mass, mono_constraint, mono_constraint_is_ppm, wt_avg_mass,
                avg_constraint, avg_constraint_is_ppm, wt_log_abundance, wt_scan, wt_net, wt_fit, max_dist, use_net,
                wt_ims_drift_time, use_cs);
        }

        public clsUMC[] GetUMCs()
        {
            int numUmcs = mobj_umc_creator.GetNumUmcs();
            clsUMC[] arr_umcs = new clsUMC[numUmcs];

            for (int umcNum = 0; umcNum < numUmcs; umcNum++)
            {
                UMC umc = mobj_umc_creator.mvect_umcs[umcNum];
                clsUMC newUmc = new clsUMC();
                newUmc.mdbl_abundance = umc.mdbl_sum_abundance;

                newUmc.mdbl_class_rep_mz = (double) umc.mdbl_class_rep_mz;
                newUmc.mdbl_mono_mass = (double) umc.mdbl_median_mono_mass;
                newUmc.mdbl_mono_mass_calibrated = newUmc.mdbl_mono_mass;

                newUmc.mint_scan = umc.mint_max_abundance_scan;
                newUmc.mint_start_scan = umc.mint_start_scan;
                newUmc.mint_end_scan = umc.mint_stop_scan;

                newUmc.mdbl_net = (double) umc.mint_max_abundance_scan;
                if (mobj_umc_creator.mint_lc_max_scan > mobj_umc_creator.mint_lc_min_scan)
                {
                    // Compute Generic NET value
                    newUmc.mdbl_net = (double) (umc.mint_max_abundance_scan - mobj_umc_creator.mint_lc_min_scan) * 1.0 /
                                      (mobj_umc_creator.mint_lc_max_scan - mobj_umc_creator.mint_lc_min_scan);
                }

                newUmc.mint_scan_aligned = umc.mint_max_abundance_scan;

                newUmc.mint_umc_index = umc.mint_umc_index;
                newUmc.mint_class_rep_charge = (int) umc.mshort_class_rep_charge;
                arr_umcs[umcNum] = newUmc;
            }

            return arr_umcs;
        }

    }
}
