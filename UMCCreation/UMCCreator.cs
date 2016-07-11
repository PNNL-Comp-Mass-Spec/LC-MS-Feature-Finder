
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace UMCCreation
{
    /// <summary>
    /// Backing class for clsUMCCreator
    /// </summary>
    public class UMCCreator
    {
        private string mstr_inputFile;

        //data filters when loading data isotopic_fit, int min_intensity, int mono_mass_start, int mono_mass_end, bool process_mass_seg, int maxDataPoints, int monoMassSegOverlap
        private float mflt_isotopic_fit_filter;
        private int mint_min_intensity;
        private float mflt_mono_mass_start;
        private float mflt_mono_mass_end;
        private bool mbln_process_mass_seg;
        private int mint_max_data_points;
        // Unused: private int mint_mono_mass_seg_overlap;
        private int mint_ims_min_scan_filter;
        private int mint_ims_max_scan_filter;
        private int mint_lc_min_scan_filter;
        private int mint_lc_max_scan_filter;

        // Weights to use when clustering points
        // Set a weight to 0 to effectively disable that weight
        private float mflt_wt_mono_mass;
        private float mflt_wt_average_mass;
        private float mflt_wt_log_abundance;
        private float mflt_wt_scan;
        private float mflt_wt_net;
        private float mflt_wt_fit;
        private float mflt_wt_ims_drift_time;

        private float mflt_constraint_mono_mass;
        private bool mbln_constraint_mono_mass_is_ppm;

        private float mflt_constraint_average_mass;
        private bool mbln_constraint_average_mass_is_ppm;

        private bool mbln_constraint_charge_state;

        private double mdbl_max_distance;
        private short mshort_percent_complete;

        private bool mbln_use_net; // When True, then uses NET and not Scan
        private bool mbln_is_ims_data;
        //private bool mbln_is_weighted_euc;

        private float mflt_segment_size;

        public int mint_lc_min_scan;
        public int mint_lc_max_scan;
        public int mint_ims_min_scan;
        public int mint_ims_max_scan;

        //public std::multimap<int, int> mmultimap_umc_2_peak_index ; 
        public readonly SortedDictionary<int, List<int>> mmultimap_umc_2_peak_index = new SortedDictionary<int, List<int>>();
        public readonly List<IsotopePeak> mvect_isotope_peaks = new List<IsotopePeak>();
        public readonly List<int> mvect_umc_num_members = new List<int>();
        public readonly List<UMC> mvect_umcs = new List<UMC>();

        public UMCCreator()
        {
            mflt_wt_mono_mass = 0.01F; // using ppms 10 ppm = length of 0.1 ppm
            mflt_wt_average_mass = 0.01F; // using ppms 10 ppm = length of 0.1 ppm
            mflt_wt_log_abundance = 0.1F;
            mflt_wt_scan = 0.01F;
            mflt_wt_fit = 0.1F;
            mflt_wt_ims_drift_time = 0.1F;

            mflt_constraint_mono_mass = 10.0F; // is in ppm
            mflt_constraint_average_mass = 10.0F; // is in ppm. 
            mdbl_max_distance = 0.1;

            mbln_use_net = true;
            mbln_constraint_mono_mass_is_ppm = true;
            mbln_constraint_average_mass_is_ppm = true;

            mshort_percent_complete = 0;

            mint_lc_min_scan = int.MaxValue;
            mint_lc_max_scan = 0;
            mint_ims_min_scan = int.MaxValue;
            mint_ims_max_scan = 0;
        }

        public short GetPercentComplete()
        {
            return mshort_percent_complete;
        }

        public double PeakDistance(IsotopePeak a, IsotopePeak b)
        {
            if (mbln_constraint_mono_mass_is_ppm)
            {
                if (a.mdbl_mono_mass > 0 &&
                    (Math.Abs((a.mdbl_mono_mass - b.mdbl_mono_mass) * mflt_wt_mono_mass / a.mdbl_mono_mass * 1000000) >
                     mflt_constraint_mono_mass))
                    return double.MaxValue;
            }
            else
            {
                if (Math.Abs((a.mdbl_mono_mass - b.mdbl_mono_mass)) * mflt_wt_mono_mass > mflt_constraint_mono_mass)
                    return double.MaxValue;
            }

            if (mbln_constraint_average_mass_is_ppm)
            {
                if (a.mdbl_average_mass > 0 &&
                    (Math.Abs((a.mdbl_average_mass - b.mdbl_average_mass) * mflt_wt_average_mass / a.mdbl_average_mass *
                              1000000) > mflt_constraint_average_mass))
                    return double.MaxValue;
            }
            else
            {
                if (Math.Abs((a.mdbl_average_mass - b.mdbl_average_mass)) * mflt_wt_average_mass >
                    mflt_constraint_average_mass)
                    return double.MaxValue;
            }

            var a_log_abundance = Math.Log10(a.mdbl_abundance);
            var b_log_abundance = Math.Log10(b.mdbl_abundance);

            double sqrDist = 0;

            sqrDist += (a.mdbl_mono_mass - b.mdbl_mono_mass) * (a.mdbl_mono_mass - b.mdbl_mono_mass) * mflt_wt_mono_mass *
                       mflt_wt_mono_mass;
            sqrDist += (a.mdbl_average_mass - b.mdbl_average_mass) * (a.mdbl_average_mass - b.mdbl_average_mass) *
                       mflt_wt_average_mass * mflt_wt_average_mass;
            sqrDist += (a_log_abundance - b_log_abundance) * (a_log_abundance - b_log_abundance) * mflt_wt_log_abundance *
                       mflt_wt_log_abundance;

            if (mbln_use_net)
            {
                // Convert scan difference to Generic NET
                var net_distance = (a.mint_lc_scan - b.mint_lc_scan) * 1.0 / (mint_lc_max_scan - mint_lc_min_scan);
                sqrDist += net_distance * net_distance * mflt_wt_net * mflt_wt_net;
            }
            else
            {
                sqrDist += (a.mint_lc_scan - b.mint_lc_scan) * (a.mint_lc_scan - b.mint_lc_scan) * mflt_wt_scan *
                           mflt_wt_scan;
            }

            sqrDist += (a.mflt_fit - b.mflt_fit) * (a.mflt_fit - b.mflt_fit) * mflt_wt_fit * mflt_wt_fit;

            // IMS Drift time
            sqrDist += (a.mflt_ims_drift_time - b.mflt_ims_drift_time) * (a.mflt_ims_drift_time - b.mflt_ims_drift_time) *
                       mflt_wt_ims_drift_time * mflt_wt_ims_drift_time;

            return Math.Sqrt(sqrDist);
        }

        public int GetNumUmcs()
        {
            return mvect_umcs.Count;
        }

        public int ReadCSVFile()
        {
            return ReadCSVFile(mstr_inputFile);
        }

        public int ReadCSVFile(string fileName)
        {
            //string startTag; //will be defined based on header

            var numPeaks = 0;
            Reset();
            // Using a stream reader with a large buffer rather than a MemoryMappedFile
            using (
                var stream = new StreamReader(new FileStream(fileName, FileMode.Open, FileAccess.Read, FileShare.Read),
                    Encoding.ASCII, true, 65535))
            {

                var file_len = stream.BaseStream.Length;

                //columns for IMS data 
                //'frame_num,ims_scan_num,charge,abundance,mz,fit,average_mw,monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance,orig_intensity,TIA_orig_intensity, drift_time,cumulative_drift_time\n'

                //columns for LC-MS data
                //scan_num,charge,abundance,mz,fit,average_mw,monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance/

                var origLineNumber = 0;
                mshort_percent_complete = 0;
                mint_lc_min_scan = int.MaxValue;
                mint_lc_max_scan = 0;

                long pos = 0;

                // not using BaseStream.Position because it will always jump by the size of the StreamReader buffer.
                mshort_percent_complete = (short) ((100.0 * pos) / file_len);
                if (mshort_percent_complete > 99)
                    mshort_percent_complete = 99;

                var buffer = stream.ReadLine();

                if (buffer == null)
                {
                    throw new Exception("Could not read file");
                }
                pos += buffer.Length;

                if (buffer.Trim().StartsWith("frame_num"))
                {
                    //startTag = "frame_num,ims_scan_num,charge,abundance,mz,fit,average_mw,monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance,orig_intensity,TIA_orig_intensity,drift_time";
                    mbln_is_ims_data = true;
                }
                else
                {
                    //startTag = "scan_num,charge,abundance,mz,fit,average_mw, monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance";
                    mbln_is_ims_data = false;
                }

                // ReSharper disable once NotAccessedVariable
                double fwhm = 0;

                // ReSharper disable once NotAccessedVariable
                double s2n = 0;

                while (!stream.EndOfStream)
                {
                    buffer = stream.ReadLine();
                    if (buffer == null)
                    {
                        continue;
                    }
                    pos += buffer.Length;

                    var pk = new IsotopePeak
                    {
                        mdbl_abundance = 0,
                        mdbl_i2_abundance = 0,
                        mdbl_average_mass = 0,
                        mflt_fit = 0,
                        mdbl_max_abundance_mass = 0,
                        mdbl_mono_mass = 0,
                        mdbl_mz = 0,
                        mshort_charge = 0,
                        mflt_ims_drift_time = 0
                    };

                    // not using BaseStream.Position because it will always jump by the size of the StreamReader buffer.
                    mshort_percent_complete = (short) ((100.0 * pos) / file_len);

                    if (mshort_percent_complete > 99)
                        mshort_percent_complete = 99;

                    var tokens = buffer.Split(',');
                    if (tokens.Length == 0)
                    {
                        continue;
                    }

                    var index = 0;
                    //in either case the first value is the lc_scan_num
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mint_lc_scan = int.Parse(tokens[index].Trim());
                    }
                    index++;

                    if (mbln_is_ims_data)
                    {
                        //then we need to parse out ims_scan number
                        if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                        {
                            pk.mint_ims_scan = int.Parse(tokens[index].Trim());
                        }
                        index++;
                    }

                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mshort_charge = short.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mdbl_abundance = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mdbl_mz = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mflt_fit = float.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mdbl_average_mass = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mdbl_mono_mass = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mdbl_max_abundance_mass = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        fwhm = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        s2n = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mdbl_mono_abundance = double.Parse(tokens[index].Trim());
                    }
                    index++;
                    if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                    {
                        pk.mdbl_i2_abundance = double.Parse(tokens[index].Trim());
                    }
                    index++;

                    //if it's ims data then we have to read four more columns of data
                    if (mbln_is_ims_data)
                    {
                        //orig intensity, TIA original intensity, drift time and cumulative drift time
                        if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                        {
                            pk.mflt_orig_intensity = float.Parse(tokens[index].Trim());
                        }
                        index++;
                        if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                        {
                            pk.mflt_tia_orig_intensity = float.Parse(tokens[index].Trim());
                        }
                        index++;
                        if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                        {
                            pk.mflt_ims_drift_time = float.Parse(tokens[index].Trim());
                        }
                        index++;
                        if (tokens.Length > index && !string.IsNullOrWhiteSpace(tokens[index]))
                        {
                            pk.mflt_cum_drift_time = float.Parse(tokens[index].Trim());
                        }
                    }

                    // Check here to see of MAP index is correct. Dameng
                    pk.mint_original_index = numPeaks;
                    //when reading from the file, this should be the line number in the original isos file
                    pk.mint_line_number_in_file = origLineNumber;

                    //check if the filter criteria is satisfied before adding a peak
                    if (ConsiderPeak(pk))
                    {
                        // Remove the following condition may cause one segment has more than mint_max_data_points
                        if (false)
                        {
                            if (mbln_process_mass_seg)
                            {
                                if (numPeaks > mint_max_data_points)
                                {
                                    Console.WriteLine("I know this happens");
                                    numPeaks--;
                                    break;
                                }
                            }
                        }

                        //check if the min scans and max scans for both lc and ims need to be fixed
                        if (pk.mint_lc_scan <= mint_lc_min_scan)
                        {
                            mint_lc_min_scan = pk.mint_lc_scan;
                        }

                        if (pk.mint_lc_scan >= mint_lc_max_scan)
                        {
                            mint_lc_max_scan = pk.mint_lc_scan;
                        }

                        if (mbln_is_ims_data)
                        {
                            //need to fix the ims min scans and max scans that were loaded
                            if (pk.mint_ims_scan <= mint_ims_min_scan)
                            {
                                mint_ims_min_scan = pk.mint_ims_scan;
                            }

                            if (pk.mint_ims_scan >= mint_ims_max_scan)
                            {
                                mint_ims_max_scan = pk.mint_ims_scan;
                            }
                        }

                        mvect_isotope_peaks.Add(pk);

                        //increment number of peaks read
                        numPeaks++;
                    }

                    origLineNumber++;
                }
            }

            return numPeaks;
        }

        /// <summary>
        /// Reads to the line in the file specified by <see cref="startString"/>
        /// </summary>
        /// <param name="stream">StreamReader</param>
        /// <param name="buffer">string buffer to store last read line</param>
        /// <param name="startString">start of line being searched for</param>
        /// <param name="pos">tracked position in the file</param>
        /// <returns>true if we found the line, false if not found (end of stream)</returns>
        private bool ReadToLine(StreamReader stream, ref string buffer, string startString, ref long pos)
        {
            while (!stream.EndOfStream && !buffer.Trim().StartsWith(startString))
            {
                buffer = stream.ReadLine() ?? "";
                pos += buffer.Length;
            }

            return !stream.EndOfStream;
        }

        public void ReadPekFileMemoryMapped(string fileName)
        {
            // Using a stream reader with a large buffer rather than a MemoryMappedFile
            Reset();
            using (
                var stream = new StreamReader(new FileStream(fileName, FileMode.Open, FileAccess.Read, FileShare.Read),
                    Encoding.ASCII, true, 65535))
            {
                var file_len = stream.BaseStream.Length;

                var fileNameTag = "Filename:";
                var startTag = "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW";
                var startTagIsotopicLabeled =
                    "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW,   Imono,   I+2";
                var stopTag = "Processing stop time:";
                var isotopically_labeled = false;
                var is_first_scan = true;
                var is_pek_file_from_wiff = false;

                var numPeaks = 0;
                mshort_percent_complete = 0;
                var mint_min_scan = int.MaxValue;
                var mint_max_scan = 0;

                long pos = 0;

                while (!stream.EndOfStream)
                {
                    // not using BaseStream.Position because it will always jump by the size of the StreamReader buffer.
                    mshort_percent_complete = (short) ((100.0 * pos) / file_len);
                    if (mshort_percent_complete > 99)
                        mshort_percent_complete = 99;

                    // Give an empty string if ReadLine returns null
                    var buffer = stream.ReadLine() ?? "";
                    pos += buffer.Length;

                    if (!ReadToLine(stream, ref buffer, fileNameTag, ref pos))
                    {
                        break;
                    }

                    var pk = new IsotopePeak
                    {
                        mdbl_abundance = 0,
                        mdbl_i2_abundance = 0,
                        mdbl_average_mass = 0,
                        mflt_fit = 0,
                        mdbl_max_abundance_mass = 0,
                        mdbl_mono_mass = 0,
                        mdbl_mz = 0,
                        mshort_charge = 0,
                        mflt_ims_drift_time = 0
                    };

                    // found file name. start at the end.
                    var index = buffer.LastIndexOf('.') + 1;

                    // wiff file pek files have a weird format. Another ICR2LS-ism.
                    if (is_first_scan)
                    {
                        if (buffer.Substring(index).ToLower().StartsWith("wiff"))
                        {
                            is_pek_file_from_wiff = true;
                            // REMEMBER TO DELETE WHEN PARAMETERS ARE SET ELSEWHERE
                            mflt_wt_mono_mass = 0.0025F; // using ppms 10 ppm = length of 0.1 ppm
                            mflt_wt_average_mass = 0.0025F; // using ppms 10 ppm = length of 0.1 ppm
                            mflt_wt_log_abundance = 0.1F;
                            mflt_wt_scan = 0.01F;
                            mflt_wt_fit = 0.1F;
                            mflt_wt_ims_drift_time = 0.1F;

                            mflt_constraint_mono_mass = 25.0F; // is in ppm
                            mflt_constraint_average_mass = 25.0F; // is in ppm. 

                            mbln_use_net = true;
                            mbln_constraint_mono_mass_is_ppm = true;
                            mbln_constraint_average_mass_is_ppm = true;

                            mdbl_max_distance = 0.1;

                            mshort_percent_complete = 0;
                        }
                    }

                    if (is_pek_file_from_wiff)
                        index += 5;
                    pk.mint_lc_scan = int.Parse(buffer.Substring(index).Trim());
                    if (pk.mint_lc_scan > mint_max_scan)
                        mint_max_scan = pk.mint_lc_scan;
                    if (pk.mint_lc_scan < mint_min_scan)
                        mint_min_scan = pk.mint_lc_scan;

                    if (!ReadToLine(stream, ref buffer, startTag, ref pos))
                    {
                        break;
                    }
                    if (is_first_scan)
                    {
                        if (buffer.StartsWith(startTagIsotopicLabeled)) // TODO: Check this for working...
                        {
                            isotopically_labeled = true;
                        }
                        is_first_scan = false;
                    }

                    while (!stream.EndOfStream)
                    {
                        buffer = stream.ReadLine() ?? "";
                        pos += buffer.Length;
                        if (string.IsNullOrWhiteSpace(buffer))
                        {
                            continue;
                        }

                        var tokens = buffer.Split('\t');
                        if (buffer.Trim().StartsWith(stopTag) || !char.IsDigit(buffer.Trim()[0]) || tokens.Length == 0)
                        {
                            break;
                        }

                        pk.mshort_charge = short.Parse(tokens[0].Trim());
                        if (tokens.Length > 1 && !string.IsNullOrWhiteSpace(tokens[1]))
                        {
                            pk.mdbl_abundance = double.Parse(tokens[1].Trim());
                        }
                        if (tokens.Length > 2 && !string.IsNullOrWhiteSpace(tokens[2]))
                        {
                            pk.mdbl_mz = double.Parse(tokens[2].Trim());
                        }
                        if (tokens.Length > 3 && !string.IsNullOrWhiteSpace(tokens[3]))
                        {
                            pk.mflt_fit = float.Parse(tokens[3].Trim());
                        }
                        if (tokens.Length > 4 && !string.IsNullOrWhiteSpace(tokens[4]))
                        {
                            pk.mdbl_average_mass = double.Parse(tokens[4].Trim());
                        }
                        if (tokens.Length > 5 && !string.IsNullOrWhiteSpace(tokens[5]))
                        {
                            pk.mdbl_mono_mass = double.Parse(tokens[5].Trim());
                        }
                        if (tokens.Length > 6 && !string.IsNullOrWhiteSpace(tokens[6]))
                        {
                            pk.mdbl_max_abundance_mass = double.Parse(tokens[6].Trim());
                        }
                        if (isotopically_labeled)
                        {
                            if (tokens.Length > 7 && !string.IsNullOrWhiteSpace(tokens[7]))
                            {
                                pk.mdbl_mono_abundance = double.Parse(tokens[7].Trim());
                            }
                            if (tokens.Length > 8 && !string.IsNullOrWhiteSpace(tokens[8]))
                            {
                                pk.mdbl_i2_abundance = double.Parse(tokens[8].Trim());
                            }
                        }

                        pk.mint_original_index = numPeaks;
                        mvect_isotope_peaks.Add(pk);
                        numPeaks++;
                    }
                }
            }
        }

        /*/// <summary>
        /// NOT USED.....................
        /// </summary>
        /// <param name="fileName"></param>
        public void ReadPekFile(string fileName)
        {
            Reset();
            using (var fp = new StreamReader(new FileStream(fileName, FileMode.Open, FileAccess.Read, FileShare.Read)))
            {
                long file_len = fp.BaseStream.Length;

                string buffer;
                string fileNameTag = "Filename:";
                string startTag = "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW";
                string startTagIsotopicLabeled =
                    "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW,   Imono,   I+2";
                string stopTag = "Processing stop time:";
                bool reading = false;
                bool isotopically_labeled = false;
                bool is_first_scan = true;

                int numPeaks = 0;
                mshort_percent_complete = 0;
                mint_lc_min_scan = int.MaxValue;
                mint_lc_max_scan = 0;

                int pos = 0;
                while (!fp.EndOfStream)
                {
                    IsotopePeak pk = new IsotopePeak();

                    mshort_percent_complete = (short) ((100.0 * pos) / file_len);
                    buffer = fp.ReadLine();
                    if (string.IsNullOrWhiteSpace(buffer))
                    {
                        continue;
                    }
                    pos += buffer.Length;
                    if (!reading)
                    {
                        if (buffer.StartsWith(startTag))
                        {
                            reading = true;
                            if (is_first_scan)
                            {
                                if (buffer.StartsWith(startTagIsotopicLabeled))
                                {
                                    isotopically_labeled = true;
                                }
                                is_first_scan = false;
                            }
                        }
                        else if (buffer.StartsWith(fileNameTag))
                        {
                            pk = new IsotopePeak();
                            pk.mdbl_abundance = 0;
                            pk.mdbl_i2_abundance = 0;
                            pk.mdbl_average_mass = 0;
                            pk.mflt_fit = 0;
                            pk.mdbl_max_abundance_mass = 0;
                            pk.mdbl_mono_mass = 0;
                            pk.mdbl_mz = 0;
                            pk.mshort_charge = 0;
                            pk.mflt_ims_drift_time = 0;

                            // found file name. start at the end. 
                            int index = buffer.LastIndexOf('.') + 1;
                            pk.mint_lc_scan = int.Parse(buffer.Substring(index).Trim());
                            if (pk.mint_lc_scan > mint_lc_max_scan)
                                mint_lc_max_scan = pk.mint_lc_scan;
                            if (pk.mint_lc_scan < mint_lc_min_scan)
                                mint_lc_min_scan = pk.mint_lc_scan;
                        }
                    }
                    else
                    {
                        //int numRead = fscanf(fp, "%hd\t%lg\t%g\t%g\t%g\t%g\t%g", &charge, &abundance, &mz, &fit, 
                        //    &averageMass, &monoMass, &maxMass) ; 
                        var tokens = buffer.Split('\t');
                        if (buffer.Trim().StartsWith(stopTag) || !char.IsDigit(buffer.Trim()[0]) || tokens.Length == 0)
                        {
                            // at the end of the reading sections
                            reading = false;
                        }

                        pk.mshort_charge = short.Parse(tokens[0].Trim());
                        if (tokens.Length > 1 && !string.IsNullOrWhiteSpace(tokens[1]))
                        {
                            pk.mdbl_abundance = double.Parse(tokens[1].Trim());
                        }
                        if (tokens.Length > 2 && !string.IsNullOrWhiteSpace(tokens[2]))
                        {
                            pk.mdbl_mz = double.Parse(tokens[2].Trim());
                        }
                        if (tokens.Length > 3 && !string.IsNullOrWhiteSpace(tokens[3]))
                        {
                            pk.mflt_fit = float.Parse(tokens[3].Trim());
                        }
                        if (tokens.Length > 4 && !string.IsNullOrWhiteSpace(tokens[4]))
                        {
                            pk.mdbl_average_mass = double.Parse(tokens[4].Trim());
                        }
                        if (tokens.Length > 5 && !string.IsNullOrWhiteSpace(tokens[5]))
                        {
                            pk.mdbl_mono_mass = double.Parse(tokens[5].Trim());
                        }
                        if (tokens.Length > 6 && !string.IsNullOrWhiteSpace(tokens[6]))
                        {
                            pk.mdbl_max_abundance_mass = double.Parse(tokens[6].Trim());
                        }
                        if (isotopically_labeled)
                        {
                            if (tokens.Length > 7 && !string.IsNullOrWhiteSpace(tokens[7]))
                            {
                                pk.mdbl_mono_abundance = double.Parse(tokens[7].Trim());
                            }
                            if (tokens.Length > 8 && !string.IsNullOrWhiteSpace(tokens[8]))
                            {
                                pk.mdbl_i2_abundance = double.Parse(tokens[8].Trim());
                            }
                        }

                        pk.mint_original_index = numPeaks;
                        mvect_isotope_peaks.Add(pk);
                        numPeaks++;
                    }
                }
            }
        }*/

        private void AddUMCToMap(int key, int value)
        {
            if (!mmultimap_umc_2_peak_index.ContainsKey(key))
            {
                mmultimap_umc_2_peak_index.Add(key, new List<int>());
            }
            mmultimap_umc_2_peak_index[key].Add(value);
        }

        public void CreateUMCsSinglyLinkedWithAll()
        {
            var chargeStateMatch = true;
            mshort_percent_complete = 0;
            mmultimap_umc_2_peak_index.Clear();
            var numPeaks = mvect_isotope_peaks.Count;
            for (var pkNum = 0; pkNum < numPeaks; pkNum++)
            {
                mvect_isotope_peaks[pkNum].mint_umc_index = -1;
            }

            var vectTempPeaks = new List<IsotopePeak>();
            //vectTempPeaks.AddRange(mvect_isotope_peaks.Select(p => p.Clone()).OrderBy(x => x.mdbl_mono_mass)); // Performs a stable sort (item ordering maintained for equivalent objects)
            vectTempPeaks.AddRange(mvect_isotope_peaks.Select(p => p.Clone()));

            // basically take all umcs sorted in mass and perform single linkage clustering.
            //vectTempPeaks.Sort((x, y) => x.mdbl_mono_mass.CompareTo(y.mdbl_mono_mass)); // Unstable sort - item ordering not guaranteed
            vectTempPeaks.Sort(); // use the item comparator (currently causes a stable sort); probably faster in release mode than OrderBy.

            // now we are sorted. Start with the first index and move rightwards.
            // For each index, 
            // 1. If it already belongs to a UMC, move right comparing to all umcs in mass tolerance. For each match
            //        a. If match is not in a UMC, calculate distance. If within tolerance add it to current umc.
            //        b. If match is in the same UMC skip.
            //        c. If match is in a different UMC, calculate distance. If within tolerance, merge umcs, update current guys UMC number.
            // 2. If it does not belong to a UMC, move right comparing to all umcs in mass tolerance. For each match, calculate distance.
            //        a. Create new UMC and follow 1.

            var currentIndex = 0;

            var numUmcsSoFar = 0;
            var tempIndices = new List<int>(128);
                // used to store indices of isotope peaks that are moved from one umc to another.

            mshort_percent_complete = 0;
            while (currentIndex < numPeaks)
            {
                mshort_percent_complete = (short) ((100.0 * currentIndex) / numPeaks);
                var currentPeak = vectTempPeaks[currentIndex];
                if (currentPeak.mint_umc_index == -1)
                {
                    // create UMC
                    AddUMCToMap(numUmcsSoFar, currentIndex);
                    currentPeak.mint_umc_index = numUmcsSoFar;
                    vectTempPeaks[currentIndex].mint_umc_index = numUmcsSoFar;
                    numUmcsSoFar++;
                }
                var matchIndex = currentIndex + 1;
                if (matchIndex == numPeaks)
                    break;

                double massTolerance = mflt_constraint_mono_mass;
                if (mbln_constraint_mono_mass_is_ppm)
                    massTolerance *= currentPeak.mdbl_mono_mass / 1000000.0; // Convert from ppm to Da tolerance

                var maxMass = currentPeak.mdbl_mono_mass + massTolerance;
                var matchPeak = vectTempPeaks[matchIndex];
                while (matchPeak.mdbl_mono_mass < maxMass)
                {
                    if (matchPeak.mint_umc_index != currentPeak.mint_umc_index)
                    {
                        var currentDistance = PeakDistance(currentPeak, matchPeak);
                        if (mbln_constraint_charge_state)
                        {
                            chargeStateMatch = (currentPeak.mshort_charge == matchPeak.mshort_charge);
                        }
                        if (currentDistance < mdbl_max_distance && chargeStateMatch)
                        {
                            if (matchPeak.mint_umc_index == -1)
                            {
                                AddUMCToMap(currentPeak.mint_umc_index, matchIndex);
                                vectTempPeaks[matchIndex].mint_umc_index = currentPeak.mint_umc_index;
                            }
                            else
                            {
                                tempIndices.Clear();
                                var numPeaksMerged = 0;
                                // merging time. Merge this guy's umc into the next guys UMC.
                                var listAtIndex = mmultimap_umc_2_peak_index[currentPeak.mint_umc_index];
                                for (var i = 0; i < listAtIndex.Count;)
                                {
                                    var deletePeakIndex = listAtIndex[i];
                                    tempIndices.Add(deletePeakIndex);
                                    vectTempPeaks[deletePeakIndex].mint_umc_index = matchPeak.mint_umc_index;
                                    listAtIndex.RemoveAt(i); // remove item at i and don't increment i
                                    numPeaksMerged++;
                                }
                                for (var mergedPeakNum = 0; mergedPeakNum < numPeaksMerged; mergedPeakNum++)
                                {
                                    AddUMCToMap(matchPeak.mint_umc_index, tempIndices[mergedPeakNum]);
                                }
                                currentPeak.mint_umc_index = matchPeak.mint_umc_index;
                            }
                        }
                    }
                    matchIndex++;
                    if (matchIndex < numPeaks)
                    {
                        matchPeak = vectTempPeaks[matchIndex];
                    }
                    else
                        break;
                }
                currentIndex++;
            }

            // At the end of all of this. The mapping from mmultimap_umc_2_peak_index is from umc_index to index in sorted stuff. 
            // Also, several of the umc indices are no longer valid. So lets step through the map, get new umc indices, renumber them,
            // and set the umc indices in the original vectors.
            numUmcsSoFar = 0;
            foreach (var item in mmultimap_umc_2_peak_index)
            {
                var numMembers = 0;
                foreach (var p in item.Value)
                {
                    var pk = vectTempPeaks[p];
                    mvect_isotope_peaks[pk.mint_original_index].mint_umc_index = numUmcsSoFar;
                    numMembers++;
                }
                mvect_umc_num_members.Add(numMembers);
                numUmcsSoFar++;
            }
            // now set the map object. 
            mmultimap_umc_2_peak_index.Clear();
            for (var pkNum = 0; pkNum < numPeaks; pkNum++)
            {
                var pk = mvect_isotope_peaks[pkNum];
                AddUMCToMap(pk.mint_umc_index, pkNum);
            }
            // DONE!! 
        }

        // Will map be affected by chunking?
        public bool PrintMapping(StreamWriter stream, int featureStartIndex)
        {

            stream.WriteLine("Feature_Index\tPeak_Index");

            foreach (var item in mmultimap_umc_2_peak_index)
            {
                var currentUmcNum = item.Key;
                foreach (var p in item.Value)
                {
                    var pk = mvect_isotope_peaks[p];
                    stream.WriteLine("{0}\t{1}", currentUmcNum + featureStartIndex, pk.mint_line_number_in_file);
                }

            }

            return true;
        }

        //method can be called with either stdout or an output file to write to 
        public bool PrintUMCs(StreamWriter stream, bool print_members, int featureStartIndex)
        {
            var success = true;
            stream.Write("Feature_Index\tMonoisotopic_Mass\tAverage_Mono_Mass\tUMC_MW_Min\tUMC_MW_Max\tScan_Start" +
                         "\tScan_End\tScan\tUMC_Member_Count\tMax_Abundance\tAbundance\tClass_Rep_MZ\tClass_Rep_Charge");
            if (print_members)
            {
                stream.Write("\tData");
            }
            stream.WriteLine();

            var numPrinted = 1;
            foreach (var item in mmultimap_umc_2_peak_index)
            {
                var currentUmcNum = item.Key;
                var current_umc = mvect_umcs[currentUmcNum];
                stream.Write(
                    "{0}\t{1,4:F4}\t{2,4:F4}\t{3,4:F4}\t{4,4:F4}\t{5}\t{6}\t{7}\t{8}\t{9,4:F4}\t{10,4:F4}\t{11,4:F4}\t{12}\t",
                    current_umc.mint_umc_index + featureStartIndex, current_umc.mdbl_median_mono_mass,
                    current_umc.mdbl_average_mono_mass, current_umc.mdbl_min_mono_mass, current_umc.mdbl_max_mono_mass,
                    current_umc.mint_start_scan, current_umc.mint_stop_scan, current_umc.mint_max_abundance_scan,
                    current_umc.min_num_members, current_umc.mdbl_max_abundance, current_umc.mdbl_sum_abundance,
                    current_umc.mdbl_class_rep_mz, current_umc.mshort_class_rep_charge);

                foreach (var p in item.Value)
                {
                    if (print_members)
                    {

                        var pk = mvect_isotope_peaks[p];
                        stream.Write("{0,4:F4}\t{1}\t{2,4:F4}", pk.mdbl_mono_mass, pk.mint_lc_scan, pk.mdbl_abundance);
                    }
                }

                numPrinted++;
                stream.WriteLine();
                stream.Flush();
            }
            stream.Close();

            if (numPrinted < 1)
            {
                success = false;
            }

            return success;
        }

        public void PrintUMCs(bool print_members)
        {
            Console.Write(
                "UMCIndex\tUMCMonoMW\tAverageMonoMass\tUMCMWMin\tUMCMWMax\tScanStart\tScanEnd\tScanClassRep\tUMCMemberCount\tMaxAbundance\tUMCAbundance");
            if (print_members)
            {
                Console.Write("\tData");
            }
            Console.WriteLine();

            foreach (var item in mmultimap_umc_2_peak_index)
            {
                var currentUmcNum = item.Key;
                var current_umc = mvect_umcs[currentUmcNum];
                Console.Write("{0}\t{1:F4}\t{2:F4}\t{3:F4}\t{4:F4}\t{5:F4}\t{6}\t{7}\t{8}\t{9:F0}\t{10:F0}",
                    current_umc.mint_umc_index, current_umc.mdbl_median_mono_mass, current_umc.mdbl_average_mono_mass,
                    current_umc.mdbl_min_mono_mass, current_umc.mdbl_max_mono_mass, current_umc.mint_start_scan,
                    current_umc.mint_stop_scan, current_umc.mint_max_abundance_scan, current_umc.min_num_members,
                    current_umc.mdbl_max_abundance, current_umc.mdbl_sum_abundance);
                foreach (var p in item.Value)
                {
                    if (print_members)
                    {
                        var pk = mvect_isotope_peaks[p];
                        Console.Write("\t{0:F4}\t{1}\t{2:F0}", pk.mdbl_mono_mass, pk.mint_lc_scan,
                            pk.mdbl_abundance);
                    }
                }                
                Console.WriteLine();
            }
        }

        public void PrintPeaks()
        {
            var numPeaks = mvect_isotope_peaks.Count;
            Console.WriteLine("Scan\tCharge\tAbundance\tMZ\tFit\tAverageMass\tMonoMass\tMaxMass\n");

            for (var i = 0; i < numPeaks; i++)
            {
                var pk = mvect_isotope_peaks[i];
                Console.WriteLine("{0}\t{1}\t{2:F0}\t{3:F4}\t{4:F3}\t{5:F4}\t{6:F4}\t{7:F4}", pk.mint_lc_scan,
                    pk.mshort_charge, pk.mdbl_abundance, pk.mdbl_mz, pk.mflt_fit, pk.mdbl_average_mass,
                    pk.mdbl_mono_mass,
                    pk.mdbl_max_abundance_mass);
            }
        }

        /// <summary>
        /// Calls the worker functions for creating and outputting the Features.
        /// </summary>
        /// <param name="baseFileName"></param>
        /// <param name="featureStartIndex"></param>
        /// <returns></returns>
        public bool CreateFeatureFiles(string baseFileName, int featureStartIndex = 0)
        {
            bool success;

            // Create the file where the LCMS Features will be written
            var completeFileName = baseFileName + "_LCMSFeatures.txt";
            using (
                var file =
                    new StreamWriter(new FileStream(completeFileName, FileMode.Create, FileAccess.Write, FileShare.Read))
                )
            {
                success = PrintUMCs(file, false, featureStartIndex);
            }
            if (!success)
            {
                return false;
            }

            // Create the file where the "Features to Peak Map" will be written
            completeFileName = baseFileName + "_LCMSFeatureToPeakMap.txt";
            using (
                var file =
                    new StreamWriter(new FileStream(completeFileName, FileMode.Create, FileAccess.Write, FileShare.Read))
                )
            {
                success = PrintMapping(file, featureStartIndex);
            }

            return success;
        }

        public void Reset()
        {
            mvect_isotope_peaks.Clear();
            mvect_umcs.Clear();
            mvect_umc_num_members.Clear();
            mmultimap_umc_2_peak_index.Clear();
            mshort_percent_complete = 0;
        }

        public void SetUseNet(bool use)
        {
            mbln_use_net = use;
        }

        /*public float GetLastMonoMassLoaded()
        {
            float mass = float.MaxValue;
            try
            {
                IsotopePeak peak = (IsotopePeak) mvect_isotope_peaks[mvect_isotope_peaks.Count];
                mass = (float) peak.mdbl_mono_mass;
            }
            catch (Exception e)
            {
            }
            return mass;
        }*/

        //could all start going into a single branch but have kept it such so that it's a little easier to read.
        //Method to determine whether the filter 
        public bool ConsiderPeak(IsotopePeak peak)
        {
            return (peak.mdbl_abundance >= mint_min_intensity && peak.mflt_fit <= mflt_isotopic_fit_filter) &&
                   //if we're processing a particular mass range then we have to check for the mono mass to be within bounds
                   (mflt_mono_mass_start <= peak.mdbl_mono_mass && peak.mdbl_mono_mass <= mflt_mono_mass_end) &&
                   //check if data is within LC scan range
                   (mint_lc_min_scan_filter <= peak.mint_lc_scan && peak.mint_lc_scan <= mint_lc_max_scan_filter) &&
                   // not ims data or...
                   (!mbln_is_ims_data ||
                    //check if data is within ims scans
                    (mbln_is_ims_data && (mint_ims_min_scan_filter <= peak.mint_ims_scan &&
                                          peak.mint_ims_scan <= mint_ims_max_scan_filter)));
        }

        public void CalculateUMCs()
        {
            mvect_umcs.Clear();
            mvect_umcs.Capacity = mvect_umc_num_members.Count;

            var vect_mass = new List<double>();

            mshort_percent_complete = 0;
            var num_umcs = mvect_umcs.Count;

            foreach (var item in mmultimap_umc_2_peak_index)
            {
                vect_mass.Clear();
                var umc_index = item.Key;
                mshort_percent_complete = (short) ((100.0 * umc_index) / num_umcs);
                var numMembers = mvect_umc_num_members[umc_index];
                var minScan = int.MaxValue;
                var maxScan = int.MinValue;
                var minMass = double.MaxValue;
                var maxMass = -1 * double.MaxValue;
                var maxAbundance = -1 * double.MaxValue;
                double sumAbundance = 0;
                double sumMonoMass = 0;
                var maxAbundanceScan = 0;
                short classRepCharge = 0;
                double classRepMz = 0;

                foreach (var p in item.Value)
                {
                    var pk = mvect_isotope_peaks[p];
                    vect_mass.Add(pk.mdbl_mono_mass);

                    if (pk.mint_lc_scan > maxScan)
                        maxScan = pk.mint_lc_scan;
                    if (pk.mint_lc_scan < minScan)
                        minScan = pk.mint_lc_scan;

                    if (pk.mdbl_mono_mass > maxMass)
                        maxMass = pk.mdbl_mono_mass;
                    if (pk.mdbl_mono_mass < minMass)
                        minMass = pk.mdbl_mono_mass;

                    if (pk.mdbl_abundance > maxAbundance)
                    {
                        maxAbundance = pk.mdbl_abundance;
                        maxAbundanceScan = pk.mint_lc_scan;
                        classRepCharge = pk.mshort_charge;
                        classRepMz = pk.mdbl_mz;
                    }
                    sumAbundance += pk.mdbl_abundance;

                    sumMonoMass += pk.mdbl_mono_mass;
                }

                vect_mass.Sort();

                var new_umc = new UMC
                {
                    mint_umc_index = umc_index,
                    min_num_members = numMembers,
                    mint_start_scan = minScan,
                    mint_stop_scan = maxScan,
                    mint_max_abundance_scan = maxAbundanceScan,
                    mdbl_max_abundance = maxAbundance,
                    mdbl_sum_abundance = sumAbundance,
                    mdbl_min_mono_mass = minMass,
                    mdbl_max_mono_mass = maxMass,
                    mdbl_average_mono_mass = sumMonoMass / numMembers,
                    mdbl_class_rep_mz = classRepMz,
                    mshort_class_rep_charge = classRepCharge
                };




                if (numMembers % 2 == 1)
                {
                    new_umc.mdbl_median_mono_mass = vect_mass[numMembers / 2];
                }
                else
                {
                    new_umc.mdbl_median_mono_mass = 0.5 * (vect_mass[numMembers / 2 - 1] + vect_mass[numMembers / 2]);
                }
                mvect_umcs.Add(new_umc);
            }
        }

        public void RemoveShortUMCs(int min_length)
        {
            // first reset all isotope peak umc indices to -1. 
            var numIsotopePeaks = mvect_isotope_peaks.Count;
            for (var peakNum = 0; peakNum < numIsotopePeaks; peakNum++)
            {
                mvect_isotope_peaks[peakNum].mint_umc_index = -1;
            }

            var numUmcsSoFar = 0;
            mshort_percent_complete = 0;

            var num_umcs = mmultimap_umc_2_peak_index.Sum(i => i.Value.Count);
            foreach (var item in mmultimap_umc_2_peak_index)
            {
                var currentOldUmcNum = item.Key;
                mshort_percent_complete = (short) ((100.0 * currentOldUmcNum) / num_umcs);

                var numMembers = mvect_umc_num_members[currentOldUmcNum];
                foreach (var p in item.Value)
                {
                    if (numMembers >= min_length)
                    {
                        mvect_isotope_peaks[p].mint_umc_index = numUmcsSoFar;
                    }
                }
                if (numMembers >= min_length)
                {
                    mvect_umc_num_members[numUmcsSoFar] = numMembers;
                    numUmcsSoFar++;
                }
            }
            mvect_umc_num_members.RemoveRange(numUmcsSoFar, mvect_umc_num_members.Count - numUmcsSoFar);
            // now set the map object. 
            mmultimap_umc_2_peak_index.Clear();
            for (var pkNum = 0; pkNum < numIsotopePeaks; pkNum++)
            {
                var pk = mvect_isotope_peaks[pkNum];
                if (pk.mint_umc_index != -1)
                {
                    AddUMCToMap(pk.mint_umc_index, pkNum);
                }
            }
            // DONE!! 
        }

        /* // Block of never used code
        public int LoadPeaksFromDatabase()
        {
            //this is to read a sql lite database and retrieve the peaks from it into the mvect_isotope_peaks array.

            return 0;
        }

        //Method not used,
        public void SerializeObjects()
        {
            //First write out all loaded isotopic peaks
            //we should be smarter about this file writing since there's going to be sufficient 
            //number of peaks that don't get used into any UMC or get removed as part of short UMCs
            //maybe each peak needs to have a bit whether to used or not
            ofstream ofs("mvect_isotope_peaks.ros", ios::binary);
    
            //this should write only peaks that are being used to the file
            int count = 0;
            for ( int i = 0; i < mvect_isotope_peaks.size(); i++){
    
                IsotopePeak pk = (IsotopePeak ) mvect_isotope_peaks[i];
                if ( pk.mint_umc_index != -1){
                       ofs.write( (char*)&mvect_isotope_peaks[i], sizeof(IsotopePeak));
                    count++;
                }
            }
    
            std::cout<<"Number of peaks written to file is " << count << " while total = " << mvect_isotope_peaks.size() << std::endl;
    
            ofs.close();
    
            //write out the umc classes
            ofstream ofs1("mvect_umcs", ios::binary);
            for ( int  i=0; i < mvect_umcs.size(); i++){
                ofs1.write( (char*)&mvect_umcs[i], sizeof(UMC));
            }
    
            ofs1.close();
    
    
            //write out the multimap that maps the indices on the isotopic peaks
            //to the indices on the umcs
            fstream fs("mmultimap_umc_2_peak_index", ios::out);
                
            for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
            {
                    
                int umc_index = (*iter).first ; 
                int peak_index = (*iter).second;
                fs << umc_index << "\t" << peak_index <<std::endl;
            }
    
            fs.close();
        }

        public void DeserializeObjects()
        {
        
        }*/


        //Functions added by Anuj Shah
        //public void CreateUMCsSingleLinkedWithAllOnline();
        //public void CreateUMCFromIsotopePeak(IsotopePeak startPeak, UMC &firstUMC);
        //public void AddPeakToUMC (IsotopePeak peak, UMC &umc);
        //public int findCandidateUMCsForPeak(IsotopePeak peak, std::vector<UMC> &umcVector, std::vector<UMC> &candidateUMCs);

        public bool withinMassTolerance(double observedMass, double realMass)
        {
            var massDifferenceInPPM = Math.Abs(realMass - observedMass) * 1000000 / realMass;
            return (massDifferenceInPPM <= mflt_constraint_mono_mass);
        }

        public void SetFilterOptions(float isotopic_fit, int min_intensity, int min_lc_scan, int max_lc_scan,
            int min_ims_scan, int max_ims_scan, float mono_mass_start, float mono_mass_end, bool process_mass_seg,
            int max_data_points, float mono_mass_seg_size)
        {
            mflt_isotopic_fit_filter = isotopic_fit;
            mint_min_intensity = min_intensity;
            mflt_mono_mass_start = mono_mass_start;
            mflt_mono_mass_end = mono_mass_end;
            mbln_process_mass_seg = process_mass_seg;
            mint_max_data_points = max_data_points;
            // Unused: mint_mono_mass_seg_overlap = mono_mass_seg_overlap;
            mint_lc_min_scan_filter = min_lc_scan;
            mint_lc_max_scan_filter = max_lc_scan;
            mint_ims_min_scan_filter = min_ims_scan;
            mint_ims_max_scan_filter = max_ims_scan;
            mflt_segment_size = mono_mass_seg_size;
        }

        public void SetMassRange(float min_mono_mass, float max_mono_mass)
        {
            mflt_mono_mass_start = min_mono_mass;
            mflt_mono_mass_end = max_mono_mass;
        }

        // This function enables the default constraints and assumes ppm units for the mass constraints
        public void SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan,
            float wt_fit, float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net,
            float wt_ims_drift_time, bool use_cs)
        {
            mflt_wt_mono_mass = wt_mono_mass;
            mflt_constraint_mono_mass = mono_constraint;
            mbln_constraint_mono_mass_is_ppm = true;

            mflt_wt_average_mass = wt_avg_mass;
            mflt_constraint_average_mass = avg_constraint;
            mbln_constraint_average_mass_is_ppm = true;

            mflt_wt_log_abundance = wt_log_abundance;
            mflt_wt_scan = wt_scan;
            mflt_wt_net = wt_net;
            mflt_wt_fit = wt_fit;
            mflt_wt_ims_drift_time = wt_ims_drift_time;

            mdbl_max_distance = max_dist;
            mbln_use_net = use_net;
            mbln_constraint_charge_state = use_cs;
        }

        // This function allows one to specify the units for the constraints
        public void SetOptionsEx(float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
            float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm, float wt_log_abundance, float wt_scan,
            float wt_net, float wt_fit, double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs)
        {
            mflt_wt_mono_mass = wt_mono_mass;
            mflt_constraint_mono_mass = mono_constraint;
            mbln_constraint_mono_mass_is_ppm = mono_constraint_is_ppm;

            mflt_wt_average_mass = wt_avg_mass;
            mflt_constraint_average_mass = avg_constraint;
            mbln_constraint_average_mass_is_ppm = avg_constraint_is_ppm;

            mflt_wt_log_abundance = wt_log_abundance;
            mflt_wt_scan = wt_scan;
            mflt_wt_net = wt_net;
            mflt_wt_fit = wt_fit;
            mflt_wt_ims_drift_time = wt_ims_drift_time;

            mdbl_max_distance = max_dist;
            mbln_use_net = use_net;

            mbln_constraint_charge_state = use_cs;
            // mbln_is_weighted_euc = use_weighted_euc;
        }

        public void SetLCMinMaxScan(int minScan, int maxScan)
        {
            mint_lc_min_scan = minScan;
            mint_lc_max_scan = maxScan;

            // Do not allow the minimum and maximum scans to be the same number (would lead to divide by zero errors)
            if (mint_lc_min_scan == mint_lc_max_scan)
                mint_lc_max_scan = mint_lc_min_scan + 1;
        }

        public void SetPeks(ref List<IsotopePeak> vectPks)
        {
            mvect_isotope_peaks.Clear();
            mvect_isotope_peaks.AddRange(vectPks);

            var numPeaks = mvect_isotope_peaks.Count;
            for (var i = 0; i < numPeaks; i++)
            {
                var pk = mvect_isotope_peaks[i];
                if (pk.mint_lc_scan > mint_lc_max_scan)
                    mint_lc_max_scan = pk.mint_lc_scan;
                if (pk.mint_lc_scan < mint_lc_min_scan)
                    mint_lc_min_scan = pk.mint_lc_scan;
            }
        }

        public string InputFileName
        {
            get { return mstr_inputFile; }
            set { mstr_inputFile = value; }
        }

        public string OutputDirectory { get; set; }

        public float GetSegmentSize()
        {
            return mflt_segment_size;
        }

        public bool SortUMCsByMonoMassAndScan(UMC a, UMC b)
        {
            if (a.mdbl_average_mono_mass < b.mdbl_average_mono_mass)
            {
                return true;
            }

            return false;
        }
    }
}