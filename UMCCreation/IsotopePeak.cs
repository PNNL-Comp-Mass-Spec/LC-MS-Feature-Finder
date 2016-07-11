using System;

namespace UMCCreation
{
    /// <summary>
    /// Original c++ native version of IsotopePeak
    /// </summary>
    public class IsotopePeak : IComparable<IsotopePeak>
    {
        public int mint_original_index;
        public int mint_line_number_in_file;
        public int mint_umc_index;

        //Anuj: This should now be used to represent the LC frame number
        public int mint_lc_scan;
        public short mshort_charge;
        public double mdbl_abundance;
        public double mdbl_mz;
        public float mflt_fit;
        public double mdbl_average_mass;
        public double mdbl_mono_mass;
        public double mdbl_max_abundance_mass;
        public double mdbl_i2_abundance;
        public double mdbl_mono_abundance;

        //Anuj: Parameters to account for added IMS data dimensions
        public int mint_ims_scan;
        public float mflt_ims_drift_time;
        public float mflt_orig_intensity;
        public float mflt_tia_orig_intensity;
        public float mflt_cum_drift_time;

        public IsotopePeak Clone()
        {
            return (IsotopePeak) this.MemberwiseClone();
        }

        public void printPeak()
        {
            System.Console.WriteLine("{0}\t{1}\t{2:F0}\t{3:F4}\t{4:F3}\t{5:F4}\t{6:F4}\t{7:F4}", mint_lc_scan,
                mshort_charge, mdbl_abundance, mdbl_mz, mflt_fit, mdbl_average_mass, mdbl_mono_mass,
                mdbl_max_abundance_mass);
        }

        public override string ToString()
        {
            return mint_lc_scan + " " + mdbl_mz + " " + mdbl_abundance;
        }

        public int CompareTo(IsotopePeak other)
        {
            var result = this.mdbl_mono_mass.CompareTo(other.mdbl_mono_mass);
            if (result == 0)
            {
                result = this.mint_lc_scan.CompareTo(other.mint_lc_scan); // This results in a stable sort, if the input is sorted by scan num (and is not IMS data)
            }
            return result;
        }
    }
}
