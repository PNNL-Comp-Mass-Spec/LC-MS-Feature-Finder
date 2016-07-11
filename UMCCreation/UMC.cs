using System;

namespace UMCCreation
{
    public struct UMC : IComparable<UMC>
    {
        public double mdbl_min_mono_mass;
        public double mdbl_max_mono_mass;
        public double mdbl_average_mono_mass;
        public double mdbl_median_mono_mass;

        public double mdbl_max_abundance;
        public double mdbl_sum_abundance;

        public int min_num_members;
        public int mint_start_scan;
        public int mint_stop_scan;
        public int mint_max_abundance_scan;

        public double mdbl_class_rep_mz;
        public short mshort_class_rep_charge;

        public int mint_umc_index;

        //Added by Anuj Shah to keep track of the last peak that belongs
        //to this UMC
        public IsotopePeak lastPeak;

        public int CompareTo(UMC other)
        {
            return this.mdbl_average_mono_mass.CompareTo(other.mdbl_average_mono_mass);
        }
    }
}
