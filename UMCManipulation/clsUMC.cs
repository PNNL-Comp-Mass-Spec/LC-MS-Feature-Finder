using System;

namespace UMCManipulation
{
	/// <summary>
	/// Summary description for clsUMC.
	/// </summary>
	public class clsUMC
	{		
		public int mint_umc_index ; 
		public double mdbl_mono_mass ; 
		public double mdbl_class_rep_mz ; 
		public int mint_class_rep_charge ; 
		public int mint_scan ;
		public int mint_start_scan ; 
		public int mint_end_scan ; 
		public double mdbl_net ; 
		public double mdbl_mono_mass_calibrated ; 
		public int mint_scan_aligned ; 
		public double mdbl_abundance ; 

		public clsUMC()
		{
			//
			// TODO: Add constructor logic here
			//
			mint_umc_index = -1 ; 
			mdbl_mono_mass = 0 ; 
			mdbl_class_rep_mz = 0 ; 
			mint_class_rep_charge = -1 ; 
			mint_scan = -1 ; 
			mint_start_scan = -1 ; 
			mint_end_scan = -1 ; 
			mdbl_net = 0 ; 
			mdbl_mono_mass_calibrated = 0 ; 
			mint_scan_aligned = - 1 ; 
			mdbl_abundance = 0 ; 
		}
	}

}
