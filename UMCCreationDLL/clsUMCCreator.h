// clsUMCCreator.h
#include "UMCCreator.h"
#pragma once

using namespace System;
namespace UMCCreation
{

	public __gc class clsIsotopePeak
	{
	public:
		int mint_original_index; 
		int mint_umc_index ; 
		int mint_scan ; 
		short mshort_charge ; 
		double mdbl_abundance ; 
		double mdbl_mz ; 
		float mflt_fit ; 
		double mdbl_average_mass ; 
		double mdbl_mono_mass ; 
		double mdbl_max_abundance_mass ;
		double mdbl_i2_abundance ; 
	} ; 

	public __value enum enmStatus { IDLE = 0, LOADING, CLUSTERING, SUMMARIZING, FAILED, COMPLETE } ; 
	public __gc class clsUMCCreator
	{
		int mint_min_umc_length ; 
		int mint_percent_done ; 
		System::String *mstr_message ; 
		System::String *mstr_file_name ; 
		enmStatus menm_status ; 

		// TODO: Add your methods for this class here.
		UMCCreator __nogc *mobj_umc_creator ; 
		void LoadFindUMCs(bool is_pek_file) ; 

	public:
		clsUMCCreator() ; 
		~clsUMCCreator() ; 

		void FindUMCs() ;
		void LoadFindUMCsPEK() ; 
		void LoadFindUMCsCSV() ; 
		void ResetStatus() ; 

		void SetIsotopePeaks(clsIsotopePeak* (&isotope_peaks) __gc[]) ; 
		int GetUmcMapping(int (&isotope_peaks_index) __gc[], int (&umc_index) __gc[]) ; 

		void SetMinMaxScans(int min, int max) ; 
		void SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan, float wt_fit,
			float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net) ; 

		void SetOptionsEx(
				float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
				float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm,
				float wt_log_abundance, float wt_scan, float wt_net, float wt_fit,
				double max_dist, bool use_net) ;

		UMCManipulation::clsUMC* GetUMCs()[] ; 

		__property short get_PercentComplete()
		{
			switch (menm_status)
			{
				case enmStatus::LOADING:
					if (mobj_umc_creator != NULL)
						return mobj_umc_creator->GetPercentComplete() ; 
					return 0 ; 
					break ; 
				case enmStatus::CLUSTERING:
					if (mobj_umc_creator != NULL)
						return mobj_umc_creator->GetPercentComplete() ; 
					return 0 ; 
					break ; 
				case enmStatus::SUMMARIZING:
					if (mobj_umc_creator != NULL)
						return mobj_umc_creator->GetPercentComplete() ; 
					return 0 ; 
					break ; 
				case enmStatus::COMPLETE:
					return 100 ; 
					break ; 
			}
			return 0 ; 
		}

		__property enmStatus get_Status()
		{
			return menm_status ; 
		}

		__property System::String* get_Message()
		{
			return mstr_message ; 
		} 

		__property System::String* get_FileName()
		{
			return mstr_file_name ; 
		} 

		__property void set_FileName(System::String *fileName)
		{
			mstr_file_name = fileName ; 
		} 

		__property int get_MinUMCLength()
		{
			return mint_min_umc_length ; 
		}

		__property void set_MinUMCLength(int len)
		{
			mint_min_umc_length = len ; 
		}

		__property int get_MinScan()
		{
			return mobj_umc_creator->mint_min_scan ; 
		}

		__property int get_MaxScan()
		{
			return mobj_umc_creator->mint_max_scan ; 
		}



	};
}
