// clsUMCCreator.h
#include "UMCCreator.h"
#include <stdlib.h>
#include <iostream> 
#include <fstream>
#include <stddef.h>
#include <ctime>
#pragma once

using namespace System;
namespace UMCCreation
{

	public __gc class clsIsotopePeak
	{
	public:
		int mint_original_index; 
		int mint_umc_index ; 
		int mint_lc_scan ; 
		short mshort_charge ; 
		double mdbl_abundance ; 
		double mdbl_mz ; 
		float mflt_fit ; 
		double mdbl_average_mass ; 
		double mdbl_mono_mass ; 
		double mdbl_max_abundance_mass ;
		double mdbl_i2_abundance ; 
		float mflt_ims_drift_time ;
	} ; 

	public __value enum enmStatus { IDLE = 0, LOADING, CLUSTERING, SUMMARIZING, FAILED, COMPLETE, CHUNKING } ; 
	public __gc class clsUMCCreator
	{
		int mint_min_umc_length ; 
		int mint_percent_done ; 
		bool mbln_process_chunks;
		float mflt_mono_mass_start;
		float mflt_mono_mass_end;
		int mint_mono_mass_overlap;

		System::String *mstr_message ; 
		System::String *mstr_file_name ; 
		System::String *mstr_options_name;
		System::String *mstr_output_direc;
		enmStatus menm_status ; 

		// TODO: Add your methods for this class here.
		UMCCreator __nogc *mobj_umc_creator ; 
		void LoadFindUMCs(bool is_pek_file) ; 

	private:
		FILE *mfile_logFile;
		System::String *mstr_baseFileName;

		char* CreateBaseFileName(char* directoryName, char* inputFileName);
		void createLogFile();
		void log(char* textToLog);
		void log(char* textToLog, int numToLog);
		
	public:
		clsUMCCreator() ; 
		~clsUMCCreator() ; 

		void LoadFindUMCs();
		void FindUMCs() ;
		void LoadFindUMCsPEK() ; 
		void LoadFindUMCsCSV() ; 
		void ResetStatus() ; 

		void SetIsotopePeaks(clsIsotopePeak* (&isotope_peaks) __gc[]) ; 
		bool LoadProgramOptions(); 
		int GetUmcMapping(int (&isotope_peaks_index) __gc[], int (&umc_index) __gc[]) ; 

		void SetLCMinMaxScans(int min, int max) ; 
		void SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan, float wt_fit,
			float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs) ; 


		bool PrintUMCsToFile();
		bool PrintUMCsToFile(int chunkIndex, int featureStartIndex);

		//MaxIsotopicFit=0.15
		//MinimumIntensity=0
		//MonoMassStart=0
		//MonoMassEnd=0
		//ProcessDataInMonoMassSegments=False
		//MaxDataPointsPerMonoMassSegment=1000000
		//MonoMassSegmentOverlapDa=2
		void SetFilterOptions(float isotopic_fit, int min_intensity, int min_lc, int max_lc, int min_ims, int max_ims, float mono_mass_start, float mono_mass_end, bool process_mass_seg, int maxDataPoints, int monoMassSegOverlap, float segmentSize);

		void SetOptionsEx(
				float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
				float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm,
				float wt_log_abundance, float wt_scan, float wt_net, float wt_fit,
				double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs, bool use_wt_euc) ;

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

		__property System::String* get_OutputFileName(){
			return mstr_output_direc;
		}

		
		__property void set_OutputFileName(System::String *output){
			mstr_output_direc = output;
		}

		__property System::String* get_OptionsFileName(){
			return mstr_options_name;
		}

		__property void set_OptionsFileName(System::String *options){
			mstr_options_name = options;
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
			return mobj_umc_creator->mint_lc_min_scan ; 
		}

		__property int get_MaxScan()
		{
			return mobj_umc_creator->mint_lc_max_scan ; 
		}



	};
}
