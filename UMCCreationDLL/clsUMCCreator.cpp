// This is the main DLL file.

#include "clsUMCCreator.h"
#using <mscorlib.dll>

namespace UMCCreation
{
	void GetStr(System::String *src, char *dest)
	{
		if (src == 0 || src == "" || src->get_Length() == 0)
		{
			dest[0] = '\0' ; 
			return ; 
		}

		int len = src->get_Length() ; 
		for (int i = 0 ; i < len ; i++)
		{
			dest[i] = (char) src->Chars[i] ; 
		}
		dest[len] = '\0' ; 
	}

	clsUMCCreator::clsUMCCreator()
	{
		mobj_umc_creator = new UMCCreator() ; 
	}

	clsUMCCreator::~clsUMCCreator()
	{
		if (mobj_umc_creator != NULL)
			delete mobj_umc_creator ; 
	}

	void clsUMCCreator::ResetStatus()
	{
		menm_status = IDLE ; 
		mobj_umc_creator->Reset() ; 
	}

	void clsUMCCreator::FindUMCs()
	{
		menm_status = CLUSTERING ; 
		mstr_message = new System::String("Clustering Isotope Peaks") ; 
		mobj_umc_creator->CreateUMCsSinglyLinkedWithAll() ; 
		menm_status = SUMMARIZING ; 
		mstr_message = new System::String("Filtering out short clusters") ; 
		mobj_umc_creator->RemoveShortUMCs(mint_min_umc_length) ;
		mstr_message = new System::String("Calculating UMC statistics") ; 
		mobj_umc_creator->CalculateUMCs() ;
		menm_status = COMPLETE ; 
	}

	void clsUMCCreator::LoadFindUMCs(bool is_pek_file)
	{
		Console::WriteLine(S"Loading UMCs") ; 
		char file_name[512] ; 
		GetStr(mstr_file_name, file_name) ; 
		menm_status = LOADING ;

		if (is_pek_file)
		{
			mstr_message = new System::String("Loading PEK file") ; 
			mobj_umc_creator->ReadPekFileMemoryMapped(file_name) ; 
		}
		else
		{
			mstr_message = new System::String("Loading CSV file") ; 
			mobj_umc_creator->ReadCSVFile(file_name) ; 
		}

		menm_status = CLUSTERING ; 
		mstr_message = new System::String("Clustering Isotope Peaks") ; 
		mobj_umc_creator->CreateUMCsSinglyLinkedWithAll() ; 
		menm_status = SUMMARIZING ; 
		mstr_message = new System::String("Filtering out short clusters") ; 
		mobj_umc_creator->RemoveShortUMCs(mint_min_umc_length) ;
		mstr_message = new System::String("Calculating UMC statistics") ; 
		mobj_umc_creator->CalculateUMCs() ;
	}

	void clsUMCCreator::LoadFindUMCsPEK()
	{
		LoadFindUMCs(true) ; 
		menm_status = COMPLETE ; 
	}

	void clsUMCCreator::LoadFindUMCsCSV()
	{
		LoadFindUMCs(false) ; 
		menm_status = COMPLETE ; 
	}

	int clsUMCCreator::GetUmcMapping(int (&isotope_peaks_index) __gc[], int (&umc_index) __gc[])
	{
		int numMappings = mobj_umc_creator->mmultimap_umc_2_peak_index.size() ; 
		isotope_peaks_index = new int __gc [numMappings] ; 
		umc_index = new int __gc [numMappings] ; 

		int mappingNum = 0 ; 
		for (std::multimap<int,int>::iterator iter = mobj_umc_creator->mmultimap_umc_2_peak_index.begin() ; iter != mobj_umc_creator->mmultimap_umc_2_peak_index.end() ; iter++)
		{
			int currentUmcNum = (*iter).first ; 
			int pkIndex = (*iter).second ; 
			isotope_peaks_index[mappingNum] = pkIndex ; 
			umc_index[mappingNum] = currentUmcNum ; 
			mappingNum++ ; 
		}
		return numMappings ; 
	}

	UMCManipulation::clsUMC* clsUMCCreator::GetUMCs()[]
	{
		int numUmcs = mobj_umc_creator->GetNumUmcs(); 
		UMCManipulation::clsUMC* arr_umcs __gc[] ; 
		arr_umcs = new UMCManipulation::clsUMC* __gc [numUmcs] ; 
		
		for (int umcNum = 0 ; umcNum < numUmcs ; umcNum++)
		{
			UMC umc = mobj_umc_creator->mvect_umcs[umcNum] ; 
			UMCManipulation::clsUMC *newUmc = __gc new UMCManipulation::clsUMC() ;
			newUmc->mdbl_abundance = umc.mdbl_sum_abundance ; 

			newUmc->mdbl_class_rep_mz = (double) umc.mdbl_class_rep_mz ; 
			newUmc->mdbl_mono_mass = (double) umc.mdbl_median_mono_mass ; 
			newUmc->mdbl_mono_mass_calibrated = newUmc->mdbl_mono_mass ; 

			newUmc->mint_scan = umc.mint_max_abundance_scan ; 
			newUmc->mint_start_scan = umc.mint_start_scan ; 
			newUmc->mint_end_scan = umc.mint_stop_scan ; 

			newUmc->mdbl_net = (double) umc.mint_max_abundance_scan ; 
            if (mobj_umc_creator->mint_max_scan > mobj_umc_creator->mint_min_scan) {
			     // Compute Generic NET value
			     newUmc->mdbl_net = (double) (umc.mint_max_abundance_scan - mobj_umc_creator->mint_min_scan) * 1.0 / (mobj_umc_creator->mint_max_scan - mobj_umc_creator->mint_min_scan) ; 
            }

			newUmc->mint_scan_aligned = umc.mint_max_abundance_scan ; 

			newUmc->mint_umc_index = umc.mint_umc_index ; 
			newUmc->mint_class_rep_charge = (int) umc.mshort_class_rep_charge ; 
			arr_umcs[umcNum] = newUmc ; 
		}

		return arr_umcs ; 
	}

	void clsUMCCreator::SetMinMaxScans(int min, int max)
	{
		mobj_umc_creator->SetMinMaxScan(min, max) ; 

		if (max <= min)
			throw new exception("Max scan must be greater than min scan");
	}

	void clsUMCCreator::SetIsotopePeaks(clsIsotopePeak* (&isotope_peaks) __gc[])
	{
		int numPeaks = isotope_peaks->Length ; 
		std::vector<IsotopePeak> vectPeaks ;
		vectPeaks.reserve(numPeaks) ; 

		clsIsotopePeak *isoPk ; 
		for (int pkNum = 0 ; pkNum < numPeaks ; pkNum++)
		{
			IsotopePeak pk ; 
			isoPk = isotope_peaks[pkNum] ; 

			pk.mint_original_index = isoPk->mint_original_index ; 
			pk.mint_umc_index = isoPk->mint_umc_index ; 
			pk.mint_scan = isoPk->mint_scan ; 
			pk.mshort_charge = isoPk->mshort_charge ; 
			pk.mdbl_abundance = isoPk->mdbl_abundance ; 
			pk.mdbl_mz = isoPk->mdbl_mz ; 
			pk.mflt_fit = isoPk->mflt_fit ; 
			pk.mdbl_average_mass = isoPk->mdbl_average_mass ; 
			pk.mdbl_mono_mass = isoPk->mdbl_mono_mass ; 
			pk.mdbl_max_abundance_mass = isoPk->mdbl_max_abundance_mass ;
			pk.mdbl_i2_abundance = isoPk->mdbl_i2_abundance ; 

			vectPeaks.push_back(pk) ; 
		}
		mobj_umc_creator->SetPeks(vectPeaks) ; 
	}
	void clsUMCCreator::SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan, float wt_fit,
			float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net)
	{
		mobj_umc_creator->SetOptions(wt_mono_mass, wt_avg_mass, wt_log_abundance, wt_scan, wt_fit, wt_net, 
			mono_constraint, avg_constraint, max_dist, use_net) ; 
	}

	void clsUMCCreator::SetOptionsEx(
					float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
					float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm,
					float wt_log_abundance, float wt_scan, float wt_net, float wt_fit,
					double max_dist, bool use_net)
	{
		mobj_umc_creator->SetOptionsEx(
					wt_mono_mass, mono_constraint, mono_constraint_is_ppm,
					wt_avg_mass, avg_constraint, avg_constraint_is_ppm,
					wt_log_abundance, wt_scan, wt_net, wt_fit, 
					max_dist, use_net) ;
	}

}
