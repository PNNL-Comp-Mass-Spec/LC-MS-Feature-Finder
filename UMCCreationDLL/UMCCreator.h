#pragma once
#include "IsotopePeak.h" 
#include <vector> 
#include <map> 
#include <math.h> 
#include <float.h> 
#include "UMC.h" 

class UMCCreator
{
	// Weights to use when clustering points
	// Set a weight to 0 to effectively disable that weight
	float mflt_wt_mono_mass ;
	float mflt_wt_average_mass ; 
	float mflt_wt_log_abundance ; 
	float mflt_wt_scan ; 
	float mflt_wt_net ; 
	float mflt_wt_fit ; 
	float mflt_wt_ims_drift_time ;

	float mflt_constraint_mono_mass ; 
	bool mbln_constraint_mono_mass_is_ppm ;

	float mflt_constraint_average_mass ; 
	bool mbln_constraint_average_mass_is_ppm ;

	double mdbl_max_distance ; 
	short mshort_percent_complete ; 
	
	bool mbln_use_net ;		// When True, then uses NET and not Scan

public:
	int mint_min_scan ; 
	int mint_max_scan ; 

	std::multimap<int, int> mmultimap_umc_2_peak_index ; 
	std::vector<IsotopePeak> mvect_isotope_peaks ; 
	std::vector<int> mvect_umc_num_members ; 
	std::vector<UMC> mvect_umcs ; 

	UMCCreator(void);
	~UMCCreator(void);

	short GetPercentComplete() { return mshort_percent_complete ; } ; 
	inline double PeakDistance(IsotopePeak &a, IsotopePeak &b) 
	{
		if (mbln_constraint_mono_mass_is_ppm) {
			if (a.mdbl_mono_mass > 0 && (abs((a.mdbl_mono_mass - b.mdbl_mono_mass) * mflt_wt_mono_mass / a.mdbl_mono_mass * 1000000) > mflt_constraint_mono_mass))
				return DBL_MAX ; 
		} else {
			if (abs((a.mdbl_mono_mass - b.mdbl_mono_mass)) * mflt_wt_mono_mass > mflt_constraint_mono_mass)
				return DBL_MAX ; 
		}

		if (mbln_constraint_average_mass_is_ppm) {
			if (a.mdbl_average_mass > 0 && (abs((a.mdbl_average_mass - b.mdbl_average_mass) * mflt_wt_average_mass / a.mdbl_average_mass * 1000000) > mflt_constraint_average_mass))
				return DBL_MAX ; 
		} else {
			if (abs((a.mdbl_average_mass - b.mdbl_average_mass)) * mflt_wt_average_mass > mflt_constraint_average_mass)
				return DBL_MAX ; 
		}
			
		double a_log_abundance = log10(a.mdbl_abundance) ; 
		double b_log_abundance = log10(b.mdbl_abundance) ; 

		double sqrDist = 0 ;
		
		sqrDist += (a.mdbl_mono_mass - b.mdbl_mono_mass) * (a.mdbl_mono_mass - b.mdbl_mono_mass) * mflt_wt_mono_mass * mflt_wt_mono_mass  ; 
		sqrDist += (a.mdbl_average_mass - b.mdbl_average_mass) * (a.mdbl_average_mass - b.mdbl_average_mass) * mflt_wt_average_mass * mflt_wt_average_mass ; 
		sqrDist += (a_log_abundance - b_log_abundance) * (a_log_abundance - b_log_abundance) * mflt_wt_log_abundance * mflt_wt_log_abundance; 
		
		if (mbln_use_net)
		{
			// Convert scan difference to Generic NET
			double net_distance = (a.mint_scan - b.mint_scan) * 1.0 / (mint_max_scan - mint_min_scan) ; 
			sqrDist += net_distance * net_distance * mflt_wt_net * mflt_wt_net ; 
		} else {
			sqrDist += (a.mint_scan - b.mint_scan) * (a.mint_scan - b.mint_scan) * mflt_wt_scan * mflt_wt_scan ; 
		}

		sqrDist += (a.mflt_fit - b.mflt_fit) * (a.mflt_fit - b.mflt_fit) * mflt_wt_fit * mflt_wt_fit ; 

		// IMS Drift time
		sqrDist += (a.mflt_ims_drift_time - b.mflt_ims_drift_time) * (a.mflt_ims_drift_time - b.mflt_ims_drift_time) * mflt_wt_ims_drift_time * mflt_wt_ims_drift_time ; 

		return sqrt(sqrDist) ; 

	}

	int GetNumUmcs() { return mvect_umcs.size() ; } ; 
	void ReadCSVFile(char *fileName) ; 
	void ReadPekFileMemoryMapped(char *fileName) ; 
	void ReadPekFile(char *fileName) ; 
	void CreateUMCsSinglyLinkedWithAll() ;
	void RemoveShortUMCs(int min_length) ; 
	void CalculateUMCs() ; 
	void PrintPeaks() ; 
	void PrintUMCs(bool print_members) ; 
	void Reset() ; 
	void SetUseNet(bool use) { mbln_use_net = use ; } ; 

	// This function enables the default constraints and assumes ppm units for the mass constraints
	void SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan, float wt_fit,
		float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net, float wt_ims_drift_time)
	{
		mflt_wt_mono_mass = wt_mono_mass; 
		mflt_constraint_mono_mass = mono_constraint ; 
		mbln_constraint_mono_mass_is_ppm = true ;

		mflt_wt_average_mass = wt_avg_mass ; 
		mflt_constraint_average_mass = avg_constraint ; 
		mbln_constraint_average_mass_is_ppm = true ;

		mflt_wt_log_abundance = wt_log_abundance ; 
		mflt_wt_scan = wt_scan ; 
		mflt_wt_net = wt_net ; 
		mflt_wt_fit = wt_fit ;
		mflt_wt_ims_drift_time = wt_ims_drift_time ;

		mdbl_max_distance = max_dist ; 
		mbln_use_net = use_net ;
	}

	// This function allows one to specify the units for the constraints
	void SetOptionsEx(
			float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
			float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm,
			float wt_log_abundance, float wt_scan, float wt_net, float wt_fit,
			double max_dist, bool use_net, float wt_ims_drift_time)
	{
		mflt_wt_mono_mass = wt_mono_mass; 
		mflt_constraint_mono_mass = mono_constraint ; 
		mbln_constraint_mono_mass_is_ppm = mono_constraint_is_ppm ;

		mflt_wt_average_mass = wt_avg_mass ; 
		mflt_constraint_average_mass = avg_constraint ; 
		mbln_constraint_average_mass_is_ppm = avg_constraint_is_ppm ;

		mflt_wt_log_abundance = wt_log_abundance ; 
		mflt_wt_scan = wt_scan ; 
		mflt_wt_net = wt_net ; 
		mflt_wt_fit = wt_fit ;
		mflt_wt_ims_drift_time = wt_ims_drift_time ;

		mdbl_max_distance = max_dist ; 
		mbln_use_net = use_net ;
	}

	void SetMinMaxScan(int minScan, int maxScan) { 
		mint_min_scan = minScan ;
		mint_max_scan = maxScan ; 

		// Do not allow the minimum and maximum scans to be the same number (would lead to divide by zero errors)
		if (mint_min_scan == mint_max_scan)
				mint_max_scan = mint_min_scan + 1 ;
	} ; 

	void SetPeks(std::vector<IsotopePeak> &vectPks) ; 

};
