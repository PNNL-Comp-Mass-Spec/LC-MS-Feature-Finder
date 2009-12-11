#pragma once
#include "IsotopePeak.h" 
#include <vector> 
#include <map> 
#include <math.h> 
#include <float.h> 
#include "UMC.h" 

class UMCCreator
{

	char mstr_inputFile[512];
	char outputDir[512];

	//data filters when loading data isotopic_fit, int min_intensity, int mono_mass_start, int mono_mass_end, bool process_mass_seg, int maxDataPoints, int monoMassSegOverlap
	float mflt_isotopic_fit_filter;
	int mint_min_intensity;
	float mflt_mono_mass_start;
	float mflt_mono_mass_end;
	bool mbln_process_mass_seg;
	int mint_max_data_points;
	int mint_mono_mass_seg_overlap;
	int mint_ims_min_scan_filter;
	int mint_ims_max_scan_filter;
	int mint_lc_min_scan_filter;
	int mint_lc_max_scan_filter;
	


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

	bool mbln_constraint_charge_state;

	double mdbl_max_distance ; 
	short mshort_percent_complete ; 
	
	bool mbln_use_net ;		// When True, then uses NET and not Scan
	bool mbln_is_ims_data;
	bool mbln_is_weighted_euc;

	float mflt_segment_size;

public:
	int mint_lc_min_scan ; 
	int mint_lc_max_scan ; 
	int mint_ims_min_scan;
	int mint_ims_max_scan;


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
			double net_distance = (a.mint_lc_scan - b.mint_lc_scan) * 1.0 / (mint_lc_max_scan - mint_lc_min_scan) ; 
			sqrDist += net_distance * net_distance * mflt_wt_net * mflt_wt_net ; 
		} else {
			sqrDist += (a.mint_lc_scan - b.mint_lc_scan) * (a.mint_lc_scan - b.mint_lc_scan) * mflt_wt_scan * mflt_wt_scan ; 
		}

		sqrDist += (a.mflt_fit - b.mflt_fit) * (a.mflt_fit - b.mflt_fit) * mflt_wt_fit * mflt_wt_fit ; 

		// IMS Drift time
		sqrDist += (a.mflt_ims_drift_time - b.mflt_ims_drift_time) * (a.mflt_ims_drift_time - b.mflt_ims_drift_time) * mflt_wt_ims_drift_time * mflt_wt_ims_drift_time ; 

		return sqrt(sqrDist) ; 

	}

	int GetNumUmcs() { return mvect_umcs.size() ; } ; 
	int ReadCSVFile(char *fileName) ; 
	int ReadCSVFile();
	void ReadPekFileMemoryMapped(char *fileName) ; 
	void ReadPekFile(char *fileName) ; 
	void CreateUMCsSinglyLinkedWithAll() ;
	void RemoveShortUMCs(int min_length) ; 
	void CalculateUMCs() ; 
	void PrintPeaks() ; 
	void PrintUMCs(bool print_members) ; 
	bool PrintUMCs(FILE *stream, bool print_members);
	bool UMCCreator::PrintUMCs(FILE *stream, bool print_members, int featureStartIndex);
	bool PrintMapping(FILE *stream);
	bool PrintMapping(FILE *stream, int featureStartIndex);

	void Reset() ; 
	void SetUseNet(bool use) { mbln_use_net = use ; } ; 
	bool ConsiderPeak(IsotopePeak pk);
	float GetLastMonoMassLoaded();
	void SerializeObjects();
	void DeserializeObjects();
	int LoadPeaksFromDatabase();


	//Functions added by Anuj Shah
	void CreateUMCsSingleLinkedWithAllOnline();
	
	void CreateUMCFromIsotopePeak(IsotopePeak startPeak, UMC &firstUMC);
	void AddPeakToUMC (IsotopePeak peak, UMC &umc);
	bool withinMassTolerance(double observedMass, double realMass);
	int findCandidateUMCsForPeak(IsotopePeak peak, std::vector<UMC> &umcVector, std::vector<UMC> &candidateUMCs);

	void SetFilterOptions(float isotopic_fit, int min_intensity, int min_lc_scan, int max_lc_scan, int min_ims_scan, int max_ims_scan, float mono_mass_start, float mono_mass_end, bool process_mass_seg, int max_data_points, int mono_mass_seg_overlap, float mono_mass_seg_size){
		mflt_isotopic_fit_filter = isotopic_fit;
		mint_min_intensity = min_intensity;
		mflt_mono_mass_start = mono_mass_start;
		mflt_mono_mass_end = mono_mass_end;
		mbln_process_mass_seg = process_mass_seg;
		mint_max_data_points = max_data_points;
		mint_mono_mass_seg_overlap = mono_mass_seg_overlap;
		mint_lc_min_scan_filter = min_lc_scan;
		mint_lc_max_scan_filter = max_lc_scan;
		mint_ims_min_scan_filter = min_ims_scan;
		mint_ims_max_scan_filter = max_ims_scan;
		mflt_segment_size = mono_mass_seg_size;
	}

	void SetMassRange ( float min_mono_mass, float max_mono_mass ){
		mflt_mono_mass_start = min_mono_mass;
		mflt_mono_mass_end = max_mono_mass;
	}



	// This function enables the default constraints and assumes ppm units for the mass constraints
	void SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan, float wt_fit,
		float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs)
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
		mbln_constraint_charge_state = use_cs;
	}

	// This function allows one to specify the units for the constraints
	void SetOptionsEx(
			float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
			float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm,
			float wt_log_abundance, float wt_scan, float wt_net, float wt_fit,
			double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs, bool use_weighted_euc)
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

		mbln_constraint_charge_state = use_cs;
		mbln_is_weighted_euc = use_weighted_euc;
	}

	void SetLCMinMaxScan(int minScan, int maxScan) { 
		mint_lc_min_scan = minScan ;
		mint_lc_max_scan = maxScan ; 

		// Do not allow the minimum and maximum scans to be the same number (would lead to divide by zero errors)
		if (mint_lc_min_scan == mint_lc_max_scan)
				mint_lc_max_scan = mint_lc_min_scan + 1 ;
	} 

	void SetPeks(std::vector<IsotopePeak> &vectPks) ; 


	void SetInputFileName ( char * filename);
	void SetOutputDiretory( char * dir);

	char * GetOutputDirectory(){
		return outputDir;
	}

	float GetSegmentSize(){
		return mflt_segment_size;
	}

};
