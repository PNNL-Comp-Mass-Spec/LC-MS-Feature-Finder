#pragma once

class IsotopePeak
{
public:
	int mint_original_index; 
	int mint_line_number_in_file;
	int mint_umc_index ; 

	//Anuj: This should now be used to represent the LC frame number
	int mint_lc_scan ; 
	short mshort_charge ; 
	double mdbl_abundance ; 
	double mdbl_mz ; 
	float mflt_fit ; 
	double mdbl_average_mass ; 
	double mdbl_mono_mass ; 
	double mdbl_max_abundance_mass ;
	double mdbl_i2_abundance ; 
	double mdbl_mono_abundance ; 
	
	//Anuj: Parameters to account for added IMS data dimensions
	int mint_ims_scan;
	float mflt_ims_drift_time ;
	float mflt_orig_intensity;
	float mflt_tia_orig_intensity;
	float mflt_cum_drift_time;


	IsotopePeak(void);
	~IsotopePeak(void);


	void IsotopePeak::printPeak();
	
	
};
