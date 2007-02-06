#pragma once

class IsotopePeak
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
	double mdbl_mono_abundance ; 
	IsotopePeak(void);
	~IsotopePeak(void);
};
