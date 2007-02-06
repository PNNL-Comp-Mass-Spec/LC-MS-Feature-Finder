#pragma once

class UMC
{
public:
	double mdbl_min_mono_mass ; 
	double mdbl_max_mono_mass ; 
	double mdbl_average_mono_mass ; 
	double mdbl_median_mono_mass ; 

	double mdbl_max_abundance ;
	double mdbl_sum_abundance ; 

	int min_num_members ; 
	int mint_start_scan ; 
	int mint_stop_scan ; 
	int mint_max_abundance_scan ; 

	double mdbl_class_rep_mz ;
	short mshort_class_rep_charge ; 

	int mint_umc_index ; 

	UMC(void);
	~UMC(void);
};
