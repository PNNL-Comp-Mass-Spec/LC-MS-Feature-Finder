#include <iostream>
#include ".\isotopepeak.h"

IsotopePeak::IsotopePeak(void)
{
}

IsotopePeak::~IsotopePeak(void)
{
}


void IsotopePeak::printPeak(){
	
		std::cout.precision(0) ; 
		std::cout<<mint_lc_scan<<"\t"<<mshort_charge<<"\t"<<mdbl_abundance<<"\t" ; 
		std::cout.precision(4) ; 
		std::cout<<mdbl_mz<<"\t" ; 
		std::cout.precision(3) ; 
		std::cout<<mflt_fit<<"\t" ; 
		std::cout.precision(4) ; 

		std::cout<<mdbl_average_mass ; 
		std::cout<<"\t"<<mdbl_mono_mass ; 
		std::cout<<"\t"<<mdbl_max_abundance_mass<<"\n" ; 
}