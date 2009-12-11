#include ".\umccreator.h"
#include "MemMappedReader.h"
#include <stdlib.h> 
#include <algorithm>
#include <iostream> 
#include <fstream>
#include <cmath>
#include <stddef.h>

#define DEBUG
#define DATAFILTERS
using namespace std;

bool SortUMCsByMonoMassAndScan(UMC &a, UMC &b)
{
	if ( a.mdbl_average_mono_mass < b.mdbl_average_mono_mass){
		return true;
	}

	return false;
}


bool SortIsotopesByMonoMassAndScan(IsotopePeak &a, IsotopePeak &b) 
{
	if (a.mdbl_mono_mass < b.mdbl_mono_mass)
		return true ; 

	return false;
}

UMCCreator::UMCCreator(void)
{
	mflt_wt_mono_mass = 0.01F ; // using ppms 10 ppm = length of 0.1 ppm
	mflt_wt_average_mass = 0.01F ; // using ppms 10 ppm = length of 0.1 ppm
	mflt_wt_log_abundance = 0.1F ; 
	mflt_wt_scan = 0.01F ; 
	mflt_wt_fit = 0.1F ; 
	mflt_wt_ims_drift_time = 0.1F ;

	mflt_constraint_mono_mass = 10.0F ; // is in ppm
	mflt_constraint_average_mass = 10.0F ; // is in ppm. 
	mdbl_max_distance = 0.1 ; 
	
	mbln_use_net = true ;
	mbln_constraint_mono_mass_is_ppm = true ;
	mbln_constraint_average_mass_is_ppm = true ;

	mshort_percent_complete = 0 ;

	mint_lc_min_scan = INT_MAX ; 
	mint_lc_max_scan = 0 ;
	mint_ims_min_scan = INT_MAX;
	mint_ims_max_scan = 0;
}

UMCCreator::~UMCCreator(void)
{
}

int UMCCreator::ReadCSVFile(){
	return ReadCSVFile(mstr_inputFile);
}

int UMCCreator::ReadCSVFile(char *fileName)
{
	char startTag[1024] ; //will be defined based on header
	int startTagLength;
	char *stopTag = "Blah" ; 
	int stopTagLen = (int)strlen(stopTag) ; 

	
	Reset() ; 
	MemMappedReader mappedReader ; 
	mappedReader.Load(fileName) ; 
	__int64 file_len = mappedReader.FileLength() ; 


	//columns for IMS data 
	//'frame_num,ims_scan_num,charge,abundance,mz,fit,average_mw,monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance,orig_intensity,TIA_orig_intensity, drift_time,cumulative_drift_time\n'

	//columns for LC-MS data
	//scan_num,charge,abundance,mz,fit,average_mw,monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance/


	bool reading = false ; 
	bool is_first_scan = true ; 

	IsotopePeak pk ; 
	short charge = 0 ; 
	double abundance = 0 ;
	float fit = 0, mz = 0 , averageMass = 0 , monoMass = 0 , maxMass = 0 ;  

	int numPeaks = 0 ; 
	int origLineNumber = 0;

	mshort_percent_complete = 0 ;

	pk.mdbl_abundance = 0 ; 
	pk.mdbl_i2_abundance = 0 ; 
	pk.mdbl_average_mass = 0 ; 
	pk.mflt_fit = 0 ; 
	pk.mdbl_max_abundance_mass = 0 ; 
	pk.mdbl_mono_mass = 0 ; 
	pk.mdbl_mz = 0 ; 
	pk.mshort_charge = 0 ; 
	pk.mflt_ims_drift_time = 0 ;

	mint_lc_min_scan = INT_MAX ; 
	mint_lc_max_scan = 0 ; 

	int pos = 0 ; 
	const int MAX_HEADER_BUFFER_LEN = 1024 ; 
	const int MAX_BUFFER_LEN = 1024 ; 
	char buffer[MAX_BUFFER_LEN] ;

	mshort_percent_complete = (short)((100.0 * mappedReader.CurrentPosition()) / file_len) ; 
	if (mshort_percent_complete > 99)
		mshort_percent_complete = 99 ; 

	bool success = mappedReader.GetNextLine(buffer, MAX_BUFFER_LEN, "\n", MAX_BUFFER_LEN);

#ifdef DBUG
	std::cout << "Success = " << success << "\n";
#endif

	if (success){
		if (buffer[0] == 'f'){
			strcpy(startTag, "frame_num,ims_scan_num,charge,abundance,mz,fit,average_mw,monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance,orig_intensity,TIA_orig_intensity,drift_time");
			mbln_is_ims_data = true;
		}
		else{
			strcpy(startTag,"scan_num,charge,abundance,mz,fit,average_mw, monoisotopic_mw,mostabundant_mw,fwhm,signal_noise,mono_abundance,mono_plus2_abundance") ; 
			mbln_is_ims_data = false;
		}
	}
	
	startTagLength = (int)strlen(startTag) ; 

#ifdef DBUG
	std::cout << "Data file is IMS? (1=Yes, 0=No):: " << mbln_is_ims_data << "\n";
#endif

	//success = mappedReader.SkipToAfterLine(startTag, buffer, startTagLength, MAX_BUFFER_LEN) ; 
	if (!success)
	{
		throw "Incorrect header for file" ; 
	}
	double fwhm = 0, s2n = 0 ;

#ifdef DATAFILTERS
	std::cout<< "***************************************\n";
	std::cout<< "Data filters \n";
	std::cout<< " Minimum LC scan = " << mint_lc_min_scan_filter;
	std::cout<< "\n Maximum LC scan = " << mint_lc_max_scan_filter;
	std::cout<< "\n Minimum IMS scan = " << mint_ims_min_scan_filter;
	std::cout<< "\n Maximum IMS scan = " << mint_ims_max_scan_filter;
	std::cout<< "\n Maximum fit = " << mflt_isotopic_fit_filter;
	std::cout<< "\n Minimum intensity = " << mint_min_intensity;
	std::cout<<"\n Mono mass start = " << mflt_mono_mass_start;
	std::cout<<"\n Mono mass end = " << mflt_mono_mass_end;

	std::cout<< "\n***************************************\n";
#endif

	while(!mappedReader.eof() && mappedReader.GetNextLine(buffer, MAX_BUFFER_LEN, stopTag, stopTagLen))
	{

		mshort_percent_complete = (short)((100.0 * mappedReader.CurrentPosition()) / file_len) ; 

		if (mshort_percent_complete > 99)
			mshort_percent_complete = 99 ; 
		char *stopPtr ; 
		char *stopPtrNext ; 

		//in either case the first value is the lc_scan_num
		pk.mint_lc_scan = strtol(buffer, &stopPtr, 10) ; 


		stopPtr++ ; 

		if (mbln_is_ims_data){
			//then we need to parse out ims_scan number
			pk.mint_ims_scan =strtol(stopPtr, &stopPtrNext, 10) ; 
			stopPtr = ++stopPtrNext;
		}

		pk.mshort_charge = (short) strtol(stopPtr, &stopPtrNext, 10) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mdbl_abundance = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mdbl_mz = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mflt_fit = (float)strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mdbl_average_mass = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mdbl_mono_mass = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mdbl_max_abundance_mass = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		fwhm = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		s2n = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mdbl_mono_abundance = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 
		pk.mdbl_i2_abundance = strtod(stopPtr, &stopPtrNext) ; 
		stopPtr = ++stopPtrNext ; 

		//if it's ims data then we have to read four more columns of data
		if (mbln_is_ims_data){
			//orig intensity, TIA original intensity, drift time and cumulative drift time
			pk.mflt_orig_intensity = strtod(stopPtr, &stopPtrNext);
			stopPtr = ++stopPtrNext;

			pk.mflt_tia_orig_intensity = strtod(stopPtr, &stopPtrNext);
			stopPtr = ++stopPtrNext;
			pk.mflt_ims_drift_time = strtod(stopPtr, &stopPtrNext);
			stopPtr = ++stopPtrNext;
			pk.mflt_cum_drift_time = strtod(stopPtr, &stopPtrNext);
			stopPtr = ++stopPtrNext;

		}

		// Check here to see of MAP index is correct. Dameng
		pk.mint_original_index = numPeaks ; //when reading from the file, this should be the line number in the original isos file
		pk.mint_line_number_in_file = origLineNumber;

		//check if the filter criteria is satisfied before adding a peak
		if ( ConsiderPeak(pk)){

			// Remove the following condition may cause one segment has more than mint_max_data_points
			if ( false ){
				if ( mbln_process_mass_seg ){
					if ( numPeaks > mint_max_data_points ){
						std::cout << "I know this happens" << "\n";
						numPeaks--;
						break;
					}
				}
			}
				#ifdef DBUG
					std::cout << "Adding peak ... "  << ConsiderPeak(pk) << "\n";
				#endif
            
				//check if the min scans and max scans for both lc and ims need to be fixed
				if (pk.mint_lc_scan <= mint_lc_min_scan ){
					mint_lc_min_scan = pk.mint_lc_scan;
				}
	
				if (pk.mint_lc_scan >= mint_lc_max_scan){
					mint_lc_max_scan = pk.mint_lc_scan;
				}

				if ( mbln_is_ims_data ){
					//need to fix the ims min scans and max scans that were loaded
					if ( pk.mint_ims_scan <= mint_ims_min_scan ){
						mint_ims_min_scan = pk.mint_ims_scan;
					}

					if ( pk.mint_ims_scan >= mint_ims_max_scan ){
						mint_ims_max_scan = pk.mint_ims_scan;
					}
				}
	
				mvect_isotope_peaks.push_back(pk) ;

				#ifdef DBUG
					pk.printPeak();
				#endif
				//increment number of peaks read
				numPeaks++ ; 

		}

		origLineNumber++;
		
	}

	
	std::cout << " Total number of peaks we'll consider = " << numPeaks << "\n";
	mappedReader.Close();
	return numPeaks;
}

int UMCCreator::LoadPeaksFromDatabase(){
	//this is to read a sql lite database and retrieve the peaks from it into the mvect_isotope_peaks array.

	return 0;

}

//Method not used, 
void UMCCreator::SerializeObjects(){


	//First write out all loaded isotopic peaks
	//we should be smarter about this file writing since there's going to be sufficient 
	//number of peaks that don't get used into any UMC or get removed as part of short UMCs
	//maybe each peak needs to have a bit whether to used or not
	ofstream ofs("mvect_isotope_peaks.ros", ios::binary);

	//this should write only peaks that are being used to the file
	int count = 0;
	for ( int i = 0; i < mvect_isotope_peaks.size(); i++){

		IsotopePeak pk = (IsotopePeak ) mvect_isotope_peaks[i];
		if ( pk.mint_umc_index != -1){
       		ofs.write( (char*)&mvect_isotope_peaks[i], sizeof(IsotopePeak));
			count++;
		}
	}

	std::cout<<"Number of peaks written to file is " << count << " while total = " << mvect_isotope_peaks.size() << std::endl;

	ofs.close();

	//write out the umc classes
	ofstream ofs1("mvect_umcs", ios::binary);
	for ( int  i=0; i < mvect_umcs.size(); i++){
		ofs1.write( (char*)&mvect_umcs[i], sizeof(UMC));
	}

	ofs1.close();


	//write out the multimap that maps the indices on the isotopic peaks
	//to the indices on the umcs
	fstream fs("mmultimap_umc_2_peak_index", ios::out);
		
	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
	{
			
		int umc_index = (*iter).first ; 
		int peak_index = (*iter).second;
		fs << umc_index << "\t" << peak_index <<std::endl;
	}

	fs.close();
}


void UMCCreator::DeserializeObjects(){

}


float UMCCreator::GetLastMonoMassLoaded(){
	float mass = FLT_MAX;
	try{
		IsotopePeak peak = (IsotopePeak) mvect_isotope_peaks[mvect_isotope_peaks.size()];
		mass = peak.mdbl_mono_mass;
	}
	catch(exception e){
	}
	return mass;

}
void UMCCreator::SetInputFileName(char *filename){
	strcpy(mstr_inputFile, filename);
}

void UMCCreator::SetOutputDiretory(char * dir){
	strcpy(outputDir, dir);
}

//could all start going into a single branch but have kept it such so that it's a little easier to read.
//Method to determine whether the filter 
bool UMCCreator::ConsiderPeak(IsotopePeak peak ){

	bool consider = false;

	if ( peak.mdbl_abundance >= mint_min_intensity && peak.mflt_fit <= mflt_isotopic_fit_filter){
			consider = true;
	}
		//if we're processing a particular mass range then we have to check for the mono mass to be within bounds
	if ( consider && mflt_mono_mass_start <= peak.mdbl_mono_mass && peak.mdbl_mono_mass <= mflt_mono_mass_end ){
		consider = true;
	}
	else{
		consider = false;
	}
		
	
	//check if data is within LC scan range
	if ( consider && mint_lc_min_scan_filter <=  peak.mint_lc_scan && peak.mint_lc_scan <= mint_lc_max_scan_filter  ){
		consider = true;
	}
	else{
		consider = false;
	}

	if ( mbln_is_ims_data ){
		//check fi data is within ims scans
		if ( consider && mint_ims_min_scan_filter <= peak.mint_ims_scan && peak.mint_ims_scan <= mint_ims_max_scan_filter ){
			consider = true;
		}
		else{
			consider = false;
		}
	}

	return consider;
}

void UMCCreator::CalculateUMCs()
{
	mvect_umcs.clear() ; 
	mvect_umcs.reserve(mvect_umc_num_members.size()) ; 

	std::vector<double> vect_mass ; 

	UMC new_umc ; 
	mshort_percent_complete = 0 ; 
	int num_umcs = mvect_umcs.size() ; 

	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
	{
		vect_mass.clear() ; 
		int umc_index = (*iter).first ; 
		mshort_percent_complete = (short) ((100.0 * umc_index) / num_umcs) ; 
		int numMembers = mvect_umc_num_members[umc_index] ; 
		int minScan = INT_MAX ; 
		int maxScan = INT_MIN ; 
		double minMass = DBL_MAX ; 
		double maxMass = -1 * DBL_MAX ; 
		double maxAbundance = -1 * DBL_MAX ; 
		double sumAbundance = 0 ; 
		double sumMonoMass = 0 ; 
		int maxAbundanceScan = 0 ; 
		short classRepCharge = 0 ; 
		double classRepMz = 0 ; 

		while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == umc_index)
		{
			IsotopePeak pk = mvect_isotope_peaks[(*iter).second] ; 
			vect_mass.push_back(pk.mdbl_mono_mass) ; 

			if (pk.mint_lc_scan > maxScan)
				maxScan = pk.mint_lc_scan ; 
			if (pk.mint_lc_scan < minScan )
				minScan = pk.mint_lc_scan ; 

			if (pk.mdbl_mono_mass > maxMass)
				maxMass = pk.mdbl_mono_mass ; 
			if (pk.mdbl_mono_mass < minMass )
				minMass = pk.mdbl_mono_mass ; 

			if (pk.mdbl_abundance > maxAbundance)
			{
				maxAbundance = pk.mdbl_abundance ; 
				maxAbundanceScan = pk.mint_lc_scan ; 
				classRepCharge = pk.mshort_charge ; 
				classRepMz = pk.mdbl_mz ; 
			}
			sumAbundance += pk.mdbl_abundance ; 

			sumMonoMass += pk.mdbl_mono_mass ; 

			iter++ ; 
		}

		sort(vect_mass.begin(), vect_mass.end()) ; 

		new_umc.mint_umc_index = umc_index ; 
		new_umc.min_num_members = numMembers ; 
		new_umc.mint_start_scan = minScan ; 
		new_umc.mint_stop_scan = maxScan ; 
		new_umc.mint_max_abundance_scan = maxAbundanceScan ; 

		new_umc.mdbl_max_abundance = maxAbundance ; 
		new_umc.mdbl_sum_abundance = sumAbundance ;

		new_umc.mdbl_min_mono_mass = minMass ; 
		new_umc.mdbl_max_mono_mass = maxMass ; 
		new_umc.mdbl_average_mono_mass = sumMonoMass/numMembers ; 

		new_umc.mdbl_class_rep_mz = classRepMz ; 
		new_umc.mshort_class_rep_charge = classRepCharge ; 

		if (numMembers % 2 == 1)
		{
			new_umc.mdbl_median_mono_mass = vect_mass[numMembers/2] ; 
		}
		else
		{
			new_umc.mdbl_median_mono_mass = 0.5 * (vect_mass[numMembers/2-1] + vect_mass[numMembers/2]) ; 
		}
		mvect_umcs.push_back(new_umc) ; 
	}
}

void UMCCreator::RemoveShortUMCs(int min_length)
{
	std::vector<int> vectNewUMCLength ; 
	// first reset all isotope peak umc indices to -1. 
	int numIsotopePeaks = mvect_isotope_peaks.size() ;
	for (int peakNum = 0 ; peakNum < numIsotopePeaks ; peakNum++)
	{
		mvect_isotope_peaks[peakNum].mint_umc_index = -1 ; 
	}

	int numUmcsSoFar = 0 ; 

	mshort_percent_complete = 0 ; 

	int num_umcs = mmultimap_umc_2_peak_index.size() ; 
	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
	{
		int currentOldUmcNum = (*iter).first ; 
		mshort_percent_complete = (short) ((100.0 * currentOldUmcNum) / num_umcs) ; 

		int numMembers = mvect_umc_num_members[currentOldUmcNum] ; 
		while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == currentOldUmcNum)
		{
			if (numMembers >= min_length)
			{
				mvect_isotope_peaks[(*iter).second].mint_umc_index = numUmcsSoFar ; 
			}
			iter++ ; 
		}
		if (numMembers >= min_length)
		{
			mvect_umc_num_members[numUmcsSoFar] = numMembers ; 
			numUmcsSoFar++ ; 
		}
	}
	mvect_umc_num_members.resize(numUmcsSoFar) ; 
	// now set the map object. 
	mmultimap_umc_2_peak_index.clear() ; 
	for (int pkNum = 0 ; pkNum < numIsotopePeaks ; pkNum++)
	{
		IsotopePeak pk = mvect_isotope_peaks[pkNum] ; 
		if (pk.mint_umc_index != -1)
		{
			mmultimap_umc_2_peak_index.insert(std::pair<int,int>(pk.mint_umc_index, pkNum)) ; 
		}
	}
	// DONE!! 
}


bool UMCCreator::withinMassTolerance(double observedMass, double realMass){
	double massDifferenceInPPM = abs(realMass - observedMass)* 1000000/realMass;
	return ( massDifferenceInPPM <= mflt_constraint_mono_mass);
}



void UMCCreator::CreateUMCsSinglyLinkedWithAll()
{
	bool chargeStateMatch = true;
	mshort_percent_complete = 0 ; 
	mmultimap_umc_2_peak_index.clear() ; 
	int numPeaks = mvect_isotope_peaks.size() ; 
	for (int pkNum = 0 ; pkNum < numPeaks ; pkNum++)
	{
		mvect_isotope_peaks[pkNum].mint_umc_index = -1 ; 
	}

	std::vector<IsotopePeak> vectTempPeaks ; 

	vectTempPeaks.insert(vectTempPeaks.begin(), mvect_isotope_peaks.begin(), mvect_isotope_peaks.end()) ; 

	// basically take all umcs sorted in mass and perform single linkage clustering. 
	sort(vectTempPeaks.begin(), vectTempPeaks.end(), &SortIsotopesByMonoMassAndScan) ; 


	// now we are sorted. Start with the first index and move rightwards.
	// For each index, 
	// 1. If it already belongs to a UMC, move right comparing to all umcs in mass tolerance. For each match
	//		a. If match is not in a UMC, calculate distance. If within tolerance add it to current umc.
	//		b. If match is in the same UMC skip.
	//		c. If match is in a different UMC, calculate distance. If within tolerance, merge umcs, update current guys UMC number.
	// 2. If it does not belong to a UMC, move right comparing to all umcs in mass tolerance. For each match, calculate distance.
	//		a. Create new UMC and follow 1.

	int currentIndex = 0 ; 

	IsotopePeak currentPeak ; 
	IsotopePeak matchPeak ; 
	int numUmcsSoFar = 0 ; 
	double currentDistance = 0 ; 
	std::multimap<int, int>::iterator iter ; 
	std::multimap<int, int>::iterator deleteIter ; 
	std::vector<int> tempIndices ; // used to store indices of isotope peaks that are moved from one umc to another. 
	tempIndices.reserve(128) ; 


	mshort_percent_complete = 0 ; 
	while(currentIndex < numPeaks)
	{
		mshort_percent_complete = (short)((100.0 * currentIndex)/numPeaks) ; 
		currentPeak = vectTempPeaks[currentIndex] ; 
		if (currentPeak.mint_umc_index == -1)
		{
			// create UMC
			mmultimap_umc_2_peak_index.insert(std::pair<int,int>(numUmcsSoFar, currentIndex)) ; 
			currentPeak.mint_umc_index = numUmcsSoFar ; 
			vectTempPeaks[currentIndex].mint_umc_index = numUmcsSoFar ; 
			numUmcsSoFar++ ; 
		}
		int matchIndex = currentIndex + 1; 
		if (matchIndex == numPeaks)
			break ;

		double massTolerance = mflt_constraint_mono_mass ; 
		if (mbln_constraint_mono_mass_is_ppm)
			massTolerance *= currentPeak.mdbl_mono_mass / 1000000.0 ;		// Convert from ppm to Da tolerance

		double maxMass = currentPeak.mdbl_mono_mass + massTolerance ; 
		matchPeak = vectTempPeaks[matchIndex] ; 
		while (matchPeak.mdbl_mono_mass < maxMass)
		{
			if (matchPeak.mint_umc_index != currentPeak.mint_umc_index)
			{		
				currentDistance = PeakDistance(currentPeak, matchPeak) ; 
				if (mbln_constraint_charge_state){
						chargeStateMatch = (currentPeak.mshort_charge == matchPeak.mshort_charge);
				}
				if (currentDistance < mdbl_max_distance && chargeStateMatch)
				{
					if (matchPeak.mint_umc_index == -1)
					{
						mmultimap_umc_2_peak_index.insert(std::pair<int,int>(currentPeak.mint_umc_index, matchIndex)) ; 
						vectTempPeaks[matchIndex].mint_umc_index = currentPeak.mint_umc_index ; 
					}
					else
					{
						tempIndices.clear() ; 
						int numPeaksMerged = 0 ; 
						// merging time. Merge this guy's umc into the next guys UMC.
						for (iter = mmultimap_umc_2_peak_index.find(currentPeak.mint_umc_index) ; iter != mmultimap_umc_2_peak_index.end()
							&& (*iter).first == currentPeak.mint_umc_index; )
						{
							deleteIter = iter ; 
							int deletePeakIndex = (*iter).second ; 
							tempIndices.push_back(deletePeakIndex) ; 
							vectTempPeaks[deletePeakIndex].mint_umc_index = matchPeak.mint_umc_index ; 
							iter++ ; 
							mmultimap_umc_2_peak_index.erase(deleteIter) ; 
							numPeaksMerged++ ; 
						}
						for (int mergedPeakNum = 0 ; mergedPeakNum < numPeaksMerged ; mergedPeakNum++)
						{
							mmultimap_umc_2_peak_index.insert(std::pair<int,int>(matchPeak.mint_umc_index, tempIndices[mergedPeakNum])) ; 
						}
						currentPeak.mint_umc_index = matchPeak.mint_umc_index ; 
					}
				}
			}
			matchIndex++ ;
			if (matchIndex < numPeaks)
			{
				matchPeak = vectTempPeaks[matchIndex] ; 
			}
			else
				break ; 
		}
		currentIndex++ ; 
	}

	// At the end of all of this. The mapping from mmultimap_umc_2_peak_index is from umc_index to index in sorted stuff. 
	// Also, several of the umc indices are no longer valid. So lets step through the map, get new umc indices, renumber them,
	// and set the umc indices in the original vectors.
	numUmcsSoFar = 0 ; 
	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
	{
		int currentOldUmcNum = (*iter).first ; 
		int numMembers = 0 ; 
		while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == currentOldUmcNum)
		{
			IsotopePeak pk = vectTempPeaks[(*iter).second] ; 
			mvect_isotope_peaks[pk.mint_original_index].mint_umc_index = numUmcsSoFar ; 
			iter++ ; 
			numMembers++ ; 
		}
		mvect_umc_num_members.push_back(numMembers) ; 
		numUmcsSoFar++ ; 
	}
	// now set the map object. 
	mmultimap_umc_2_peak_index.clear() ; 
	for (int pkNum = 0 ; pkNum < numPeaks ; pkNum++)
	{
		IsotopePeak pk = mvect_isotope_peaks[pkNum] ; 
		mmultimap_umc_2_peak_index.insert(std::pair<int,int>(pk.mint_umc_index, pkNum)) ; 
	}
	// DONE!! 
}
	
void UMCCreator::SetPeks(std::vector<IsotopePeak> &vectPks)
{
	mvect_isotope_peaks.clear() ; 
	mvect_isotope_peaks.insert(mvect_isotope_peaks.begin(), vectPks.begin(), vectPks.end()) ; 

	int numPeaks = mvect_isotope_peaks.size() ; 
	for (int i = 0 ; i < numPeaks  ; i++)
	{
		IsotopePeak pk = mvect_isotope_peaks[i] ; 
		if (pk.mint_lc_scan > mint_lc_max_scan)
			mint_lc_max_scan = pk.mint_lc_scan ; 
		if (pk.mint_lc_scan < mint_lc_min_scan)
			mint_lc_min_scan = pk.mint_lc_scan ; 
	}
}
void UMCCreator::ReadPekFileMemoryMapped(char *fileName)
{
	Reset() ; 
	MemMappedReader mappedReader ; 
	mappedReader.Load(fileName) ; 
	__int64 file_len = mappedReader.FileLength() ; 

	char *fileNameTag = "Filename:" ; 
	int fileNameTagLength = (int)strlen(fileNameTag) ; 
	char *startTag = "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW" ; 
	char *startTagIsotopicLabeled = "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW,   Imono,   I+2" ; 
	int startTagLength = (int)strlen(startTag) ; 
	char *stopTag = "Processing stop time:" ; 
	int stopTagLen = (int)strlen(stopTag) ; 
	bool reading = false ; 
	bool isotopically_labeled = false ;
	bool is_first_scan = true ; 
	bool is_pek_file_from_wiff = false ; 

	IsotopePeak pk ; 

	short charge = 0 ; 
	double abundance = 0 ;
	float fit = 0, mz = 0 , averageMass = 0 , monoMass = 0 , maxMass = 0 ;  

	int numPeaks = 0 ; 

	mshort_percent_complete = 0 ;

	pk.mdbl_abundance = 0 ; 
	pk.mdbl_i2_abundance = 0 ; 
	pk.mdbl_average_mass = 0 ; 
	pk.mflt_fit = 0 ; 
	pk.mdbl_max_abundance_mass = 0 ; 
	pk.mdbl_mono_mass = 0 ; 
	pk.mdbl_mz = 0 ; 
	pk.mshort_charge = 0 ; 
	pk.mflt_ims_drift_time = 0 ;

	int mint_min_scan = INT_MAX ; 
	int mint_max_scan = 0 ; 

	int pos = 0 ; 
	const int MAX_HEADER_BUFFER_LEN = 1024 ; 
	char headerBuffer[MAX_HEADER_BUFFER_LEN] ;
	const int MAX_BUFFER_LEN = 1024 ; 
	char buffer[MAX_BUFFER_LEN] ;

	while (!mappedReader.eof())
	{
		mshort_percent_complete = (short)((100.0 * mappedReader.CurrentPosition()) / file_len) ; 
		if (mshort_percent_complete > 99)
			mshort_percent_complete = 99 ; 

		bool success = mappedReader.SkipToAfterLine(fileNameTag, headerBuffer, fileNameTagLength, MAX_HEADER_BUFFER_LEN) ; 
		if (!success)
		{
			break ; 
		}
		else
		{
			// found file name. start at the end. 
			int index = (int)strlen(headerBuffer) ; 
			while (index > 0 && headerBuffer[index] != '.')
			{
				index-- ; 
			}
			index++ ; 
			// wiff file pek files have a weird format. Another ICR2LS-ism.
			if (is_first_scan)
			{
				if(_strnicmp(&headerBuffer[index], "wiff", 4) ==0)
				{
					is_pek_file_from_wiff = true ; 
					// REMEMBER TO DELETE WHEN PARAMETERS ARE SET ELSEWHERE
					mflt_wt_mono_mass = 0.0025F ; // using ppms 10 ppm = length of 0.1 ppm
					mflt_wt_average_mass = 0.0025F ; // using ppms 10 ppm = length of 0.1 ppm
					mflt_wt_log_abundance = 0.1F ; 
					mflt_wt_scan = 0.01F ; 
					mflt_wt_fit = 0.1F ; 
					mflt_wt_ims_drift_time = 0.1F ;

					mflt_constraint_mono_mass = 25.0F ; // is in ppm
					mflt_constraint_average_mass = 25.0F ; // is in ppm. 

					mbln_use_net = true ;
					mbln_constraint_mono_mass_is_ppm = true ;
					mbln_constraint_average_mass_is_ppm = true ;

					mdbl_max_distance = 0.1 ; 

					mshort_percent_complete = 0 ;
				}
			}
			if (is_pek_file_from_wiff)
				index += 5 ; 
			pk.mint_lc_scan = atoi(&headerBuffer[index]) ; 
			if (pk.mint_lc_scan > mint_max_scan)
				mint_max_scan = pk.mint_lc_scan ; 
			if (pk.mint_lc_scan < mint_min_scan)
				mint_min_scan = pk.mint_lc_scan ; 
		}

		if (is_first_scan)
		{
			bool success = mappedReader.SkipToAfterLine(startTag, headerBuffer, startTagLength, MAX_HEADER_BUFFER_LEN) ; 
			if (!success)
				break ; 
			if (strncmp(headerBuffer, startTagIsotopicLabeled, strlen(startTagIsotopicLabeled)) == 0)
			{
				isotopically_labeled = true ; 
			}
			is_first_scan = false ; 
		}
		else
		{
			bool success = mappedReader.SkipToAfterLine(startTag, startTagLength) ; 
			if (!success)
				break ; 
		}
		while(!mappedReader.eof() && mappedReader.GetNextLine(buffer, MAX_BUFFER_LEN, stopTag, stopTagLen))
		{
			if (isotopically_labeled)
			{
				char *stopPtr ; 
				char *stopPtrNext ; 
				pk.mshort_charge = (short) strtol(buffer, &stopPtr, 10) ; 
				pk.mdbl_abundance = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_mz = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mflt_fit = (float)strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_average_mass = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_mono_mass = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_max_abundance_mass = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_mono_abundance = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_i2_abundance = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 

				pk.mflt_ims_drift_time = 0 ;
			}
			else
			{
				char *stopPtr ; 
				char *stopPtrNext ; 
				pk.mshort_charge = (short) strtol(buffer, &stopPtr, 10) ; 
				pk.mdbl_abundance = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_mz = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mflt_fit = (float)strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_average_mass = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_mono_mass = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 
				pk.mdbl_max_abundance_mass = strtod(stopPtr, &stopPtrNext) ; 
				stopPtr = stopPtrNext ; 

				pk.mflt_ims_drift_time = 0 ;
			}
			pk.mint_original_index = numPeaks ; 
			mvect_isotope_peaks.push_back(pk) ; 
			numPeaks++ ; 
		}
	}
	mappedReader.Close() ; 	
}


void UMCCreator::ReadPekFile(char *fileName)
{
	Reset() ; 
	FILE *fp = fopen(fileName, "r") ;
	int success = fseek(fp, 0, SEEK_END) ; 
	int file_len = ftell(fp) ; 
	success = fseek(fp, 0, SEEK_SET) ; 

	const int MAX_LINE_LENGTH = 512 ; 
	char buffer[MAX_LINE_LENGTH] ; 
	char *fileNameTag = "Filename:" ; 
	int fileNameTagLength = strlen(fileNameTag) ; 
	char *startTag = "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW" ; 
	char *startTagIsotopicLabeled = "CS,  Abundance,   m/z,   Fit,    Average MW, Monoisotopic MW,    Most abundant MW,   Imono,   I+2" ; 
	int startTagLength = strlen(startTag) ; 
	char *stopTag = "Processing stop time:" ; 
	bool reading = false ; 
	bool isotopically_labeled = false ;
	bool is_first_scan = true ; 

	IsotopePeak pk ; 
	int stopTagLen = (int)strlen(stopTag) ; 

	short charge = 0 ; 
	double abundance = 0 ;
	float fit = 0, mz = 0 , averageMass = 0 , monoMass = 0 , maxMass = 0 ;  

	int numPeaks = 0 ; 

	mshort_percent_complete = 0 ;

	pk.mdbl_abundance = 0 ; 
	pk.mdbl_i2_abundance = 0 ; 
	pk.mdbl_average_mass = 0 ; 
	pk.mflt_fit = 0 ; 
	pk.mdbl_max_abundance_mass = 0 ; 
	pk.mdbl_mono_mass = 0 ; 
	pk.mdbl_mz = 0 ; 
	pk.mshort_charge = 0 ; 
	pk.mflt_ims_drift_time = 0 ;

	mint_lc_min_scan = INT_MAX ; 
	mint_lc_max_scan = 0 ; 


	int pos = 0 ; 
	while (!feof(fp))
	{
		mshort_percent_complete = (short)((100.0 * pos) / file_len) ; 
		if (!reading)
		{
			fgets(buffer, MAX_LINE_LENGTH, fp) ; 
			pos += (int) strlen(buffer) ; 
			if (strncmp(buffer, startTag, startTagLength) == 0)
			{
				reading = true ; 
				if (is_first_scan)
				{
					if (strncmp(buffer, startTagIsotopicLabeled, strlen(startTagIsotopicLabeled)) == 0)
					{
						isotopically_labeled = true ; 
					}
					is_first_scan = false ; 
				}
			}
			else if (strncmp(buffer, fileNameTag, fileNameTagLength) == 0)
			{
				// found file name. start at the end. 
				int index = strlen(buffer) ; 
				while (index > 0 && buffer[index] != '.')
				{
					index-- ; 
				}
				index++ ; 
				pk.mint_lc_scan = atoi(&buffer[index]) ; 
				if (pk.mint_lc_scan > mint_lc_max_scan)
					mint_lc_max_scan = pk.mint_lc_scan ; 
				if (pk.mint_lc_scan < mint_lc_min_scan)
					mint_lc_min_scan = pk.mint_lc_scan ; 
			}
		}
		else
		{
			//int numRead = fscanf(fp, "%hd\t%lg\t%g\t%g\t%g\t%g\t%g", &charge, &abundance, &mz, &fit, 
			//	&averageMass, &monoMass, &maxMass) ; 
			int numRead = 0 ; 
			if (isotopically_labeled)
			{
				numRead = fscanf(fp, "%hd\t%lg\t%lg\t%g\t%lg\t%lg\t%lg\t%lg\t%lg", &pk.mshort_charge, &pk.mdbl_abundance, &pk.mdbl_mz, &pk.mflt_fit, 
					&pk.mdbl_average_mass, &pk.mdbl_mono_mass, &pk.mdbl_max_abundance_mass, &pk.mdbl_mono_abundance, &pk.mdbl_i2_abundance) ; 
				pos += 71 ; // this is really approximate, but its better than using ftellp
			}
			else
			{
				numRead = fscanf(fp, "%hd\t%lg\t%lg\t%g\t%lg\t%lg\t%lg", &pk.mshort_charge, &pk.mdbl_abundance, &pk.mdbl_mz, &pk.mflt_fit, 
					&pk.mdbl_average_mass, &pk.mdbl_mono_mass, &pk.mdbl_max_abundance_mass) ; 
				pos += 51 ; // this is really approximate, but its better than using ftellp
			}
			if (numRead == 0)
			{
				// at the end of the reading sections
				reading = false ; 
				pos += stopTagLen ; 
			}
			else
			{
				pk.mint_original_index = numPeaks ; 
				mvect_isotope_peaks.push_back(pk) ; 
				numPeaks++ ; 
			}
		}
	}
}

// Will map be affected by chunking?
bool UMCCreator::PrintMapping(FILE *stream, int featureStartIndex){

	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; ){
		int currentUmcNum = (*iter).first;
		while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == currentUmcNum)
		{
				IsotopePeak pk = mvect_isotope_peaks[(*iter).second] ; 
				
				fprintf(stream, "%d\t",currentUmcNum + featureStartIndex) ; 
				fprintf(stream, "%d\n",pk.mint_line_number_in_file);
				iter++ ;
		}
			 
	}

	return true;

}

bool UMCCreator::PrintMapping(FILE *stream){

	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; ){
		int currentUmcNum = (*iter).first;
		while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == currentUmcNum)
		{
				IsotopePeak pk = mvect_isotope_peaks[(*iter).second] ; 
				
				fprintf(stream, "%d\t",currentUmcNum) ; 
				fprintf(stream, "%d\n",pk.mint_line_number_in_file);
				iter++ ;
		}
			 
	}

	return true;

}

//method can be called with either stdout or an output file to write to 
bool UMCCreator::PrintUMCs(FILE *stream, bool print_members, int featureStartIndex){
	bool success = true;
	fprintf(stream, "Feature_index\tmonoisotopic_mass\tAverageMonoMass\tUMCMWMin\tUMCMWMax\tScanStart\tScanEnd\tScan\tUMCMemberCount\tMaxAbundance\tUMCAbundance") ; 
	if (print_members){
		fprintf(stream, "\tData") ; 
	}
	
	fprintf(stream, "\n") ; 

	
	int numPrinted = 1 ; 
	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
	{
		int currentUmcNum = (*iter).first; 
		UMC current_umc = mvect_umcs[currentUmcNum] ; 
		fprintf(stream, "%d\t", current_umc.mint_umc_index + featureStartIndex) ; 		
		fprintf(stream, "%4.4f\t", current_umc.mdbl_median_mono_mass);
		fprintf(stream, "%4.4f\t", current_umc.mdbl_average_mono_mass) ; 
		fprintf(stream, "%4.4f\t", current_umc.mdbl_min_mono_mass);
		fprintf(stream, "%4.4f\t", current_umc.mdbl_max_mono_mass);
		fprintf(stream, "%d\t", current_umc.mint_start_scan) ; 
		fprintf(stream, "%d\t", current_umc.mint_stop_scan);
		fprintf(stream, "%d\t",current_umc.mint_max_abundance_scan);
		fprintf(stream, "%d\t",current_umc.min_num_members) ; 
		fprintf(stream, "%4.4f\t", current_umc.mdbl_max_abundance);
		fprintf(stream, "%4.4f\t", current_umc.mdbl_sum_abundance) ; 

		
			while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == currentUmcNum)
			{
				if (print_members){	
			
				IsotopePeak pk = mvect_isotope_peaks[(*iter).second] ; 
				
				fprintf(stream, "%4.4f\t",pk.mdbl_mono_mass) ; 
				fprintf(stream, "%d\t",pk.mint_lc_scan);
				fprintf(stream, "%4.4\t", pk.mdbl_abundance) ; 
				}
				iter++ ; 
			}
		

		numPrinted++ ; 
		fprintf(stream, "\n") ; 
		fflush(stream);
	}

	fclose(stream);

	if ( numPrinted < 1 ){
		success = false;
	}

	return success;

}

//method can be called with either stdout or an output file to write to 
bool UMCCreator::PrintUMCs(FILE *stream, bool print_members){
	bool success = true;
	fprintf(stream, "Feature_index\tmonoisotopic_mass\tAverageMonoMass\tUMCMWMin\tUMCMWMax\tScanStart\tScanEnd\tScan\tUMCMemberCount\tMaxAbundance\tUMCAbundance") ; 
	if (print_members){
		fprintf(stream, "\tData") ; 
	}
	
	fprintf(stream, "\n") ; 

	
	int numPrinted = 1 ; 
	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
	{
		int currentUmcNum = (*iter).first; 
		UMC current_umc = mvect_umcs[currentUmcNum] ; 
		fprintf(stream, "%d\t", current_umc.mint_umc_index) ; 		
		fprintf(stream, "%4.4f\t", current_umc.mdbl_median_mono_mass);
		fprintf(stream, "%4.4f\t", current_umc.mdbl_average_mono_mass) ; 
		fprintf(stream, "%4.4f\t", current_umc.mdbl_min_mono_mass);
		fprintf(stream, "%4.4f\t", current_umc.mdbl_max_mono_mass);
		fprintf(stream, "%d\t", current_umc.mint_start_scan) ; 
		fprintf(stream, "%d\t", current_umc.mint_stop_scan);
		fprintf(stream, "%d\t",current_umc.mint_max_abundance_scan);
		fprintf(stream, "%d\t",current_umc.min_num_members) ; 
		fprintf(stream, "%4.4f\t", current_umc.mdbl_max_abundance);
		fprintf(stream, "%4.4f\t", current_umc.mdbl_sum_abundance) ; 

		
			while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == currentUmcNum)
			{
				if (print_members){	
			
				IsotopePeak pk = mvect_isotope_peaks[(*iter).second] ; 
				
				fprintf(stream, "%4.4f\t",pk.mdbl_mono_mass) ; 
				fprintf(stream, "%d\t",pk.mint_lc_scan);
				fprintf(stream, "%4.4\t", pk.mdbl_abundance) ; 
				}
				iter++ ; 
			}
		

		numPrinted++ ; 
		fprintf(stream, "\n") ; 
		fflush(stream);
	}

	fclose(stream);

	if ( numPrinted < 1 ){
		success = false;
	}

	return success;

}

void UMCCreator::PrintUMCs(bool print_members)
{
	std::cout<<"UMCIndex\tUMCMonoMW\tAverageMonoMass\tUMCMWMin\tUMCMWMax\tScanStart\tScanEnd\tScanClassRep\tUMCMemberCount\tMaxAbundance\tUMCAbundance" ; 
	if (print_members)
		std::cout<<"\tData" ; 
		std::cout<<"\n" ; 
	std::cout.precision(4);     
	std::cout.flags(std::ios::right + std::ios::fixed);
	int numPrinted = 1 ; 
	for (std::multimap<int,int>::iterator iter = mmultimap_umc_2_peak_index.begin() ; iter != mmultimap_umc_2_peak_index.end() ; )
	{
		int currentUmcNum = (*iter).first; 
		UMC current_umc = mvect_umcs[currentUmcNum] ; 
		std::cout.precision(0) ; 
		std::cout<<current_umc.mint_umc_index<<"\t" ; 		
		std::cout.precision(4) ; 
		std::cout<<current_umc.mdbl_median_mono_mass<<"\t"<<current_umc.mdbl_average_mono_mass<<"\t" ; 
		std::cout<<current_umc.mdbl_min_mono_mass<<"\t"<<current_umc.mdbl_max_mono_mass<<"\t"<<current_umc.mint_start_scan<<"\t" ; 
		std::cout.precision(0) ; 
		std::cout<<current_umc.mint_stop_scan<<"\t"<<current_umc.mint_max_abundance_scan<<"\t"<<current_umc.min_num_members<<"\t" ; 
		std::cout<<current_umc.mdbl_max_abundance<<"\t"<<current_umc.mdbl_sum_abundance ; 
		while(iter != mmultimap_umc_2_peak_index.end() && (*iter).first == currentUmcNum)
		{
			if (print_members)
			{
				IsotopePeak pk = mvect_isotope_peaks[(*iter).second] ; 
				std::cout.precision(4) ; 
				std::cout<<"\t"<<pk.mdbl_mono_mass<<"\t" ; 
				std::cout.precision(0) ; 
				std::cout<<pk.mint_lc_scan<<"\t"<<pk.mdbl_abundance<<"\t" ; 
			}
			iter++ ; 
		}
		numPrinted++ ; 
		std::cout<<"\n" ; 
	}
}

void UMCCreator::PrintPeaks()
{
	int numPeaks = mvect_isotope_peaks.size() ; 
	std::cout<<"Scan\tCharge\tAbundance\tMZ\tFit\tAverageMass\tMonoMass\tMaxMass\n" ; 

	std::cout.precision(4);     
	std::cout.flags(std::ios::right + std::ios::fixed);

	for (int i = 0 ; i < numPeaks  ; i++)
	{
		IsotopePeak pk = mvect_isotope_peaks[i] ; 
		std::cout.precision(0) ; 
		std::cout<<pk.mint_lc_scan<<"\t"<<pk.mshort_charge<<"\t"<<pk.mdbl_abundance<<"\t" ; 
		std::cout.precision(4) ; 
		std::cout<<pk.mdbl_mz<<"\t" ; 
		std::cout.precision(3) ; 
		std::cout<<pk.mflt_fit<<"\t" ; 
		std::cout.precision(4) ; 

		std::cout<<pk.mdbl_average_mass ; 
		std::cout<<"\t"<<pk.mdbl_mono_mass ; 
		std::cout<<"\t"<<pk.mdbl_max_abundance_mass<<"\n" ; 
	}
}

void UMCCreator::Reset()
{
	mvect_isotope_peaks.clear() ; 
	mvect_umcs.clear() ; 
	mvect_umc_num_members.clear() ;
	mmultimap_umc_2_peak_index.clear() ; 
	mshort_percent_complete = 0 ; 
}
