//This is the main DLL file.

#include "clsUMCCreator.h"
#include "IniReader.h"
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

	

	bool clsUMCCreator::LoadProgramOptions(){
		char settings_file[1024];
		char logText[1024];
		bool success = false;

		GetStr(mstr_options_name, settings_file);
		CIniReader iniReader(settings_file);

		//first load incoming and outgoing filenames and folder options
		char *isos_file = iniReader.ReadString("Files", "InputFileName", "");
		char *output_dir = iniReader.ReadString("Files", "OutputDirectory", ".");
		mobj_umc_creator->SetInputFileName(isos_file);
		mobj_umc_creator->SetOutputDiretory(output_dir);

		mstr_baseFileName = CreateBaseFileName(output_dir, isos_file);

		createLogFile();

		strcpy(logText, "Loading settings from INI file: ");
		strcat(logText, settings_file);
		log(logText);
		
		//next load data filters
		float isotopicFit = iniReader.ReadFloat("DataFilters", "MaxIsotopicFit", 1);
		if ( isotopicFit == 0 ){
			isotopicFit = 1;
		}

		int intensityFilter = iniReader.ReadInteger("DataFilters", "MinimumIntensity", 500);
		mflt_mono_mass_start = iniReader.ReadFloat("DataFilters", "MonoMassStart", 0);
		mflt_mono_mass_end = iniReader.ReadFloat("DataFilters", "MonoMassEnd", FLT_MAX);
		mbln_process_chunks = iniReader.ReadBoolean("DataFilters", "ProcessDataInChunks", false);

		//mint_mono_mass_overlap = iniReader.ReadInteger("DataFilters", "MonoMassSegmentOverlapDa", 2);
		if ( mflt_mono_mass_end == 0){
			if (mbln_process_chunks){
				mflt_mono_mass_end = mflt_mono_mass_start + 250;
			}
		}
		
		int maxPoints = iniReader.ReadInteger("DataFilters", "MaxDataPointsPerChunk", INT_MAX);
		if ( maxPoints == 0){
			maxPoints = INT_MAX;
		}

		int chunkSize = iniReader.ReadInteger("DataFilters", "ChunkSize", 3000);
		if (chunkSize == 0){
			chunkSize = 3000;
		}
		
		int imsMinScan = iniReader.ReadInteger("DataFilters", "IMSMinScan", 0);

		int imsMaxScan = iniReader.ReadInteger("DataFilters", "IMSMaxScan", 50000);
		if (imsMaxScan == 0 ){
			imsMaxScan = INT_MAX;
		}

		int lcMinScan = iniReader.ReadInteger("DataFilters", "LCMinScan", 0);
		int lcMaxScan = iniReader.ReadInteger("DataFilters", "LCMaxScan", 50000);
		if (lcMaxScan == 0){
			lcMaxScan = INT_MAX;
		}

		mobj_umc_creator->SetFilterOptions(isotopicFit, intensityFilter, lcMinScan, lcMaxScan, imsMinScan, imsMaxScan, mflt_mono_mass_start, mflt_mono_mass_end, mbln_process_chunks, maxPoints, mint_mono_mass_overlap, chunkSize);

		//next load the UMC creation options
		float monoMassWeight = iniReader.ReadFloat("UMCCreationOptions", "MonoMassWeight", 0.01);
		float monoMassConstraint = iniReader.ReadFloat("UMCCreationOptions", "MonoMassConstraint", 50);
		bool monoMassPPM = iniReader.ReadBoolean("UMCCreationOptions", "MonoMassConstraintIsPPM", true);
		float imsDriftWeight = iniReader.ReadFloat("UMCCreationOptions","IMSDriftTimeWeight", 0.1);
		float logAbundanceWeight = iniReader.ReadFloat("UMCCreationOptions", "LogAbundanceWeight", 0.1);
		float netWeight = iniReader.ReadFloat("UMCCreationOptions", "NETWeight", 0.01);
		float fitWeight = iniReader.ReadFloat("UMCCreationOptions", "FitWeight", 0.01);
		float avgMassWeight = iniReader.ReadFloat("UMCCreationOptions", "AvgMassWeight", 0.01);
		float avgMassConstr = iniReader.ReadFloat("UMCCreationOptions", "AvgMassConstraint", 10);
		bool avgMassPPM = iniReader.ReadBoolean("UMCCreationOptions", "AvgMassConstraintIsPPM", true);
		float scanWeight = iniReader.ReadFloat("UMCCreationOptions", "ScanWeight", 0);
		float maxDist = iniReader.ReadFloat("UMCCreationOptions", "MaxDistance", 0.1);
		bool useGeneric = iniReader.ReadBoolean("UMCCreationOptions", "UseGenericNET", true);
		mint_min_umc_length = iniReader.ReadInteger("UMCCreationOptions", "MinFeatureLengthPoints", 2);
		bool useCharge = iniReader.ReadBoolean("UMCCreationOptions", "UseCharge", false);

		//this one is not sent over for now
		bool useWeightedEuclidean = iniReader.ReadBoolean("UMCCreationOptions", "UseWeightedEuclidean", false);

		log("Data Filters - ");
		log(" Minimum LC scan = ", lcMinScan);
		log(" Maximum LC scan = ", lcMaxScan);
		log(" Minimum IMS scan = ", imsMinScan);
		log(" Maximum IMS scan = ", imsMaxScan);
		log(" Maximum fit = ", isotopicFit);
		log(" Minimum intensity = ", intensityFilter);
		log(" Mono mass start = ", mflt_mono_mass_start);
		log(" Mono mass end = ", mflt_mono_mass_end);
		log(" Require matching charge state = ", useCharge);

		//load all the umc creation options
		mobj_umc_creator->SetOptionsEx(monoMassWeight,monoMassConstraint, monoMassPPM, avgMassWeight,avgMassConstr, avgMassPPM, logAbundanceWeight, scanWeight, netWeight, fitWeight, maxDist, useGeneric, imsDriftWeight, useCharge);

		return success;
	}

	/**
	 * Calls methods for creating UMCs and writing them to files
	 */
	bool clsUMCCreator::PrintUMCsToFile(){
		bool success = false;

		char baseFileName[1024];

		GetStr(mstr_baseFileName, baseFileName);
		success = mobj_umc_creator->CreateFeatureFiles(baseFileName);

		return success;		
	}

	/**
	 * This overloaded version of the method is used to support chunking
	 *		- Chunking was only thought to be necessary because we ran out of memory on 32-bit machines.
	 *		- Running this program on a 64-bit machine with enough memory will allow it to process very large isos files
	 */
	bool clsUMCCreator::PrintUMCsToFile(int chunkIndex, int featureStartIndex){
		bool success = false;
		
		char baseFileName[1024];
		char chunk[1024];

		itoa(chunkIndex, chunk, 10);

		GetStr(mstr_baseFileName, baseFileName);
		strcat(baseFileName, "_chunk");
		strcat(baseFileName, chunk);

		success = mobj_umc_creator->CreateFeatureFiles(baseFileName, featureStartIndex);

		return success;		
	}

	/**
	Anuj added this method to load all the necessary details for the LC ms feature finder
	Program options, data filters and output directories are all loaded here
	The file is read and data filters applied and finally UMCs are found.
	*/
	void clsUMCCreator::LoadFindUMCs(){
		menm_status = LOADING;
		LoadProgramOptions();

		if ( mbln_process_chunks )
		{
			double mflt_mono_mass_chunk_start = mflt_mono_mass_start;
			float chunk_size = mobj_umc_creator->GetSegmentSize();
			
			bool objectsSerialized = false;
			log("Processing with Chunks ...");
			int iChunk = 0;
			int UMC_count = 0;

			while ( menm_status != COMPLETE ){
				//here's where you need to read a file process maxPoints and then continue from where you left off
				menm_status = CHUNKING;

				//here's where we have to figuure out how to laod the isos sqlite database file into the mvect_isotope_peaks
				//vector using custom reading code. Once we have that information, we can process the rest of the data as is
				//and write the UMC's out to the final output file. We can create the output file in APPEND mode if necessary

				menm_status = LOADING;
				if (mflt_mono_mass_chunk_start + chunk_size > mflt_mono_mass_end)
				{
					menm_status= COMPLETE;
					continue;
				}

				//set the range of mass in which we'll operate
				mobj_umc_creator->SetMassRange(mflt_mono_mass_chunk_start, mflt_mono_mass_chunk_start + chunk_size);
				
				//instead of READCSV, here's where we need to read from the sqLite database and get a list of peaks
				//also you need to do this read in such a manner that you get close to the maxPoints per segment range
				//int numPeaks = mobj_umc_creator->LoadPeaksFromDatabase();
				
				mflt_mono_mass_chunk_start = mflt_mono_mass_chunk_start + chunk_size;
				log("Processing one Chunk");
				int numPeaks = mobj_umc_creator->ReadCSVFile();

				log("Total number of peaks we'll consider = ", numPeaks); 
				if (numPeaks == 0){
					menm_status= COMPLETE;
					log("Nothing read in");
					continue;
				}

				//if ( false ){
				//since mbln_process is true, the file reading will stop when it has loaded maxPoints of peaks that pass filters
				//so the mono mass of the last peak loaded will be the end mass 
				menm_status = CLUSTERING;
				mobj_umc_creator->CreateUMCsSinglyLinkedWithAll();
				menm_status = SUMMARIZING;
				
				mstr_message = new System::String("Filtering out short clusters") ; 
				mobj_umc_creator->RemoveShortUMCs(mint_min_umc_length) ;
				mstr_message = new System::String("Calculating UMC statistics") ; 
				mobj_umc_creator->CalculateUMCs() ;

				//now here's where we'll have to perform some magic of writing to temporary files and reloading those files in 
				//conjunction with loading data from the original isos files. The isos file data can be directly loaded using ReadCSV
				//but it'll have to stack the peaks on top of existing UMC peaks. Maybe we could just serialize all objects
				//to a temporary file and reload from there? Something to consider
				// mobj_umc_creator->SerializeObjects();
				// objectsSerialized = true;
				
				PrintUMCsToFile(iChunk, UMC_count);
				UMC_count += mobj_umc_creator->GetNumUmcs();
				//}
				iChunk++;
			}

			// Combine several partial UMC files into one file.

		} 
		else 
		{
			log("Processing without Chunks...");
			menm_status = LOADING;
			int numPeaks = mobj_umc_creator->ReadCSVFile();
			log("Total number of peaks we'll consider = ", numPeaks);

			menm_status = CLUSTERING;
			log("Creating UMCs...");
			mobj_umc_creator->CreateUMCsSinglyLinkedWithAll();

			menm_status = SUMMARIZING;
			log("Filtering out short UMCs...");
			mobj_umc_creator->RemoveShortUMCs(mint_min_umc_length);

			log("Calculating UMC statistics...");
			mobj_umc_creator->CalculateUMCs();
			menm_status = COMPLETE;

			log("Total number of UMCs = ", mobj_umc_creator->GetNumUmcs());

			log("Writing output files...");
			PrintUMCsToFile();
		}

		fclose(mfile_logFile);
	
	}

	void clsUMCCreator::FindUMCs()
	{
		menm_status = CLUSTERING ; 
		mstr_message = new System::String("Clustering Isotope Peaks") ; 
		mobj_umc_creator->CreateUMCsSinglyLinkedWithAll(); 
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
			Console::WriteLine(mstr_message);
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
            if (mobj_umc_creator->mint_lc_max_scan > mobj_umc_creator->mint_lc_min_scan) {
			     // Compute Generic NET value
			     newUmc->mdbl_net = (double) (umc.mint_max_abundance_scan - mobj_umc_creator->mint_lc_min_scan) * 1.0 / (mobj_umc_creator->mint_lc_max_scan - mobj_umc_creator->mint_lc_min_scan) ; 
            }

			newUmc->mint_scan_aligned = umc.mint_max_abundance_scan ; 

			newUmc->mint_umc_index = umc.mint_umc_index ; 
			newUmc->mint_class_rep_charge = (int) umc.mshort_class_rep_charge ; 
			arr_umcs[umcNum] = newUmc ; 
		}

		return arr_umcs ; 
	}

	void clsUMCCreator::SetLCMinMaxScans(int min, int max)
	{
		mobj_umc_creator->SetLCMinMaxScan(min, max) ; 

		//if (max <= min)
			//cout << " Max scan must be greater than min scan" << endl;
			//throw new exception("Max scan must be greater than min scan");
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
			pk.mint_lc_scan = isoPk->mint_lc_scan ; 
			pk.mshort_charge = isoPk->mshort_charge ; 
			pk.mdbl_abundance = isoPk->mdbl_abundance ; 
			pk.mdbl_mz = isoPk->mdbl_mz ; 
			pk.mflt_fit = isoPk->mflt_fit ; 
			pk.mflt_ims_drift_time = isoPk->mflt_ims_drift_time ;
			pk.mdbl_average_mass = isoPk->mdbl_average_mass ; 
			pk.mdbl_mono_mass = isoPk->mdbl_mono_mass ; 
			pk.mdbl_max_abundance_mass = isoPk->mdbl_max_abundance_mass ;
			pk.mdbl_i2_abundance = isoPk->mdbl_i2_abundance ; 

			vectPeaks.push_back(pk) ; 
		}
		mobj_umc_creator->SetPeks(vectPeaks) ; 
	}

	void clsUMCCreator::SetFilterOptions(float isotopic_fit, int min_intensity, int min_lc_scan, int max_lc_scan, int min_ims_scan, int max_ims_scan, float mono_mass_start, float mono_mass_end, bool process_mass_seg, int maxDataPoints, int monoMassSegOverlap, float segmentSize){
		mobj_umc_creator->SetFilterOptions(isotopic_fit, min_intensity, min_lc_scan, max_lc_scan, min_ims_scan, max_ims_scan, mono_mass_start, mono_mass_end, process_mass_seg, maxDataPoints, monoMassSegOverlap, segmentSize);
	}
	void clsUMCCreator::SetOptions(float wt_mono_mass, float wt_avg_mass, float wt_log_abundance, float wt_scan, float wt_fit,
			float wt_net, float mono_constraint, float avg_constraint, double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs)
	{
		mobj_umc_creator->SetOptions(wt_mono_mass, wt_avg_mass, wt_log_abundance, wt_scan, wt_fit, wt_net, 
			mono_constraint, avg_constraint, max_dist, use_net, wt_ims_drift_time, use_cs) ; 
	}

	void clsUMCCreator::SetOptionsEx(
					float wt_mono_mass, float mono_constraint, bool mono_constraint_is_ppm,
					float wt_avg_mass, float avg_constraint, bool avg_constraint_is_ppm,
					float wt_log_abundance, float wt_scan, float wt_net, float wt_fit,
					double max_dist, bool use_net, float wt_ims_drift_time, bool use_cs)
	{
		mobj_umc_creator->SetOptionsEx(
					wt_mono_mass, mono_constraint, mono_constraint_is_ppm,
					wt_avg_mass, avg_constraint, avg_constraint_is_ppm,
					wt_log_abundance, wt_scan, wt_net, wt_fit, 
					max_dist, use_net, wt_ims_drift_time, use_cs) ;
	}

	/*
	 * Creates the baseFileName that will be used for properly naming the output files.
	 *		- Desired Output format: baseFileName_LCMSFeatures.txt and baseFileName_LCMSFeatureToPeakMap.txt
	 */
	char* clsUMCCreator::CreateBaseFileName(char* directoryName, char* inputFileName){

		char baseFileName[1024];
		char* inputFileNameTok;
		char* inputFileNameTokHelper;

		// Remove any directory information from the filename
		inputFileNameTok = strtok(inputFileName, "\\/");
		inputFileNameTokHelper = inputFileNameTok;
		while(inputFileNameTokHelper != NULL){
			inputFileNameTokHelper = strtok(NULL, "\\/");
			if(inputFileNameTokHelper != NULL){
				inputFileNameTok = inputFileNameTokHelper;
			}
		}

		// Remove the "_isos.csv" file extension
		char* fileExtension = strstr(inputFileNameTok, "_isos.csv");
		strcpy(fileExtension, "");
		
		// Append the Input File Name to the Output Directory to create the Base File Name
		strcpy(baseFileName, directoryName);
		strcat(baseFileName, "\\");
		strcat(baseFileName, inputFileNameTok);

		return baseFileName;

	}

	void clsUMCCreator::createLogFile(){
		char logFileName[1024];

		GetStr(mstr_baseFileName, logFileName);
		strcat(logFileName, "_FeatureFinder_Log.txt");

		mfile_logFile = fopen(logFileName, "w");
	}

	void clsUMCCreator::log(char* textToLog){
		time_t now = time(NULL);
		struct tm *localTime = localtime(&now);

		fprintf(mfile_logFile, "%.2d/%.2d/%.2d %.2d:%.2d:%.2d\t%s\n", localTime->tm_mon+1, localTime->tm_mday, localTime->tm_year+1900, localTime->tm_hour, localTime->tm_min, localTime->tm_sec, textToLog);
		fflush(mfile_logFile);
		Console::WriteLine(textToLog);
	}

	void clsUMCCreator::log(char* textToLog, int numToLog){
		time_t now = time(NULL);
		struct tm *localTime = localtime(&now);

		fprintf(mfile_logFile, "%.2d/%.2d/%.2d %.2d:%.2d:%.2d\t%s%d\n", localTime->tm_mon+1, localTime->tm_mday, localTime->tm_year+1900, localTime->tm_hour, localTime->tm_min, localTime->tm_sec, textToLog, numToLog);
		fflush(mfile_logFile);
		Console::Write(textToLog);
		Console::WriteLine(numToLog);
	}

	void clsUMCCreator::log(char* textToLog, float numToLog){
		time_t now = time(NULL);
		struct tm *localTime = localtime(&now);

		fprintf(mfile_logFile, "%.2d/%.2d/%.2d %.2d:%.2d:%.2d\t%s%4.4f\n", localTime->tm_mon+1, localTime->tm_mday, localTime->tm_year+1900, localTime->tm_hour, localTime->tm_min, localTime->tm_sec, textToLog, numToLog);
		fflush(mfile_logFile);
		Console::Write(textToLog);
		Console::WriteLine(numToLog);
	}

}
