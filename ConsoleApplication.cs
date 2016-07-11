using System;
using System.IO;
using UMCCreation;

namespace FeatureFinder
{
	class FeatureFinderConsoleApplication
	{
		IniFile mIniFile_ControlFile;
		string mstr_controlFile;

		//input output file names
		string mstr_InputFileName;
		DirectoryInfo mDirInfo_OutputDir;

		//data filters for incoming peaks
		float mflt_isotopicFit = 0.15F;
		int mint_minimumIntensity = 500;
		float mflt_monoMassStart = Single.MinValue;
		float mflt_monoMassEnd = Single.MaxValue;
		bool mbln_processDataInSegments = false;
		int mint_pointsToLoadPerSegment = Int32.MaxValue;
		byte mbyt_monoMassSegmentOverlapDa = 1;
		int mint_min_lc_scan = 0;
		int mint_max_lc_scan = Int32.MaxValue;

		int mint_min_ims_scan = 0;
		int mint_max_ims_scan = Int32.MaxValue;


		//options for UMC creation dll
		float mflt_avgMassContraint = 10;
		float mflt_monoMassConstraint = 50;
		
		bool mbln_monoMassConstraintIsPPM = false;
		bool mbln_avgMassConstraintIsPPM = false;
		bool mbln_useGenericNet = false;


		//These are totally arbitrarily determined for now and that's not the correct way of doing it,
		//These should be defined on a 0-1 or a percentage scale for better operations
		float mflt_imsDriftTimeWeight = 0.01f;
		float mflt_logAbundanceWeight = 0.01f;
		float mflt_monoMassWeight = 0.01f;
		float mflt_netWeight = 15;
		float mflt_fitWeight = 0.1f;
		float mflt_avgMassWeight = 0.1f;
		float mflt_scanWeight = 0.1f;


		//This is the distance between two points in Euclidean distance space to be considered the same
		//after accounting for all factors. Again an abitrarily chosen value. There has to be some method 
		//to this madness.
		float mflt_MaxDistance = 0.1f;

		
		ushort mushrt_MinScan = 0;
		ushort mushrt_MaxScan = 0;
		byte mbyt_MinFeatureLengthPoints = 0;


		/**
		 * This method is used to check if the input file exists and whether we have access to create the output
		 * directories. It returns a boolean value indicating whether we can proceed with the rest of the program.
		 */
		private bool loadFiles()
		{
			bool success = false;
			try
			{
				mstr_InputFileName = mIniFile_ControlFile.GetStringValue("InputFileName");
				
				success = File.Exists(mstr_InputFileName);
				if (!success)
				{
					Console.WriteLine("Input file does not exist. Please make sure you've provided the correct path");
				}
			
				string outputName = mIniFile_ControlFile.GetStringValue("OutputDirectory");
				if ( !Directory.Exists(outputName))
				{
					Console.WriteLine("Output directory not present. Creating " + outputName);
					mDirInfo_OutputDir = Directory.CreateDirectory(outputName);
					success = mDirInfo_OutputDir.Exists;
				}
				
			}
			catch(Exception e )
			{
				Console.WriteLine("Exception caused when loading input files ");
				Console.WriteLine(e.Message);
			}
			return success;

		}


		/**
		 * This method is used to load in all data filters and assign values depending on control file. If values are
		 * not present then we leave it at reasonable defaults.
		 * Returns a bool indicating whether we can proceed with the rest of the program
		 */
		private bool loadDataFilters()
		{
			bool success = false;
			try
			{
				mflt_isotopicFit = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("MaxIsotopicFit"));
				mint_minimumIntensity = Convert.ToInt32(mIniFile_ControlFile.GetStringValue("MinimumIntensity"));

				mbln_processDataInSegments = Convert.ToBoolean(mIniFile_ControlFile.GetStringValue("ProcessDataInChunks"));
				if (mbln_processDataInSegments)
				{
					mint_pointsToLoadPerSegment = Convert.ToInt32(mIniFile_ControlFile.GetStringValue("MaxDataPointsPerChunk"));
					mbyt_monoMassSegmentOverlapDa = Convert.ToByte(mIniFile_ControlFile.GetStringValue("MonoMassSegmentOverlapDa"));
				}

				mflt_monoMassStart = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("MonoMassStart"));
				mflt_monoMassEnd = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("MonoMassEnd"));
				if (mflt_monoMassStart == mflt_monoMassEnd)
				{
					mflt_monoMassStart = Single.MinValue;
					mflt_monoMassEnd = Single.MaxValue;
				}
				
				mint_min_lc_scan = Convert.ToInt32(mIniFile_ControlFile.GetStringValue("LCMinScan"));
				mint_max_lc_scan = Convert.ToInt32(mIniFile_ControlFile.GetStringValue("LCMaxScan"));

				if ( mint_min_lc_scan < 0 )
				{
					mint_min_lc_scan = 0;
				}
				if (mint_max_lc_scan == 0)
				{
					mint_max_lc_scan = Int32.MaxValue;
				}
				
				mint_min_ims_scan = Convert.ToInt32(mIniFile_ControlFile.GetStringValue("IMSMinScan"));
				mint_max_ims_scan = Convert.ToInt32(mIniFile_ControlFile.GetStringValue("IMSMaxScan"));

				if (mint_max_ims_scan == 0)
				{
					mint_max_ims_scan = Int32.MaxValue;
				}

				success = true;
			}
			catch(Exception e)
			{
				Console.WriteLine("Exception while reading data filtering parameters. Check control file.\n");
				Console.WriteLine(e.Message);
			}

			return success;

		}

		private bool loadUMCCreationOptions()
		{
			bool success = false;
			try
			{
				mflt_avgMassContraint = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("AvgMassConstraint"));
				mflt_monoMassConstraint = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("MonoMassConstraint"));
				
				if ( mIniFile_ControlFile.GetStringValue("MonoMassConstraintIsPPM").Equals("True"))
				{
					mbln_monoMassConstraintIsPPM = true;
				}
				else
				{
					mbln_monoMassConstraintIsPPM = false;
				}

				if (mIniFile_ControlFile.GetStringValue("AvgMassConstraintIsPPM").Equals("True"))
				{
					mbln_avgMassConstraintIsPPM = true;
				}
				else
				{
					mbln_avgMassConstraintIsPPM=false;
				}

				if (mIniFile_ControlFile.GetStringValue("UseGenericNet").Equals("True"))
				{
					mbln_useGenericNet = true;
				}
				else
				{
					mbln_useGenericNet = false;
				}

				//These are totally arbitrarily determined for now and that's not the correct way of doing it,
				//These should be defined on a 0-1 or a percentage scale for better operations
				mflt_imsDriftTimeWeight = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("IMSDriftTimeWeight"));
				mflt_logAbundanceWeight = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("LogAbundanceWeight"));
				mflt_monoMassWeight = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("MonoMassWeight"));
				mflt_netWeight = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("NETWeight"));
				mflt_fitWeight = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("FitWeight"));
				mflt_avgMassWeight = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("AvgMassWeight"));
				mflt_scanWeight = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("ScanWeight"));

				Console.WriteLine("Finished loading all weights");

				//This is the distance between two points in Euclidean distance space to be considered the same
				//after accounting for all factors. Again an abitrarily chosen value. There has to be some method 
				//to this madness.
				mflt_MaxDistance = Convert.ToSingle(mIniFile_ControlFile.GetStringValue("MaxDistance"));
		
				mushrt_MinScan = Convert.ToUInt16(mIniFile_ControlFile.GetStringValue("IMSMinScan"));
				if (mushrt_MinScan <= 0)
				{
					mushrt_MinScan = UInt16.MinValue;
				}

				mushrt_MaxScan = Convert.ToUInt16(mIniFile_ControlFile.GetStringValue("IMSMaxScan"));
				if (mushrt_MaxScan <= 0)
				{
					mushrt_MaxScan = UInt16.MaxValue;
				}
				
				mbyt_MinFeatureLengthPoints = Convert.ToByte(mIniFile_ControlFile.GetStringValue("MinFeatureLengthPoints"));
				success = true;
				
			}
			catch(Exception e)
			{
				Console.WriteLine("Exception while reading UMC creation parameters. Check control file.\n");
				Console.WriteLine(e.StackTrace);
			}
			return success;
		}

		private string ProcessFilename(string filename)
		{
			// Replace all slashes to backslashes since we are working with a Windows directory 
			filename = filename.Replace("/", "\\");

			// If the string does not contain ":\" or "\\", move on.
			if (!filename.Contains(":\\") && !filename.StartsWith("\\\\"))
			{
				// Append "." to the front of the string if in the form of "\blabla"
				if (filename.StartsWith("\\"))
				{
					return "." + filename;
				}
				// Append ".\" to the front of the string if in the form of "blabla"
				else
				{
					return ".\\" + filename;
				}
			}

			// filename is in the form of "C:\blabla" or "\\blabla"
			return filename;
		}


		public FeatureFinderConsoleApplication(string filename)
		{
			string newFileName = ProcessFilename(filename);

			mstr_controlFile = newFileName;
			mIniFile_ControlFile = new IniFile(mstr_controlFile);
		}


		public bool processControlFile()
		{
			bool success = false;
			try
			{
				mIniFile_ControlFile = new IniFile(mstr_controlFile);
				Console.WriteLine("Reading file " + mstr_controlFile + " ...");

				success = loadFiles();

				//move on only if successful on loading input files
				if (success)
				{
					Console.WriteLine("Now loading filter options");
					success = loadDataFilters();
				}

				//move on and load the UMCCreation options 
				if (success)
				{
					Console.WriteLine("Now loading umc creation options");
					success = loadUMCCreationOptions();
				}
			}
			catch(Exception e)
			{
				Console.WriteLine(e.Message);
				PrintUsage();
			}

			return success;
		}

		public void invokeUMCCreator()
		{
#if !(DEBUG)
			try
			{
#endif
				UMCCreation.clsUMCCreator c = new clsUMCCreator();
				c.OptionsFileName = mstr_controlFile;

                // This is where UMC creator code is called
				c.LoadFindUMCs();
#if !(DEBUG)
			}
			catch(Exception c)
			{
				Console.WriteLine(c.Message);
				Console.WriteLine(c.StackTrace);
			}
#endif
		}
	
		/// <summary>
		/// The main entry point for the application.
		/// </summary>
		[STAThread]
		static void Main(string[] args)
		{
			if ( args.Length < 1)
			{
				PrintUsage();
				Console.ReadLine();
			}
			else
			{
				FeatureFinderConsoleApplication cs = new FeatureFinderConsoleApplication(args[0]);
				//cs.processControlFile();
				if (cs.loadFiles())
				{
					cs.invokeUMCCreator();
				}
			}
		}


		static void PrintUsage()
		{
			Console.WriteLine("***************************************************");
			Console.WriteLine("FeatureFinder usage : FeatureFinder controlfile ");
			Console.WriteLine("Press enter to terminate program");
			Console.WriteLine("***************************************************");
			Console.WriteLine("***************************************************");
		}

	}
}
