using System;
using System.Collections ; 

namespace UMCCreation
{

	/// <summary>
	/// Summary description for Class1.
	/// </summary>
	public class clsUMCData
	{
		// the umcs are sorted in terms of the datasets and 
		// indices are kept for start and end of each dataset.
	    private Hashtable mhash_umc_2_dataset_num ;
	    private ArrayList marr_file_names ;
	    private ArrayList marr_dataset_start_index ;
	    private ArrayList marr_dataset_stop_index ;
	    private ArrayList marr_dataset_min_scan ;
	    private ArrayList marr_dataset_max_scan ;

	    private int mint_num_datasets ;
	    private int mint_num_umcs ;
	    private const int mint_default_umc_allocation = 100000 ; 
		
		public clsUMC []  marr_umcs ; 
		public clsUMCData()
		{
			//
			// TODO: Add constructor logic here
			//
			Allocate() ; 
		}

		public void Clear()
		{
			Allocate() ; 
		}

		public int NumUMCS => mint_num_umcs;

	    public void Allocate() 
		{
			mhash_umc_2_dataset_num = new Hashtable() ;
			marr_file_names = new ArrayList() ; 
			mint_num_datasets = 0 ; 
			mint_num_umcs = 0 ; 
			marr_umcs = new clsUMC[mint_default_umc_allocation] ; 
			marr_dataset_start_index = new ArrayList() ; 
			marr_dataset_stop_index = new ArrayList() ; 
			marr_dataset_min_scan = new ArrayList() ; 
			marr_dataset_max_scan = new ArrayList() ; 
		}

		public string [] DatasetName
		{
			get
			{
				var datasets = new string [mint_num_datasets] ; 
				marr_file_names.CopyTo(datasets) ; 
				return datasets ; 
			}
		}

		public string DatasetNameAt(int dataset_num)
		{
			if (dataset_num >= marr_file_names.Count || dataset_num < 0)
				return null ; 
			return (string) marr_file_names[dataset_num] ; 
		}

		public clsUMC[] GetUMCS(int dataset_num)
		{
			int start_index = 0 , stop_index = 0 ; 
			GetDataIndex(dataset_num, ref start_index, ref stop_index) ; 
			var num_pts = stop_index - start_index ;
			var umc_arr = new clsUMC[num_pts] ;
		    try
			{
			    int index ;
			    for (index = start_index ; index < stop_index ; index++)
				{
					umc_arr[index-start_index] = marr_umcs[index] ; 
				}
			}
			catch (Exception e)
			{
				Console.WriteLine(e.Message + e.StackTrace) ; 
			}
			return umc_arr ; 
		}

		public int SetUMCS(string file_path, clsUMC [] umcs, int min_scan, int max_scan)
		{
			try
			{
				var index = file_path.LastIndexOf("\\", StringComparison.Ordinal) ;
				var file = file_path.Substring(index+1) ; 

				marr_file_names.Add(file) ; 
				mhash_umc_2_dataset_num [file] = mint_num_datasets ; 
				marr_dataset_start_index.Add(mint_num_umcs) ; 

				if (marr_umcs.Length < mint_num_umcs + umcs.Length)
				{
					var temp = new clsUMC[2*marr_umcs.Length]; 
					for (var i = 0 ; i < mint_num_umcs ; i++)
					{
						temp[i] = marr_umcs[i] ; 
					}
					marr_umcs = temp ; 
				}
				var num_to_add = umcs.Length ; 

				for (var i = 0 ; i < num_to_add ; i++)
				{
					marr_umcs[i+mint_num_umcs] = umcs[i] ; 
				}
				mint_num_umcs += num_to_add ; 
				marr_dataset_stop_index.Add(mint_num_umcs) ; 
				marr_dataset_min_scan.Add(min_scan) ; 
				marr_dataset_max_scan.Add(max_scan) ; 
			}
			catch (Exception ex)
			{
				Console.WriteLine(ex.Message + ex.StackTrace) ; 
			}

			mint_num_datasets++ ; 
			return mint_num_datasets - 1; 
		}

		public int SetUMCS(string file_path, clsUMC [] umcs)
		{
			try
			{
				var index = file_path.LastIndexOf("\\", StringComparison.Ordinal) ;
				var file = file_path.Substring(index+1) ; 

				marr_file_names.Add(file) ; 
				mhash_umc_2_dataset_num [file] = mint_num_datasets ; 
				marr_dataset_start_index.Add(mint_num_umcs) ; 

				if (marr_umcs.Length < mint_num_umcs + umcs.Length)
				{
					var temp = new clsUMC[2*marr_umcs.Length]; 
					for (var i = 0 ; i < mint_num_umcs ; i++)
					{
						temp[i] = marr_umcs[i] ; 
					}
					marr_umcs = temp ; 
				}
				var num_to_add = umcs.Length ; 
				var min_scan = int.MaxValue ; 
				var max_scan = int.MinValue ; 

				for (var i = 0 ; i < num_to_add ; i++)
				{
					marr_umcs[i+mint_num_umcs] = umcs[i] ; 
					if (umcs[i].mint_scan < min_scan)
						min_scan = umcs[i].mint_scan ; 
					if (umcs[i].mint_scan > max_scan)
						max_scan = umcs[i].mint_scan ; 
				}
				mint_num_umcs += num_to_add ; 
				marr_dataset_stop_index.Add(mint_num_umcs) ; 
				marr_dataset_min_scan.Add(min_scan) ; 
				marr_dataset_max_scan.Add(max_scan) ; 
			}
			catch (Exception ex)
			{
				Console.WriteLine(ex.Message + ex.StackTrace) ; 
			}

			mint_num_datasets++ ; 
			return mint_num_datasets - 1; 
		}

		public void GetMassesAndScans(int dataset_num, bool aligned, ref float [] masses, ref float [] scans)
		{
			var start_index = (int) marr_dataset_start_index[dataset_num] ; 
			var stop_index =(int) marr_dataset_stop_index[dataset_num] ; 
			var num_pts = stop_index - start_index ;
		    try
			{

				masses = new float[num_pts] ; 
				scans = new float[num_pts] ;
			    int index ;
			    for (index = start_index ; index < stop_index ; index++)
				{
					if (!aligned)
					{
						masses[index-start_index] = Convert.ToSingle(marr_umcs[index].mdbl_mono_mass) ; 
						scans[index-start_index] =  Convert.ToSingle(marr_umcs[index].mint_scan) ; 
					}
					else
					{
						masses[index-start_index] = Convert.ToSingle(marr_umcs[index].mdbl_mono_mass_calibrated) ; 
						scans[index-start_index] =  Convert.ToSingle(marr_umcs[index].mint_scan_aligned) ; 
					}
				}
			}
			catch (Exception e)
			{
				Console.WriteLine(e.Message + e.StackTrace) ; 
			}
		}

		// returns the masses and scans of the features in the dataset. Also returns the starting 
		// index for the umcs in the array marr_umcs.
		public int GetMassesAndScans(int dataset_num, ref double [] masses, ref int [] scans)
		{
			var start_index = (int) marr_dataset_start_index[dataset_num] ; 
			var stop_index =(int) marr_dataset_stop_index[dataset_num] ; 
			var num_pts = stop_index - start_index ;
		    try
			{

				masses = new double[num_pts] ; 
				scans = new int[num_pts] ;
			    int index ;
			    for (index = start_index ; index < stop_index ; index++)
				{
					masses[index-start_index] = marr_umcs[index].mdbl_mono_mass ; 
					scans[index-start_index] =  marr_umcs[index].mint_scan ; 
				}
			}
			catch (Exception e)
			{
				Console.WriteLine(e.Message + e.StackTrace) ; 
			}
			return start_index ; 
		}

		public void GetDataIndex(int dataset_num, ref int start_index, ref int stop_index)
		{
			start_index = (int) marr_dataset_start_index[dataset_num] ; 
			stop_index =(int) marr_dataset_stop_index[dataset_num] ; 
			return ; 
		}
		
		public void GetMinMaxScan(int dataset_num, ref int min_scan, ref int max_scan)
		{
			try
			{
				min_scan = (int) marr_dataset_min_scan[dataset_num] ; 
				max_scan = (int) marr_dataset_max_scan[dataset_num] ; 
			}
			catch (Exception e)
			{
				Console.WriteLine(e.Message + e.StackTrace) ; 
			}
		}


		public int NumDatasets
		{
			get
			{
				return mint_num_datasets ; 
			}
		}

	}
}
