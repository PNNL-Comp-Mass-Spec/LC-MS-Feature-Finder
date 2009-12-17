using System;
using System.Collections ; 

namespace UMCManipulation
{
	public enum enmClusterRepresentative { MEAN = 0, MEDIAN} ; 
	public enum enmClusterIntensity {MAX_PER_DATASET = 0, SUM_PER_DATASET} ; 

	/// <summary>
	/// Summary description for clsClusters.
	/// </summary>
	public class clsClusters
	{
		public clsUMCData mobj_umc_data ; 
		private int mint_num_clusters ; 
		private double mdbl_mass_tolerance ;
		private double mdbl_net_tolerance ; 
		private int [] marr_dataset_num ; 
		public double [] marr_mass ;
		public double [] marr_mass_calibrated ;
		public int [] marr_scans ; 
		public double [] marr_nets ; 
		public bool mbln_msms_aligned  ; 
		// the following nets are the nets after alignment to a mass tag database. 
		public double [] marr_nets_transformed ; 
		public int [] marr_charges ; 
		public int [] marr_umc_2_cluster ; 
		public ArrayList marr_clusters_2_umcs ; 
		public double [] marr_cluster_intensity ; 
		public double [] marr_cluster_intensity_normalized ; 
		public int [] marr_cluster_main_member ; 

		enmClusterRepresentative menm_cluster_representative ; 
		enmClusterIntensity menm_cluster_intensity ;

		public clsClusters()
		{
			//
			// TODO: Add constructor logic here
			//
			marr_clusters_2_umcs = new ArrayList() ; 
			mobj_umc_data = new clsUMCData() ; 
			mbln_msms_aligned = false ; 
		}

		public void Clear()
		{
			mint_num_clusters = 0 ; 
			marr_mass = null ; 
			marr_scans = null ; 
			marr_nets = null ; 
			marr_charges = null ; 
			marr_umc_2_cluster = null ; 
			marr_cluster_intensity = null ; 
			marr_dataset_num = null ; 
			mobj_umc_data.Clear() ; 
			marr_clusters_2_umcs.Clear() ; 
			menm_cluster_representative = enmClusterRepresentative.MEAN ; 
			menm_cluster_intensity = enmClusterIntensity.SUM_PER_DATASET ; 
		}


		public clsClusters(clsUMCData umc_data)
		{
			mobj_umc_data = umc_data ; 
		}

		public double MassTolerance
		{
			get
			{
				return mdbl_mass_tolerance ; 
			}
			set
			{
				mdbl_mass_tolerance = value ; 
			}
		}

		public double NETTolerance
		{
			get
			{
				return mdbl_net_tolerance ; 
			}
			set
			{
				mdbl_net_tolerance = value ; 
			}
		}
		
		public int NumClusters
		{
			get
			{
				return mint_num_clusters ; 
			}
			set
			{
				mint_num_clusters = value ; 
			}
		}

		public void CalculateClusterToUMCMaps()
		{
			int num_datasets = mobj_umc_data.NumDatasets ; 


			marr_mass = new double [mint_num_clusters] ; 
			marr_mass_calibrated = new double [mint_num_clusters] ; 
			marr_scans = new int [mint_num_clusters] ; 
			marr_nets = new double [mint_num_clusters] ; 
			marr_charges = new int [mint_num_clusters] ; 


			marr_dataset_num = new int [mobj_umc_data.NumUMCS] ; 

			for (int dataset_num = 0 ; dataset_num < num_datasets ; dataset_num++)
			{
				int start_index = 0 , stop_index = 0 ; 
				mobj_umc_data.GetDataIndex(dataset_num, ref start_index, ref stop_index) ; 
				for (int index = start_index ; index < stop_index ; index++)
				{
					clsPair pair = new clsPair() ; 
					int cluster_num = marr_umc_2_cluster[index] ; 
					pair.Set(cluster_num, index) ; 
					marr_clusters_2_umcs.Add(pair) ; 
					marr_dataset_num[index] = dataset_num ; 
				}	
			}
			// Sort by cluster num.
			marr_clusters_2_umcs.Sort() ; 
		}

		public void CalculateStatistics()
		{
			int num_datasets = mobj_umc_data.NumDatasets ; 
			marr_cluster_intensity = new double [mint_num_clusters*num_datasets] ; 
			marr_cluster_main_member = new int [mint_num_clusters*num_datasets] ; 

			for (int index = 0 ; index < mint_num_clusters*num_datasets ; index++)
			{
				marr_cluster_intensity[index] = 0 ; 
				marr_cluster_main_member[index] = -1 ; 
			}

			ArrayList mzs = new ArrayList() ; 
			ArrayList scans = new ArrayList() ; 
			ArrayList nets = new ArrayList() ; 
			ArrayList charges = new ArrayList() ; 

			// Go through each cluster. 
			int num_elements = marr_clusters_2_umcs.Count ; 
			if (num_elements == 0)
				return ;

			clsUMC umc ; 
			clsPair pair = (clsPair) marr_clusters_2_umcs[0] ; 
			int last_cluster_num = (int) pair.mint_key ;
			int num_nonzero ; 
			double mz ; 
			int charge, avg_charge ; 
			int scan, mid_index ; 
			int avg_scan ;
			double avg_mz ; 
			double net, avg_net ;  
			int pt_num ; 

 
			for (int element_num = 0 ; element_num < num_elements ; element_num++)
			{

				pair = (clsPair) marr_clusters_2_umcs[element_num] ;
				if ((int) pair.mint_key != last_cluster_num)
				{
					num_nonzero = mzs.Count ; 

					switch (menm_cluster_representative)
					{
						case enmClusterRepresentative.MEDIAN:
							if (num_nonzero == 1)
							{
								mz = (double) mzs[0] ; 
								scan = (int) scans[0] ; 
								net = (double) nets[0] ; 
								charge = (int) charges[0] ; 
							}
							else
							{
								mzs.Sort() ; 
								scans.Sort() ; 
								nets.Sort() ; 
								charges.Sort() ; 

								if (num_nonzero%2 == 0)
								{
									mid_index = num_nonzero/2 - 1 ; 
									mz = ((double)mzs[mid_index] + (double)mzs[mid_index+1])/2 ; 
									scan = ((int)scans[mid_index] + (int)scans[mid_index+1])/2 ; 
									net = ((double)nets[mid_index] + (double)nets[mid_index+1])/2 ; 
									charge = (int) charges[mid_index] ; 
								}
								else
								{
									mid_index = num_nonzero/2 ; 
									mz = (double)mzs[mid_index] ; 
									scan = (int)scans[mid_index] ; 
									net = (double)nets[mid_index] ;
									charge = (int) charges[mid_index];
								}
							}
							marr_mass[last_cluster_num] = mz ; 
							marr_mass_calibrated[last_cluster_num] = mz ; 
							marr_scans[last_cluster_num] = scan ; 
							marr_nets[last_cluster_num] = net ; 
							marr_charges[last_cluster_num] = charge ; 
							break ; 
						case enmClusterRepresentative.MEAN:
							avg_mz = 0 ; 
							avg_scan = 0 ; 
							avg_net = 0 ; 
							avg_charge = 0 ; 
							for (pt_num = 0 ; pt_num < num_nonzero ; pt_num++)
							{
								avg_mz += (double)mzs[pt_num] ; 
								avg_scan += (int)scans[pt_num] ; 
								avg_net += (double)nets[pt_num] ; 
								avg_charge += (int) charges[pt_num] ; 
							}
							avg_mz /= num_nonzero ;
							avg_scan /= num_nonzero ; 
							avg_net /= num_nonzero ; 
							avg_charge /= num_nonzero ; 
							marr_mass[last_cluster_num] = avg_mz ; 
							marr_mass_calibrated[last_cluster_num] = avg_mz ; 
							marr_scans[last_cluster_num] = avg_scan ; 
							marr_nets[last_cluster_num] = avg_net ; 
							marr_charges[last_cluster_num] = avg_charge  ; 
							break ; 
					}

					last_cluster_num = (int) pair.mint_key ;
					mzs.Clear() ; 
					scans.Clear() ;
					nets.Clear() ; 
					charges.Clear() ; 
				}
				int index = (int)pair.mobj_val ; 
				umc = mobj_umc_data.marr_umcs[index] ; 

				double rep_mass = umc.mdbl_mono_mass_calibrated ; 
				int rep_scan = umc.mint_scan_aligned ; 
				double rep_net = umc.mdbl_net ; 
				double intensity = umc.mdbl_abundance ; 

				int pt_index = last_cluster_num * num_datasets + marr_dataset_num[index] ; 

				if(marr_cluster_intensity[pt_index] < intensity)
				{
					marr_cluster_main_member[pt_index] = index ;
				}

				if (menm_cluster_intensity == enmClusterIntensity.SUM_PER_DATASET)
				{
					marr_cluster_intensity[pt_index] += intensity ;
				}
				else
				{
					if(marr_cluster_intensity[pt_index] < intensity)
					{
						marr_cluster_intensity[pt_index] = intensity ;
					}
				}


				mzs.Add(rep_mass) ; 
				scans.Add(rep_scan) ; 
				nets.Add(rep_net) ; 
				charges.Add(umc.mint_class_rep_charge) ; 
			}

			// the last cluster.
			num_nonzero = mzs.Count ; 

			switch (menm_cluster_representative)
			{
				case enmClusterRepresentative.MEDIAN:
					if (num_nonzero == 1)
					{
						mz = (double) mzs[0] ; 
						scan = (int) scans[0] ; 
						net = (double) nets[0] ; 
						charge = (int) charges[0] ; 
					}
					else
					{
						mzs.Sort() ; 
						scans.Sort() ; 
						nets.Sort() ; 
						charges.Sort() ; 

						if (num_nonzero%2 == 0)
						{
							mid_index = num_nonzero/2 - 1 ; 
							mz = ((double)mzs[mid_index] + (double)mzs[mid_index+1])/2 ; 
							scan = ((int)scans[mid_index] + (int)scans[mid_index+1])/2 ; 
							net = ((double)nets[mid_index] + (double)nets[mid_index+1])/2 ; 
							charge = (int) charges[mid_index] ; 
						}
						else
						{
							mid_index = num_nonzero/2 ; 
							mz = (double)mzs[mid_index] ; 
							scan = (int)scans[mid_index] ; 
							net = (double)nets[mid_index] ; 
							charge = (int) charges[mid_index] ; 
						}
					}
					marr_mass[last_cluster_num] = mz ; 
					marr_mass_calibrated[last_cluster_num] = mz ; 
					marr_scans[last_cluster_num] = scan ; 
					marr_nets[last_cluster_num] = net ; 
					marr_charges[last_cluster_num] = charge ; 
					break ; 
				case enmClusterRepresentative.MEAN:
					avg_mz = 0 ; 
					avg_scan = 0 ; 
					avg_net = 0 ; 
					avg_charge = 0 ; 
					for (pt_num = 0 ; pt_num < num_nonzero ; pt_num++)
					{
						avg_mz += (double)mzs[pt_num] ; 
						avg_scan += (int)scans[pt_num] ; 
						avg_net += (double)nets[pt_num] ; 
						avg_charge += (int) charges[pt_num] ; 
					}
					avg_mz /= num_nonzero ;
					avg_scan /= num_nonzero ; 
					avg_net /= num_nonzero ; 
					avg_charge /= num_nonzero ; 

					marr_mass[last_cluster_num] = avg_mz ; 
					marr_mass_calibrated[last_cluster_num] = avg_mz ; 
					marr_scans[last_cluster_num] = avg_scan ; 
					marr_nets[last_cluster_num] = avg_net ; 
					marr_charges[last_cluster_num] = (int) (avg_charge +0.5) ; 
					break ; 
			}

			last_cluster_num = (int) pair.mint_key ;
			mzs.Clear() ; 
			scans.Clear() ;
			nets.Clear() ; 
			charges.Clear() ; 
		}
	}
}
