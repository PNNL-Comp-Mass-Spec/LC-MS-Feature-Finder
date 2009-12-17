using System;
using System.Collections ; 

namespace UMCManipulation
{
	/// <summary>
	/// Summary description for clsClusterMatches.
	/// </summary>
	public class clsClusterMatches
	{
		public clsClusterMassTagMatch [] marr_matches ;
		public clsMassTags [] marr_mass_tags ; 
		public clsClusters mobj_clusters ; 
		public clsClusterMatches()
		{
			//
			// TODO: Add constructor logic here
			//
		}
		public void Set(clsClusters clusters, clsMassTags [] mass_tags, clsClusterMassTagMatch [] matches)
		{
			try
			{
				mobj_clusters = clusters ; 
				int num_matches = matches.Length ; 
				marr_mass_tags = new clsMassTags [num_matches] ;
				marr_matches = new clsClusterMassTagMatch [num_matches] ;
				for (int match_num = 0 ; match_num < num_matches ; match_num++)
				{
					clsClusterMassTagMatch match = matches[match_num] ; 
					marr_mass_tags[match_num] = mass_tags[match.TagIndex] ; 
					marr_matches[match_num] = new clsClusterMassTagMatch(match_num, match.ClusterId) ; 
				}
			}
			catch (Exception e)
			{
				Console.WriteLine(e.Message + e.StackTrace) ; 
			}

		}
	}
}
