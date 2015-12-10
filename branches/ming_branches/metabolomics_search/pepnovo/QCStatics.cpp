#include "QuickClustering.h"
#include "auxfun.h"





// initializes the minimal number of occurences necessary so a 
// peak doesn't get filtered out of the cluster spectrum
// assumes approximately 100 peaks per 100 Daltons...
void ClusterSpectrum::init_min_num_occurences(mass_t tolerance)
{
	static const double target_prob = 0.975;
	mass_t p;

	if (tolerance<0.01)
	{
		p = tolerance * 1.5;
	}
	else if (tolerance < 0.1)
	{
		p = tolerance * 0.75;
	}
	else
		p = tolerance * 0.2;


	int i;

	int max_min_num = cluster_reset_values[num_cluster_reset_values-1]*2+2;
	min_num_occurences.resize(max_min_num);

	min_num_occurences[0]=1;
	min_num_occurences[1]=1;
	for (i=2; i<max_min_num; i++)
	{
		min_num_occurences[i] = get_min_number_from_binomial_prob(i,p,target_prob);
		if (min_num_occurences[i]<2)
			min_num_occurences[i]=2;

	//	cout << i << " " << min_num_occurences[i] << endl;
	}

	// for each final peak mark how many spectra it can appear in (max)
	// and use the binomoal min_occurrences accordingly...
	// for the peaks that don't get in, choose the strongest one in every group of n

	
	large_cluster_size = min_num_occurences.size();
	large_num_occurence_ratios = (float)min_num_occurences[min_num_occurences.size()-1]/
								 (float)(min_num_occurences.size());
}



void ClusterSpectrum::increase_tmp_storage_size(int num_peaks)
{
	if (tmp_peak_area1.size()>=num_peaks)
		return;

	tmp_peak_area1.resize((int)(num_peaks*1.5));
	tmp_peak_area2.resize((int)(num_peaks*1.5));
}


