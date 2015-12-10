#include "QuickClustering.h"
#include "auxfun.h"
#include "AnnotatedSpectrum.h"


// static member dclr

vector<QCPeak> ClusterSpectrum::tmp_peak_area1;
vector<QCPeak> ClusterSpectrum::tmp_peak_area2;

vector<int> ClusterSpectrum::min_num_occurences;
int   ClusterSpectrum::num_top_peaks_per_1000_da;
float ClusterSpectrum::large_num_occurence_ratios;
int   ClusterSpectrum::large_cluster_size;


void ClusterSpectrum::filter_peaks_with_slidinig_window()
{
	vector<bool> inds;
	vector<QCPeak> new_peaks;

	mark_top_peaks_with_sliding_window(&peaks[0],
									   peaks.size(), 
									   config->get_local_window_size(),
									   config->get_max_number_peaks_per_local_window(),
									   inds);

	int i;
	for (i=0; i<inds.size(); i++)
	{
		if (inds[i])
			new_peaks.push_back(peaks[i]);
	}

	peaks=new_peaks;

}


struct peak_idx_pair
{
	bool operator< (const peak_idx_pair& other) const
	{
		return (intensity>other.intensity);
	}
	int idx;
	intensity_t intensity;
};
/***********************************************************************
Uses a heuristic approach jumps every half window
************************************************************************/
bool mark_top_peaks_with_sliding_window(const QCPeak *peaks, 
										int num_peaks, 
										mass_t window_size, 
										int num_peaks_per_window, 
										vector<bool>& indicators)
{
	// filter low intensity noise
	// and mark those that are good peaks
	const mass_t half_window_size = 0.5 * window_size;
	const int max_peak_idx = num_peaks -1;

	if (num_peaks<=5)
	{
		indicators.resize(num_peaks,true);
		return false;
	}
	int i;
	for (i=0; i<5; i++)
	{
		if (peaks[i].scaled_intensity<=0)
			break;
	}

	const bool use_scaled_intensity = (i==5);
	int start_window_idx =0;

	indicators.resize(num_peaks,false);
	indicators[0]=true;
	indicators[max_peak_idx]=true;

	while (start_window_idx<max_peak_idx)
	{
		const mass_t max_window_mass = peaks[start_window_idx].mass + window_size;

		int end_window_idx=start_window_idx;
		while (end_window_idx<max_peak_idx && peaks[end_window_idx].mass<max_window_mass)
			end_window_idx++;


		if (end_window_idx - start_window_idx>num_peaks_per_window)
		{
			const int num_peaks_in_window = end_window_idx - start_window_idx+1;
			vector<peak_idx_pair> pairs;
			pairs.resize(num_peaks_in_window);

			if (use_scaled_intensity)
			{
				int i;
				for (i=0; i<num_peaks_in_window; i++)
				{
					const int peak_idx = i+start_window_idx;
					peak_idx_pair& pair = pairs[i];
					pair.idx = peak_idx ;
					pair.intensity = peaks[peak_idx].scaled_intensity;
				}
			}
			else
			{
				int i;
				for (i=0; i<num_peaks_in_window; i++)
				{
					const int peak_idx = i+start_window_idx;
					peak_idx_pair& pair = pairs[i];
					pair.idx = peak_idx ;
					pair.intensity = peaks[peak_idx].intensity;
				}
			}

			sort(pairs.begin(),pairs.end());

			if (pairs[0].intensity<pairs[1].intensity)
			{
				printf("Error: with peak intensity order (possible corruption in the files)!\n");
			//	int i;
			//	for (i=0; i<pairs.size(); i++)
			//		cout << i << " " << pairs[i].intensity << endl;
			//	exit(1);
				return false;
			}

			int i;
			for (i=0; i<num_peaks_per_window; i++)
				indicators[pairs[i].idx]=true;	
		}
		else 
		{
			int i;
			for (i=start_window_idx; i<=end_window_idx; i++)
				indicators[i]=true;
		}

		// advance half a window
		const mass_t mid_mass = peaks[start_window_idx].mass + half_window_size;
		start_window_idx++;
		while (start_window_idx<max_peak_idx && peaks[start_window_idx].mass<mid_mass)
			start_window_idx++;

	}

	return true;
}






/****************************************************************
*****************************************************************/
void ClusterSpectrum::create_new_cluster(Config *config,
										 BasicSpectrum& bs,
										 int cluster_idx)
{
	int i;
	this->tmp_cluster_idx = cluster_idx;
	this->config = config;
	this->tolerance = config->get_tolerance();
	this->m_over_z = bs.ssf->m_over_z;
	this->peptide_str = bs.ssf->peptide.as_string(config);

	this->num_spectra_in_cluster =1;

	if (bs.ssf->sqs >=0)
		this->best_sqs_spec_idx = 0;

	retention_time = bs.ssf->retention_time;

	maximum_good_peaks_to_output = (int)( (config->get_number_of_strong_peaks_per_local_window()/
									 config->get_local_window_size()) 
									 * bs.peaks[bs.num_peaks-1].mass);

	maximum_peaks_vector_size = (int)(2.5*maximum_good_peaks_to_output);

	if (bs.num_peaks>maximum_peaks_vector_size)
	{
		maximum_peaks_vector_size = bs.num_peaks;
	}

	peaks.reserve(maximum_peaks_vector_size);
	peaks.resize(bs.num_peaks);
	for (i=0; i<bs.num_peaks; i++)
		peaks[i]=bs.peaks[i];

	bs.ssf->assigned_cluster = cluster_idx;
	basic_spectra.push_back(bs);	
}






/*************************************************************************
// adds the spectrum to this cluster.
**************************************************************************/
void ClusterSpectrum::add_spectrum_to_cluster(BasicSpectrum& bs, 
											  const vector<int>& spec_top_idxs,
											  float top_x_masses[NUM_TOP_CLUSTER_PEAKS])
{
	int i;

	bs.ssf->assigned_cluster = tmp_cluster_idx;
	basic_spectra.push_back(bs);
	const mass_t tolerance = (config->get_tolerance()>0.1? config->get_tolerance()*0.8 : config->get_tolerance());

	if (bs.ssf->sqs<0 || basic_spectra.size()>MAX_SIZE_FOR_SQS_REP)
	{
		add_peak_list(bs.peaks,bs.num_peaks,tolerance,1);
		static vector<int> cluster_top_idxs;
		set_adjusted_inten(&peaks[0],peaks.size());
		set_cluster_m_over_z();
		select_top_peak_idxs(&peaks[0],peaks.size(),m_over_z,tolerance,
				cluster_top_idxs, top_peak_masses, num_top_peaks_per_1000_da, config);
		set_top_ranked_idxs(cluster_top_idxs);
	}
	else
	{
		if (basic_spectra.size() == MAX_SIZE_FOR_SQS_REP)
		{
			create_consensus_by_binning_basic_spectra();
		}
		else
		{
			if (best_sqs_spec_idx<0)
			{
				cout << "Error: best_sqs_idx<0 ! " << endl;
				best_sqs_spec_idx=0;
			}

			if (bs.ssf->sqs>basic_spectra[best_sqs_spec_idx].ssf->sqs)
			{
				best_sqs_spec_idx = basic_spectra.size()-1;
				peaks.resize(bs.num_peaks);
				int i;
				for (i=0; i<bs.num_peaks; i++)
					peaks[i]=bs.peaks[i];
				
				m_over_z = bs.ssf->m_over_z;
				set_top_ranked_idxs(spec_top_idxs);
				set_top_masses(top_x_masses);
			}
		}
	}



	// update retention time
	if (retention_time>0)
	{
		retention_time =0;
		int count=0;
		for (i=0; i<basic_spectra.size(); i++)
		{
			const float &rt = basic_spectra[i].ssf->retention_time;
			if (rt>=0)
			{
				count++;
				retention_time += rt;
			}
		}
		retention_time /= (float)count;
	}

	// look for a consensus string
	vector<string> peps;
	vector<int> counts;

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->peptide.get_num_aas()>3)
		{
			string pep_str = basic_spectra[i].ssf->peptide.as_string(config);
			int j;
			for (j=0; j<peps.size(); j++)
				if (! strcmp(pep_str.c_str(),peps[j].c_str()) )
					break;
				
			if (j==peps.size())
			{
				peps.push_back(pep_str);
				counts.push_back(1);
			}
			else
				counts[j]++;
		}
	}

	if (peps.size() == 1)
	{
		this->peptide_str = peps[0];
	}
}


/*************************************************************************
 tries to add the cluster
 succeeds only if the similarity of the two originals to the new consensus
 is above the sim_tresh (returns true if it made the addition, false otherwise)
**************************************************************************/
bool ClusterSpectrum::add_cluster(ClusterSpectrum& cs, float sim_thresh)
{
	const int size_before = basic_spectra.size();

	int i;
	for (i=0; i<cs.basic_spectra.size(); i++)
	{
		BasicSpectrum& bs = cs.basic_spectra[i];

		bs.ssf->assigned_cluster = tmp_cluster_idx;
		basic_spectra.push_back(bs);
	}

	mass_t tolerance = (config->get_tolerance()>0.1? config->get_tolerance()*0.8 : config->get_tolerance());


	if (basic_spectra.size()<=MAX_SIZE_FOR_SQS_REP)
	{
		if (best_sqs_spec_idx<0)
		{
		//	cout << "Error: found best_sqs_idx<0 !!" << endl;
			best_sqs_spec_idx = 0;
		}

		float best_sqs = basic_spectra[best_sqs_spec_idx].ssf->sqs;
		int new_idx=-1;
		int i;
		for (i=size_before; i<basic_spectra.size(); i++)
		{
			if (basic_spectra[i].ssf->sqs > best_sqs)
			{
				new_idx = i;
				best_sqs = basic_spectra[i].ssf->sqs;
			}
		}

		if (new_idx>0)
		{
			best_sqs_spec_idx = new_idx;
			peaks.resize(cs.peaks.size());
			int i;
			for (i=0; i<cs.peaks.size(); i++)
				peaks[i]=cs.peaks[i];

			m_over_z = cs.get_m_over_z();
			set_top_masses(cs.get_top_peak_masses());
			set_top_ranked_idxs(cs.get_top_ranked_idxs());
		}
	}
	else
	{
		// avoid creating consensus for very large clusters
		if (size_before <= MAX_SIZE_FOR_SQS_REP && 
			cs.basic_spectra.size() <= MAX_SIZE_FOR_SQS_REP)
		{
			create_consensus_by_binning_basic_spectra(); 
		}
		else
		{
			add_peak_list(cs.get_peaks_pointer(),
						  cs.get_num_peaks(),tolerance, 
					      cs.num_spectra_in_cluster);

			static vector<int> cluster_top_idxs;
			set_adjusted_inten(&peaks[0],peaks.size());

			set_cluster_m_over_z();

			select_top_peak_idxs(&peaks[0],peaks.size(),m_over_z,tolerance,
				cluster_top_idxs, top_peak_masses, num_top_peaks_per_1000_da, config);

			set_top_ranked_idxs(cluster_top_idxs);
		}
	}



	// look for a consensus string
	vector<string> peps;
	vector<int> counts;

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->peptide.get_num_aas()>3)
		{
			string pep_str = basic_spectra[i].ssf->peptide.as_string(config);
			int j;
			for (j=0; j<peps.size(); j++)
				if (! strcmp(pep_str.c_str(),peps[j].c_str()) )
					break;
			
			if (j==peps.size())
			{
				peps.push_back(pep_str);
				counts.push_back(1);
			}
			else
				counts[j]++;

		}
	}

	if (peps.size() == 1)
	{
		this->peptide_str = peps[0];
	}

	return true;
}

/*************************************************************************
// adds the given list to the clusters current list.
// if the new cluster is larger than X, peak scores are weighted according
// to the percentage in which the peaks appear
// list of peaks is then filtered to remove excess peaks
**************************************************************************/
bool ClusterSpectrum::add_peak_list(
		const QCPeak *second_peaks, 
		int num_second_peaks, 
		mass_t tolerance, 
		int num_basic_spectra_added,
		bool need_to_scale)
{
	static vector<QCPeak> tmp_peaks1, tmp_peaks2;
	static vector<float> scaling_factors;

	if (scaling_factors.size()<100)
	{
		scaling_factors.resize(101);
		int i;
		for (i=0; i<=100; i++)
			scaling_factors[i]=0.2 + 0.2 *pow(1+i*0.01,5);
	}


	num_spectra_in_cluster+=num_basic_spectra_added;

	const int    num_joined = this->peaks.size() + num_second_peaks;	

	if (tmp_peaks1.size() < num_joined)
	{
		tmp_peaks1.resize(2*num_joined);
		tmp_peaks2.resize(2*num_joined);
	}

	if (num_joined<=0)
		return true;

	// place a merged list of peaks in tmp_peaks1
	const int num_peaks = peaks.size();
	int a_idx=0, b_idx=0, p_idx=0;

	while (a_idx<num_peaks && b_idx<num_second_peaks)
	{
		if (peaks[a_idx].mass<= second_peaks[b_idx].mass)
		{
			tmp_peaks1[p_idx++]=peaks[a_idx++];
		}
		else
			tmp_peaks1[p_idx++]=second_peaks[b_idx++];
	}

	while (a_idx<num_peaks)
		tmp_peaks1[p_idx++]=peaks[a_idx++];

	while (b_idx<num_second_peaks)
		tmp_peaks1[p_idx++]=second_peaks[b_idx++];

	// use 3 rounds with increasing tolerance to join peaks
	// each time put them in the new list in area2
	vector<mass_t> tolerances;
	tolerances.resize(3,0);
	tolerances[0]=tolerance*0.25;
	tolerances[1]=tolerance*0.5;
	tolerances[2]=tolerance;
	
	int round;
	int num_area1 = p_idx;
	for (round=0; round<3; round++)
	{
		const mass_t join_tolerance = tolerances[round];
		vector<QCPeak>& area1 = ( (round == 1) ? tmp_peaks2 : tmp_peaks1);
		vector<QCPeak>& area2 = ( (round == 1) ? tmp_peaks1 : tmp_peaks2);
		int j=0;
		int i;

		area2[0]=area1[0];
		for (i=1; i<num_area1; i++)
		{
			if (area1[i].mass - area2[j].mass<join_tolerance)
			{
				const QCPeak& peak1 = area1[i];
				QCPeak& peak2 = area2[j];

				const int num_occurences = peak2.num_occurences + peak1.num_occurences;

				intensity_t sum_intens = peak1.intensity + peak2.intensity;
				mass_t weight = peak1.intensity / sum_intens;
				peak2.mass = (peak1.mass * weight) + (1.0-weight)* peak2.mass;
				peak2.intensity = sum_intens;
			
				peak2.num_occurences = num_occurences;
			}
			else
				area2[++j]=area1[i];
		}
		num_area1 = j+1;
	}

	const int merged_num_peaks = num_area1;

	// scale intensity according to peak probability
	if (need_to_scale && num_spectra_in_cluster>2)
	{
		const float idx_mult = 100.0 / num_spectra_in_cluster;
		int i;

		for (i=0; i<merged_num_peaks; i++)
		{
			QCPeak& peak = tmp_peaks2[i];
			int idx = (int)(idx_mult*peak.num_occurences);
			if (idx>100)
				idx=100;

			peak.scaled_intensity = peak.intensity * scaling_factors[idx];
		}
	}

	// check if we can just copy the peaks
	if (merged_num_peaks <= maximum_peaks_vector_size)
	{
		// copy over peaks
		peaks.clear();
		int i;
		for (i=0; i<merged_num_peaks; i++)
			peaks.push_back(tmp_peaks2[i]);

		return true;
	}

	vector<bool> indicators;
	mark_top_peaks_with_sliding_window(
		&tmp_peaks2[0],
		merged_num_peaks,
		config->get_local_window_size(),
		(int)(config->get_max_number_peaks_per_local_window()*2.5),
		indicators);

//	cout << config->get_local_window_size() << "\t" <<  config->get_max_number_peaks_per_local_window() << endl;

	peaks.clear();
	int i;
	for (i=0; i<merged_num_peaks; i++)
		if (indicators[i])
			peaks.push_back(tmp_peaks2[i]);

	return true;
}


/***********************************************************************
// recursively merges the peak lists from the various spectra
// performs a merge until the pointers list has only one entry
************************************************************************/
void ClusterSpectrum::merge_peak_lists(vector<QCPeak>& org_peaks,
									   vector<QCPeak>& new_peaks,
									   vector<PeakListPointer>& pointers)
{
	QCPeak *org_peak_area = &org_peaks[0];
	QCPeak *new_peak_area = &new_peaks[0];
	while (pointers.size()>1)
	{
		vector<PeakListPointer> new_pointers;
		QCPeak *n_pos = new_peak_area;
		new_pointers.clear();
		int i;
		for (i=0; i<pointers.size(); i+=2)
		{
			// merge the two lists
			if (i<pointers.size()-1)
			{
				QCPeak* first_list = pointers[i].peaks;
				QCPeak* second_list = pointers[i+1].peaks;
				QCPeak* merge_start = n_pos;

				const int n1 = pointers[i].num_peaks;
				const int n2 = pointers[i+1].num_peaks;
			
				// merge
				int i1=0,i2=0;
				while (i1<n1 && i2<n2)
				{
					if (first_list[i1].mass < second_list[i2].mass)
					{
						*n_pos++=first_list[i1++];
					}
					else
						*n_pos++=second_list[i2++];
				}

				while (i1<n1)
					*n_pos++=first_list[i1++];

				while (i2<n2)
					*n_pos++=second_list[i2++];

				// add pointer for new merged list
				PeakListPointer new_pointer;
				new_pointer.peaks = merge_start;
				new_pointer.num_peaks = n1+n2;
				new_pointers.push_back(new_pointer);
			}
			else // write peaks directly to new area
			{	
				QCPeak* first_list = pointers[i].peaks;
				QCPeak* start_copy = n_pos;
				const int n1 = pointers[i].num_peaks;
				int i1=0;
				while (i1<n1)
					*n_pos++=first_list[i1++];
		
				PeakListPointer new_pointer;
				new_pointer.peaks = start_copy;
				new_pointer.num_peaks = n1;
				new_pointers.push_back(new_pointer);
				n_pos+= n1;
			}
		}

		pointers = new_pointers;

		// switch between the pointers of the peak storage areas
		QCPeak *tmp = org_peak_area;
		org_peak_area = new_peak_area;
		new_peak_area = tmp;
	}
}



/************************************************************************
// joins adjacent peaks, marks invaldiated peaks by assigning their mass to -1
// then condences list to contain only good peaks.
// also calculates for each peak mass what is the maximal number of peaks
// that could be detected at that range (based on all the spectra's min/max
// peak values).
*************************************************************************/
void ClusterSpectrum::join_merged_peak_lists(
										PeakListPointer& plp,
										PeakListPointer& alt_plp,
										int num_merged_spectra,
										mass_t tolerance)
{
	vector<mass_t> join_tolerances;
	int i,t;
	int mid_idx = plp.num_peaks/2;

	// set maximal number of detected peaks
	vector<bool> used_ind;
	used_ind.resize(num_merged_spectra,false);

	// fill from left to middle
	int p_idx,max_peaks_at_idx=0;
	for (p_idx=0; p_idx<=mid_idx; p_idx++)
	{
		QCPeak& peak = plp.peaks[p_idx];
		const int& spec_idx = peak.source_spec_idx;
		if (! used_ind[spec_idx])
		{
			used_ind[spec_idx]=true;
			max_peaks_at_idx++;
		}
		peak.max_num_occurences = max_peaks_at_idx;
	}

	// fill from right to middle
	for (i=0; i<num_merged_spectra; i++)
		used_ind[i]=false;

	max_peaks_at_idx=0;
	for (p_idx = plp.num_peaks-1; p_idx>mid_idx; p_idx--)
	{
		QCPeak& peak = plp.peaks[p_idx];
		const int& spec_idx = peak.source_spec_idx;
		if (! used_ind[spec_idx])
		{
			used_ind[spec_idx]=true;
			max_peaks_at_idx++;
		}
		peak.max_num_occurences = max_peaks_at_idx;
	}

	// join peaks

	join_tolerances.push_back(tolerance*0.15);
	join_tolerances.push_back(tolerance*0.3);
	join_tolerances.push_back(tolerance*0.5);

	for (t=0; t<join_tolerances.size();t++)
	{
		const mass_t join_tolerance = join_tolerances[t];
		const int max_p_idx = plp.num_peaks;

		int p_idx=0, np_idx=0;

		while (p_idx<max_p_idx)
		{
			QCPeak& curr_peak = plp.peaks[p_idx];
			if (curr_peak.mass<0)
				continue;

			int next_idx = p_idx+1;
			while (next_idx<max_p_idx && plp.peaks[next_idx].mass<0)
				next_idx++;
			
			// try joining the peak with the next one ahead
			while (next_idx < max_p_idx &&
				   plp.peaks[next_idx].mass - curr_peak.mass < join_tolerance)
			{
				QCPeak& next_peak = plp.peaks[next_idx];

				// don't try and join peaks that should stay apart
				if (t>=2 &&
					(curr_peak.num_occurences >= curr_peak.max_num_occurences ||
					 next_peak.num_occurences >= next_peak.max_num_occurences ) )
					break;
					 
				const intensity_t sum_inten = curr_peak.intensity + next_peak.intensity;

				mass_t new_mass = (curr_peak.mass * curr_peak.intensity + next_peak.mass * next_peak.intensity) / 
								  sum_inten;

				curr_peak.mass      = new_mass;
				curr_peak.intensity = sum_inten;
				curr_peak.num_occurences += next_peak.num_occurences;

				if (curr_peak.max_num_occurences<next_peak.max_num_occurences)
					curr_peak.max_num_occurences = next_peak.max_num_occurences;
				
				next_peak.mass = -1;
				next_peak.intensity = -100000000;
				next_peak.num_occurences=0;
				
				while (next_idx<max_p_idx && plp.peaks[next_idx].mass<0)
					next_idx++;
			}

			// copy peak to new area
			alt_plp.peaks[np_idx++] = plp.peaks[p_idx];

			// advance to next peak
			p_idx = next_idx;
		}

		alt_plp.num_peaks = np_idx;

		// switch the plps
		PeakListPointer tmp;
		tmp=plp;
		plp=alt_plp;
		alt_plp=plp;
	}
}



/************************************************************************
// selects the consensus peaks - those that appear more than the expected cutoff
// also takes some of the stronger peaks that didn't make the cutoff
// writes the selected peaks into the peaks of the cluster spectrum
*************************************************************************/
void ClusterSpectrum::select_consensus_peaks(PeakListPointer& plp, 
											 PeakListPointer& alt_plp,
											 int num_org_spectra)
{
	vector<bool> keep_indicators;
	keep_indicators.clear();

	keep_indicators.resize(plp.num_peaks,false);


	if (num_org_spectra<this->large_cluster_size)
	{
		int i;
		int num_saved_peaks=0;
		for (i=0; i<plp.num_peaks; i++)
		{
			QCPeak& peak = plp.peaks[i];

			if (peak.num_occurences>= min_num_occurences[peak.max_num_occurences])
				keep_indicators[i]=true;
		}
			
	}
	else // for large clusters might need to use the ratio to calc min num_occurrences
	{
		int i;
		int num_saved_peaks=0;
		for (i=0; i<plp.num_peaks; i++)
		{
			QCPeak& peak = plp.peaks[i];

			if (peak.max_num_occurences < large_cluster_size)
			{
				if (peak.num_occurences>= min_num_occurences[peak.max_num_occurences])
					keep_indicators[i]=true;
			}
			else
			{
				int min_num_occurences = (int)(peak.max_num_occurences * this->large_num_occurence_ratios +0.5);
				if (peak.num_occurences>= min_num_occurences)
					keep_indicators[i]=true;
			}
		}
	}


	// copy peaks
	int i;
	int alt_idx=0;
	QCPeak* peak_list = alt_plp.peaks;
	
	for (i=0; i<plp.num_peaks; i++)
		if (keep_indicators[i])
			peak_list[alt_idx++]=plp.peaks[i];
			
	// filter low intensity peaks

	const mass_t half_window_size = 0.5 * config->get_local_window_size();
	const int num_peaks_in_window = config->get_max_number_peaks_per_local_window();
	int max_peak_idx = alt_idx -1;
	int min_idx=1;
	int max_idx=1;
	
	peaks.clear();
	peaks.reserve(max_peak_idx);

	peaks.push_back(peak_list[0]);

	// check the rest of the peaks
	for (i=1; i<max_peak_idx; i++)
	{
		const mass_t& peak_mass=peak_list[i].mass;
		mass_t min_mass = peak_list[min_idx].mass;
		mass_t max_mass = peak_list[max_idx].mass;

	
		// advance min/max pointers
		while (peak_mass-min_mass > half_window_size)
			min_mass=peak_list[++min_idx].mass;

		while (max_idx < max_peak_idx && max_mass - peak_mass <= half_window_size)
			max_mass=peak_list[++max_idx].mass;

		if (max_mass - peak_mass > half_window_size)
			max_idx--;

		// if there are less than the maximum number of peaks in the window, keep it.
		if (max_idx-min_idx < num_peaks_in_window)
		{
			peaks.push_back(peak_list[i]);
			continue;
		}

		// check if this is one of the top peaks in the window
		int higher_count=0;
		for (int j=min_idx; j<=max_idx; j++)
			if (peak_list[j].intensity > peak_list[i].intensity)
				higher_count++;

		if (higher_count < num_peaks_in_window)
			peaks.push_back(peak_list[i]);
	}
	peaks.push_back(peak_list[max_peak_idx]);
}




void ClusterSpectrum::create_consensus_sepctrum_from_peak_list_pointers(
					  vector<PeakListPointer>& plp, int total_num_peaks)
{
	int i;

	// the merged list pointer is in plp[0]
	merge_peak_lists(tmp_peak_area1,tmp_peak_area2,plp);

	QCPeak *peak_list = plp[0].peaks;
	if (plp[0].num_peaks != total_num_peaks)
	{
		cout << "Error: mismatch in peak numbers of merged lists: " <<
			total_num_peaks << " vs. " << plp[0].num_peaks << endl;
		exit(1);
	}

	
	// assign the appropriate plp pointers
	PeakListPointer merged_peaks_plp = plp[0];
	PeakListPointer alt_plp;
	QCPeak *peak_area1 = &tmp_peak_area1[0];
	QCPeak *peak_area2 = &tmp_peak_area2[0];

	if (merged_peaks_plp.peaks == peak_area1)
	{
		alt_plp.peaks = peak_area2;
	}
	else if (merged_peaks_plp.peaks == peak_area2)
	{
		alt_plp.peaks = peak_area1;
	}
	else
	{
		cout << "Error: mismatch in peak areas!" << endl;
		exit(0);
	}


	join_merged_peak_lists(merged_peaks_plp, alt_plp, 
						   basic_spectra.size(), config->get_tolerance());


	select_consensus_peaks(merged_peaks_plp, alt_plp, basic_spectra.size());


	static vector<int> cluster_top_idxs;
	set_adjusted_inten(&peaks[0],peaks.size());

	set_cluster_m_over_z();

	select_top_peak_idxs(&peaks[0],peaks.size(),m_over_z,tolerance,
						cluster_top_idxs, top_peak_masses, 
						num_top_peaks_per_1000_da, config);

	set_top_ranked_idxs(cluster_top_idxs);


	// update retention time
	if (retention_time>0)
	{
		retention_time =0;
		int count=0;
		for (i=0; i<basic_spectra.size(); i++)
		{
			const float &rt = basic_spectra[i].ssf->retention_time;
			if (rt>=0)
			{
				count++;
				retention_time += rt;
			}
		}
		retention_time /= (float)count;
	}

	// look for a consensus string
	vector<string> peps;
	vector<int> counts;

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->peptide.get_num_aas()>3)
		{
			string pep_str = basic_spectra[i].ssf->peptide.as_string(config);
			int j;
			for (j=0; j<peps.size(); j++)
				if (! strcmp(pep_str.c_str(),peps[j].c_str()) )
					break;
			
			if (j==peps.size())
			{
				peps.push_back(pep_str);
				counts.push_back(1);
			}
			else
				counts[j]++;

		}
	}

	if (peps.size() == 1)
	{
		this->peptide_str = peps[0];
	}
}




/************************************************************************
	Creates a single consensus spectrum from the basic spectra.
	First merges all the peak lists into a single (sorted) list.
*************************************************************************/
void ClusterSpectrum::create_cluster_by_binning_basic_spectra()
{
	int i;

	create_new_cluster(config,basic_spectra[0],0);
	
	for (i=1; i<basic_spectra.size(); i++)
	{

		add_peak_list(basic_spectra[i].peaks,
			basic_spectra[i].num_peaks,
			config->get_tolerance()*0.8,1);
	}


	static vector<int> cluster_top_idxs;
	set_adjusted_inten(&peaks[0],peaks.size());

	set_cluster_m_over_z();

	select_top_peak_idxs(&peaks[0],peaks.size(),m_over_z,tolerance,
						cluster_top_idxs, top_peak_masses, 
						num_top_peaks_per_1000_da, config);

	set_top_ranked_idxs(cluster_top_idxs);


	// update retention time
	if (retention_time>0)
	{
		retention_time =0;
		int count=0;
		for (i=0; i<basic_spectra.size(); i++)
		{
			const float &rt = basic_spectra[i].ssf->retention_time;
			if (rt>=0)
			{
				count++;
				retention_time += rt;
			}
		}
		retention_time /= (float)count;
	}

	// look for a consensus string
	vector<string> peps;
	vector<int> counts;

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->peptide.get_num_aas()>3)
		{
			string pep_str = basic_spectra[i].ssf->peptide.as_string(config);
			int j;
			for (j=0; j<peps.size(); j++)
				if (! strcmp(pep_str.c_str(),peps[j].c_str()) )
					break;
			
			if (j==peps.size())
			{
				peps.push_back(pep_str);
				counts.push_back(1);
			}
			else
				counts[j]++;

		}
	}

	if (peps.size() == 1)
	{
		this->peptide_str = peps[0];
	}

/*	increase_tmp_storage_size(total_num_peaks);
	
	// write all peaks from the spectra into a single area
	// and create peak list pointer
	int num_peaks_written=0;
	vector<PeakListPointer> plp;
	plp.resize(basic_spectra.size());
	for (i=0; i<basic_spectra.size(); i++)
	{
		QCPeak *list_start  = &tmp_peak_area1[0] + num_peaks_written;
		const int num_peaks_in_list = basic_spectra[i].num_peaks;
		const QCPeak *org_peaks  = basic_spectra[i].peaks;
		plp[i].peaks = list_start;
		plp[i].num_peaks = num_peaks_in_list;

		int j;
		for (j=0; j<num_peaks_in_list; j++)
		{
			QCPeak& tcp = list_start[j];
			tcp.mass = org_peaks[j].mass;
			tcp.intensity = org_peaks[j].intensity;
			tcp.num_occurences=1;
			tcp.source_spec_idx=i;
		}
		num_peaks_written+=num_peaks_in_list;
	}

	create_consensus_sepctrum_from_peak_list_pointers(plp,total_num_peaks);*/

}


/************************************************************************
Clears the current consensus peak list and creates a consensus spectrum
from all cluster members.
*************************************************************************/
void ClusterSpectrum::create_consensus_by_binning_basic_spectra()
{
	int i;
	for (i=0; i<basic_spectra.size(); i++)
	{
		if (i == best_sqs_spec_idx)
			continue;

		add_peak_list(basic_spectra[i].peaks,basic_spectra[i].num_peaks,tolerance,1);
	}

	static vector<int> cluster_top_idxs;
	set_adjusted_inten(&peaks[0],peaks.size());
	set_cluster_m_over_z();
	select_top_peak_idxs(&peaks[0],peaks.size(),m_over_z,tolerance,
			cluster_top_idxs, top_peak_masses, num_top_peaks_per_1000_da, config);
	set_top_ranked_idxs(cluster_top_idxs);
	
	best_sqs_spec_idx = -1;
}


/************************************************************************
// sets the consensus spectrum to be the basic spectrum with the maximal 
// similarity to other spectra.
*************************************************************************/
int ClusterSpectrum::select_max_similarity_spectrum_as_consensus()
{	
	mass_t tolerance = config->get_tolerance();
	int i;
	vector<float> sim_sums;
	vector< vector<int> > top_idxs;

	top_idxs.resize(basic_spectra.size());
	sim_sums.resize(basic_spectra.size(),0);

	for (i=0; i<basic_spectra.size(); i++)
	{

		BasicSpectrum& spec = basic_spectra[i];
		float top_x_masses[NUM_TOP_CLUSTER_PEAKS];

		

		set_adjusted_inten(spec.peaks,spec.num_peaks);
		select_top_peak_idxs(spec.peaks,spec.num_peaks,spec.ssf->m_over_z,
			tolerance, top_idxs[i], top_x_masses, 20);
	}


	
	for (i=0; i<basic_spectra.size()-1; i++)
	{
		int j;
		for (j=i+1; j<basic_spectra.size(); j++)
		{
			float sim = calc_selected_dot_prod(tolerance,
				basic_spectra[i].peaks,basic_spectra[i].num_peaks, top_idxs[i],
 				basic_spectra[j].peaks,basic_spectra[j].num_peaks, top_idxs[j]);

			sim_sums[i]+=sim;
			sim_sums[j]+=sim;
		}
	}

	float max_sim=-1;
	int max_sim_idx=-1;
	for (i=0; i<basic_spectra.size(); i++)
	{
		if (sim_sums[i]>max_sim)
		{
			max_sim=sim_sums[i];
			max_sim_idx=i;
		}
	}

	
	peaks.resize(basic_spectra[max_sim_idx].num_peaks);
	for (i=0; i<basic_spectra[max_sim_idx].num_peaks; i++)
		peaks[i]=basic_spectra[max_sim_idx].peaks[i];

	set_top_ranked_idxs(top_idxs[max_sim_idx]);
	set_cluster_m_over_z();

	return max_sim_idx;
}



int ClusterSpectrum::select_max_sqs_spectrum_as_consensus()
{
	int max_idx=0;
	float max_sqs=0;
	int i;

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->sqs>max_sqs)
		{
			max_idx=i;
			max_sqs = basic_spectra[i].ssf->sqs;
		}
	}

	BasicSpectrum& spec = basic_spectra[max_idx];
	float top_x_masses[NUM_TOP_CLUSTER_PEAKS];
	vector<int> top_idxs;

	set_adjusted_inten(spec.peaks,spec.num_peaks);
	select_top_peak_idxs(spec.peaks,spec.num_peaks,spec.ssf->m_over_z,
			tolerance, top_idxs, top_x_masses, 20);

	peaks.resize(basic_spectra[max_idx].num_peaks);
	for (i=0; i<basic_spectra[max_idx].num_peaks; i++)
		peaks[i]=basic_spectra[max_idx].peaks[i];


	set_top_ranked_idxs(top_idxs);
	set_cluster_m_over_z();

	return max_idx;
}

/************************************************************************
// sets the m_over_z as the average of the m_over_z of the basic_spectra
*************************************************************************/
void ClusterSpectrum::set_cluster_m_over_z()
{
	int i;

	m_over_z=0;

/*	int c = basic_spectra[0].ssf->charge;
	for (i=1; i<basic_spectra.size(); i++)
		if (basic_spectra[i].ssf->charge != c)
			break;

	// not all charges are the same, use the first m_over_z 
	if (i<basic_spectra.size())
	{
		m_over_z = basic_spectra[0].ssf->m_over_z;
		return;
	} */

	for (i=0; i<basic_spectra.size(); i++)
		m_over_z += basic_spectra[i].ssf->m_over_z;

	if (basic_spectra.size() == 0)
	{
		cout << "Error: cluster has no basic spectra!" << endl;
		exit(1);
	}
	m_over_z /= basic_spectra.size();
}


/************************************************************************
Choose Majority charge
*************************************************************************/
void ClusterSpectrum::set_charge()
{

	vector<int> charge_counts;
	charge_counts.resize(8,0);

	if (basic_spectra.size() == 1)
	{
		charge = basic_spectra[0].ssf->charge;
		return;
	}

	int i;
	for (i=0; i<basic_spectra.size(); i++)
		charge_counts[basic_spectra[i].ssf->charge]++;
	
	int max_charge=0;
	for (i=0; i<8; i++)
		if (charge_counts[i]>=charge_counts[max_charge])
			max_charge=i;

	charge = max_charge;

	if (charge ==0)
		return;

	const float ratio_max = (float)charge_counts[max_charge]/(float)basic_spectra.size();

	if (charge>1 && ratio_max<=0.5)
		charge=0;

	if (charge == 2 || charge == 3)
	{
		float ratio2 = (float)charge_counts[2]/(float)basic_spectra.size();
		float ratio3 = (float)charge_counts[3]/(float)basic_spectra.size();

		if (ratio2>=0.33 && ratio3>=0.33)
			charge=0;
	}
}



// creates the title (file name) for the cluster
void ClusterSpectrum::make_title(string& name, int batch_idx, int cluster_idx)
{
	ostringstream os1;
	ostringstream os2;

	os1 << batch_idx;
	os2 << cluster_idx;

	title = name + "." + os1.str() + "." + os2.str();
}


// finds how many spectra have a peptide sequence that 
// doesn't match the majority assignment
// adds spectra to the mismatch category only if they contain 
// masses that are not within tolerance of each other (modulo -3,-2,-1,+1,+2,+3)
int  ClusterSpectrum::get_num_misassigned_spectra(mass_t pm_tolerance, int* mismatched_with_pep) const
{
	vector<string> peptide_strings;
	vector<int> counts;
	int i;

	if (mismatched_with_pep)
		*mismatched_with_pep=0;

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->peptide.get_num_aas()<2)
			continue;

		const string& pep_str = basic_spectra[i].ssf->peptide.as_string(config);
	
		int j;
		for (j=0; j<peptide_strings.size(); j++)
		{
			if (! strcmp(peptide_strings[j].c_str(),pep_str.c_str()))
			{
				counts[j]++;
				break;
			}
		}

		if (j<peptide_strings.size())
			continue;

		peptide_strings.push_back(pep_str);
		counts.push_back(1);
	}

	if (counts.size()<=1)
		return 0;

	int max_idx=0;
	int max_count=counts[0];

	for (i=1; i<peptide_strings.size(); i++)
	{
		if (counts[i]>max_count)
		{
			max_idx=i;
			max_count=counts[i];
		}
	}

	Peptide max_pep;
	string max_pep_str = peptide_strings[max_idx];
	max_pep.parse_from_string(config,max_pep_str);
	mass_t max_pep_mass_with_19 = max_pep.get_mass() + MASS_OHHH;

	int num_missmatched=0;
	int num_missmatched_with_peptide=0;
	for (i=0; i<basic_spectra.size(); i++)
	{
		string pep_str = basic_spectra[i].ssf->peptide.as_string(config);
		if (pep_str.length()>0 && ! strcmp(pep_str.c_str(),max_pep_str.c_str()) )
			continue;

		// check that if parent mass is wrong
		mass_t pep_m_over_z = basic_spectra[i].ssf->m_over_z;
		mass_t best_off = 9999;
		mass_t best_pep_mass = -1;
		int c;
		for (c=1; c<=4; c++) // assume the spectrum's charge is between 1-4
		{
			mass_t pep_mass = pep_m_over_z * c - (c-1)* MASS_PROTON;
			mass_t pm_offset;
			for (pm_offset=-4.0; pm_offset<=4.0; pm_offset+=1.0)
			{
				mass_t offset = fabs(max_pep_mass_with_19 - pep_mass - pm_offset);
				if (offset<best_off)
				{
					best_off = offset;
					best_pep_mass = pep_mass;
				}
			}
		}

		if (best_off>pm_tolerance)
		{
		//	cout << setprecision(8) << max_pep_mass_with_19 << " " << best_pep_mass << "  " << 
		//		max_pep_mass_with_19 - best_pep_mass << endl;
			num_missmatched++;

			if (pep_str.length()>0)
			{
			//	cout << pep_str << " " << num_missmatched_with_peptide << endl;
				num_missmatched_with_peptide++;
			}
		}
	}

	if (mismatched_with_pep)
		*mismatched_with_pep=num_missmatched_with_peptide;

	return num_missmatched;
}



// returns true if there is a mjority annotation (with 75% of the annotated spectra)
// if so, returns its string and mass in the variables
bool ClusterSpectrum::has_majority_annotation(string& consensus_pep_str, mass_t& consensus_pep_mass) const
{
	vector<string> peptide_strings;
	vector<int> counts;
	int num_annoatated_spectra=0;
	int i;

	peptide_strings.clear();
	counts.clear();
	

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->peptide.get_num_aas()<2)
			continue;

		const string& pep_str = basic_spectra[i].ssf->peptide.as_string(config);
		num_annoatated_spectra++;
	
		int j;
		for (j=0; j<peptide_strings.size(); j++)
		{
			if (! strcmp(peptide_strings[j].c_str(),pep_str.c_str()))
			{
				counts[j]++;
				break;
			}
		}

		if (j<peptide_strings.size())
			continue;

		peptide_strings.push_back(pep_str);
		counts.push_back(1);
	}

	if (counts.size()==0)
		return false;

//	cout << ">>> " << basic_spectra.size() << " <" << counts.size() << "> " <<
//		num_annoatated_spectra << endl;

	

	int max_idx=0;
	int max_count=counts[0];

	for (i=1; i<counts.size(); i++)
	{
		if (counts[i]>max_count)
		{
			max_idx=i;
			max_count=counts[i];
		}
	}

//	cout <<  max_count << " " << num_annoatated_spectra << endl;
	if ( ((double)max_count/(double)num_annoatated_spectra) <0.75)
		return false;

	Peptide max_pep;

	consensus_pep_str = peptide_strings[max_idx];
	max_pep.parse_from_string(config,consensus_pep_str);
	consensus_pep_mass = max_pep.get_mass() + MASS_OHHH;
	
	return true;
}


int ClusterSpectrum::get_num_basic_spectra_with_peptide() const
{
	int num_peps=0;
	int i;

	for (i=0; i<basic_spectra.size(); i++)
		if (basic_spectra[i].ssf->peptide.get_num_aas()>1)
			num_peps++;
	
	return num_peps;
}

	
void ClusterSpectrum::print_cluster_peptides() const
{
	vector<string> peptide_strings;
	vector<int> counts;
	int i;

	for (i=0; i<basic_spectra.size(); i++)
	{
		if (basic_spectra[i].ssf->peptide.get_num_aas()<3)
			continue;

		const string& pep_str = basic_spectra[i].ssf->peptide.as_string(config);

		int j;
		for (j=0; j<peptide_strings.size(); j++)
		{
			if (! strcmp(peptide_strings[j].c_str(),pep_str.c_str()))
			{
				counts[j]++;
				break;
			}
		}

		if (j<peptide_strings.size())
			continue;

		peptide_strings.push_back(pep_str);
		counts.push_back(1);
	}

	
	int num_with_peptide=0;
	for (i=0; i<counts.size(); i++)
		num_with_peptide+=counts[i];

	cout << "Cluster: " << this->tmp_cluster_idx << "  " << basic_spectra.size() << " (";
	cout << num_with_peptide << ")" << endl;
	for (i=0; i<counts.size(); i++)
		cout << setw(4) << left << counts[i] << peptide_strings[i].c_str() << endl;
}



void ClusterSpectrum::print_cluster_similarities()
{
	int i;
	mass_t tolerance = config->get_tolerance();

	for (i=0; i<basic_spectra.size(); i++)
	{
		BasicSpectrum& spec = basic_spectra[i];
		float top_x_masses[NUM_TOP_CLUSTER_PEAKS];
		vector<int> spec_top_idxs;
		set_adjusted_inten(spec.peaks,spec.num_peaks);
		select_top_peak_idxs(spec.peaks,spec.num_peaks,spec.ssf->m_over_z,
			tolerance,spec_top_idxs, top_x_masses, num_top_peaks_per_1000_da);


		float sim = calc_selected_dot_prod(tolerance,
					&peaks[0], peaks.size(), top_ranked_peak_idxs,
					spec.peaks,spec.num_peaks, spec_top_idxs);

		cout << i << " " << sim << " " << spec.ssf->single_name;
	}
}



void BasicSpectrum::output_to_mgf(ostream& mgf, Config *config, const char *seq) const
{
	mgf << "BEGIN IONS" << endl;
	mgf << "TITLE=" <<  ssf->single_name << endl;
	
	if (ssf->peptide.get_num_aas()>0)
	{
		mgf << "SEQ=" << ssf->peptide.as_string(config) << endl;
	}
	else if (seq && strlen(seq)>2)
		mgf << "SEQ=" << seq << endl;
	
	if (ssf->type == MZXML)
	{
		MZXML_single *mzxml_single = (MZXML_single *)ssf;
		if (mzxml_single->scan_number>=0)
			mgf << "SCANS=" <<mzxml_single->scan_number << endl;

//		if (mzxml_single->retention_time>=0)
//			mgf << "RT=" << mzxml_single->retention_time << endl;
	}

	if (ssf->type == MGF)
	{
		MGF_single *mgf_single = (MGF_single *)ssf;
		if (mgf_single->scan_number >=0)
			mgf << "SCANS=" << mgf_single->scan_number << endl;
	}

	if (this->ssf->retention_time>0)
		mgf << "RTINSECONDS="<< ssf->retention_time << endl;

	mgf << "CHARGE=+" << ssf->charge << endl;
		
	mgf << "PEPMASS=" << ssf->m_over_z << endl;
	
	int i;
	for (i=0; i<this->num_peaks; i++)
		mgf << fixed << setprecision(3) << peaks[i].mass << " " << peaks[i].intensity << endl;

	mgf << "END IONS" << endl << endl;
}



void BasicSpectrum::output_to_mgf(FILE* mgf_stream, Config *config, const char *seq) const
{
	fprintf(mgf_stream,"BEGIN IONS\n");
	fprintf(mgf_stream,"TITLE=%s\n",ssf->single_name.c_str());
	
	if (ssf->peptide.get_num_aas()>0)
	{
		fprintf(mgf_stream,"SEQ=%s\n",ssf->peptide.as_string(config).c_str());
	}
	else if (seq && strlen(seq)>2)
		fprintf(mgf_stream,"SEQ=%s\n", seq);
	
	if (ssf->type == MZXML)
	{
		MZXML_single *mzxml_single = (MZXML_single *)ssf;
		if (mzxml_single->scan_number>=0)
			fprintf(mgf_stream,"SCANS=%d\n",mzxml_single->scan_number);
	}

	if (ssf->type == MGF)
	{
		MGF_single *mgf_single = (MGF_single *)ssf;
		if (mgf_single->scan_number >=0)
			fprintf(mgf_stream,"SCANS=%d\n",mgf_single->scan_number);
	}

	fprintf(mgf_stream,"CHARGE=+%d\n",ssf->charge);
		
	fprintf(mgf_stream,"PEPMASS=%.3f\n",ssf->m_over_z);
	
	int i;
	for (i=0; i<this->num_peaks; i++)
		fprintf(mgf_stream,"%.3f %.3f\n",peaks[i].mass, peaks[i].intensity);

	fprintf(mgf_stream,"END IONS\n\n");
}




void ClusterSpectrum::write_cluster_basic_spectra_to_mgf(char *mgf_file) const
{
	ofstream mgf(mgf_file);
	int i;
	for (i=0; i<basic_spectra.size(); i++)
		basic_spectra[i].output_to_mgf(mgf,config);
	mgf.close();
}

// writes the spectrum to the output file in the mgf format
void ClusterSpectrum::write_spectrum_to_mgf(ostream& mgf,bool write_peak_count, bool ind_write_charge) const
{
	const float min_peak_intensity = 0.0005;

	mgf << "BEGIN IONS" << endl;
	mgf << "TITLE=" <<  title << endl;

	if (peptide_str.length()>0)
		mgf << "SEQ=" << peptide_str << endl;

	if (this->charge>0 && (config->get_use_spectrum_charge() || ind_write_charge))
		mgf << "CHARGE=" << this->charge << "+" << endl;
		
	mgf << "PEPMASS=" << this->m_over_z << endl;

	// filter peaks if there are too many

	if (this->basic_spectra.size()>1)
	{
		if (peaks.size() <= maximum_good_peaks_to_output)
		{	
			int i;
			for (i=0; i<peaks.size(); i++)
			{
				const float intensity = (peaks[i].scaled_intensity>0 ? peaks[i].scaled_intensity : peaks[i].intensity);
				if (intensity<min_peak_intensity)
					continue;

				mgf << fixed << setprecision(3) << peaks[i].mass << " " << setprecision(2) << intensity;
					
				if (write_peak_count)
					mgf << "   " << peaks[i].num_occurences << "/" << peaks[i].max_num_occurences;
				mgf << endl;
			}
		}
		else
		{
			vector<bool> good_peak_indicators;
			mark_top_peaks_with_sliding_window(&peaks[0], peaks.size(),
				config->get_local_window_size(),
				config->get_max_number_peaks_per_local_window(),
				good_peak_indicators);

			int i;
			for (i=0; i<peaks.size(); i++)
			{
				if (! good_peak_indicators[i])
					continue;

				const float intensity = (peaks[i].scaled_intensity>0 ? peaks[i].scaled_intensity : peaks[i].intensity);
				if (intensity<min_peak_intensity)
					continue;

				mgf << fixed << setprecision(3) << peaks[i].mass << " " << setprecision(2) << intensity;

				if (write_peak_count)
					mgf << "   " << peaks[i].num_occurences << "/" << peaks[i].max_num_occurences;
				mgf << endl;
			}
		}
	}
	else // output all peaks in the single basic spectrum - no filtering
	{
		const QCPeak *peaks = basic_spectra[0].peaks;
		const int num_peaks = basic_spectra[0].num_peaks;
		int i;
		for (i=0; i<num_peaks; i++)
		{
			if (peaks[i].intensity<min_peak_intensity)
				continue;

			mgf << fixed << setprecision(3) << peaks[i].mass << " " << setprecision(3) << 
					peaks[i].intensity << endl;
		}
	}

	mgf << "END IONS" << endl << endl;
	
}


// writes the spectrum to the output file in the mgf format
void ClusterSpectrum::write_spectrum_to_pkl_single(string& file_name) const
{
	ofstream pkl(file_name.c_str());
	if (! pkl.is_open())
	{
		cout << "Error: couldn't open pkl for writing " << file_name << endl;
		exit(1);
	}

	double total_ion_current =0;
	int i;
	for (i=0; i<basic_spectra.size(); i++)
		total_ion_current += basic_spectra[i].ssf->precursor_intensity;

//	cout << "TIC: " << total_ion_current << endl;

	pkl << setprecision(4) << fixed << m_over_z << "\t";
	pkl << setprecision(0) << total_ion_current << "\t";
	pkl << charge << endl;


	bool use_scaled=false;
	double scale_factor=1.0;
	double peaks_inten=0;
	if (peaks[0].scaled_intensity>0)
	{
		int i;
		for (i=0; i<peaks.size(); i++)
			peaks_inten += peaks[i].scaled_intensity;

		use_scaled = true;
	}
	else
	{
		int i;
		for (i=0; i<peaks.size(); i++)
			peaks_inten += peaks[i].intensity;
	}

	// avoid very large peak intensity numbers
//	if (total_ion_current > 1000000000)
//		total_ion_current = 1000000000; 

	if (peaks_inten>0 && total_ion_current>0)
	{
		scale_factor = total_ion_current / peaks_inten;
	}

	// filter peaks if there are too many
	if (peaks.size() <= maximum_good_peaks_to_output)
	{	
		int i;
		for (i=0; i<peaks.size(); i++)
		{
			pkl << fixed << setprecision(2) << peaks[i].mass << "\t" << setprecision(1) << 
				scale_factor * (peaks[i].scaled_intensity>0 ? peaks[i].scaled_intensity : peaks[i].intensity);
		//	if (write_peak_count)
		//		pkl << "   " << peaks[i].num_occurences << "/" << peaks[i].max_num_occurences;
			pkl << endl;
		}
	}
	else
	{
		vector<bool> good_peak_indicators;
		mark_top_peaks_with_sliding_window(&peaks[0], peaks.size(),
			config->get_local_window_size(),
			config->get_max_number_peaks_per_local_window(),
			good_peak_indicators);

		int i;
		for (i=0; i<peaks.size(); i++)
		{
			if (! good_peak_indicators[i])
				continue;

			pkl << fixed << setprecision(2) << peaks[i].mass << "\t" << setprecision(1) << 
				scale_factor * (peaks[i].scaled_intensity>0 ? peaks[i].scaled_intensity : peaks[i].intensity);

		//	if (write_peak_count)
		//		pkl << "   " << peaks[i].num_occurences << "/" << peaks[i].max_num_occurences;
			pkl << endl;
		}

	}

	pkl.close();

//	cout << "Wrote: " << file_name << endl;
	
}

void ClusterSpectrum::print_explained_intensity_stats(Peptide& pep) const
{
	int i;
	float total_inten=0;
	vector<mass_t> masses;
	vector<intensity_t> intensities;

	masses.resize(peaks.size());
	intensities.resize(peaks.size());
	for (i=0; i<peaks.size(); i++)
	{
		masses[i]=peaks[i].mass;
		intensities[i]=peaks[i].intensity;
	}

	mass_t pm_with_19 = pep.get_mass()+MASS_OHHH;
	AnnotatedSpectrum as;
	as.read_from_peak_arrays(config,pep,pm_with_19,charge,peaks.size(),
		&masses[0],&intensities[0]);
	as.init_spectrum();

	as.annotate_spectrum(pm_with_19);
	as.print_expected_by();
	
	int len = pep.get_num_aas();
			
	// count b,y
	const int b_frag_idx = config->get_frag_idx_from_label("b");
	const int y_frag_idx = config->get_frag_idx_from_label("y");
	int num_b = as.get_num_observed_frags(b_frag_idx);
	int num_y = as.get_num_observed_frags(y_frag_idx);
	float exp_int = as.get_explianed_intensity();

	cout << "exp: " << exp_int << "  b:" << num_b << "  y:" << num_y <<endl;


}


void BasicSpectrum::print_peaks() const
{
	int i;
	for (i=0; i<this->num_peaks; i++)
		cout << left << setw(5) << i << this->peaks[i].mass << "  " << peaks[i].intensity << endl;
}


// Assumes PTMs are already set
void extractAnnoatedScansFromFiles(Config *config, char *file_list, char *anns, 
								   char *output_file, bool extract_no_process)
{
	vector< vector<int> >    annotation_idxs;
	vector<mzXML_annotation> annotations;
	read_mzXML_annotations(file_list,anns,annotation_idxs,annotations,35000);
	
	cout << "Read annotations: " << annotations.size() << endl;

	FileManager fm;
	fm.init_from_list_file_and_add_annotations(config,file_list,
			annotation_idxs, annotations,true);

	FileSet fs;
	fs.select_all_files(fm,true);
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	ofstream mgf_stream(output_file,ios::out);
	
	BasicSpecReader bsr;
	QCPeak peaks[5000];
	int i;
	for (i=0; i<all_ssf.size(); i++)
	{
		BasicSpectrum bs;
		MZXML_single *ssf = (MZXML_single *)all_ssf[i];
	
		bs.peaks = peaks;
		bs.ssf = ssf;
		bs.ssf->charge   = ssf->charge;

		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks,extract_no_process);
		
		if (ssf->scan_number<0)
		{
			cout << "Error: no scan number read from MGF!!!" << endl;
			exit(1);
		}

		bs.output_to_mgf(mgf_stream,config);
	}
	mgf_stream.close();
}





void read_mzXML_annotations_to_map(char *ann_file, 
								   map<mzXML_annotation,int>& ann_map)
{

	int i;
	char buff[256];


	FILE *ann_stream = fopen(ann_file,"r");
	if (! ann_stream)
	{
		cout << "Error: couldn't open annotation files for run!: " << ann_file << endl;
		exit(1);
	}



	int anns=0;
	i=0;
	while (fgets(buff,256,ann_stream))
	{
		int file_idx=-1, mzXML_idx=-1; 
		char only_peptide[128];
		int scan=-1,charge=0;

		if (sscanf(buff,"%d %d %d %d %s",&file_idx,&mzXML_idx,&scan,&charge,only_peptide)<5)
			continue;

		mzXML_annotation ann;
		ann.charge =charge;
		ann.scan = scan;
		ann.mzXML_file_idx = file_idx;
		ann.pep = only_peptide;


		ann_map.insert(make_pair(ann,1));
		anns++;
	}
	cout << "Read " << anns << " annotations..." << endl;
}


// reads an annotation file. mzXML
// saves the annotation for each file number and scan number
// assume max scan number is 30000
// file_idx  mzXML_file_idx  scan peptide
// 169 01 8854 2 VAQGVSGAVQDK
void read_mzXML_annotations(char *mzXML_list, 
							char *ann_file, 
							vector< vector<int> >& annotation_idxs, 
							vector<mzXML_annotation>& annotations,
							int max_ann_size) 
{
	int i;
	char buff[256];
	FILE *mzxml_stream = fopen(mzXML_list,"r");
	if (! mzxml_stream)
	{
		cout << "Error: couldn't open annotation file for mzXML run!: " << mzXML_list << endl;
		exit(1);
	}

	int count=0;
	while (fgets(buff,256,mzxml_stream))
	{
		count++;
	}
	fclose(mzxml_stream);

	annotation_idxs.resize(count+1);
	for (i=0; i<=count; i++)
		annotation_idxs[i].resize(max_ann_size,-1);

	FILE *ann_stream = fopen(ann_file,"r");
	if (! ann_stream)
	{
		cout << "Error: couldn't open annotation files for run!: " << ann_file << endl;
		exit(1);
	}

/*	b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2-19.mzXML'...
20989 spectra...
255 Parse spectra from 'C:/Work/Data/Briggs\H293b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2\H293
b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2-20.mzXML'...
20712 spectra...
256 Parse spectra from 'C:/Work/Data/Briggs\H293b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2\H293
b-total-try-2nd-digest-b-400ug-2D34-121505-LTQ2-21.mzXML'...
20603 spectra...*/

	

	i=0;
	while (fgets(buff,256,ann_stream))
	{
		int file_idx=-1, mzXML_idx=-1; 
		char only_peptide[128];
		int scan=-1,charge=0;

		if (sscanf(buff,"%d %d %d %d %s",&file_idx,&mzXML_idx,&scan,&charge,only_peptide)<5)
			continue;

		mzXML_annotation ann;
		ann.charge =charge;
		ann.mzXML_file_idx = file_idx;
		ann.pep = only_peptide;


		annotation_idxs[file_idx][scan]=annotations.size();
		annotations.push_back(ann);
		
	//	cout << file_idx << " : >> " << scan << "   " << only_peptide << " " << charge << endl;

	//	if (++i>=200)
	//		break;
	
	}
}






void FileManager::init_from_dat_list_extract_only_annotated(Config *config, 
						char* dat_list_file, char *ann_file)
{

	int i;

	// read annotations
	map<mzXML_annotation,int> ann_map;
	read_mzXML_annotations_to_map(ann_file,ann_map);

	vector<string> list;
	read_paths_into_list(dat_list_file,list);
	dat_files.clear();
	for (i=0; i<list.size(); i++)
    {
		if (list[i][0] == '#')
			continue;

		DAT_file dat;
		dat.dat_name =list[i];

		dat.initial_read(config,dat_files.size());

		cout << dat.dat_name << " .. ";

		// change the single spectrum pointers in the mgf file record
		// to include only those that have a mass that is in the permitted range

		vector<DAT_single> good_singles;
		int j;
		for (j=0; j<dat.single_spectra.size(); j++)
		{
			int mzxml_file_idx = dat.single_spectra[j].mzxml_file_idx;
			int scan_number = dat.single_spectra[j].scan_number;

			map<mzXML_annotation,int>::const_iterator it;
			mzXML_annotation ann_pos;

			ann_pos.mzXML_file_idx = mzxml_file_idx;
			ann_pos.scan = scan_number;

			it = ann_map.find(ann_pos);

			if (it == ann_map.end())
				continue;

			Peptide pep;
			int charge = it->first.charge;
			const string& pep_str = it->first.pep;

			pep.parse_from_string(config,pep_str);
			pep.calc_mass(config);
			mass_t m_over_z = (pep.get_mass()+19.0183 + (mass_t)charge)/charge;

			if (fabs(m_over_z - dat.single_spectra[j].m_over_z)>7.0)
			{
				cout << "Error: mismatch between ann " << pep_str << " " << m_over_z << endl;
				cout << "       and dat  " << mzxml_file_idx  << ", " << scan_number << "  " << dat.single_spectra[j].m_over_z << endl;

				cout << "Mass Cys: " << config->get_aa2mass()[Cys] << endl;
				continue;
			}
					
			dat.single_spectra[j].peptide.parse_from_string(config,
					it->first.pep);
		
			good_singles.push_back(dat.single_spectra[j]);
		}

		dat.single_spectra = good_singles;
		dat_files.push_back(dat);

		cout << good_singles.size() << " ..." << endl;

	}

}



void extract_annotate_scans_from_dat(Config *config, char *dat_list, char *anns_file, 
									 char *out_name)
{
	FileManager fm;
	fm.init_from_dat_list_extract_only_annotated(config,dat_list,anns_file);

	FileSet fs;
	fs.select_all_files(fm,true);
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	cout << "Processing " << all_ssf.size() << " spectra..." << endl;


	vector<FILE *> len_streams,size_streams;
	vector<int> len_counts, size_counts;

	const int max_len = 60;
	const int max_size_idx = 65;

	len_streams.resize(max_len+1,NULL);
	size_streams.resize(max_size_idx+1,NULL);

	len_counts.resize(max_len+1,0);
	size_counts.resize(max_size_idx+1,0);

	
	BasicSpecReader bsr;
	QCPeak peaks[5000];
	int i;
	for (i=0; i<all_ssf.size(); i++)
	{
		BasicSpectrum bs;
		DAT_single *ssf = (DAT_single *)all_ssf[i];
	
		bs.peaks = peaks;
		bs.ssf = ssf;
		bs.ssf->charge   = ssf->charge;

		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks,false,true);

		const int len_idx = bs.ssf->peptide.get_num_aas();

		bs.ssf->peptide.calc_mass(config);
		const int size_idx = (int)((bs.ssf->peptide.get_mass()+MASS_OHHH)/100);

		if (len_idx>max_len || size_idx>max_size_idx)
			continue;
		
		if (ssf->scan_number<0)
		{
			cout << "Error: no scan number read from MGF!!!" << endl;
			exit(1);
		}

		char title[32];

		sprintf(title,"%d_%d",ssf->mzxml_file_idx,ssf->scan_number);
		ssf->single_name = string(title);



		ofstream len_stream,size_stream;

		if (! len_streams[len_idx])
		{
			char len_name[256];
			sprintf(len_name,"%s_%d.mgf",out_name,len_idx);
			len_streams[len_idx] = fopen(len_name,"w");
			cout << "open: " << len_name << endl;
		}
	

		if (! size_streams[size_idx])
		{
			char size_name[256];
			sprintf(size_name,"%s_%d.mgf",out_name,size_idx*100);
		
			size_streams[size_idx] = fopen(size_name,"w");

			cout << "open: " << size_name << endl;
		}
	
		bs.output_to_mgf(len_streams[len_idx],config);
		bs.output_to_mgf(size_streams[size_idx],config);

		len_counts[len_idx]++;
		size_counts[size_idx]++;
	}


	for (i=0; i<len_streams.size(); i++)
		if (len_streams[i])
			fclose(len_streams[i]);

	for (i=0; i<size_streams.size(); i++)
		if (size_streams[i])
			fclose(size_streams[i]);

	char len_sum_name[256],size_sum_name[256];
	sprintf(len_sum_name,"%s_len_sum.txt",out_name);
	sprintf(size_sum_name,"%s_size_sum.txt",out_name);

	ofstream len_sum(len_sum_name,ios::out);
	for (i=0; i<=max_len; i++)
		if (len_counts[i]>0)
			len_sum << i << "\t" << len_counts[i] << endl;
	len_sum.close();

	ofstream size_sum(size_sum_name,ios::out);
	for (i=0; i<=max_size_idx; i++)
		if (size_counts[i]>0)
			size_sum << i*100 << "\t" << size_counts[i] << endl;
	size_sum.close();
}



void convert_dat_to_mgf(Config *config, 
						char *dat_list, 
						char *out_name, 
						char *out_dir,
						char *anns_file)
{
	QCOutputter qco_all,qco_labeled;

	qco_all.init(string(out_name) + "_all",out_dir);
	qco_labeled.init(string(out_name) + "_labeled",out_dir);

	map<mzXML_annotation,int> ann_map;
	bool use_map = false;
	if (anns_file)
	{
		read_mzXML_annotations_to_map(anns_file,ann_map);
		if (ann_map.size()>1)
			use_map=true;
	}

	FileManager fm;
	FileSet fs;
	fm.init_from_list_file(config,dat_list);
	fs.select_all_files(fm);
	const vector<SingleSpectrumFile *>& all_ssf = fs.get_ssf_pointers();

	BasicSpecReader bsr;
	QCPeak peaks[5000];

	cout << "Processing " << all_ssf.size() << " spectra..." << endl;
	int num_wrote = 0;
	int num_annotated = 0;
	int bad_anns= 0;
	int i;
	for (i=0; i<all_ssf.size(); i++)
    {
		DAT_single *ssf = (DAT_single *)all_ssf[i];
		const int mzxml_file_idx = ssf->mzxml_file_idx;
		const int scan_number = ssf->scan_number;

		map<mzXML_annotation,int>::const_iterator it;
		mzXML_annotation ann_pos;

		ann_pos.mzXML_file_idx = mzxml_file_idx;
		ann_pos.scan = scan_number;

		it = ann_map.find(ann_pos);

		bool has_ann = false;
		if (it != ann_map.end())
		{
			ssf->peptide.parse_from_string(config,it->first.pep);
			ssf->peptide.calc_mass(config);
			const int charge = it->first.charge;
			mass_t exp_m_over_z = (ssf->peptide.get_mass() + 18.0 + charge)/charge;
			if (fabs(ssf->m_over_z-exp_m_over_z)>6.0)
			{
				bad_anns++;
				has_ann=false;
				ssf->peptide.clear();
				continue;
			}
			{
				num_annotated++;
				has_ann=true;
			}
		}
	
		BasicSpectrum bs;
		
		bs.peaks = peaks;
		bs.ssf = ssf;
		bs.ssf->charge   = ssf->charge;

		bs.num_peaks = bsr.read_basic_spec(config,fm,ssf,peaks,false,true);

		if (bs.num_peaks<10)
			continue;

		char title[32];

		sprintf(title,"%d_%d",ssf->mzxml_file_idx,ssf->scan_number);
		ssf->single_name = string(title);

		qco_all.output_basic_spectrum_to_mgf(bs,config);
		if (it != ann_map.end())
			qco_labeled.output_basic_spectrum_to_mgf(bs,config);
		num_wrote++;
	}

	cout << "Wrote " << num_wrote << " spectra to MGF" << endl;
	cout << num_annotated << " were annotated. " << endl;
	cout << "Found " << bad_anns << " bad anns..." << endl;
}
