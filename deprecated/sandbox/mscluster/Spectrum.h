#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__


#include "Config.h"
#include "BasicDataStructs.h"
#include "includes.h"


// Peak structure
struct Peak {
	Peak() : mass(-1),  intensity(0), log_intensity(NEG_INF), rank(-1), iso_level(0), 
			 log_local_rank(NEG_INF), log_rand_prob(NEG_INF), multiplicity(0) {};

	bool operator< ( const Peak& other) const
	{
		return mass<other.mass;
	}

	mass_t       mass;      // in Daltons

	intensity_t  intensity;

	float		 log_intensity;

	int    rank;		   // amongst all peaks (including iso)	

	float  iso_level;     // if >0, means this might be isotopic peak
					      // the higher the value, the more likely this is an isotopic peak

	float  log_local_rank;

	float  log_rand_prob; // the porbability of observing this peak at random (computed according to the neighbor's density)

	float  multiplicity;  // how many spectra have this peak (for clusters)
};



/////////////////////////////////////////////////////////////////
// This class holds the common characteristics of a spectrum
// A spectrum is an entity whose charge is known!
class Spectrum {
public:
	Spectrum() : org_pm_with_19(-1),       m_over_z(-1),
				 corrected_pm_with_19(-1), 
			     secondary_pm_with_19(-1), corrected_pm_score(-1), 
				 secondary_pm_score(-1),   total_original_intensity(-1),
				 max_mass(-1),			   min_peak_mass(-1), 
				 max_peak_mass(-1),		   charge(-1),
				 scan_number(-1),		   retention_time(-1),
				 cluster_size(1),		   size_idx(-1),             
				 config(NULL),			   file_name("") {};


	// This function will be re-implemented for cluster and paired spectra
	// were there will be data structures to support 
	virtual float       get_peak_percent(int p_idx =-1, int charge = -1, int orientation = -1) const
	{
		return 1.0;
	}

	mass_t      get_true_mass() const         { return  peptide.get_mass(); }
	mass_t      get_true_mass_with_19() const { return  peptide.get_mass() + MASS_OHHH; }
	mass_t      get_org_pm_with_19() const { return  org_pm_with_19; }
	void        set_org_pm_with_19( mass_t pm) { org_pm_with_19 = pm; }

	mass_t	    get_m_over_z() const { return m_over_z; }
	void        set_m_over_z(mass_t mz) { m_over_z = mz; }

	void        calc_original_pm_with_19_for_different_charges(vector<mass_t>& pms_with_19) const;
	mass_t      get_corrected_pm_with_19() const { return corrected_pm_with_19; }
	mass_t      get_secondary_pm_with_19() const { return secondary_pm_with_19; }
	score_t     get_corrected_pm_score() const { return corrected_pm_score; }
	score_t     get_secondary_pm_score() const { return secondary_pm_score; }

	void        set_corrected_pm_with_19(mass_t pm) { corrected_pm_with_19 = pm; }
	void        set_secondary_pm_with_19(mass_t pm) { secondary_pm_with_19 = pm; }
	void        set_corrected_pm_score(score_t score) { corrected_pm_score = score; }
	void        set_secondary_pm_score(score_t score) { secondary_pm_score = score; }

	intensity_t get_total_original_intensity() const { return total_original_intensity; }
	int         get_size_idx() const { return size_idx; }
	void        set_size_idx(int s) { size_idx = s; }
	int         get_charge() const { return charge; }
	void        set_charge(int c) { charge = c; }
	int         get_scan_number() const { return scan_number; }
	int         get_cluster_size() const { return cluster_size; }
	Config*     get_config() const { return config; }
	const string&  get_file_name() const { return file_name; }
	void set_file_name(const string& fname) { file_name = fname; }
	const Peptide& get_peptide()   const { return peptide; }
	void  set_peptide(const Peptide& p) { peptide = p; }
	void  set_peptide(const string& pep_string) { peptide.parse_from_string(config,pep_string); }
	void  convert_peptide_ILQK() { peptide.convert_ILQK(config); }

	mass_t      get_min_peak_mass() const { return min_peak_mass; }
	mass_t      get_max_peak_mass() const { return max_peak_mass; }
	mass_t      get_peak_mass(int peak_idx) const { return peaks[peak_idx].mass; }
	intensity_t get_peak_intensity(int peak_idx) const { return peaks[peak_idx].intensity; }
	int         get_peak_rank(int peak_idx) const { return peaks[peak_idx].rank; }
	float       get_peak_iso_level(int peak_idx) const { return peaks[peak_idx].iso_level; }

	int         get_num_peaks()     const { return peaks.size(); }
	const Peak& get_peak(int p_idx) const { return peaks[p_idx]; }
	const vector<Peak>& get_peaks() const { return peaks; }
	const vector<int>& get_strong_peak_idxs() const { return strong_peak_idxs; }
	const vector<int>& get_ranked_mono_peaks() const { return ranked_mono_peaks; }

	PeakRange get_peaks_in_range(mass_t min_r, mass_t max_r) const;
	int get_min_dis_peak(mass_t mass, mass_t tolerance) const;
	int get_max_inten_peak(mass_t mass, mass_t tolerance) const;

	void get_iso_intens(int p_idx, vector<float>& iso_intens, mass_t iso_tolerance, int charge) const;

	bool read_from_dta(Config* _config, const char *dta_name, 
					   int charge = -1, bool verbose =false);

	bool read_from_MGF_stream(Config*_config, FILE *stream, int charge = 0, 
					   bool verbose = false);

	bool read_from_peak_arrays(Config* _config, Peptide& _pep, mass_t _pm_with_19, 
		int _charge, int _num_peaks, mass_t *_masses, intensity_t *_intensities);

	bool init_from_QCPeaks(Config* _config, void *QCPeaks, int num_peaks, void *ssf,
			bool ind_perfrom_filtering = true);

	void output_as_MGF(ostream& os) const;

	
	void init_spectrum(bool perform_filtering = true);

	void create_random_spec_with_1_good_peak();

	void print_spectrum(ostream& os = cout) const;

	void print_expected_by(ostream& os = cout) const;

	void print_expected_fragment_peaks(vector<string>& frag_labels, ostream& os = cout) const;

	bool check_m_over_z_and_sequence(ostream& os = cout);

	void mark_all_non_iso_peaks_as_strong();

protected:
	mass_t  org_pm_with_19, m_over_z;
	mass_t  corrected_pm_with_19, secondary_pm_with_19;
	score_t corrected_pm_score, secondary_pm_score;
	intensity_t  total_original_intensity; // sum of all peak intensities
	mass_t  max_mass;   // the maximal mass that can be aniticipated
	mass_t  min_peak_mass, max_peak_mass;
	int charge;
	int scan_number;        // if originated from mzXML
	float  retention_time;  // if originated from mzXML
	int cluster_size;      // if originated from a cluster
	int size_idx;
	Config *config; // The config used to read this Spectrum
	string file_name;
	Peptide peptide;
	vector<int> index_array;  // used to quickly locate peaks

	vector<int> strong_peak_idxs; // holds the indexes of peaks that should be considered strong
							     // these generally are the top peaks in the spectrum that are
							    // not isotopic peaks, and some peaks that are locally strong

	vector<int> ranked_mono_peaks; // holds the idxs of peaks sorted according to rank
								   // peaks that have an isotopic score appear after the ones
								   // without such a score. This vector is filled automatically
								   // after isotopic levels are calculated

	vector<Peak> peaks;


	// functions

	void init_index_array();     // creates the fast lookup array for the peaks
	void join_adjacent_peaks();  // joins peaks that are close to each other
	void filter_peaks();         // removes the locally weak peaks
	void normalize_intensities();
	void calc_ranks();
	void calc_log_local_ranks(); // calcs the logarithm of the local rank for each
	void calc_isotope_levels();  // calcs for each peak the level it resembles an isotopic peak
	void select_strong_peaks();
	void set_log_random_probs(); // calcs for each peak the probability of observing it at random (based on
								 // the neighbor's distribution
};


ostream& operator << (ostream& os, const Spectrum& spec);

/*****************************************************************
Returns the indices of all peaks that are within the mass range
******************************************************************/
inline PeakRange  Spectrum::get_peaks_in_range(mass_t min_r, mass_t max_r) const
{
	PeakRange pr;
	const int max_peak_idx = peaks.size()-1;

	if (max_r>= max_mass+0.1)
		max_r = max_mass+0.1;

	if (min_r<0)
		min_r=0;

	if (max_r<=min_r)
		return pr;

	int i_min=index_array[(int)(min_r)];
	int i_max=index_array[(int)(max_r)];
		
	if (i_max<max_peak_idx)
		i_max++;

	while (peaks[i_min].mass<min_r && i_min<i_max)
		i_min++;

	while (peaks[i_max].mass>max_r && i_max>0)
		i_max--;

	if (peaks[i_min].mass > max_r || peaks[i_max].mass<min_r)
		return pr;

	pr.num_peaks = i_max - i_min+1;
	pr.low_idx = i_min;
	pr.high_idx = i_max;

	return pr;
}


/***********************************************************
returns the idx of the peak that is closes to the mass
(-1 is returned if no peak is found)
************************************************************/
inline int Spectrum::get_min_dis_peak(mass_t mass, mass_t tolerance) const
{
	PeakRange pr = get_peaks_in_range(mass-tolerance,mass+tolerance);
	if (pr.num_peaks==0)
		return -1;
	mass_t min_dis=fabs(mass - peaks[pr.low_idx].mass);
	int min_idx=pr.low_idx;
	int j;

	// find closest peak to exp_pos
	for (j=1; j<pr.num_peaks; j++)
	{
		mass_t dis = fabs(mass - peaks[pr.low_idx+j].mass);
		if (dis<min_dis)
		{
			min_dis=dis;
			min_idx = pr.low_idx+j;
		}
	}
	return min_idx;
}

/***********************************************************
returns the idx of the peak that has the highest intensity
in the range (-1 is returned if no peak is found)
Balances between intensity and proximity to the expected position
A peak at the edge needs to be at least 2 times stronger than
a peak exactly at the middle to be chosen.
************************************************************/
inline int Spectrum::get_max_inten_peak(mass_t exp_mass, mass_t tolerance) const
{
	PeakRange pr = get_peaks_in_range(exp_mass-tolerance,exp_mass+tolerance);
	if (pr.num_peaks==0)
		return -1;

	if (pr.num_peaks ==1)
		return pr.low_idx;

	mass_t max_val = (2*tolerance - fabs(peaks[pr.low_idx].mass - exp_mass)) *
					  peaks[pr.low_idx].intensity;

	int max_idx=pr.low_idx;
	int j;
	// find closest peak to exp_pos
	for (j=1; j<pr.num_peaks; j++)
	{
		mass_t peak_val = (2*tolerance - fabs(peaks[pr.low_idx].mass - exp_mass)) *
					  peaks[pr.low_idx].intensity;
		if (peak_val>max_val)
		{
			max_val=peak_val;
			max_idx = pr.low_idx+j;
		}
		else
			break;
	}
	return max_idx;
}



// returns intensities at positions -2,-1,+1,+2
inline void Spectrum::get_iso_intens(int p_idx, vector<float>& iso_intens, mass_t iso_tolerance, int charge) const
{
	const mass_t one_over = 1.0/charge;
	const mass_t one_proton = MASS_PROTON*one_over;
	const mass_t two_proton = 2*one_proton;
	const mass_t p_mass=peaks[p_idx].mass;
	const mass_t min_minus2_mass =  p_mass - two_proton - iso_tolerance;
	const mass_t min_minus1_mass =  p_mass - one_proton - iso_tolerance;
	const mass_t min_plus1_mass  =  p_mass + one_proton - iso_tolerance;
	const mass_t min_plus2_mass  =  p_mass + two_proton - iso_tolerance;
	const int num_peaks = peaks.size();

	iso_intens.clear();
	iso_intens.resize(4,NEG_INF);

	int p=p_idx;
	while (p>=0 && peaks[p].mass>min_minus2_mass)
		p--;

	p++;

	if (fabs(peaks[p].mass- p_mass + two_proton)<=iso_tolerance)
		iso_intens[0]=peaks[p].log_intensity;

	while (p<p_idx && peaks[p].mass<min_minus1_mass)
		p++;

	if (fabs(peaks[p].mass - p_mass + one_proton)<=iso_tolerance)
		iso_intens[1]=peaks[p].log_intensity;

	p=p_idx+1;

	while (p<num_peaks && peaks[p].mass<min_plus1_mass)
		p++;

	if (p<num_peaks && fabs(peaks[p].mass - p_mass - one_proton)<=iso_tolerance)
		iso_intens[2]=peaks[p].log_intensity;

	while (p<num_peaks && peaks[p].mass<min_plus2_mass)
		p++;

	if (p<num_peaks && fabs(peaks[p].mass - p_mass - two_proton)<=iso_tolerance)
		iso_intens[3]=peaks[p].log_intensity;
}


#endif




