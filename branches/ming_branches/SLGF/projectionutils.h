#ifndef PROJECTIONUTIL_H
#define PROJECTIONUTIL_H

#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include <tr1/tuple>

#include "spectrum.h"
#include "ms1.h"
#include "SpectralLibrary.h"
#include "utils.h"


using namespace specnets;

string getAdditionalMassString(int parent_mass_difference);
string cleanAnnotationEnds(string annotation);
string create_annotation_ends(string annotation);
string remove_annotation_ends(string annotation);
int getpeptideLength(string annotation);

map<string, float> getAminoAcidSubstitutionLookup();
string reverse_string(string input);

map<char, float> getAminoAcidLookup();
float getSubstitutionCosineDepression(char a, char b, map<string, float> amino_acid_transform_map);
float getMass(string peptide, map<char, float> amino_acid_map);
int getDifferenceIndex(string in1, string in2);
int getStringDifference(string in1, string in2);
float cosine(vector<float> &u, vector<float> &v);
void preprocess_spectrum_intensities_max_intensity(Spectrum * s, float max);
void preprocess_spectrum_intensities(Spectrum * s, int do_early_normalize, int do_filter);
void normalize_extracted_ions(vector<float> &v);
void extractIons(psmPtr psm, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, vector<float> &ions, int do_early_normalize, int do_filter);


void extractIons(psmPtr psm, int peptideLength, MS2ScoringModel &model,
        vector<string> &ionsToExtract, vector<pair <float, float> > &ions, int do_early_normalize, int do_filter);
        
void deleteIons(psmPtr psm, int peptideLength, MS2ScoringModel &model, vector<string> &ionsToExtract);
        
void norm_vector(vector<float> & input);
        
float spectrum_similarity(psmPtr psm1, 
                          string annotation1, 
                          psmPtr psm2, 
                          string annotation2, 
                          int peptideLength, 
                          MS2ScoringModel &model, 
                          vector<string> &ionsToExtract, 
                          string allIons);
                          
float spectrum_similarity(psmPtr library, 
                          psmPtr query, 
                          int peptideLength, 
                          MS2ScoringModel &model, 
                          vector<string> &ionsToExtract, 
                          string allIons);
                          
//Subtracting off average spectrum
float spectrum_similarity(psmPtr psm1,
                          psmPtr psm2,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons, 
                          vector<vector<float> > average_ions);
                          
float spectrum_similarity_sqrt_librarypeaks_isocombine(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons, 
                          float &percent_explained_intensity);
                          
float spectrum_similarity_sqrt_librarypeaks(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons,
                          float &percent_explained_intensity);
                          
float spectrum_similarity_sqrt_librarypeaks_speed_optimization(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons,
                          float &percent_explained_intensity);

float spectrum_similarity_sqrt_librarypeaks_isocombine_speed_optimization(psmPtr library,
                          psmPtr query,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons, 
                          float &percent_explained_intensity);
                          

                          
float spectrum_similarity_nosqrt(psmPtr psm1, 
                          psmPtr psm2, 
                          int peptideLength, 
                          MS2ScoringModel &model, 
                          vector<string> &ionsToExtract, 
                          string allIons);
                          
                          
float preprocess_library_ion_extraction(psmPtr library,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons);

float preprocess_library_ion_extraction_isocombine(psmPtr library,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons);
                          
float full_spectrum_similarity(Spectrum spec1, Spectrum spec2);
float full_spectrum_similarity_nosqrt(Spectrum spec1, Spectrum spec2);
bool search_results_comparator (score_results_tuple i, score_results_tuple j);

float full_spectrum_dotbias(Spectrum spec1, Spectrum spec2, float spec_sim);

void sorted_vector_library(vector<Spectrum *> &library_ptr, SpectralLibrary & real_library);

bool spectrum_ptr_compare (Spectrum* first, Spectrum* second);

int spectrum_ptr_startend(vector<Spectrum *> &library_ptr, float parentMZ, float tolerance, int &start_idx, int &end_idx);

int hashstring(string input);
int hashpeptide(string input, int charge);

/**
 * Model = 1 (Gaussian with zero'd means)
 * Model = 2 (Gaussian with real means)
 * Model = 3 (Empirical Ion Variance Distributions with mean offset)
 * Model = 4 (Empirical Ion Variance Distributions without mean offset)
**/
void dynamicprogramming_spectrum(Spectrum &spectrum, MS2ScoringModel model, vector<string> ions_to_extract, string allIons, vector<float> ion_mean, vector<float> ion_variance, vector< pair < float, float > > insertion_probability , vector< pair < float, float > > deletion_probability, vector< pair< float, float> > &output_distribution, vector<vector<pair<float, float> > > histograms, int training_model);
void dynamicprogramming_ionvector(vector<float> ion_intensities,
                                  vector<float> ion_mean, 
                                  vector<float> ion_variance, 
                                  vector< pair < float, float > > insertion_probability , 
                                  vector< pair < float, float > > deletion_probability, 
                                  vector< pair< float, float> > &output_distribution, 
                                  vector<vector<pair<float, float> > > histograms, 
                                  int training_model,
                                  int cosine_dimension,
                                  int intensity_dimension);
                                 
void montecarlo_spectrum(Spectrum &spectrum, MS2ScoringModel model, vector<string> ions_to_extract, string allIons, vector<float> ion_mean, vector<float> ion_variance, vector< pair < float, float > > insertion_probability , vector< pair < float, float > > deletion_probability, vector< pair< float, float> > &output_distribution, vector<vector<pair<float, float> > > histograms, int training_model);
float max_ion_intensity(vector<float> ion_intensity);
vector<pair<float, float> > create_histogram(int buckets, float start_range, float end_range, vector<float> input, bool normalize);

vector<pair<float, float> > create_histogram2(int buckets, float start_range, float end_range, vector<float> input, bool normalize);

inline double get_3d_array(int x_dim, int y_dim, int z_dim, int x_addr, int y_addr, int z_addr, double * ptr);
inline int set_3d_array(int x_dim, int y_dim, int z_dim, int x_addr, int y_addr, int z_addr, double * ptr, double value);


//SLGF Training
void save_cosine_distribution(vector<vector<pair<float, float> > > histograms, string filename);
void load_cosine_distribution(vector<vector<pair<float, float> > > & histograms, string filename);
void save_deletion_distribution(vector < pair < float, float > > &deletion_prob, string filename);
void load_deletion_distribution(vector < pair < float, float > > &deletion_prob, string filename);

void sqrt_vector(vector<float> & input);

float SLGF_rescore(vector<pair < float, float> > &cosine_distribution, float cosine);


vector<string> create_deliminated_aminoacids(string input);


float percent_explained_intensity(psmPtr psm,
                                    int peptideLength, 
                                    MS2ScoringModel &model, 
                                    vector<string> &ionsToExtract, 
                                    string allIons);
                                    
float percent_sequence_similarity(string peptide1, string peptide2);

double ion_variance_cdf(vector<vector<pair<float, float> > > &histograms, int bandidx, float x, float mean_adjustment);
double ion_variance_integration(vector<vector<pair<float, float> > > &histograms, int bandidx, double left, double right, double mean_adjustment);

double random_ion_variance(vector<vector<pair<float, float> > > &histograms, int bandidx);

vector<float> sort_descending(vector<float> input);
vector<float> reverse_vector(vector<float> input);

void get_band_meanvariance(vector<float> &band_mean, vector<float> &band_variance, vector<vector<pair<float, float> > > histograms);

void combine_isotopic(vector<string> ions, vector<float> &extracted_ions, int length);

int spectrum_signal_peak_count(MS2ScoringModel &model, vector<string> &ionsToExtract, string allIons, psmPtr psm, float percent_minimum);

void save_SLGF(SpectralLibrary lib, string filename);
void load_SLGF(SpectralLibrary &lib, string filename);

void printIon(vector<string> &ionsToExtract, vector<float> ions);
void printDoubleIon(vector<string> &ionsToExtract, vector<float> library, vector<float> query);

#endif
