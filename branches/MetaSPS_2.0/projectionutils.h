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
        
void norm_vector(vector<float> & input);
        
float spectrum_similarity(psmPtr psm1, 
                          string annotation1, 
                          psmPtr psm2, 
                          string annotation2, 
                          int peptideLength, 
                          MS2ScoringModel &model, 
                          vector<string> &ionsToExtract, 
                          string allIons);
                          
float spectrum_similarity(psmPtr psm1, 
                          psmPtr psm2, 
                          int peptideLength, 
                          MS2ScoringModel &model, 
                          vector<string> &ionsToExtract, 
                          string allIons);
                          
float spectrum_similarity(psmPtr psm1,
                          psmPtr psm2,
                          int peptideLength, 
                          MS2ScoringModel &model,
                          vector<string> &ionsToExtract, 
                          string allIons, 
                          vector<vector<float> > average_ions);
                          
float full_spectrum_similarity(Spectrum spec1, Spectrum spec2);

bool search_results_comparator (score_results_tuple i, score_results_tuple j);

float full_spectrum_dotbias(Spectrum spec1, Spectrum spec2, float spec_sim);

void sorted_vector_library(vector<Spectrum *> &library_ptr, SpectralLibrary & real_library);

bool spectrum_ptr_compare (Spectrum* first, Spectrum* second);

int spectrum_ptr_startend(vector<Spectrum *> &library_ptr, float parentMZ, float tolerance, int &start_idx, int &end_idx);

#endif
