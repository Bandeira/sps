#ifndef SPECTRAL_LIBRARY_H
#define SPECTRAL_LIBRARY_H

#include "SpecSet.h"
#include "PeptideSpectrumMatch.h"
#include "mzxml.h"

#include <vector>
#include <tr1/tuple>
#include <string>

namespace specnets
{
    typedef std::tr1::shared_ptr<PeptideSpectrumMatch> psmPtr;

    typedef tr1::tuple<Spectrum *, float, int, std::string, float> score_results_tuple;

    class SpectralLibrary : public SpecSet
    {
        public:
            enum MatchScoreType {
                MatchScoreType_DotProduct = 0, MatchScoreType_DotProduct_DeltaD = 1, MatchScoreType_DotProduct_DeltaD_DotBias = 2
            };
            
            SpectralLibrary(unsigned int sz = 0){
                specs.resize(sz);
            }

            int createlibrary(  float envelope_score_filter,
                                float pvalue_filter,
                                MS2ScoringModel &model,
                                vector<string> &ionsToExtract,
                                string allIons,
                                string aminoacidexclusions,
                                vector<int> charge_filter,
                                bool filter_best_rep);

            int projection( string target_annotation,
                            MS2ScoringModel model,
                            vector<string> ions_to_extract,
                            string allIons,
                            Spectrum & outputspectrum);
            virtual int search(Spectrum query_spec, psmPtr output_psm, float parentmass_tolerance);
            virtual int search_target_decoy(SpectralLibrary &decoy, 
                                            Spectrum query_spec, 
                                            psmPtr output_psm, 
                                            float parentmz_tolerance, 
                                            vector<Spectrum *> target_library_ptr, 
                                            vector<Spectrum *> decoy_library_ptr,
                                            int scoring_method);
                                            
            unsigned int LoadSpecSet_additionalmgf(const char * filename);

            vector<vector<float> > find_global_average_spectrum(MS2ScoringModel model,
                                                                vector<string> ions_to_extract,
                                                                string allIons);

            SpectralLibrary create_decoy_spectral_library(MS2ScoringModel model, vector<string> ions_to_extract, string allIons);

        protected:
            void filter_no_psm();
            void filter_multiple_interpretations();
            void filter_envelope_filter(float envelope_score_filter);
            void filter_pvalue(float pvalue_filter);
            void filter_best_representative(MS2ScoringModel &model,
                                            vector<string> &ionsToExtract,
                                            string allIons);
            void filter_aminoacids(string aminoacidexclusions);
            void filter_forcharge(vector<int> charge);
            void filter_forfragmentation(vector<int> accepted_fragmentation);


            void get_consensus( vector<Spectrum *> spectrum_group,
                                MS2ScoringModel model,
                                vector<string> ionsToExtract,
                                string allIons,
                                Spectrum &consensus_spectrum);

            int get_possible_projections(string target_annotation,
                                         MS2ScoringModel model,
                                         vector<string> ions_to_extract,
                                         string allIons, SpecSet &specs,
                                         vector<float> &projections_set_cosine_depression);
                                         
            string create_decoy_peptide(string peptide);
            
            void generate_prefix_suffix_peptide(string peptide, vector<string> &prefix, vector<string> &suffix);

    };
}

#endif
