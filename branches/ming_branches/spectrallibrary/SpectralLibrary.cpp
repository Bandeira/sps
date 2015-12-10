#include "PeptideSpectrumMatchSet.h"
#include "SpectralLibrary.h"
#include "projectionutils.h"
#include "spectrum_window_filter.h"
#include "spectrum.h"

//System Includes
#include <map>
#include <vector>
#include <algorithm>

namespace specnets
{
    int SpectralLibrary::createlibrary( float envelope_score_filter, 
                                        float pvalue_filter,
                                        MS2ScoringModel &model, 
                                        vector<string> &ionsToExtract, 
                                        string allIons,
                                        string aminoacidexclusions, 
                                        vector<int> charge_filter, 
                                        bool filter_best_rep){
        
        cout<<specs.size()<<endl;
        
        //Creating annotation ends
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                string annot = (*it)->m_annotation;
                string new_annot = create_annotation_ends(annot);
                (*it)->m_annotation = new_annot;
            }
        }
        
        cout<<"After Creating Annotation Ends: "<<specs.size()<<endl;
        
        filter_no_psm();
        cout<<"Filtering No PSM"<<specs.size()<<endl;
        
        if(aminoacidexclusions.length() > 0){    
            filter_aminoacids(aminoacidexclusions);
            filter_no_psm();
        }
        cout<<"Filter Amino Acids: "<<specs.size()<<endl;
        
        if(envelope_score_filter > -1){
            filter_envelope_filter(envelope_score_filter);
            filter_no_psm();
        }
        cout<<"Envelope Filter: "<<specs.size()<<endl;
        
        if(pvalue_filter > -1){
            filter_pvalue(pvalue_filter);
            filter_no_psm();
        }
        cout<<"P Value Filter: "<<specs.size()<<endl;
        
        filter_multiple_interpretations();
        filter_no_psm();
        cout<<"Filtering Multiple Interpretations: "<<specs.size()<<endl;
        
        filter_forcharge(charge_filter);
        cout<<"Charge Filter: "<<specs.size()<<endl;
        
        vector<int> accepted_fragmentation;
        accepted_fragmentation.push_back(Spectrum::FragType_CID);
        filter_forfragmentation(accepted_fragmentation);
        cout<<"Fragmentation Filter: "<<specs.size()<<endl;
        
        if(filter_best_rep)
            filter_best_representative(model, ionsToExtract, allIons);
        cout<<"After Best Rep: "<<specs.size()<<endl;
        
        return 0;
    }
    
    int SpectralLibrary::projection(string target_annotation, 
                                    MS2ScoringModel model, 
                                    vector<string> ions_to_extract, 
                                    string allIons, 
                                    Spectrum & outputspectrum){
        
        //Making sure the annotation looks nice
        target_annotation = create_annotation_ends(target_annotation);
        
        SpecSet projection_spectra;
        vector<float> projections_set_cosine_depression;
        get_possible_projections(target_annotation, 
                                 model, 
                                 ions_to_extract, 
                                 allIons, 
                                 projection_spectra, 
                                 projections_set_cosine_depression);
                                 
        if(projection_spectra.size() == 0) return -1;
        
        vector<Spectrum *> spectrum_group;
        for(int spectrum_idx = 0; spectrum_idx < projection_spectra.size(); spectrum_idx++) 
            spectrum_group.push_back(&projection_spectra[spectrum_idx]);
        
        get_consensus(spectrum_group, model, ions_to_extract, allIons, outputspectrum);
        
        outputspectrum.psmList.clear();
        psmPtr psm(new PeptideSpectrumMatch());
        outputspectrum.psmList.push_back(psm);
        psm->m_annotation = target_annotation;
        psm->m_spectrum = &outputspectrum;
        
        
        return 0;
    }
    
    
    int SpectralLibrary::search_target_decoy(SpectralLibrary &decoy, 
                                             Spectrum query_spec, 
                                             psmPtr output_psm, 
                                             float parentmz_tolerance,
                                             vector<Spectrum *> target_library_ptr, 
                                             vector<Spectrum *> decoy_library_ptr, 
                                             int scoring_method){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;
        
        
        PeptideSpectrumMatchSet search_results_decoy;
        vector<score_results_tuple> scores_tuple_decoy;
        //psmPtr psm(new PeptideSpectrumMatch);
        
        int decoy_start_search_idx;
        int decoy_end_search_idx;
        spectrum_ptr_startend(decoy_library_ptr, query_spec.parentMZ, parentmz_tolerance, decoy_start_search_idx, decoy_end_search_idx);

        cout<<decoy_start_search_idx<<"\t"<<decoy_end_search_idx<<endl;
        for(int library_idx = decoy_start_search_idx; library_idx <= decoy_end_search_idx; library_idx++){

            float library_mass = decoy_library_ptr[library_idx]->parentMZ;
            charge = decoy_library_ptr[library_idx]->parentCharge;
            
            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;
            
            float sim = full_spectrum_similarity(*decoy_library_ptr[library_idx], query_spec);
            
            float dot_bias = full_spectrum_dotbias(*decoy_library_ptr[library_idx], query_spec, sim);
            
            score_results_tuple similarity_tuple;
            decoy.specs[library_idx].scan = decoy_library_ptr[library_idx]->scan;
            tr1::get<0>(similarity_tuple) = decoy_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = library_idx;
            tr1::get<3>(similarity_tuple) = (string)decoy_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            scores_tuple_decoy.push_back(similarity_tuple);
        }

        
        
        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple_decoy.begin(), scores_tuple_decoy.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple_decoy.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum->scan = psm->m_spectrum->scan;
            psm->m_spectrum = tr1::get<0>(scores_tuple_decoy[i]);
            psm->m_score = tr1::get<1>(scores_tuple_decoy[i]);
            psm->m_scanNum = tr1::get<2>(scores_tuple_decoy[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple_decoy[i]);
            search_results_decoy.push_back(psm);
        }
        

        
        
        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;
        
        int target_start_search_idx;
        int target_end_search_idx;
        spectrum_ptr_startend(target_library_ptr, query_spec.parentMZ, parentmz_tolerance, target_start_search_idx, target_end_search_idx);
        
        
        for(int library_idx = target_start_search_idx; library_idx < target_end_search_idx; library_idx++){
            
            float library_mass = target_library_ptr[library_idx]->parentMZ;
            charge = target_library_ptr[library_idx]->parentCharge;
            
            if(abs(query_spec.parentMZ - library_mass) > parentmz_tolerance) continue;
            if(query_spec.parentCharge != 0 && query_spec.parentCharge != charge) continue;
            float sim = full_spectrum_similarity(*target_library_ptr[library_idx], query_spec);
            
            float dot_bias = full_spectrum_dotbias(*target_library_ptr[library_idx], query_spec, sim);
            
            score_results_tuple similarity_tuple;
            
            tr1::get<0>(similarity_tuple) = target_library_ptr[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = library_idx;
            tr1::get<3>(similarity_tuple) = (string)target_library_ptr[library_idx]->psmList.front()->m_annotation;
            tr1::get<4>(similarity_tuple) = dot_bias;
            scores_tuple.push_back(similarity_tuple);
        }
        
        
        //Finding the best
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_scanNum = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            search_results.push_back(psm);
        }
        
        float target_top_scoring = 0.f;
        float target_second_scoring = 0.f;
        float decoy_top_scoring = 0.f;
        float decoy_second_scoring = 0.f;
        
        if(search_results.size() > 0){
            target_top_scoring = search_results[0]->m_score;
        }
        
        if(search_results.size() > 1){
            target_second_scoring = search_results[1]->m_score;
        }
        
        if(search_results_decoy.size() > 0){
            decoy_top_scoring = search_results_decoy[0]->m_score;
        }
        
        if(search_results_decoy.size() > 1){
            decoy_second_scoring = search_results_decoy[1]->m_score;
        }
        
        
        if(target_top_scoring > decoy_top_scoring && search_results.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple[0]);
            float dot_product = target_top_scoring;
            float deltaD = target_top_scoring - max(target_second_scoring, decoy_top_scoring);
            float match_score = target_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            
            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                default:
                    break;
            }
            
            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_scanNum       = search_results[0]->m_scanNum;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_isDecoy = false;
            
            return 0;
        }
        if(target_top_scoring < decoy_top_scoring && search_results_decoy.size() > 0){
            float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
            float dot_product = decoy_top_scoring;
            float deltaD = decoy_top_scoring - max(target_top_scoring, decoy_second_scoring);
            float match_score = decoy_top_scoring * 0.6 + deltaD * 0.4 - dot_bias;
            
            switch(scoring_method){
                case SpectralLibrary::MatchScoreType_DotProduct:
                    match_score = dot_product;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD:
                    match_score = dot_product * 0.6 + deltaD * 0.4;
                    break;
                case SpectralLibrary::MatchScoreType_DotProduct_DeltaD_DotBias:
                    match_score = dot_product * 0.6 + deltaD * 0.4 - dot_bias;
                    break;
                default:
                    break;
            }
            
            output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
            output_psm->m_score         = search_results_decoy[0]->m_score;
            output_psm->m_scanNum       = search_results_decoy[0]->m_scanNum;
            output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
            output_psm->m_score         = match_score;
            output_psm->m_isDecoy = true;
            
            return 0;
        }
        
        return -1;
        
        
        
        
        
        /*
        if(search_results.size() > 1){
            
            if(search_results[0]->m_score > search_results_decoy[0]->m_score){
                
                //float match_score = search_results[0]->m_score;
                float dot_bias = tr1::get<4>(scores_tuple[0]);
                float match_score = search_results[0]->m_score;
                cout<<"Dot:\t"<<match_score<<"\t";
                match_score = match_score * 0.6;
                float delta_D = search_results[0]->m_score;
                if(search_results.size() == 1) delta_D = (delta_D-0)/search_results[0]->m_score;
                else delta_D = (delta_D-(max(search_results[1]->m_score,search_results_decoy[0]->m_score)))/search_results[0]->m_score;
                
                //checking if nan
                if((delta_D != delta_D)){
                    delta_D = 0.f;
                }

                cout<<"DeltaD:\t"<<delta_D<<"\t";
                cout<<"DotBias:\t"<<dot_bias<<"\t";
                match_score += 0.4*delta_D - dot_bias;
                cout<<"MatchScore\t"<<match_score<<"\t";
                
                
                output_psm->m_spectrum      = search_results[0]->m_spectrum;
                output_psm->m_score         = search_results[0]->m_score;
                output_psm->m_scanNum       = search_results[0]->m_scanNum;
                output_psm->m_annotation    = search_results[0]->m_annotation;
                output_psm->m_score         = match_score;
                output_psm->m_isDecoy = false;
                return 0;
            }
            else{
                //float match_score = search_results[0]->m_score;
                float dot_bias = tr1::get<4>(scores_tuple_decoy[0]);
                float match_score = search_results[0]->m_score;
                cout<<"Dot:\t"<<match_score<<"\t";
                match_score = match_score * 0.6;
                float delta_D = search_results_decoy[0]->m_score;
                if(search_results_decoy.size() == 1) delta_D = (delta_D-0)/search_results_decoy[0]->m_score;
                else delta_D = (delta_D-(max(search_results_decoy[1]->m_score, search_results[0]->m_score)))/search_results_decoy[0]->m_score;
                
                //checking if nan
                if((delta_D != delta_D)){
                    delta_D = 0.f;
                }

                cout<<"DeltaD:\t"<<delta_D<<"\t";
                cout<<"DotBias:\t"<<dot_bias<<"\t";
                match_score += 0.4*delta_D - dot_bias;
                cout<<"MatchScore\t"<<match_score<<"\t";
                
                output_psm->m_spectrum      = search_results_decoy[0]->m_spectrum;
                output_psm->m_score         = search_results_decoy[0]->m_score;
                output_psm->m_scanNum       = search_results_decoy[0]->m_scanNum;
                output_psm->m_annotation    = search_results_decoy[0]->m_annotation;
                output_psm->m_score         = match_score;
                output_psm->m_isDecoy = true;
                
                return 0;
            }
            
            
        }

        return -1;*/
    }
    
    
    /*! \brief Default Spectral Library Search

     */
    
    int SpectralLibrary::search(Spectrum query_spec, psmPtr output_psm, float parentmass_tolerance){
        AAJumps aajumps(1);
        vector<float> masses;
        int charge;
        
        PeptideSpectrumMatchSet search_results;
        vector<score_results_tuple> scores_tuple;
        
        //psmPtr psm(new PeptideSpectrumMatch);
        
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            masses.clear();
            aajumps.getPRMMasses(specs[library_idx].psmList.front()->m_annotation.c_str(), masses);
            charge = specs[library_idx].parentCharge;
            float library_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
            
            if(abs(query_spec.parentMass - library_mass) > parentmass_tolerance) continue;
            
            float sim = full_spectrum_similarity(specs[library_idx], query_spec);
            
            score_results_tuple similarity_tuple;
            
            tr1::get<0>(similarity_tuple) = &specs[library_idx];
            tr1::get<1>(similarity_tuple) = sim;
            tr1::get<2>(similarity_tuple) = library_idx;
            tr1::get<3>(similarity_tuple) = (string)specs[library_idx].psmList.front()->m_annotation;
            scores_tuple.push_back(similarity_tuple);
            
        }
        
        //Finding the best
        //sort(scores.begin(), scores.end(), search_results_comparator);
        sort(scores_tuple.begin(), scores_tuple.end(), search_results_comparator);
        for(int i = 0; i < scores_tuple.size(); i++){
            psmPtr psm(new PeptideSpectrumMatch);
            psm->m_spectrum = tr1::get<0>(scores_tuple[i]);
            psm->m_score = tr1::get<1>(scores_tuple[i]);
            psm->m_scanNum = tr1::get<2>(scores_tuple[i]) + 1;
            psm->m_annotation = tr1::get<3>(scores_tuple[i]);
            search_results.push_back(psm);
        }
        
        if(search_results.size() > 0){
            float match_score = 0.6*search_results[0]->m_score;
            if(search_results.size() == 1) match_score += 0.4 * search_results[0]->m_score;
            else match_score += 0.4 * search_results[1]->m_score;
            
            output_psm->m_spectrum      = search_results[0]->m_spectrum;
            output_psm->m_score         = search_results[0]->m_score;
            output_psm->m_scanNum       = search_results[0]->m_scanNum;
            output_psm->m_annotation    = search_results[0]->m_annotation;
            output_psm->m_score         = match_score;
            
            return 0;
        }
        
        /*for(int score_idx = 0; score_idx < min((int)scores_tuple.size(), 1); score_idx++){
            cout<<tr1::get<0>(scores_tuple[score_idx])->psmList.front()->m_annotation<<"\t"<<tr1::get<1>(scores_tuple[score_idx])<<"\t0BasedIdx\t"<<tr1::get<2>(scores_tuple[score_idx])<<endl;
        }*/
        
        return -1;
    }
    
    /*! \brief Loads an mgf file and adds it to current specset

     */
    unsigned int SpectralLibrary::LoadSpecSet_additionalmgf(const char * filename){
        SpecSet tempspecs;
        tempspecs.LoadSpecSet_mgf(filename);
        
        for(int specidx = 0; specidx < tempspecs.size(); specidx++){
            specs.push_back(tempspecs[specidx]);
            list<psmPtr>::iterator it;
            for ( it=specs[specs.size()-1].psmList.begin() ; it != specs[specs.size()-1].psmList.end(); it++ ){
                (*it)->m_spectrum = &(specs[specs.size()-1]);
            }
        }
        
        return 0;
    }
    
    void SpectralLibrary::filter_no_psm(){
        SpectralLibrary temp_lib;
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            if(specs[spectrum_idx].psmList.size() >= 1){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }
    
    /*! \brief Filters out PSMs that are not the best scoring PValue

     */
    void SpectralLibrary::filter_multiple_interpretations(){
        SpectralLibrary temp_lib;
        //Filtering out spectra with multiple interpretations, or 0 interpretations
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            /*
            if(specs[spectrum_idx].psmList.size() >= 1){
                float best_pvalue = 1.f;
                list<psmPtr>::iterator it;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if((*it)->m_pValue < best_pvalue){
                        best_pvalue = (*it)->m_pValue;
                    }
                }
                
                psmPtr best_psm;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if((*it)->m_pValue == best_pvalue){
                        best_psm = (*it);
                    }
                }
                
                specs[spectrum_idx].psmList.clear();
                specs[spectrum_idx].psmList.push_back(best_psm);
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }*/
            if(specs[spectrum_idx].psmList.size() > 1){
                specs[spectrum_idx].psmList.clear();
            }
            
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }
    
    /*! \brief Filters out spectra that exceed the envelope score

     */
    void SpectralLibrary::filter_envelope_filter(float envelope_score_filter){
        SpectralLibrary temp_lib;
        //Filtering based on envelope scores
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            //cout<<"Envelope: "<<specs[spectrum_idx].psmList.front()->m_strict_envelope_score<<endl;
            
            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_strict_envelope_score > envelope_score_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
            
            /*
            
            if(specs[spectrum_idx].psmList.front()->m_strict_envelope_score < envelope_score_filter){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }*/
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }
    
    /*! \brief Filters out spectra that exceed the pvalue threshold

     */
    void SpectralLibrary::filter_pvalue(float pvalue_filter){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues
        
        
        
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_pValue > pvalue_filter){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
            /*if(specs[spectrum_idx].psmList.front()->m_pValue < pvalue_filter){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }*/
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }
    
    
    void SpectralLibrary::filter_forcharge(vector<int> charge){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            
            if(specs[spectrum_idx].parentCharge == 0){
                specs[spectrum_idx].parentCharge = specs[spectrum_idx].psmList.front()->m_charge;
            }
            
            int spec_charge = specs[spectrum_idx].psmList.front()->m_charge;
            if(spec_charge == -1) spec_charge = specs[spectrum_idx].parentCharge;
            bool valid_charge = false;
            for(int charge_idx = 0; charge_idx < charge.size(); charge_idx++){
                if(spec_charge == charge[charge_idx]){
                    valid_charge = true;
                    break;
                }
            }
            if(valid_charge){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }
    
    void SpectralLibrary::filter_forfragmentation(vector<int> accepted_fragmentation){
        SpectralLibrary temp_lib;
        //Filtering based on fragmentation
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            bool valid_fragmentation = false;
            int spec_fragmentation = specs[spectrum_idx].msFragType;
            for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
                if(spec_fragmentation == accepted_fragmentation[fragmentation_idx]){
                    valid_fragmentation = true;
                    break;
                }
            }
            if(valid_fragmentation){
                temp_lib.specs.push_back(specs[spectrum_idx]);
            }
        }
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }
    
    
    void SpectralLibrary::filter_aminoacids(string aminoacidexclusions){
        SpectralLibrary temp_lib;
        //Filtering based on pvalues
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            //string annotation = specs[spectrum_idx].psmList.front()->m_annotation;
            
            list<psmPtr>::iterator it;
            while(1){
                bool rerun = false;
                for ( it=specs[spectrum_idx].psmList.begin() ; it != specs[spectrum_idx].psmList.end(); it++ ){
                    if( (*it)->m_annotation.find_first_of(aminoacidexclusions) != string::npos){
                        specs[spectrum_idx].psmList.erase(it);
                        rerun = true;
                        break;
                    }
                }
                if(!rerun){
                    break;
                }
            }
            //if(annotation.find_first_of(aminoacidexclusions) == string::npos){
            //    temp_lib.specs.push_back(specs[spectrum_idx]);
            //}
        }
        /*
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }*/
    }
    
    
    /*! \brief Finds spectrum that best represents each peptide

     */
    void SpectralLibrary::filter_best_representative(   MS2ScoringModel &model,
                                                        vector<string> &ionsToExtract,
                                                        string allIons){
        //Grouping up peptides with the same annotation
        map<string, vector<Spectrum *> > same_annotation_clusters;
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            string annotation = specs[spectrum_idx].psmList.front()->m_annotation;
            same_annotation_clusters[annotation].push_back(&specs[spectrum_idx]);
        }
        
        cout<<"Unique Annotations:"<<same_annotation_clusters.size()<<endl;
        
        //Choosing representative peptide
        SpectralLibrary temp_lib;
        map<string, vector<Spectrum *> >::iterator it;
        for ( it=same_annotation_clusters.begin() ; it != same_annotation_clusters.end(); it++ ){
            Spectrum consensus_spectrum;
            get_consensus((*it).second, model, ionsToExtract, allIons, consensus_spectrum);
            temp_lib.specs.push_back(consensus_spectrum);
        }
        
        specs.clear();
        specs.insert(specs.begin(), temp_lib.specs.begin(), temp_lib.specs.end());
        temp_lib.specs.clear();
        
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            specs[spectrum_idx].psmList.front()->m_spectrum = &(specs[spectrum_idx]);
        }
    }
    
    void SpectralLibrary::get_consensus( vector<Spectrum *> spectrum_group, 
                                MS2ScoringModel model, 
                                vector<string> ionsToExtract, 
                                string allIons,
                                Spectrum &consensus_spectrum){
        //Setting psms to be correct
        for(int spec_idx = 0; spec_idx < spectrum_group.size(); spec_idx++){
            spectrum_group[spec_idx]->psmList.front()->m_spectrum = spectrum_group[spec_idx];
        }
        
        
        float best_score = 0.f;
        int best_idx = 0;
        for(int cluster_idx1 = 0; cluster_idx1 < spectrum_group.size(); cluster_idx1++){   
            float spec_score = 0.f;
            for(int cluster_idx2 = 0; cluster_idx2 < spectrum_group.size(); cluster_idx2++){
                if(cluster_idx1 == cluster_idx2) continue;
                spec_score += spectrum_similarity(  spectrum_group[cluster_idx1]->psmList.front(), 
                                                    spectrum_group[cluster_idx2]->psmList.front(),
                                                    spectrum_group[cluster_idx2]->psmList.front()->m_annotation.length() - 4,
                                                    model,
                                                    ionsToExtract, 
                                                    allIons);
            }
            if (spec_score > best_score){
                best_idx = cluster_idx1;
                best_score = spec_score;
            }
        }
        consensus_spectrum = *(spectrum_group[best_idx]);
        
    }
    
    int SpectralLibrary::get_possible_projections(string target_annotation, MS2ScoringModel model, vector<string> ions_to_extract, string allIons, SpecSet &projection_specs, vector<float> &projections_set_cosine_depression){
        for(int spectrum_idx = 0; spectrum_idx < specs.size(); spectrum_idx++){
            string library_annotation = specs[spectrum_idx].psmList.front()->m_annotation;
            if(getStringDifference(target_annotation, library_annotation) == 1){
                Spectrum new_synth_spec;
                
                psmPtr psm = specs[spectrum_idx].psmList.front();
                string source_annotation = psm->m_annotation;
                
                vector<pair<float, float> > ion_mass_intensity_pair_vector;
        
                string annotation = cleanAnnotationEnds(source_annotation);
                psm->annotate(annotation,allIons,model,0,0,.45);
                
                extractIons(psm,annotation.length()-4,model,ions_to_extract,ion_mass_intensity_pair_vector, 0, 0);
                
                int diff_location = getDifferenceIndex(annotation, target_annotation) - 2;
                
                //cout<<"Dif Location: "<<diff_location<<endl;
                //cout<<"Source Annotation: "<<annotation<<endl;
                
                //Calculating the expected cosine depression 
                map<string, float> aa_substitution_map = getAminoAcidSubstitutionLookup();
                //cout<<annotation[diff_location+2]<<"\t"<<target_annotation[diff_location+2]<<endl;
                float cosine_depression = getSubstitutionCosineDepression(annotation[diff_location+2], target_annotation[diff_location+2], aa_substitution_map);
                //cout<<"Cosine Depression: "<<cosine_depression<<endl;
                
                //If we have a cosine depression of less than a value, we skip this projection
                if(cosine_depression < 0.5){
                    //cout<<"Transformation Not Supported"<<endl;
                    //continue;
                }
                
                AAJumps aajumps(1);
                vector<float> masses;
                //cout<<"AAJUMPS annot: "<<annotation<<"\t"<<target_annotation<<endl;
                aajumps.getPRMMasses(annotation.c_str(), masses);
               
                int charge = 2;
                map<char, float> aamasses = getAminoAcidLookup();
                float library_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                //float library_mass = (getMass(annotation, aamasses)+ AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                aajumps.getPRMMasses(target_annotation.c_str(), masses);
                float target_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                //float target_mass = (getMass(target_annotation, aamasses)+ AAJumps::massH2O + AAJumps::massHion*charge)/charge;

                float mass_difference = library_mass - target_mass;
                
                //cout<<"Library Mass: "<<library_mass<<"\t"<<"Target mass: "<<target_mass<<"\t"<<mass_difference<<endl;
                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                    //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                    //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                }
                
                
                for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                    if(ions_to_extract[i/(annotation.length()-4)].find('b') != -1){
                        if( (i%(annotation.length()-4) + 1) < diff_location + 1) //For B ions don't need to edit
                            continue;

                        if(ions_to_extract[i/(annotation.length()-4)].find("++") != -1){
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(annotation.length()-4)].find('y') != -1){
                        if( (i%(annotation.length()-4) + 1) < ((annotation.length()-4) - diff_location) ){ //For Y ions don't need to edit
                            //cout<<"Contineu: "<<(i%(annotation.length()-4) + 1)<<endl;
                            continue;
                        }

                        if(ions_to_extract[i/(annotation.length()-4)].find("++") != -1){
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                        }
                        else{
                            //cout<<ions_to_extract[i/(annotation.length()-4)]<<"\t"<<i%(annotation.length()-4) + 1<<"\t";
                            //cout<<ion_mass_intensity_pair_vector[i].first<<"\t"<<ion_mass_intensity_pair_vector[i].second<<endl;
                            ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                        }
                    }

                    if(ions_to_extract[i/(annotation.length()-4)].find('P') != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                }
                
                new_synth_spec.resize(0);
                sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);
                int new_peaklist_size = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].first > 0.f){
                        new_peaklist_size++;
                        //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                    }
                }

                new_synth_spec.resize(new_peaklist_size);
                int new_peaklist_idx = 0;
                for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                    if(ion_mass_intensity_pair_vector[p].second > 0.f){
                        new_synth_spec[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                        new_synth_spec[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                        new_peaklist_idx++;
                    }
                }
                
                //adding projection to this
                projections_set_cosine_depression.push_back(cosine_depression);
                psmPtr target_psm(new PeptideSpectrumMatch);
                target_psm->m_annotation = target_annotation;
                new_synth_spec.psmList.push_back(target_psm);
                new_synth_spec.parentMZ = target_mass;
                new_synth_spec.parentCharge = specs[spectrum_idx].parentCharge;
                projection_specs.push_back(new_synth_spec);
                target_psm->m_spectrum = &projection_specs[projection_specs.size()-1];
            }
        }
        return 0;
    }
    
    
    vector<vector<float> > SpectralLibrary::find_global_average_spectrum(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        psmPtr psm(new PeptideSpectrumMatch);
        vector<vector<float> > averages;
        
        for(int sizes = 0; sizes < 100; sizes++){
            vector<float> average;
            int psm_count = 0;
            for(int i = 0; i < specs.size(); i++){
                if(specs[i].psmList.front()->m_annotation.length()-4 == sizes){
                    psm_count++;
                    //cout<<"I: "<<i<<endl;
                    //cout<<m_psmSet[i]->m_annotation<<endl;
                    Spectrum tempspec = specs[i];
                    //Making Max to 1000
                    preprocess_spectrum_intensities_max_intensity(&tempspec, 1000.f);
                    //Applying SQRT
                    preprocess_spectrum_intensities(&tempspec, 0, 1);
                    
                    psm->m_annotation = specs[i].psmList.front()->m_annotation;
                    psm->m_spectrum = &tempspec;
                    
                    psm->annotate(psm->m_annotation,allIons,model,0,0,.45);
                    
                    vector<float> ion_mass_intensity_pair_vector;
                    
                    extractIons(psm,psm->m_annotation.length()-4,model,ions_to_extract,ion_mass_intensity_pair_vector, 0, 0);
                    norm_vector(ion_mass_intensity_pair_vector);
                    
                    if(average.size() == 0){
                        for(int j = 0; j < ion_mass_intensity_pair_vector.size(); j++){
                            average.push_back(ion_mass_intensity_pair_vector[j]);
                        }
                        continue;
                    }
                    else{
                        for(int j = 0; j < ion_mass_intensity_pair_vector.size(); j++){
                            average[j] += (ion_mass_intensity_pair_vector[j]);
                        }
                    }
                }
            }
            cout<<sizes<<"\t";
            for(int j = 0; j < average.size(); j++){
                average[j] = average[j]/psm_count;
                cout<<average[j]<<"\t";
            }
            cout<<endl;
            averages.push_back(average);
        }
        
        /*
        average_spectrum->peakList.clear();
        sort( average.begin(), average.end(), mass_intensity_pair_mass_comp);
        int new_peaklist_size = 0;
        for(int p = 0; p < average.size(); p++){
            if(average[p].first > 0.f){
                new_peaklist_size++;
                //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
            }
        }

        average_spectrum->peakList.resize(new_peaklist_size);
        int new_peaklist_idx = 0;
        for(int p = 0; p < average.size(); p++){
            if(average[p].second > 0.f){
                average_spectrum->peakList[new_peaklist_idx][0] = average[p].first;
                average_spectrum->peakList[new_peaklist_idx][1] = average[p].second;
                new_peaklist_idx++;
            }
        }*/
        
        return averages;
    }
    
    SpectralLibrary SpectralLibrary::create_decoy_spectral_library(MS2ScoringModel model, vector<string> ions_to_extract, string allIons){
        unsigned int seed = 0;
        srand(seed);
        
        vector<string> library_peptides;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;
            annotation = create_annotation_ends(annotation);
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            stripped_annotation += ('0' + specs[library_idx].parentCharge);
            library_peptides.push_back(stripped_annotation);
            
        }
        
        sort(library_peptides.begin(), library_peptides.end());
        
        SpectralLibrary decoy_library;
        for(int library_idx = 0; library_idx < specs.size(); library_idx++){
            string annotation = specs[library_idx].psmList.front()->m_annotation;
            
            //Randomize the annotation
            string stripped_annotation = annotation.substr(2, annotation.length() - 4);
            string annotation_orig_stripped = stripped_annotation;
            
            string random_annotation;
            
            bool valid_decoy_found = false;
            
            for(int random_retries = 0 ; random_retries < 50; random_retries++){
                random_annotation = create_decoy_peptide(annotation);
                
                string search_random_annotation = random_annotation;
                search_random_annotation += ('0'  + specs[library_idx].parentCharge);
                
                if(binary_search(library_peptides.begin(), library_peptides.end(), search_random_annotation)){
                    continue;   //try again, already in library
                }
                else{
                    valid_decoy_found = true;
                    break;
                }
            }
            
            if(!valid_decoy_found){
                cout<<"No Valid Decoy Found"<<endl;
                continue;
            }
            
            cout<<decoy_library.size()<<"\t"<<annotation_orig_stripped<<"\t"<<random_annotation<<endl;
            
            int peptide_length = getpeptideLength(annotation_orig_stripped);
            
            //Extracting the ions
            vector<pair <float, float> > ion_mass_intensity_pair_vector;
            specs[library_idx].psmList.front()->annotate(annotation,allIons,model,0,0,0.45);
            extractIons(specs[library_idx].psmList.front(), peptide_length, model, ions_to_extract, ion_mass_intensity_pair_vector, 0, 0);
            
            
            //Collecting the peaks that are not annotated
            vector<pair <float, float> > unannotated_peaks;
            for(int peak_idx = 0; peak_idx < specs[library_idx].size(); peak_idx++){
                bool annotated = false;
                for(int annotated_idx = 0; annotated_idx < ion_mass_intensity_pair_vector.size(); annotated_idx++){
                    if(specs[library_idx][peak_idx][0] == ion_mass_intensity_pair_vector[annotated_idx].first &&
                        specs[library_idx][peak_idx][1] == ion_mass_intensity_pair_vector[annotated_idx].second){
                        annotated = true;
                        break;
                    }
                }
                
                if(!annotated){
                    pair<float, float> unannotated_peak;
                    unannotated_peak.first = specs[library_idx][peak_idx][0];
                    unannotated_peak.second = specs[library_idx][peak_idx][1];
                    unannotated_peaks.push_back(unannotated_peak);
                }
            }
            
            
            vector<string> original_prefix_array;
            vector<string> original_suffix_array;
            vector<string> random_prefix_array;
            vector<string> random_suffix_array;
            generate_prefix_suffix_peptide(annotation_orig_stripped, original_prefix_array, original_suffix_array);
            generate_prefix_suffix_peptide(random_annotation, random_prefix_array, random_suffix_array);
            
            
            AAJumps aajumps(1);
            
            for(int i = 0; i < ion_mass_intensity_pair_vector.size(); i++){
                if( ions_to_extract[i/(peptide_length)].find('b') != -1 || 
                    ions_to_extract[i/(peptide_length)].find('a') != -1){
                    
                    int charge = 2;
                    vector<float> masses;
                    

                    string orig_prefix = original_prefix_array[i % (peptide_length)];
                    string random_prefix = random_prefix_array[i % (peptide_length)];
                    
                    aajumps.getPRMMasses(create_annotation_ends(orig_prefix).c_str(), masses);
                    float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    aajumps.getPRMMasses(create_annotation_ends(random_prefix).c_str(), masses);
                    float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    float mass_difference = original_mass - random_mass;
                    
                    
                    if(ions_to_extract[i/(peptide_length)].find("++") != -1){
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                    else{
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                    }
                }
                
                

                if(ions_to_extract[i/(peptide_length)].find('y') != -1){
                    int charge = 2;
                    vector<float> masses;
                    
                    
                    string orig_suffix = original_suffix_array[i % (peptide_length)];
                    string random_suffix = random_suffix_array[i % (peptide_length)];
                    
                    aajumps.getPRMMasses(create_annotation_ends(orig_suffix).c_str(), masses);
                    float original_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    aajumps.getPRMMasses(create_annotation_ends(random_suffix).c_str(), masses);
                    float random_mass = (masses[masses.size() - 1] + AAJumps::massH2O + AAJumps::massHion*charge)/charge;
                    float mass_difference = original_mass - random_mass;
                    

                    if(ions_to_extract[i/(peptide_length)].find("++") != -1){
                        
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference;
                    }
                    else{
                        ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first - mass_difference*2;
                    }
                }

                if(ions_to_extract[i/(peptide_length)].find('P') != -1){
                    ion_mass_intensity_pair_vector[i].first = ion_mass_intensity_pair_vector[i].first;
                }
            }
            
            
            
            Spectrum new_synth_spec;
        
            //Adding back in the unannotated peaks
            //ion_mass_intensity_pair_vector.clear();
            ion_mass_intensity_pair_vector.insert(ion_mass_intensity_pair_vector.end(), unannotated_peaks.begin(), unannotated_peaks.end());
            
            //Adding these peaks to a spectrum
            new_synth_spec.resize(0);
            sort( ion_mass_intensity_pair_vector.begin(), ion_mass_intensity_pair_vector.end(), mass_intensity_pair_mass_comp);
            int new_peaklist_size = 0;
            for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                if(ion_mass_intensity_pair_vector[p].second > 0.f){
                    new_peaklist_size++;
                    //cout<<ion_mass_intensity_pair_vector[p].first<<" "<<ion_mass_intensity_pair_vector[p].second<<endl;
                }
            }
            

            new_synth_spec.resize(new_peaklist_size);
            int new_peaklist_idx = 0;
            for(int p = 0; p < ion_mass_intensity_pair_vector.size(); p++){
                if(ion_mass_intensity_pair_vector[p].second > 0.f){
                    new_synth_spec[new_peaklist_idx][0] = ion_mass_intensity_pair_vector[p].first;
                    new_synth_spec[new_peaklist_idx][1] = ion_mass_intensity_pair_vector[p].second;
                    new_peaklist_idx++;
                }
            }
            
            
            //adding projection to this
            psmPtr decoy_psm(new PeptideSpectrumMatch);
            decoy_psm->m_annotation = create_annotation_ends(random_annotation);
            decoy_psm->m_charge = specs[library_idx].psmList.front()->m_charge;
            //decoy_psm->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];
            
            new_synth_spec.psmList.push_back(decoy_psm);
            new_synth_spec.parentCharge = specs[library_idx].parentCharge;
            new_synth_spec.parentMass = specs[library_idx].parentMass;
            new_synth_spec.parentMZ = specs[library_idx].parentMZ;
            decoy_library.specs.push_back(new_synth_spec);
            decoy_library.specs[decoy_library.specs.size()-1].psmList.front()->m_spectrum = &decoy_library.specs[decoy_library.specs.size()-1];
            
        }
        
        return decoy_library;
    }
    
    
    string SpectralLibrary::create_decoy_peptide(string peptide){
        vector<string> deliminated_aminoacids;
        string stripped_peptide = remove_annotation_ends(peptide);
        
        //cout<<stripped_peptide<<endl;
        
        for(int pepidx = 0; pepidx < stripped_peptide.length(); pepidx++){
            if(stripped_peptide[pepidx] != '(' && stripped_peptide[pepidx] != '['){
                string temp_str = "";
                temp_str += stripped_peptide[pepidx];
                deliminated_aminoacids.push_back(temp_str);
                continue;
            }
            if(stripped_peptide[pepidx] == '('){
                int pepidx_right = pepidx+1;
                bool found_parenthesis = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ')'){
                        found_parenthesis = true;
                        break;
                    }
                    pepidx_right++;
                }
                
                if(!found_parenthesis){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }
                
                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;
                
                continue;
            }
            
            if(stripped_peptide[pepidx] == '['){
                int pepidx_right = pepidx+1;
                bool found_bracket = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ']'){
                        found_bracket = true;
                        break;
                    }
                    pepidx_right++;
                }
                
                if(!found_bracket){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }
                
                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;
                
                continue;
            }
        }

        vector<string> randomized_decoy_array;
        vector<string> remaining_aa_to_assign;
        vector<int> empty_indices;
        
        //Transferring over immobile amino acids
        for(int i = 0 ; i < deliminated_aminoacids.size(); i++){
            if( deliminated_aminoacids[i].find("R") != string::npos || 
                deliminated_aminoacids[i].find("K") != string::npos || 
                deliminated_aminoacids[i].find("P") != string::npos){
                
                randomized_decoy_array.push_back(deliminated_aminoacids[i]);
                continue;
            }
            else{
                string temp = "";
                randomized_decoy_array.push_back(temp);
            }
        }
        
        
        //Now we clear out R, K, P from original list
        for(int i = 0; i < deliminated_aminoacids.size(); i++){
            if(!( deliminated_aminoacids[i].find("R") != string::npos || 
                deliminated_aminoacids[i].find("K") != string::npos || 
                deliminated_aminoacids[i].find("P") != string::npos)){
                remaining_aa_to_assign.push_back(deliminated_aminoacids[i]);
                empty_indices.push_back(i);
            }
        }
        
        while(empty_indices.size() > 0){
            int rand_idx = rand()%remaining_aa_to_assign.size();
            randomized_decoy_array[empty_indices[0]] = remaining_aa_to_assign[rand_idx];
            empty_indices.erase(empty_indices.begin());
            remaining_aa_to_assign.erase(remaining_aa_to_assign.begin() + rand_idx);
        }
        
        string randomized = "";
        for(int i = 0; i < randomized_decoy_array.size(); i++){
            randomized +=randomized_decoy_array[i];
        }
        
        return randomized;
    }
    
    void SpectralLibrary::generate_prefix_suffix_peptide(string peptide, vector<string> &prefix, vector<string> &suffix){
        vector<string> deliminated_aminoacids;
        string stripped_peptide = remove_annotation_ends(peptide);
        
        //cout<<stripped_peptide<<endl;
        
        for(int pepidx = 0; pepidx < stripped_peptide.length(); pepidx++){
            if(stripped_peptide[pepidx] != '(' && stripped_peptide[pepidx] != '['){
                string temp_str = "";
                temp_str += stripped_peptide[pepidx];
                deliminated_aminoacids.push_back(temp_str);
                continue;
            }
            if(stripped_peptide[pepidx] == '('){
                int pepidx_right = pepidx+1;
                bool found_parenthesis = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ')'){
                        found_parenthesis = true;
                        break;
                    }
                    pepidx_right++;
                }
                
                if(!found_parenthesis){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }
                
                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;
                
                continue;
            }
            
            if(stripped_peptide[pepidx] == '['){
                int pepidx_right = pepidx+1;
                bool found_bracket = false;
                while(pepidx_right < stripped_peptide.length() ){
                    if(stripped_peptide[pepidx_right] == ']'){
                        found_bracket = true;
                        break;
                    }
                    pepidx_right++;
                }
                
                if(!found_bracket){
                    cout<<"Bad Annotation"<<endl;
                    cout<<peptide<<endl;
                    exit(1);
                }
                
                string temp_str = "";
                for(int i = pepidx; i <= pepidx_right; i++){
                    temp_str += stripped_peptide[i];
                }
                deliminated_aminoacids.push_back(temp_str);
                pepidx = pepidx_right;
                
                continue;
            }
        }
        
        prefix.clear();
        suffix.clear();
        for(int i = 0; i < deliminated_aminoacids.size(); i++){
            string prefix_annotation = "";
            string suffix_annotation = "";
            for(int j = 0; j <= i; j++){
                prefix_annotation += deliminated_aminoacids[j];
                suffix_annotation += deliminated_aminoacids[deliminated_aminoacids.size() - i + j - 1];
            }
            prefix.push_back(prefix_annotation);
            suffix.push_back(suffix_annotation);
        }
    }
}

