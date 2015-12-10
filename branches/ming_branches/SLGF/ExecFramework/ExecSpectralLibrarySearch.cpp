//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecSpectralLibrarySearch.h"
#include "ExecFdrPeptide.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include <omp.h>

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecSpectralLibrarySearch::ExecSpectralLibrarySearch(void){
    m_name = "ExecSpectralLibrarySearch";
    m_type = "ExecSpectralLibrarySearch";
  }


  ExecSpectralLibrarySearch::ExecSpectralLibrarySearch(const ParameterList & inputParams){
    m_name = "ExecSpectralLibrarySearch";
    m_type = "ExecSpectralLibrarySearch";
    inputParams.print(cout);
    
    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file");    
    
    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("annotatedmgf");
    string input_annotated_mgf_decoy_file = inputParams.getValue("annotatedmgf_decoy");

    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");
    
    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
    
    //Search Spectra
    m_input_projection_target_string = inputParams.getValue("searchspectra");
    
    //Mass Tolerance when searching
    m_search_parentmass_tolerance = inputParams.getValueFloat("search_parentmass_tolerance", 20.f);
    
    //Specfiying amino acid masses file
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile", "");
    
    //Specifying output psm file
    m_output_psm_filename = inputParams.getValue("output_psm", "");
    
    //Specify whether to search the decoy
    m_do_decoy_search = inputParams.getValueInt("search_decoy", 0);
    
    m_score_threshold = inputParams.getValueFloat("score_threshold", 0.7);
    
    m_do_score_threshold = inputParams.getValueInt("use_score_threshold", 0);
    
    m_score_calculation_method = inputParams.getValueInt("search_scoring_method", 2);
    
    //Splitting inputs
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
    stringSplit(input_annotated_mgf_decoy_file, m_mgf_decoy_file_names);
    stringSplit(m_input_projection_target_string, m_input_projection_target_spectra);
    
    stringSplit(inputParams.getValue("ANNOTATEDMGF_TARGET_ISOCOMBINED"), m_mgf_isocombined_file_names);
    stringSplit(inputParams.getValue("ANNOTATEDMGF_DECOY_ISOCOMBINED"), m_mgf_decoy_isocombined_file_names);
  }


  ExecSpectralLibrarySearch::~ExecSpectralLibrarySearch(void){
  }


  ExecBase * ExecSpectralLibrarySearch::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibrarySearch(inputParams);
  }

  bool ExecSpectralLibrarySearch::invoke(void){
    //Saving SLGF as a binary
    if(true){
        if(false){
            save_SLGF(m_library, m_mgf_file_names[0] + ".SLGFDIST");
            //save_SLGF(m_library_isocombined, m_mgf_isocombined_file_names[0] + ".SLGFDIST");
            save_SLGF(m_decoy_library, m_mgf_decoy_file_names[0] + ".SLGFDIST");
            //save_SLGF(m_library_isocombined_decoy, m_mgf_decoy_isocombined_file_names[0] + ".SLGFDIST");
            
            exit(0);
        }
            
        load_SLGF(m_library, m_mgf_file_names[0] + ".SLGFDIST");
        load_SLGF(m_decoy_library, m_mgf_decoy_file_names[0] + ".SLGFDIST");
        
        if(m_mgf_isocombined_file_names.size() > 0){
            load_SLGF(m_library_isocombined, m_mgf_isocombined_file_names[0] + ".SLGFDIST");
            load_SLGF(m_library_isocombined_decoy, m_mgf_decoy_isocombined_file_names[0] + ".SLGFDIST");
        }
        
    }
      
      
    cout<<"Invoke"<<endl;
    vector<int> charge_filter;
    
    m_library.createlibrary(-1, 
                          -1, 
                          model, 
                          ionsToExtract, 
                          allIons, 
                          aminoacidexclusions, 
                          charge_filter,
                          false);
                          
    m_library_isocombined.createlibrary(-1, 
                          -1, 
                          model, 
                          ionsToExtract, 
                          allIons, 
                          aminoacidexclusions, 
                          charge_filter,
                          false);
                          
    cout<<"Done Creating Library"<<endl;
    if(m_do_decoy_search){
        //Creating a decoy spectral library
        if(m_mgf_decoy_file_names.size() == 0){
            cout<<"Generating Decoy"<<endl;
            m_decoy_library = m_library.create_decoy_spectral_library(model, ionsToExtract, allIons);
        }
        else{
            m_decoy_library.createlibrary(-1, 
                            -1, 
                            model, 
                            ionsToExtract, 
                            allIons, 
                            aminoacidexclusions, 
                            charge_filter,
                            false);
            m_library_isocombined_decoy.createlibrary(-1, 
                            -1, 
                            model, 
                            ionsToExtract, 
                            allIons, 
                            aminoacidexclusions, 
                            charge_filter,
                            false);
        }
        cout<<"Done Creating Decoy Library"<<endl;
    }
    
    //Removing Low complexity Spectra
    //m_library.remove_low_complexity_spectra(m_library_isocombined, model, ionsToExtract, allIons);
    //m_decoy_library.remove_low_complexity_spectra(m_library_isocombined_decoy, model, ionsToExtract, allIons);
    
    
    //Doing the search
    cout<<"============================================="<<endl;
    cout<<"SEARCH"<<endl;
    vector<psmPtr> all_target_search_results;
    vector<psmPtr> all_decoy_search_results;
    
    cout<<searchable_spectra.size()<<endl;
    
    for(int i = 0; i < m_library.size(); i++){
        m_library[i].psmList.front()->m_dbIndex = i;
    }
    
    for(int i = 0; i < m_decoy_library.size(); i++){
        m_decoy_library[i].psmList.front()->m_dbIndex = i;
    }
    
    //Creating Sorted library spectra pointers
    vector<Spectrum *> target_library_ptr;
    vector<Spectrum *> decoy_library_ptr;
    
    sorted_vector_library(target_library_ptr, m_library);
    sorted_vector_library(decoy_library_ptr, m_decoy_library);

    //High Abundance Iso Combined
    vector<Spectrum *> target_library_isocombined_ptr;
    vector<Spectrum *> decoy_library_isocombined_ptr;
    
    sorted_vector_library(target_library_isocombined_ptr, m_library_isocombined);
    sorted_vector_library(decoy_library_isocombined_ptr, m_library_isocombined_decoy);
    
    
    
    
    vector<int> accepted_fragmentation;
    accepted_fragmentation.push_back(Spectrum::FragType_CID);
        
        
    PeptideSpectrumMatchSet all_search_results;
    
    
    /*

    for(int query_idx = 0; query_idx < searchable_spectra.size(); query_idx++){
        //Filtering in acceptable fragmentation types
        bool valid_fragmentation = false;
        for(int fragmentation_idx = 0; fragmentation_idx < accepted_fragmentation.size(); fragmentation_idx++){
            if(searchable_spectra[query_idx].msFragType == accepted_fragmentation[fragmentation_idx]){
                valid_fragmentation = true;
                break;
            }
        }
        
        if(!valid_fragmentation) continue;
        
        int thread_id = omp_get_thread_num();
            
        if(thread_id == 0){
            cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<query_idx*6<<"\tof\t"<<searchable_spectra.size()<<endl;
        }
        else{
            cout<<"Searching Scan OTHER THREAD:\t"<<searchable_spectra[query_idx].scan<<"\tAND Index\t"<<query_idx<<"\tof\t"<<searchable_spectra.size()<<endl;
        }
        //cout<<"Searching Scan:\t"<<searchable_spectra[query_idx].scan<<"\t";
        //cout<<"mslevel\t"<<searchable_spectra[query_idx].msLevel
        //cout<<searchable_spectra[query_idx].parentMass<<"\t"<<searchable_spectra[query_idx].parentMZ<<"\t"<<searchable_spectra[query_idx].parentCharge<<"\t";
        psmPtr targetdecoy_psm(new PeptideSpectrumMatch);
        int target_decoy_search = m_library.search_target_decoy(m_decoy_library, 
                                                                searchable_spectra[query_idx], 
                                                                targetdecoy_psm, 
                                                                m_search_parentmass_tolerance, 
                                                                target_library_ptr, 
                                                                decoy_library_ptr, 
                                                                m_score_calculation_method,
                                                                model, 
                                                                ionsToExtract,
                                                                allIons);
        
        if(target_decoy_search == 0){
            targetdecoy_psm->m_scanNum = searchable_spectra[query_idx].scan;
            float match_score = targetdecoy_psm->m_score;
            //cout<<"ISDECOY:\t"<<targetdecoy_psm->m_isDecoy<<"\t"<<targetdecoy_psm->m_annotation<<"\t"<<targetdecoy_psm->m_score<<"\t"<<targetdecoy_psm->m_spectrum->scan<<"\t";
            if( ((m_do_score_threshold == 0)) || match_score > m_score_threshold){
                //cout<<match_score<<"\t"<<m_do_score_threshold<<endl;
                //#pragma omp critical
                {
                    all_search_results.push_back(targetdecoy_psm);    
                }
            }
            
        }
        cout<<endl;
        
        continue;
    }*/
    
    if(m_library_isocombined.size() == 0){
    m_library.search_target_decoy_specset(                      m_decoy_library, 
                                                                searchable_spectra,
                                                                m_search_parentmass_tolerance, 
                                                                target_library_ptr, 
                                                                decoy_library_ptr,
                                                                m_score_calculation_method,
                                                                model, 
                                                                ionsToExtract, 
                                                                allIons,
                                                                m_do_score_threshold, 
                                                                m_score_threshold, all_search_results);
    }
    else{
    m_library.search_target_decoy_specset(                      m_decoy_library, 
                                                                searchable_spectra,
                                                                m_search_parentmass_tolerance, 
                                                                target_library_ptr, 
                                                                decoy_library_ptr,
                                                                m_library_isocombined, 
                                                                m_library_isocombined_decoy,
                                                                target_library_isocombined_ptr,
                                                                decoy_library_isocombined_ptr,
                                                                m_score_calculation_method,
                                                                model, 
                                                                ionsToExtract, 
                                                                allIons,
                                                                m_do_score_threshold, 
                                                                m_score_threshold, all_search_results);
    }
    
    if(m_output_psm_filename.length() > 0){
        all_search_results.saveToFile(m_output_psm_filename.c_str());
    }
    
    if(m_do_decoy_search){
        vector<string> fdr_outputs;
        for(float fdr = 0.f; fdr < 0.05; fdr += 0.001){
            ParameterList fdr_params;
            char buf[1000];
            sprintf(buf,"%f", fdr);
            
            fdr_params.setValue("PEPTIDE_FDR_CUTOFF", buf);
            fdr_params.setValue("TDA_TYPE", "concatenated");
            PeptideSpectrumMatchSet fdr_peptides;
            ExecFdrPeptide calculateFDR(fdr_params, &all_search_results, &fdr_peptides);
            if (!calculateFDR.invoke())
            {
                //DEBUG_MSG("Unable to generate fdr results!);
                return false;
            }
            
            stringstream fdr_output (stringstream::in | stringstream::out);
            
            fdr_output<<"FDR\t"<<fdr<<"\t"<<fdr_peptides.size()<<endl;
            fdr_outputs.push_back(fdr_output.str());
        }
        
        for(int i = 0; i < fdr_outputs.size();i++){
            cout<<fdr_outputs[i];
        }
    }
    
    
    
    return true;
  }

  bool ExecSpectralLibrarySearch::loadInputData(void){      
    
    load_aminoacid_masses();
      
    //Loading models and setting up ions
    model.LoadModel(m_model_file_name.c_str());
    allIons = "all";
    /*
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b-NH3");
    ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b-H2O-NH3");
    ionsToExtract.push_back("b-H2O-H2O");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("b++-H2O");
    ionsToExtract.push_back("b++-H2O-H2O");
    ionsToExtract.push_back("b++-NH3");
    ionsToExtract.push_back("b++-iso");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y-NH3");
    ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y-H2O-NH3");
    ionsToExtract.push_back("y-H2O-H2O");
    ionsToExtract.push_back("y++");
    ionsToExtract.push_back("y++-H2O");
    ionsToExtract.push_back("y++-H2O-H2O");
    ionsToExtract.push_back("y++-NH3");
    ionsToExtract.push_back("y++-iso");
    ionsToExtract.push_back("P++");
    ionsToExtract.push_back("P++-H2O");
    ionsToExtract.push_back("P++-NH3");
    ionsToExtract.push_back("P++-H2O-H2O"); 
    */
    
    //Old Set
    /*ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y++");*/
    
    //Low Abundance
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    //ionsToExtract.push_back("b-iso");
    //ionsToExtract.push_back("b-NH3");
    //ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b++");
    //ionsToExtract.push_back("b++-NH3");
    //ionsToExtract.push_back("b++-H2O");
    ionsToExtract.push_back("y");
    //ionsToExtract.push_back("y-iso");
    //ionsToExtract.push_back("y-NH3");
    //ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y++");
    //ionsToExtract.push_back("y++-NH3");
    //ionsToExtract.push_back("y++-H2O");
    
    //High Abundance
    /*
    ionsToExtract.push_back("a");
    ionsToExtract.push_back("b");
    ionsToExtract.push_back("b-iso");
    ionsToExtract.push_back("b-NH3");
    ionsToExtract.push_back("b-H2O");
    ionsToExtract.push_back("b++");
    ionsToExtract.push_back("b++-NH3");
    ionsToExtract.push_back("b++-H2O");
    ionsToExtract.push_back("y");
    ionsToExtract.push_back("y-iso");
    ionsToExtract.push_back("y-NH3");
    ionsToExtract.push_back("y-H2O");
    ionsToExtract.push_back("y++");
    ionsToExtract.push_back("y++-NH3");
    ionsToExtract.push_back("y++-H2O");*/
    
    //Loading annotated MGF
    for(int file_idx = 0; file_idx < m_mgf_file_names.size(); file_idx++){
        cout<<"LOADING TARGET MGF\n"<<endl;
        string mgf_file_name = m_mgf_file_names[file_idx];
        m_library.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }
    
    //Loading annotated MGF for isotopic combined 
    for(int file_idx = 0; file_idx < m_mgf_isocombined_file_names.size(); file_idx++){
        cout<<"LOADING TARGET MGF\n"<<endl;
        string mgf_file_name = m_mgf_isocombined_file_names[file_idx];
        m_library_isocombined.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
        cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf file"<<endl;
    }
    
    if(m_do_decoy_search){
        //Loading annotated MGF
        for(int file_idx = 0; file_idx < m_mgf_decoy_file_names.size(); file_idx++){
            string mgf_file_name = m_mgf_decoy_file_names[file_idx];
            cout<<"LOADING DECOY MGF\t"<<mgf_file_name<<endl;
            m_decoy_library.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
            cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf decoy file"<<endl;
        }
        
        for(int file_idx = 0; file_idx < m_mgf_decoy_isocombined_file_names.size(); file_idx++){
            string mgf_file_name = m_mgf_decoy_isocombined_file_names[file_idx];
            cout<<"LOADING DECOY MGF\t"<<mgf_file_name<<endl;
            m_library_isocombined_decoy.LoadSpecSet_additionalmgf(mgf_file_name.c_str());
            cout<<"LOADED\t"<<mgf_file_name<<"\t as an annotated mgf decoy file"<<endl;
        }
    }

    
    cout<<"Loading Search Spectra"<<endl;
    
    //Loading searchable spectral
    //Looping through mzxml/pklbin files to load
    for(int file_idx = 0; file_idx < m_input_projection_target_spectra.size(); file_idx++){
        string spectra_file_name = m_input_projection_target_spectra[file_idx];
        int dotpos = spectra_file_name.find_last_of('.');
        string extension = spectra_file_name.substr(dotpos+1);
        if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
            SpecSet temp_specs;
            cout<<"LOADED Search\t"<<spectra_file_name<<"\t as an mzxml file"<<endl;
            LoadMzxml(spectra_file_name.c_str(), temp_specs);
            searchable_spectra.insert(searchable_spectra.begin(), temp_specs.begin(), temp_specs.end());
            //projection_spectra.LoadSpecSet_mzxml(spectra_file_name.c_str(), 2);
            //cout<<"MZXML not Supported"<<endl;
        }
        else if(strcmp(extension.c_str(), "pklbin") == 0){
            cout<<"LOADED Search\t"<<spectra_file_name<<"\t as a pklbin file"<<endl;
            searchable_spectra.LoadSpecSet_pklbin(spectra_file_name.c_str());
        }
        else{
            cerr<<"CANNOTLOAD Search\t"<<spectra_file_name<<endl;
        }
    }
    
    
    
    
    return true;
  }


  bool ExecSpectralLibrarySearch::saveOutputData(void){
        return true;
  }


  bool ExecSpectralLibrarySearch::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectralLibrarySearch::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibrarySearch::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectralLibrarySearch::merge(void){
    return true;
  }


  bool ExecSpectralLibrarySearch::validateParams(std::string & error){
    return true;
  }
  
  void ExecSpectralLibrarySearch::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
  
}
