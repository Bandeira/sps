//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecSpectralLibrarySearchSLGFProjection.h"
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

  ExecSpectralLibrarySearchSLGFProjection::ExecSpectralLibrarySearchSLGFProjection(void){
    m_name = "ExecSpectralLibrarySearchSLGFProjection";
    m_type = "ExecSpectralLibrarySearchSLGFProjection";
  }


  ExecSpectralLibrarySearchSLGFProjection::ExecSpectralLibrarySearchSLGFProjection(const ParameterList & inputParams){
    m_name = "ExecSpectralLibrarySearchSLGFProjection";
    m_type = "ExecSpectralLibrarySearchSLGFProjection";
    inputParams.print(cout);
    
    this->m_params = inputParams;
    
    //Grabbing the location of the model file
    m_model_file_name = inputParams.getValue("model_file");    
    
    //Annotated MGF files most likely comes from some other spectral library
    //And the annotations are embedded in the file, so don't bother with other
    //annotation file
    string input_annotated_mgf_file = inputParams.getValue("annotatedmgf", inputParams.getValue("EXISTING_LIBRARY_MGF"));
    string input_annotated_mgf_decoy_file = inputParams.getValue("annotatedmgf_decoy", inputParams.getValue("EXISTING_LIBRARY_MGF_DECOY"));

    
    //Characters we are excluding
    aminoacidexclusions = inputParams.getValue("aminoacidexclusion");
    
    //File for amino acid masses
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile");
    
    //Search Spectra
    m_input_projection_target_string = inputParams.getValue("searchspectra");
    
    //Mass Tolerance when searching
    m_search_parentmass_tolerance = inputParams.getValueFloat("search_parentmass_tolerance",
                                                              inputParams.getValueFloat("tolerance.PM_tolerance", 2));
    
    //Specfiying amino acid masses file
    m_aminoacids_masses_file = inputParams.getValue("aminoacidmassfile", "");
    
    //Specifying output psm file
    m_output_psm_filename = inputParams.getValue("output_psm", inputParams.getValue("RESULTS_DIR", ""));
    
    //Specify whether to search the decoy
    m_do_decoy_search = inputParams.getValueInt("search_decoy", 0);
    
    m_score_threshold = inputParams.getValueFloat("score_threshold", 0.7);
    
    m_do_score_threshold = inputParams.getValueInt("use_score_threshold", 0);
    
    m_score_calculation_method = inputParams.getValueInt("search_scoring_method", 2);
    
    m_SLGF_load_mode = inputParams.getValueInt("SLGFLOADMODE", 1);
    
    //Default All Abundance
    m_abundance_mode = inputParams.getValueInt("SEARCHABUNDANCE" , 0);
    
    //Node Idx and count
    m_node_idx = inputParams.getValueInt("NODEIDX", 0);
    m_node_count = inputParams.getValueInt("NODECOUNT", 1);
    
    //Splitting inputs
    stringSplit(input_annotated_mgf_file, m_mgf_file_names);
    stringSplit(input_annotated_mgf_decoy_file, m_mgf_decoy_file_names);
    //stringSplit(m_input_projection_target_string, m_input_projection_target_spectra);
    
    stringSplit(inputParams.getValue("ANNOTATEDMGF_TARGET_ISOCOMBINED"), m_mgf_isocombined_file_names);
    stringSplit(inputParams.getValue("ANNOTATEDMGF_DECOY_ISOCOMBINED"), m_mgf_decoy_isocombined_file_names);
    
    if(inputParams.getValueInt("CCMS", 0) == 0){
        //Search Spectra
        m_input_projection_target_string = inputParams.getValue("searchspectra");
        stringSplit(m_input_projection_target_string, m_input_projection_target_spectra);    
    }
    
    if(inputParams.getValueInt("CCMS", 0) == 1){
        //Reading Mapping
        //Reading in file mapping
        int mapping_count = 0;
        map<string,string> file_name_mapping;
        while(1){
            char buf[100];
            sprintf(buf, "upload_file_mapping%i", mapping_count);
            if(!inputParams.exists(buf))
                break;
            std::string mapping = inputParams.getValue(buf);
            std::string mangled_name = mapping.substr(0, mapping.find("|"));
            std::string original_name = mapping.substr(mapping.find("|")+1);
            file_name_mapping[original_name] = mangled_name;
            mapping_count++;
            std::string path_to_spectra = inputParams.getValue("SPECTRA_DIR", "");
            m_input_projection_target_spectra.push_back(path_to_spectra + "/" + mangled_name);
        }
    }
    
  }


  ExecSpectralLibrarySearchSLGFProjection::~ExecSpectralLibrarySearchSLGFProjection(void){
  }


  ExecBase * ExecSpectralLibrarySearchSLGFProjection::clone(const ParameterList & inputParams) const{
    return new ExecSpectralLibrarySearchSLGFProjection(inputParams);
  }

  bool ExecSpectralLibrarySearchSLGFProjection::invoke(void){
    //Loading SLGF as a binary
    if(m_SLGF_load_mode == 2){
        //Save out file as binary
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
    
    
    int spectra_count = searchable_spectra.size();
    //Calculating ranges
    //spectra_count = 10;
    int chunk_size = spectra_count/m_node_count;
    int start_idx = chunk_size * m_node_idx;
    int end_idx = chunk_size * (m_node_idx+1);
    if(m_node_idx == (m_node_count-1)) end_idx = spectra_count;
    
    DEBUG_MSG(m_node_idx<<"\t"<<m_node_count<<"\t"<<chunk_size);
    DEBUG_MSG("SEARCHING ["<<start_idx<<" TO " <<end_idx<< ")");
    
    m_library.search_target_decoy_specset_SLGFProjection(       m_decoy_library, 
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
                                                                m_score_threshold, 
                                                                m_abundance_mode,
                                                                m_all_search_results,
                                                                start_idx,
                                                                end_idx);

    
    return true;
                                                                
    if(m_do_decoy_search){
        vector<string> fdr_outputs;
        for(float fdr = 0.f; fdr < 0.05; fdr += 0.001){
            ParameterList fdr_params;
            char buf[1000];
            sprintf(buf,"%f", fdr);
            
            fdr_params.setValue("PEPTIDE_FDR_CUTOFF", buf);
            fdr_params.setValue("TDA_TYPE", "concatenated");
            PeptideSpectrumMatchSet fdr_peptides;
            ExecFdrPeptide calculateFDR(fdr_params, &m_all_search_results, &fdr_peptides);
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

  bool ExecSpectralLibrarySearchSLGFProjection::loadInputData(void){      
    
    load_aminoacid_masses();
      
    //Loading models and setting up ions
    model.LoadModel(m_model_file_name.c_str());
    allIons = "all";
    
    
    
    if(m_abundance_mode == 1){
        //Low Abundance
        ionsToExtract.push_back("a");
        ionsToExtract.push_back("b");
        ionsToExtract.push_back("b++");
        ionsToExtract.push_back("y");
        ionsToExtract.push_back("y++");
    }
     
    if(m_abundance_mode == 2 || m_abundance_mode == 0){
        //High Abundance
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
        ionsToExtract.push_back("y++-H2O");
    }
    
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
        else if(strcmp(extension.c_str(), "mgf") == 0){
            cout<<"LOADED Search\t"<<spectra_file_name<<"\t as a mgf file"<<endl;
            
            SpecSet temp_specs;
            temp_specs.LoadSpecSet_mgf(spectra_file_name.c_str());
            for(int i = 0; i < temp_specs.size(); i++){
                temp_specs[i].scan = i+1;
                temp_specs[i].fileName = spectra_file_name;
            }
            searchable_spectra.insert(searchable_spectra.begin(), temp_specs.begin(), temp_specs.end());
            
        }
        else{
            cerr<<"CANNOTLOAD Search\t"<<spectra_file_name<<endl;
        }
    }
    
    
    
    
    return true;
  }


  bool ExecSpectralLibrarySearchSLGFProjection::saveOutputData(void){
    DEBUG_MSG("SAVE OUTPUT"<<m_output_psm_filename);
    if(m_output_psm_filename.length() > 0){
        m_all_search_results.saveToFile(m_output_psm_filename.c_str());
    }
    return true;
  }


  bool ExecSpectralLibrarySearchSLGFProjection::saveInputData(std::vector<std::string> & filenames){
    DEBUG_MSG("SAVEINPUT");
    string dataDir = m_params.getValue("GRID_DATA_DIR", ".");
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();
    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);
    DEBUG_VAR(paramFilename);
    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    
    return true;
  }

  bool ExecSpectralLibrarySearchSLGFProjection::loadOutputData(void){
    DEBUG_MSG("LOADING OUTPUT DATA\t"<<getName()<<"\t"<<m_output_psm_filename);
    m_all_search_results.loadFromFile(m_output_psm_filename.c_str());
    
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectralLibrarySearchSLGFProjection::split(int numSplit){
    DEBUG_MSG("SPLIT\t"<<numSplit<<"==================================================");
    for(int i = 0; i < numSplit; i++){
        ParameterList childParams(m_params);
        
        string dataDir = m_params.getValue("GRID_DATA_DIR", ".");
        string baseDirectory = dataDir + "/";
        string baseFilename = baseDirectory + makeName(getName(), i);
        string psmFilename = baseFilename + ".psms";
        DEBUG_VAR(psmFilename);
        childParams.setValue("output_psm", psmFilename);
        
        stringstream ss1 (stringstream::in | stringstream::out);
        stringstream ss2 (stringstream::in | stringstream::out);
        
        ss1 << i;
        childParams.setValue("NODEIDX", ss1.str());
        
        ss2<<numSplit;
        childParams.setValue("NODECOUNT", ss2.str());
        
        ExecBase * theClone = new ExecSpectralLibrarySearchSLGFProjection(childParams);
        
        theClone->setName(makeName(m_name, i));
        m_subModules.push_back(theClone);
    }
    
    //m_subModules.push_back
    return m_subModules;
  }

  bool ExecSpectralLibrarySearchSLGFProjection::merge(void){
    DEBUG_MSG("MERGING\t"<<m_subModules.size());
    for(int i = 0; i < m_subModules.size(); i++){
        for(int j = 0; j < ((ExecSpectralLibrarySearchSLGFProjection*)m_subModules[i])->m_all_search_results.size(); j++){
            m_all_search_results.push_back(((ExecSpectralLibrarySearchSLGFProjection*)m_subModules[i])->m_all_search_results[j]);
        }
    }
    DEBUG_MSG("MERGED SIZE\t"<<m_all_search_results.size());
    
    return true;
    
    if(m_do_decoy_search){
        vector<string> fdr_outputs;
        for(float fdr = 0.f; fdr < 0.05; fdr += 0.001){
            ParameterList fdr_params;
            char buf[1000];
            sprintf(buf,"%f", fdr);
            
            fdr_params.setValue("PEPTIDE_FDR_CUTOFF", buf);
            fdr_params.setValue("TDA_TYPE", "concatenated");
            PeptideSpectrumMatchSet fdr_peptides;
            ExecFdrPeptide calculateFDR(fdr_params, &m_all_search_results, &fdr_peptides);
            if (!calculateFDR.invoke())
            {
                DEBUG_MSG("Unable to generate fdr results!");
                //return false;
                break;
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


  bool ExecSpectralLibrarySearchSLGFProjection::validateParams(std::string & error){
    return true;
  }
  
  void ExecSpectralLibrarySearchSLGFProjection::load_aminoacid_masses(){
      if(m_aminoacids_masses_file.length() != 0){
          AAJumps jumps(1);
          jumps.loadJumps(m_aminoacids_masses_file.c_str(), true);  //load globally
      }
  }
  
}
