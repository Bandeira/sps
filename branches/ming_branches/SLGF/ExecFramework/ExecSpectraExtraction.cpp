//Module Includes
#include "SpectralLibrary.h"

// Header Include
#include "ExecSpectraExtraction.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include "utils.h"
#include "mzxml.h"
#include <fstream>

using namespace specnets;
using namespace std;


namespace specnets
{

  ExecSpectraExtraction::ExecSpectraExtraction(void){
    m_name = "ExecSpectraExtraction";
    m_type = "ExecSpectraExtraction";
  }


  ExecSpectraExtraction::ExecSpectraExtraction(const ParameterList & inputParams){
    m_name = "ExecSpectraExtraction";
    m_type = "ExecSpectraExtraction";
    inputParams.print(cout);
    
    input_spectra_data_file = inputParams.getValue("inputSpectraDataFile");
    
    submission_user = inputParams.getValue("submission_user");
    submission_ID = inputParams.getValue("submission_ID");
    submission_date = inputParams.getValue("submission_date", get_time_string().c_str());
    
    output_mgf_name = inputParams.getValue("output_mgf", "");
    
    existing_library_name = inputParams.getValue("existing_library", "");
    
    results_dir = inputParams.getValue("results_dir", ".");
    spectra_dir = inputParams.getValue("spectra_dir", "");
    
    //Reading in file mapping
    int mapping_count = 0;
    while(1){
      char buf[100];
      sprintf(buf, "upload_mapping%i", mapping_count);
      if(!inputParams.exists(buf))
        break;
      std::string mapping = inputParams.getValue(buf);
      std::string mangled_name = mapping.substr(0, mapping.find("|"));
      std::string original_name = mapping.substr(mapping.find("|")+1);
      file_name_mapping[original_name] = mangled_name;
      mapping_count++;
    }
  }


  ExecSpectraExtraction::~ExecSpectraExtraction(void){
  }


  ExecBase * ExecSpectraExtraction::clone(const ParameterList & inputParams) const{
    return new ExecSpectraExtraction(inputParams);
  }

  bool ExecSpectraExtraction::invoke(void){
    
    
    
    return true;
  }

  bool ExecSpectraExtraction::loadInputData(void){     
    
    ifstream input_spectra_metadata(input_spectra_data_file.c_str());
    string line;
    SpectralLibrary extracted_specs;
    
    if(existing_library_name.length() > 0){
        cout<<"Loading Existing library"<<endl;
        extracted_specs.LoadSpecSet_mgf(existing_library_name.c_str());
    }
    
    
    cout<<"Opened\t"<<input_spectra_data_file<<endl;
    
    stringstream new_annotations_string (stringstream::in | stringstream::out);
    stringstream all_annotations_string (stringstream::in | stringstream::out);
    
    //Results Headers
    new_annotations_string<<"Filename\t"<<"Scan\t"<<"Peptide\t"<<"Organism\t"<<"Smiles\t"<<"Inchi\t"<<"InchiAux\t"<<"SpectrumQuality"<<endl;
    all_annotations_string<<"Filename\t"<<"Scan\t"<<"Peptide\t"<<"Organism\t"<<"Smiles\t"<<"Inchi\t"<<"InchiAux\t"<<"SpectrumQuality"<<endl;
    
    if(input_spectra_metadata.is_open()){
        getline (input_spectra_metadata ,line);
        cout<<line<<endl;
        while(input_spectra_metadata.good()){
            getline (input_spectra_metadata ,line);
            
            cout<<line<<endl;
            
            if(line.length() < 1)
                continue;   //empty line
            
            stringstream metadata_line_stream;
            metadata_line_stream<<line;
            
            std::string input_mz_file_name;
            std::string spectrum_sequence;
            std::string compound_name;
            float molecule_mass;
            std::string protein;
            std::string organism;
            float ITOL;
            std::string ITOLU;
            float TOL;
            std::string TOLU;
            std::string source_instrument;
            int extract_scan;
            int spectrum_quality;
            std::string smiles;
            std::string InChI;
            std::string InChI_Aux;
            
            metadata_line_stream>>input_mz_file_name;
            metadata_line_stream>>spectrum_sequence;
            metadata_line_stream>>compound_name;
            metadata_line_stream>>molecule_mass;
            metadata_line_stream>>protein;
            metadata_line_stream>>organism;
            metadata_line_stream>>ITOL;
            metadata_line_stream>>ITOLU;
            metadata_line_stream>>TOL;
            metadata_line_stream>>TOLU;
            metadata_line_stream>>source_instrument;
            metadata_line_stream>>extract_scan;
            metadata_line_stream>>spectrum_quality;
            metadata_line_stream>>smiles;
            metadata_line_stream>>InChI;
            metadata_line_stream>>InChI_Aux;
            
            
            int dotpos = input_mz_file_name.find_last_of('.');
            string extension = input_mz_file_name.substr(dotpos+1);
            
            SpecSet specs;
            vector<short> mslevel;
	    
	    string mangled_name;
	    
	    if( file_name_mapping.find(input_mz_file_name) != file_name_mapping.end()){
	      mangled_name = spectra_dir + "/" + file_name_mapping[input_mz_file_name];
	      cout<<"Mangled from "<<input_mz_file_name<<" to "<<mangled_name<<endl;
	    }
	    else{
	      mangled_name = input_mz_file_name;
	    }
            
            if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
                LoadMzxml(mangled_name.c_str(), specs, &mslevel, 0);
            }
            else if(strcmp(extension.c_str(), "mgf") == 0){
                specs.LoadSpecSet_mgf(mangled_name.c_str());
            }
            
            
            Spectrum extracted_spec;
            bool scan_found = false;
            int found_index = -1;
            for(int i = 0; i < specs.size(); i++){
                if(specs[i].scan == extract_scan){
                    scan_found = true;
                    found_index = i;
                    extracted_spec = specs[i];
                    
                    extracted_spec.instrument_name = source_instrument;
                    extracted_spec.ITOL = ITOL;
                    extracted_spec.ITOLU = ITOLU;
                    extracted_spec.TOL = TOL;
                    extracted_spec.TOLU = TOLU;
                    extracted_spec.msLevel = mslevel[i];
                    extracted_spec.spectrum_quality = spectrum_quality;
                    
                    if(input_mz_file_name.find_last_of("/") != string::npos)
                        extracted_spec.fileName = input_mz_file_name.substr(input_mz_file_name.find_last_of("/")+1);
                    
                    psmPtr psm(new PeptideSpectrumMatch());
                    
                    psm->m_scanNum = extract_scan;
                    psm->m_annotation = spectrum_sequence;
                    string submission_metadata = submission_user  + ":"  + submission_ID + ":" + submission_date;
                    psm->m_submission_metadata.push_back(submission_metadata);
                    
                    psm->m_organism.push_back(organism);
                    psm->m_compound_name.push_back(compound_name);
                    psm->m_smiles.push_back(smiles);
                    psm->m_InChI.push_back(InChI);
                    psm->m_InChI_Aux.push_back(InChI_Aux);
                    
                    extracted_spec.psmList.push_back(psm);
                    
                    break;
                }
            }
            
            if(extracted_spec.psmList.size() == 0)
	      continue;
            
            new_annotations_string<<mangled_name<<"\t"<<extracted_spec.scan<<"\t"<<extracted_spec.psmList.front()->m_annotation<<"\t";
	    new_annotations_string<<extracted_spec.psmList.front()->m_organism[extracted_spec.psmList.front()->m_organism.size()]<<"\t";
	    new_annotations_string<<extracted_spec.psmList.front()->m_smiles[extracted_spec.psmList.front()->m_smiles.size()]<<"\t";
	    new_annotations_string<<extracted_spec.psmList.front()->m_InChI[extracted_spec.psmList.front()->m_InChI.size()]<<"\t";
	    new_annotations_string<<extracted_spec.psmList.front()->m_InChI_Aux[extracted_spec.psmList.front()->m_InChI_Aux.size()]<<"\t";
	    new_annotations_string<<extracted_spec.spectrum_quality<<endl;
            
            //Checking if Spectra are already present in library
            extracted_specs.add_update_spectrum_to_Library(extracted_spec);
            //extracted_specs.push_back(extracted_spec);
        }
    }
    
    for(int library_idx = 0; library_idx < extracted_specs.size(); library_idx++){
      all_annotations_string<<extracted_specs[library_idx].fileName<<"\t"<<extracted_specs[library_idx].scan<<"\t"<<extracted_specs[library_idx].psmList.front()->m_annotation<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_organism[extracted_specs[library_idx].psmList.front()->m_organism.size()]<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_smiles[extracted_specs[library_idx].psmList.front()->m_smiles.size()]<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_InChI[extracted_specs[library_idx].psmList.front()->m_InChI.size()]<<"\t";
      all_annotations_string<<extracted_specs[library_idx].psmList.front()->m_InChI_Aux[extracted_specs[library_idx].psmList.front()->m_InChI_Aux.size()]<<"\t";
      all_annotations_string<<extracted_specs[library_idx].spectrum_quality<<endl;
    }
    
    cout<<"Writing\t"<<(results_dir + "/new_annotations.tsv")<<endl;
    cout<<"Writing\t"<<(results_dir + "/all_annotations.tsv")<<endl;
    
    fstream new_file_stream( (results_dir + "/new_annotations.tsv").c_str(), fstream::out);
    fstream all_file_stream( (results_dir + "/all_annotations.tsv").c_str(), fstream::out);
    
    
    new_file_stream<<new_annotations_string.str();
    new_file_stream.close();
    
    all_file_stream<<all_annotations_string.str();
    all_file_stream.close();
    
    extracted_specs.SaveSpecSet_mgf(output_mgf_name.c_str());
    
    return true;
  }


  bool ExecSpectraExtraction::saveOutputData(void){
        return true;
  }


  bool ExecSpectraExtraction::saveInputData(std::vector<std::string> & filenames){
    return true;
  }

  bool ExecSpectraExtraction::loadOutputData(void){
    return true;
  }

  std::vector<ExecBase *> const & ExecSpectraExtraction::split(int numSplit){
    return m_subModules;
  }

  bool ExecSpectraExtraction::merge(void){
    return true;
  }


  bool ExecSpectraExtraction::validateParams(std::string & error){
    return true;
  }
  
}
