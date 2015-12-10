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
    submission_date = inputParams.getValue("submission_date");
    
    output_mgf_name = inputParams.getValue("output_mgf");
    
    existing_library_name = inputParams.getValue("existing_library", "");
    
   
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
    SpecSet extracted_specs;
    
    if(existing_library_name.length() > 0){
        cout<<"Loading Existing library"<<endl;
        extracted_specs.LoadSpecSet_mgf(existing_library_name.c_str());
    }
    
    
    cout<<"Opened\t"<<input_spectra_data_file<<endl;
    
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
            
            int dotpos = input_mz_file_name.find_last_of('.');
            string extension = input_mz_file_name.substr(dotpos+1);
            
            SpecSet specs;
            vector<short> mslevel;
            
            if(strcmp(extension.c_str(), "mzXML") == 0 || strcmp(extension.c_str(), "mzxml") == 0){
                LoadMzxml(input_mz_file_name.c_str(), specs, &mslevel, 0);
            }
            else if(strcmp(extension.c_str(), "mgf") == 0){
                specs.LoadSpecSet_mgf(input_mz_file_name.c_str());
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
                    
                    psmPtr psm(new PeptideSpectrumMatch());
                    
                    psm->m_scanNum = extract_scan;
                    psm->m_annotation = spectrum_sequence;
                    psm->m_submission_user = submission_user;
                    psm->m_submission_id = submission_ID;
                    psm->m_submission_date = submission_date;
                    psm->m_organism = organism;
                    psm->m_molecule_mass = molecule_mass;
                    psm->m_compound_name = compound_name;
                    
                    extracted_spec.psmList.push_back(psm);
                    
                    break;
                }
            }
            extracted_specs.push_back(extracted_spec);
        }
    }
    
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
