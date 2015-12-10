// Header Includes
#include "ExecPrmScoring.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"

// External Includes

// System Includes
#include <stdlib.h>

using namespace std;
using namespace specnets;

namespace specnets
{

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(void) :
    ownInput(true), m_inputSpectra(0x0),
    ownOutput(true), m_outputSpectra(0x0)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams) :
    ExecBase(inputParams), ownInput(true), m_inputSpectra(0x0),
    ownOutput(true), m_outputSpectra(0x0)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams,
                           SpecSet * inputSpectra,
                           SpecSet * outputSpectra) :
  ExecBase(inputParams), ownInput(false), m_inputSpectra(inputSpectra),
      ownOutput(false), m_outputSpectra(outputSpectra)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::~ExecPrmScoring(void)
  {
    if (ownInput)
    {
      delete m_inputSpectra;
    }
    if (ownOutput)
    {
      delete m_outputSpectra;
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecPrmScoring::clone(const ParameterList & inputParams) const
  {
    return new ExecPrmScoring(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::invoke(void)
  {
    DEBUG_TRACE;

    DEBUG_VAR(m_inputSpectra->size());

    // Fake high parent charges (4+ and higher) as 3+ for Pepnovo scoring
    vector<int> originalCharges(m_inputSpectra->size());
    for(unsigned int i = 0; i < m_inputSpectra->size(); i++)
    {
    	originalCharges[i] = (*m_inputSpectra)[i].parentCharge;
    	if((*m_inputSpectra)[i].parentCharge > 3)
    	{
    		(*m_inputSpectra)[i].parentCharge = 3;
    	}
    }

    m_inputSpectra->SaveSpecSet_mgf("spectra/specs_ms.mgf");

    // Reset original parent charges
    for(unsigned int i = 0; i < m_inputSpectra->size(); i++)
    {
    	(*m_inputSpectra)[i].parentCharge = originalCharges[i];
    }
    originalCharges.resize(0);
    m_inputSpectra->SaveSpecSet_mgf("spectra/specs_ms_z.mgf");

    string pepnovoMode;
    if (m_params.exists("MIN_SPECTRUM_QUALITY"))
    {
      pepnovoMode = "-min_filter_prob ";
      pepnovoMode += m_params.getValue("MIN_SPECTRUM_QUALITY");
    }
    else
    {
      pepnovoMode = "-no_quality_filter";
    }

    // TODO: Lower case the value?
    if (m_params.getValue("CORRECT_PM") == "yes")
    {
      pepnovoMode += " -correct_pm";
    }
    else
    {
      pepnovoMode += " -use_spectrum_mz";
    }

    // TODO: Lower case the value?
    if (m_params.getValue("GUESS_CHARGE") == "no")
    {
      pepnovoMode += " -use_spectrum_charge";
    }

    if (m_params.exists("PEPNOVO_PTMS"))
    {
      pepnovoMode += " -PTMs M+16:C+57:";
      pepnovoMode += m_params.getValue("PEPNOVO_PTMS");
    }
    else
    {
      pepnovoMode += " -PTMs M+16:C+57";
    }

    DEBUG_VAR(pepnovoMode);

    string exeDir = m_params.getValue("EXE_DIR");
    string pepnovoCmd(exeDir);
    pepnovoCmd += "/PepNovo_bin -prm_norm -model_dir ";
    pepnovoCmd += exeDir;
    pepnovoCmd += "/Models_pepnovo -model ";
    pepnovoCmd += m_params.getValue("PEPNOVO_MODEL");
    pepnovoCmd += " -file spectra/specs_ms.mgf -fragment_tolerance ";
    pepnovoCmd += m_params.getValue("TOLERANCE_PEAK");
    pepnovoCmd += " -digest NON_SPECIFIC ";
    pepnovoCmd += pepnovoMode;
    pepnovoCmd += " > spectra/specs_scored.prms";

    DEBUG_VAR(pepnovoCmd);
    int status = system(pepnovoCmd.c_str());
    if (status != 0)
    {
      string errorString = "Executing ";
      errorString += exeDir;
      errorString += "/Pepnovo_bin!";
      ERROR_MSG(errorString);
      return false;
    }

    DEBUG_TRACE;

    // Initialize output spectra parent mass/charge
    m_outputSpectra->resize(m_inputSpectra->size());
    DEBUG_VAR(m_outputSpectra->size());

    for(unsigned int i = 0; i < m_inputSpectra->size(); i++)
    {
    	(*m_outputSpectra)[i].copyNP((*m_inputSpectra)[i]);
    	(*m_outputSpectra)[i].resize(0);
    }

    DEBUG_TRACE;

    SpecSet * pepnovoSpectra; // PRM spectra scored by Pepnovo
    pepnovoSpectra = new SpecSet(m_inputSpectra->size());
    // Load Pepnovo PRM spectra
    pepnovoSpectra->LoadSpecSet_prmsv3("spectra/specs_scored.prms");

    DEBUG_VAR(pepnovoSpectra->size());

    // Indicates whether Pepnovo should guess parent charges
    bool guessCharge = (m_params.getValue("GUESS_CHARGE") == "yes");
    //bool guessCharge = true;

    for(unsigned int i = 0; i < pepnovoSpectra->size(); i++)
    {
    	unsigned int specIdx = (*pepnovoSpectra)[i].scan;
			unsigned int keptIdx;  // Index of the last kept spectrum peak

    	(*m_outputSpectra)[specIdx] = (*pepnovoSpectra)[i];
    	(*m_outputSpectra)[specIdx].scan = (*m_inputSpectra)[specIdx].scan;
    	if(!guessCharge and (*m_inputSpectra)[specIdx].parentCharge>0)
    	{
    		(*m_outputSpectra)[specIdx].parentCharge = (*m_inputSpectra)[specIdx].parentCharge;
    	}
    	(*m_inputSpectra)[specIdx].parentCharge = (*m_outputSpectra)[specIdx].parentCharge;

    	// Keep only peaks with PRM score > -1
    	keptIdx = 0;
    	for(unsigned int peakIdx=0; peakIdx<(*m_outputSpectra)[specIdx].size(); peakIdx++)
    	{
    		if((*m_outputSpectra)[specIdx][peakIdx][1]>=-1)
    		{
    			(*m_outputSpectra)[specIdx][keptIdx] = (*m_outputSpectra)[specIdx][peakIdx];
    			(*m_outputSpectra)[specIdx][keptIdx++][1] += 1.0;
    		}
    	}
    	(*m_outputSpectra)[specIdx].resize(keptIdx);
    }

    DEBUG_TRACE;

    delete pepnovoSpectra;

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmScoring::loadInputData(void)
  {
    DEBUG_TRACE;
    ownInput = true;
    ownOutput = true;

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::saveInputData(std::vector<std::string> & filenames)
  {
    //SpecSet m_inputSpectra; // the input spectra
    std::string paramFilename = getName() + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmScoring::saveOutputData(void)
  {
    if(m_outputSpectra)
    	m_outputSpectra->SaveSpecSet_pklbin("spectra/specs_scored.pklbin");
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecPrmScoring::split(int numSplit)
  {
	  m_subModules.resize(0);
	  return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("EXE_DIR");
    VALIDATE_PARAM_EXIST("CORRECT_PM");
    VALIDATE_PARAM_EXIST("GUESS_CHARGE");

    m_isValid = true;
    return true;
  }

} // namespace specnets
